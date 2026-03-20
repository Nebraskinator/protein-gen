#!/usr/bin/env python3
"""
Evolved Structural Mimic / Motif Grafting Pipeline

This script automates the computational grafting of functional motifs (e.g., active sites) 
from a target protein onto an ancestral scaffold. It utilizes ProteinMPNN for sequence 
design and ESM3 for in-silico structural validation, iteratively optimizing for the lowest 
pocket RMSD and highest binding site parity.
"""

import os
import sys
import json
import argparse
import shutil
import subprocess
import numpy as np
import pandas as pd
import torch
from Bio import Align
from Bio.SVDSuperimposer import SVDSuperimposer
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data

# ESM3 Imports
try:
    from esm.models.esm3 import ESM3
    from esm.sdk.api import ESMProtein, GenerationConfig
except ImportError:
    print("Error: ESM3 library not found. Please install via: pip install esm")
    sys.exit(1)

# --- GLOBAL CONSTANTS ---
RESTYPE_3TO1 = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E',
    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F',
    'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'MSE':'M', 'KCX':'K'
}

def get_standard_aligner():
    """Returns a configured Biopython pairwise aligner using standard penalty scores."""
    aligner = Align.PairwiseAligner(mode='global')
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    return aligner

def isolate_chain(input_pdb, output_pdb, target_chain="A"):
    """Extracts a specific chain from a PDB file and writes it to a new file."""
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[21] == target_chain:
                outfile.write(line)
        outfile.write("END\n")

def get_pdb_sequence_and_mapping(pdb_path):
    """
    Parses a PDB file to extract the full sequence and physically present CA atoms.
    Resolves gaps where ProteinMPNN expects residues but physical ATOMs are missing.
    """
    res = parse_PDB(pdb_path)
    parsed = res[0] if isinstance(res, tuple) else res
    mpnn_seq = parsed[0]['seq']
    
    ca_seq, ca_nums = "", []
    seen = set()
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[12:16].strip() == "CA":
                    res_num = int(line[22:26].strip())
                    res_id = (res_num, line[26])
                    if res_id not in seen:
                        seen.add(res_id)
                        ca_seq += RESTYPE_3TO1.get(line[17:20].strip(), 'X')
                        ca_nums.append(res_num)
                        
    if len(mpnn_seq) == len(ca_nums): 
        return mpnn_seq, ca_nums
    
    # Handle missing density/gaps
    aligner = get_standard_aligner()
    safe_mpnn = mpnn_seq.replace('-', 'X')
    m_aln, c_aln = aligner.align(safe_mpnn, ca_seq)[0]
    final_nums, c_ptr = [], 0
    for m, c in zip(m_aln, c_aln):
        if m != '-':
            if c != '-':
                final_nums.append(ca_nums[c_ptr])
                c_ptr += 1
            else:
                final_nums.append(-999) 
    return mpnn_seq, final_nums

def get_ca_coords_and_seq(pdb_path):
    """Extracts sequence and physical 3D coordinates of CA atoms from a PDB."""
    res = parse_PDB(pdb_path)
    parsed = res[0] if isinstance(res, tuple) else res
    mpnn_seq = parsed[0]['seq']
    
    ca_seq, ca_coords = "", []
    seen = set()
    with open(pdb_path, 'r') as f:
        for line in f:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[12:16].strip() == "CA":
                res_num = int(line[22:26].strip())
                res_id = (res_num, line[26])
                if res_id not in seen:
                    seen.add(res_id)
                    ca_seq += RESTYPE_3TO1.get(line[17:20].strip(), 'X')
                    ca_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        
    if len(mpnn_seq) == len(ca_coords): 
        return ca_coords, mpnn_seq
    
    # Align to resolve missing density
    aligner = get_standard_aligner()
    safe_mpnn = mpnn_seq.replace('-', 'X')
    m_aln, c_aln = aligner.align(safe_mpnn, ca_seq)[0]
    final_coords, c_ptr = [], 0
    for m, c in zip(m_aln, c_aln):
        if m != '-':
            if c != '-':
                final_coords.append(ca_coords[c_ptr])
                c_ptr += 1
            else:
                final_coords.append(None) 
    return final_coords, mpnn_seq

def create_alignment_map(target_pdb_seq, ancestor_seq):
    """Creates a dictionary mapping target indices to ancestral amino acids."""
    aligner = get_standard_aligner()    
    
    safe_t = target_pdb_seq.replace('-', 'X')
    safe_a = ancestor_seq.replace('-', 'X')
    
    alignments = aligner.align(safe_t, safe_a)
    t_aln, a_aln = alignments[0][0], alignments[0][1]
    
    mapping = {}
    t_ptr = 0
    for t_char, a_char in zip(t_aln, a_aln):
        if t_char != '-':
            if a_char != '-':
                mapping[t_ptr] = a_char
            t_ptr += 1
    return mapping

def create_bias_jsonl(target_id, target_pdb_seq, target_res_nums, ancestor_map, output_path, weight, lock_map):
    """Generates a JSONL file directing ProteinMPNN to bias towards the ancestor sequence."""
    alphabet = "ACDEFGHIKLMNPQRSTVWYX"
    seq_len = len(target_pdb_seq)
    bias_list = [[0.0] * 21 for _ in range(seq_len)]
    
    for i in range(seq_len):
        current_pos = i + 1 
        if current_pos in lock_map:
            intended_aa = lock_map[current_pos]
            bias_list[i][alphabet.index(intended_aa)] = 100.0
        elif i in ancestor_map:
            anc_aa = ancestor_map[i]
            if anc_aa in alphabet:
                bias_list[i][alphabet.index(anc_aa)] = weight
            
    with open(output_path, 'w') as f:
        json.dump({target_id: {"A": bias_list}}, f)
        f.write("\n")

def run_mpnn(pdb_path, bias_file, fixed_file, out_dir, mpnn_dir, num_seqs=4):
    """Executes the ProteinMPNN model via subprocess to generate candidate sequences."""
    cmd = [
        "python", f"{mpnn_dir}/protein_mpnn_run.py", 
        "--pdb_path", pdb_path, 
        "--pdb_path_chains", "A", 
        "--out_folder", out_dir, 
        "--num_seq_per_target", str(num_seqs), 
        "--sampling_temp", "0.2", 
        "--bias_by_res_jsonl", bias_file, 
        "--fixed_positions_jsonl", fixed_file
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    
    fasta = os.path.join(out_dir, "seqs", f"{os.path.basename(pdb_path)[:-4]}.fa")
    seqs = []
    with open(fasta, 'r') as f:
        lines = f.readlines()
        for i in range(3, len(lines), 2):
            seqs.append(lines[i].strip().replace('X','A'))
    return seqs

def sanitize_and_fold(seqs, model, steps=32):
    """Folds sequences using ESM3 and sanitizes unknown atom records."""
    results = []
    for seq in seqs:
        folded = model.generate(ESMProtein(sequence=seq), GenerationConfig(track="structure", num_steps=steps, temperature=0.7))
        pdb_str = folded.to_pdb_string()
        sanitized = "\n".join([(l[:17]+"ALA"+l[20:] if l.startswith("ATOM") and "UNK" in l[17:20] else l) for l in pdb_str.splitlines()])
        results.append(sanitized)
    return results

def calculate_pocket_rmsd(target_pdb, decoy_pdb, lock_map):
    """Calculates structural RMSD specifically over the locked active site residues."""
    tar_coords, t_seq = get_ca_coords_and_seq(target_pdb)
    decoy_coords, d_seq = get_ca_coords_and_seq(decoy_pdb)

    # Convert 1-based lock_map to 0-based dict for the mapping function
    target_dict = {idx - 1: aa for idx, aa in lock_map.items()}
    t_to_d_map = map_aligned_indices(t_seq, d_seq, target_dict)

    ref_pocket, decoy_pocket = [], []
    for t_ptr, d_ptr in t_to_d_map.items():
        t_c = tar_coords[t_ptr]
        d_c = decoy_coords[d_ptr]
        # Only use residues where CA atoms physically exist in both structures
        if t_c is not None and d_c is not None:
            ref_pocket.append(t_c)
            decoy_pocket.append(d_c)

    if len(ref_pocket) < 3: 
        return 99.9
        
    sup = SVDSuperimposer()
    sup.set(np.array(ref_pocket, dtype=np.float64), np.array(decoy_pocket, dtype=np.float64))
    sup.run()
    return sup.get_rms()

def find_optimal_bridge(tar_pdb, ref_in_silico_pdb, t_seq, t_nums, run_path, target_id, ancestor_map, lock_map, f_file, model, mpnn_dir, num_seqs=128):
    """Iteratively adjusts ancestor bias to find the optimal structural bridge."""
    current_bias = 8.0
    step = 0.5
    best_overall_rmsd = 99.9
    best_overall_seq = ""
    best_overall_bias = current_bias
    
    POCKET_RMSD_THRESHOLD = 0.4
    GLOBAL_TM_SANITY = 0.60
    
    while current_bias >= 1.0:
        print(f"    - Testing Bias: {current_bias:.1f}")
        b_file = os.path.join(run_path, f"bias_{current_bias}.jsonl")
        create_bias_jsonl(target_id, t_seq, t_nums, ancestor_map, b_file, current_bias, lock_map)
        
        cand_seqs = run_mpnn(tar_pdb, b_file, f_file, run_path, mpnn_dir, num_seqs=num_seqs)
        folded_pdbs = sanitize_and_fold(cand_seqs, model, steps=32)        
        
        best_batch_rmsd = 99.9
        best_batch_seq = ""
        
        for i, p_str in enumerate(folded_pdbs):
            tmp = os.path.join(run_path, f"temp_{i}.pdb")
            with open(tmp, 'w') as f: f.write(p_str)
            
            tar_chain = next(get_structure(tar_pdb).get_chains())
            tar_coords_tm, tm_t_seq = get_residue_data(tar_chain)
            
            tmp_chain = next(get_structure(tmp).get_chains())
            tmp_coords_tm, tm_tmp_seq = get_residue_data(tmp_chain)
            
            res = tm_align(tar_coords_tm, tmp_coords_tm, tm_t_seq, tm_tmp_seq)
            tm = res.tm_norm_chain1
            pocket_rmsd = calculate_pocket_rmsd(ref_in_silico_pdb, tmp, lock_map)
            
            if pocket_rmsd < best_batch_rmsd and tm >= GLOBAL_TM_SANITY:
                best_batch_rmsd = pocket_rmsd
                best_batch_seq = cand_seqs[i]
                os.replace(tmp, os.path.join(run_path, f"champion_{current_bias}.pdb"))
            else:
                if os.path.exists(tmp): os.remove(tmp)

        print(f"      Result: Pocket RMSD = {best_batch_rmsd:.2f}A")

        if best_batch_rmsd < best_overall_rmsd:
            best_overall_rmsd = best_batch_rmsd
            best_overall_seq = best_batch_seq
            best_overall_bias = current_bias
            shutil.copy(os.path.join(run_path, f"champion_{current_bias}.pdb"), os.path.join(run_path, "champion.pdb"))

        if best_batch_rmsd <= POCKET_RMSD_THRESHOLD:
            return best_batch_seq, current_bias, best_batch_rmsd
        
        current_bias -= step
        
    return best_overall_seq, best_overall_bias, best_overall_rmsd

def map_aligned_indices(seq_ref, seq_query, ref_dict):
    """
    General function to map 0-based indices from seq_ref to seq_query using alignment.
    ref_dict format: {0_based_index: 'Expected_AA_Letter'}
    """
    aligner = get_standard_aligner()
    
    # Sanitize literal dashes so the aligner doesn't mistake missing atoms for sequence gaps
    safe_ref = seq_ref.replace('-', 'X')
    safe_query = seq_query.replace('-', 'X')
    
    aln = aligner.align(safe_ref, safe_query)[0]
    r_aln, q_aln = aln[0], aln[1]
    
    mapping = {}
    r_ptr, q_ptr = 0, 0
    
    for r_char, q_char in zip(r_aln, q_aln):
        if r_char != '-' and q_char != '-':
            if r_ptr in ref_dict:
                mapping[r_ptr] = q_ptr
        
        # Because we sanitized, missing atoms (X) will correctly increment the pointers
        if r_char != '-': r_ptr += 1
        if q_char != '-': q_ptr += 1
        
    return mapping

def main():
    parser = argparse.ArgumentParser(description="Evolved Structural Mimic / Motif Grafting Pipeline")
    parser.add_argument("--config", type=str, required=True, help="Path to targets.json configuration file")
    parser.add_argument("--seqs", type=int, default=4, help="Number of sequences to generate per step")
    parser.add_argument("--matrix_dir", type=str, default="./superfamily_matrix", help="Output directory for matrix runs")
    parser.add_argument("--mpnn_dir", type=str, default="./ProteinMPNN", help="Path to ProteinMPNN repository")
    parser.add_argument("--esm_model", type=str, default="esm3-sm-open-v1", help="ESM3 model variant to use for folding")
    args = parser.parse_args()

    # Dynamic Module Import Check for ProteinMPNN
    sys.path.append(args.mpnn_dir)
    try:
        global parse_PDB
        from protein_mpnn_utils import parse_PDB
    except ImportError:
        print(f"Error: Could not import protein_mpnn_utils. Ensure ProteinMPNN is located at '{args.mpnn_dir}'")
        sys.exit(1)

    with open(args.config, 'r') as f:
        config = json.load(f)
    os.makedirs(args.matrix_dir, exist_ok=True)
    
    # Initialize ESM3 Engine
    print(f"Initializing ESM3 ({args.esm_model}) on CUDA...")
    model = ESM3.from_pretrained(args.esm_model).to("cuda")
    all_summary = []

    for anc_name, anc_v in config.items():
        for tar_name, tar_v in config.items():
            if anc_name == tar_name: continue
            run_id = f"{anc_name}_to_{tar_name}"
            run_path = os.path.join(args.matrix_dir, run_id)
            os.makedirs(run_path, exist_ok=True)
            
            tar_pdb = os.path.join(run_path, "target.pdb")
            anc_pdb = os.path.join(run_path, "ancestor_clean.pdb")
            isolate_chain(tar_v['pdb'], tar_pdb)
            isolate_chain(anc_v['pdb'], anc_pdb)
            t_seq, t_nums = get_pdb_sequence_and_mapping(tar_pdb)
            a_seq, _ = get_pdb_sequence_and_mapping(anc_pdb)

            ancestor_map = create_alignment_map(t_seq, a_seq)
            lock_map = {} 
            site_tags = tar_v.get('anchor', []) + tar_v.get('motif', [])
            
            print(f"\n[ALIGNMENT LOCKDOWN] Mapping {tar_name}:")
            for tag in site_tags:
                intended_aa = tag[0]
                # Defense: Extract only digits to prevent crashes on insertion codes
                res_num_str = "".join([c for c in tag[1:] if c.isdigit()])
                res_num = int(res_num_str)
                
                if res_num in t_nums:
                    list_idx = t_nums.index(res_num)
                    actual_pdb_aa = t_seq[list_idx]
                    
                    # IDENTITY GUARD: Prevent hallucinated mutations
                    if intended_aa != actual_pdb_aa:
                        print(f"  [WARNING] Mismatch: Tag says {intended_aa}{res_num}, but PDB {tar_name} has {actual_pdb_aa}! Forcing PDB identity.")
                        intended_aa = actual_pdb_aa
                        
                    lock_map[list_idx + 1] = intended_aa
                    print(f"  -> {tag} aligned to Index {list_idx+1}. Identity: {intended_aa}")

            f_list = []
            for list_idx, intended in lock_map.items():
                if t_seq[list_idx - 1] == intended:
                    f_list.append(list_idx)

            f_file = os.path.join(run_path, "fixed.jsonl")
            with open(f_file, 'w') as f:
                json.dump({"target": {"A": f_list}}, f); f.write("\n")
            
            print(f"  [FOLDING] Generating 64-step In Silico Reference for {tar_name}...")
            ref_in_silico_str = sanitize_and_fold([t_seq], model, steps=64)[0]
            ref_in_silico_path = os.path.join(run_path, "target_in_silico.pdb")
            with open(ref_in_silico_path, 'w') as f:
                f.write(ref_in_silico_str)
            
            best_seq, success_bias, final_rmsd = find_optimal_bridge(
                tar_pdb, ref_in_silico_path, t_seq, t_nums, run_path, "target", 
                ancestor_map, lock_map, f_file, model, args.mpnn_dir, num_seqs=args.seqs
            )
            
            id_aligner = get_standard_aligner()

            safe_a = a_seq.replace('-', 'X')
            safe_b = best_seq.replace('-', 'X')
            id_alns = id_aligner.align(safe_a, safe_b)
            match_count = sum(1 for a, b in zip(id_alns[0][0], id_alns[0][1]) if a == b and a != '-' and a != 'X')
            global_identity = (match_count / max(len(a_seq), len(best_seq))) * 100
            
            summary_path = os.path.join(run_path, "summary.txt")
            with open(summary_path, 'w') as f:
                f.write(f"BRIDGE RUN: {anc_name} -> {tar_name}\n")
                f.write("="*60 + "\n")
                f.write(f"Pocket RMSD        : {final_rmsd:.3f} A\n")
                f.write(f"Equilibrium Bias   : {success_bias}\n")
                f.write(f"Ancestor Identity  : {global_identity:.1f}%\n")
                f.write("="*60 + "\n\n")
                
                f.write("BINDING SITE PHYSICAL AUDIT:\n")
                f.write(f"{'Tag':<10} | {'MPNN_Idx':<10} | {'PDB_Num':<10} | {'Target':<10} | {'Graft':<10} | {'Status'}\n")
                f.write("-" * 80 + "\n")
                
                target_dict = {idx - 1: aa for idx, aa in lock_map.items()}
                t_to_b_map = map_aligned_indices(t_seq, best_seq, target_dict)
                
                site_matches = 0
                for curr_mpnn_idx, intended in lock_map.items():
                    t_ptr = curr_mpnn_idx - 1
                    pdb_num = t_nums[t_ptr]
                    t_char = t_seq[t_ptr]
                    
                    if t_ptr in t_to_b_map:
                        b_ptr = t_to_b_map[t_ptr]
                        actual_aa = best_seq[b_ptr]
                    else:
                        actual_aa = "GAP"
                        
                    if actual_aa == intended:
                        status = "MATCH"
                        site_matches += 1
                    else:
                        status = "!! FAIL !!"
                        print(f"      [!] AUDIT FAIL: Tag {intended}{pdb_num} -> Graft generated '{actual_aa}'")
                    
                    f.write(f"{intended+str(pdb_num):<10} | {curr_mpnn_idx:<10} | {pdb_num:<10} | {t_char:<10} | {actual_aa:<10} | {status}\n")
                
                site_parity = (site_matches / len(lock_map)) * 100
                f.write("-" * 80 + "\n")
                f.write(f"TOTAL BINDING SITE PARITY: {site_parity:.1f}%\n\n")
                f.write("="*40 + "\n\n")
                
                f.write(f">Target_Scaffold ({tar_name})\n{t_seq}\n\n")
                f.write(f">Ancestor_Sequence ({anc_name})\n{a_seq}\n\n")
                f.write(f">Graft_Sequence (Bias {success_bias})\n{best_seq}\n")
                
            all_summary.append({
                "Ancestor": anc_name, "Target": tar_name,
                "Final_RMSD": final_rmsd, "Equilibrium_Bias": success_bias
            })
            print(f"  [DONE] {run_id} | RMSD: {final_rmsd:.3f}")

    pd.DataFrame(all_summary).to_csv("evolution_matrix_results.csv", index=False)

if __name__ == "__main__":
    main()