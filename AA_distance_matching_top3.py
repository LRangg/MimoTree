import pandas as pd
import math
import re

#applied on the surface_residue_pdf
#distance between the N atoms of two surface residues
def surface_distance_N(residue_num1, residue_num2, pdb_df):
    residue1 = pdb_df[pdb_df['residue_number'] == residue_num1]
    residue2 = pdb_df[pdb_df['residue_number'] == residue_num2]
    n1 = residue1[residue1['atom_name'] == 'N'].iloc[0].tolist()
    n2 = residue2[residue2['atom_name'] == 'N'].iloc[0].tolist()
    x1, y1, z1 = float(n1[11]), float(n1[12]), float(n1[13])
    x2, y2, z2 = float(n2[11]), float(n2[12]), float(n2[13])
    N_distance = math.sqrt(math.pow(x1-x2, 2) + math.pow(y1-y2, 2) + math.pow(z1-z2, 2))
    return(N_distance)

'''
def remove_unavbl(seeds, mimo_seq):
    final_seeds = []
    for seed in seeds:
        seed_seq = ''
        for s in seed:
            s_char = ''.join(list(filter(str.isalpha, s)))
            seed_seq = seed_seq + s_char
        if seed_seq in mimo_seq:
            final_seeds.append(seed)
    return final_seeds
'''

#seed: one element in final_seeds
#final_seeds: the whole list
#mimo_seq: string
def find_corr_seed(seed, final_seeds, mimo_seq):
    #seed = ['A_64', 'H_63']
    seed_seq = ''
    for i in seed:
        aa_char = ''.join(list(filter(str.isalpha, i)))
        seed_seq = seed_seq + aa_char
    seq_list = mimo_seq.split(seed_seq)
    for j in range(len(seq_list)-1):
        seq_list[j] = seq_list[j] + seed_seq
    f_seq = []
    for k in range(len(seq_list)-1):
        follow_seq = ''.join(seq_list[k+1:len(seq_list)])
        f_seq.append(follow_seq)
    corr_seed = []
    for seq in f_seq:
        for s in final_seeds:
            s_seq = ''
            for ss in s:
                ss_char = ''.join(list(filter(str.isalpha, ss)))
                s_seq = s_seq + ss_char
            if s_seq in seq:
                corr_seed.append(s)
    return corr_seed
#corr_seed: a list of all corresponding seeds

def filter_seeds(seed, corr_seed, mimo_distance_param, mimo_seq, pdb_df):
    #seed = ['A_64', 'H_63']
    seed_seq = ''
    avbl_corr_seed = []
    #mimo_distance_param = 4.32
    for s in seed:
        s_char = ''.join(list(filter(str.isalpha, s)))
        seed_seq = seed_seq + s_char
    seed_end = [substr.end() for substr in re.finditer(seed_seq, mimo_seq)]
    for corseed in corr_seed:
        corseed_seq = ''
        for c in corseed:
            c_char = ''.join(list(filter(str.isalpha, c)))
            corseed_seq = corseed_seq + c_char
        corseed_start = [substr.start() for substr in re.finditer(corseed_seq, mimo_seq)]
        gap_residue = max(corseed_start) - min(seed_end) + 1
        gap_length = gap_residue * mimo_distance_param
        surface_distance = surface_distance_N(int(''.join(list(filter(str.isdigit, seed[-1])))), int(''.join(list(filter(str.isdigit, corseed[0])))), pdb_df)
        if gap_length >= surface_distance:
            avbl_corr_seed.append(corseed)
    return avbl_corr_seed

#candidate_mimo = []
def candidate(mimo_seq,final_seeds,mimo_distance_param,pdb_df,candidate_mimo):
    seed_pairs = {}
    for x in range(len(final_seeds)):
        x_corr_seed = find_corr_seed(final_seeds[x],final_seeds,mimo_seq)
        avbl_corr_seed = filter_seeds(final_seeds[x],x_corr_seed,mimo_distance_param,mimo_seq,pdb_df)
        seed_pairs[x] = avbl_corr_seed
    if any(list(seed_pairs.values())):
        candidate_mimo.append(mimo_seq)
    return candidate_mimo
#if the query mimotope is a candidate

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def epi_seq(final_seeds):
    epi_pred = []
    for x in final_seeds:
        for y in x:
            if y not in epi_pred:
                epi_pred.append(y)
    return epi_pred
#predicted sequence of one mimotope

#seed_dictionary: adjacent seeds indices
#key: index of every seed in seeds
#value: indices of adjacent seeds
def seed_dictionary(seeds, mimo_seq, mimo_distance_param, pdb_df):
    seed_dict = {k:[] for k in range(len(seeds))}
    for seed in seeds:
        corr_seed = find_corr_seed(seed, seeds, mimo_seq)
        avbl = filter_seeds(seed, corr_seed, mimo_distance_param, mimo_seq, pdb_df)
        seed_idx = seeds.index(seed)
        for c_seed in avbl:
            c_seed_idx = seeds.index(c_seed)
            seed_dict[seed_idx].append(c_seed_idx)
    return seed_dict

#recursively find the seed paths starting from one seed
#path:path of the indices of the seeds in a path
def find_seed_paths(graph, start, path):
    path = path + [start]
    if graph[start] == []:
        return [path]
    if start not in graph:
        return []
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = find_seed_paths(graph, node, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths

#apply find_seed_paths and translate the content from indices to amino aicds
def link_seeds(seeds, seed_dict):
    seed_path = []
    for i in range(len(seeds)):
        path = find_seed_paths(seed_dict, i, [])
        for p in path:
            if p not in seed_path:
                seed_path.append(p)
    seed_s = []
    for item in seed_path:
        seed_residue = []
        for idx in item:
            for i in seeds[idx]:
                seed_residue.append(i)
        seed_s.append(seed_residue)
    return seed_s, seed_path

#Return top 5 linked seeds from seed_s having maximum length
def merge_top_seeds(seed_s):
    filtered_seed_path = sorted(seed_s, key=len, reverse=True)[:3]
    epitope_prediction = []
    for seedpath in filtered_seed_path:
        for j in seedpath:
            if j not in epitope_prediction:
                epitope_prediction.append(j)
    return epitope_prediction