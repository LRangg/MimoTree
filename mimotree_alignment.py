import math
PC_properties = {'A':[0.008, 0.134, -0.475, -0.039, 0.181],
                 'R':[0.171, -0.361, 0.107, -0.258, -0.364],
                'N':[0.255, 0.038, 0.117, 0.118, -0.055],
                'D':[0.303, -0.057, -0.014, 0.225, 0.156],
                'C':[-0.132, 0.174, 0.070, 0.565, -0.374],
                'Q':[0.149, -0.184, -0.030, 0.035, -0.112],
                'E':[0.221, -0.280, -0.315, 0.157, 0.303],
                'G':[0.218, 0.562, -0.024, 0.018, 0.106],
                'H':[0.023, -0.177, 0.041, 0.280, -0.021],
                'I':[-0.353, 0.071, -0.088, -0.195, -0.107],
                'L':[-0.267, 0.018, -0.265, -0.274, 0.206],
                'K':[0.243, -0.339, -0.044, -0.325, -0.027],
                'M':[-0.239, -0.141, -0.155, 0.321, 0.077],
                'F':[-0.329, -0.023, 0.072, -0.002, 0.208],
                'P':[0.173, 0.286, 0.407, -0.215, 0.384],
                'S':[0.199, 0.238, -0.015, -0.068, -0.196],
                'T':[0.068, 0.147, -0.015, -0.132, -0.274],
                'W':[-0.296, -0.186, 0.389, 0.083, 0.297],
                'Y':[-0.141, -0.057, 0.425, -0.096, -0.091],
                'V':[-0.274, 0.136, -0.187, -0.196, -0.299],
                'lambda':[1961.504, 788.200, 539.776, 276.624, 244.106]}
AA_abbreviation = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','ASX':'B','CYS':'C','GLU':'E','GLN':'Q','GLX':'Z','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def property_distance(AA1, AA2):
    PD = 0
    AAa1 = AA_abbreviation[AA1[0:3]]
#    AAa2 = AA_abbreviation[AA2[0:3]]
    for i in range(5):
        PD = PD + PC_properties['lambda'][i]*math.pow(PC_properties[AAa1][i]-PC_properties[AA2][i],2)
    return(math.sqrt(PD))

def avg_PD(seed, pdb_df):
    pd = []
    for i in seed:
        residue_n = i[2:]
        AA = pdb_df[pdb_df['residue_number']==int(residue_n)]['residue_name'].tolist()[0]
        pd.append(property_distance(AA, i[0]))
    return sum(pd)/len(pd)

def mimo_avg(input_list, pdb_df):
    if input_list == []:
        return None
    avg_list = []
    for seed in input_list:
        pd = []
        for i in seed:
            residue_n = i[2:]
            AA = pdb_df[pdb_df['residue_number']==int(residue_n)]['residue_name'].tolist()[0]
            pd.append(property_distance(AA, i[0]))
        avg_list.append(sum(pd)/len(pd))
    return sum(avg_list)/len(avg_list)

#PD_cutoff = 8
#import all mimotopes, data includes all mimotope sequences
#mimotope sequences should be splitted by '>'
def import_mimo_data(mimotope_txt):
    with open(mimotope_txt,'r') as mimo:
        data = mimo.read().replace('\n', '').split('>')
        del data[0]
        data = [''.join(list(filter(str.isalpha, x))) for x in data]
    return(data)


#find matches of ONE MIMOTOPE on the surface map
#outputs: matches(dictionary), seed_graph(empty dictionary)
def find_matches(mimo_seq, PD_cutoff, surface_dict):
    matches = {}
    #process only one sequence each time
    for i in range(len(mimo_seq)):
        matches[i] = []
    #find all similar residues of every residue in one mimotope
        for key in surface_dict.keys():
            PD_value = property_distance(key, mimo_seq[i])
            if PD_value <= PD_cutoff:
                matches[i].append(key)
    return(matches)
#MATCHES: key: index of residues in the query mimotope
#value: similar residues on the surface map of the corresponding mimo residue

#applied on MATCHES
def check_keyname(dictionary,AA):
    key_name = []
    for key, value in dictionary.items():
        for i in range(len(value)):
            if value[i] == AA:
                key_name.append(key)
    return(key_name)

def ini_seed_graph(matches, surface_dict):
    keys = list(surface_dict.keys())
    seed_graph = {ks:[] for ks in keys}
    matches_value_list = []
    del_list = []
    for key,value in surface_dict.items():
        key_idx_list = check_keyname(matches, key)
        for i in range(len(value)):
            value_idx_set = set(check_keyname(matches, value[i]))
            key_idx_set = set([x+1 for x in key_idx_list])
            if list(value_idx_set.intersection(key_idx_set)) != []:
                seed_graph[key].append(value[i])
    for val in matches.values():
        matches_value_list = matches_value_list + [x for x in val]
    for ke in seed_graph.keys():
        if ke not in matches_value_list:
            del_list.append(ke)
    for item in del_list:
        seed_graph.pop(item)
    return(seed_graph)
#create seed_graph and remove all unmatched residues

def change_knv_name(temp_matches, seed_graph, mimo_seq):
    temp_seed_graph = seed_graph
    for key, value in temp_matches.items():
        for r in range(len(value)):
            if value[r] in temp_seed_graph.keys():
                temp_seed_graph[mimo_seq[key]+'_'+value[r][3:]] = temp_seed_graph.pop(value[r])
            if value[r] not in temp_seed_graph.keys():
                key_idx = check_keyname(temp_matches, value[r])
                mimoAA = []
                valuer_num = "".join(list(filter(str.isdigit, value[r])))
                for l in range(len(key_idx)):
                    mimoAA.append(mimo_seq[key_idx[l]])
                count_mimoAA = mimoAA.count(mimoAA[0])
                if count_mimoAA != len(mimoAA):
                    mimoAA[0:len(mimoAA)] = [''.join(mimoAA[0:len(mimoAA)])]
                    for key_item in temp_seed_graph.keys():
                        key_num = "".join(list(filter(str.isdigit, key_item)))
                        if key_num == valuer_num:
                            old_key_name = key_item
                    temp_seed_graph[mimoAA[0]+'_'+value[r][3:]] = temp_seed_graph.pop(old_key_name)
    #change key name
    for val in temp_seed_graph.values():
        for i in range(len(val)):
            value_idx = check_keyname(temp_matches, val[i])
            value_mimo = []
            vali_num = "".join(list(filter(str.isdigit, val[i])))
            for j in range(len(value_idx)):
                value_mimo.append(mimo_seq[value_idx[j]])
            count_value_mimo = value_mimo.count(value_mimo[0])
            if count_value_mimo != len(value_mimo):
                value_mimo = list(set(value_mimo))
                value_mimo[0:len(value_mimo)] =  [''.join(value_mimo[0:len(value_mimo)])]
                val[i] = value_mimo[0]+'_'+vali_num
            if count_value_mimo == len(value_mimo):
                val[i] = value_mimo[0]+'_'+vali_num
    #change value name
    return(temp_seed_graph)
#the function only runs once, initiate the seed_graph if the second run is needed
#seed_graph is also changed

def list_idx(lis, char):
    listidx = []
    for i in range(len(lis)):
        if lis[i] == char:
            listidx.append(i)
    return listidx

#separate residues similar to more than one amino acid in the mimotope
def sort_name(seed_graph, mimo_seq):
    sub_seed_graph = {}
    old_key = []
    mimo_seq_list = [aam for aam in mimo_seq]
    for key in seed_graph.keys():
        key_char = ''.join(list(filter(str.isalpha, key)))
        if len(key_char) > 1:
            for aa in key_char:
                new_key = aa+'_'+''.join(list(filter(str.isdigit, key)))
                #print(new_key)
                sub_seed_graph[new_key] = [residue for residue in seed_graph[key]]
            old_key.append(key)
    for item in old_key:
        del seed_graph[item]
    for sub_key in sub_seed_graph.keys():
        seed_graph[sub_key] = [resi for resi in sub_seed_graph[sub_key]]
    #change key
    remove_graph = {k:[] for k in seed_graph.keys()}
    append_graph = {e:[] for e in seed_graph.keys()}
    for ke, value in seed_graph.items():
        for v in value:
            value_char = ''.join(list(filter(str.isalpha, v)))
            if len(value_char) > 1:
                for a in value_char:
                    new_element = a+'_'+''.join(list(filter(str.isdigit, v)))
                    append_graph[ke].append(new_element)
                remove_graph[ke].append(v)
    #make remove_graph & append_graph
    for r_key, remove_value in remove_graph.items():
        for r_v in remove_value:
            #if r_v in seed_graph[r_key]:
            seed_graph[r_key].remove(r_v)
    for a_key, append_value in append_graph.items():
        for a_v in append_value:
            #if a_v not in seed_graph[a_key]:
            seed_graph[a_key].append(a_v)
    #change value
    remove_graph = {kii:[] for kii in seed_graph.keys()}
    for m_key, m_value in seed_graph.items():
        m_k_char = ''.join(list(filter(str.isalpha, m_key)))
        key_idx = list_idx(mimo_seq_list, m_k_char)
        key_idx_plusone = set([ki+1 for ki in key_idx])
        for mv in m_value:
            #print(mv)
            m_v_char = ''.join(list(filter(str.isalpha, mv)))
            value_idx = set(list_idx(mimo_seq_list, m_v_char))
            if list(key_idx_plusone.intersection(value_idx)) == []:
                remove_graph[m_key].append(mv)
                #print(mv)
    for rm_key, rm_value in remove_graph.items():
        for rm_v in rm_value:
            seed_graph[rm_key].remove(rm_v)
    return(seed_graph)

def pre_seed(matches, surface_dict, mimo_seq):
    graph_1 = ini_seed_graph(matches, surface_dict)
    graph_2 = change_knv_name(matches, graph_1, mimo_seq)
    graph_3 = sort_name(graph_2, mimo_seq)
    return graph_3


#search from node0 in seed_graph
def BFS(node0, seed_graph):
    queue, order = [], []
    queue.append(node0)
    order.append(node0)
    while queue:
        v = queue.pop(0)
        for w in seed_graph[v]:
            if w not in order:
                order.append(w)
                queue.append(w)
    return order

def seed_lists_sets(seed_graph):
    sets = []
    seed_lists = []
    for key in seed_graph.keys():
        seed_lists.append(BFS(key, seed_graph))
    for seed_list in seed_lists:
        sets.append(set(seed_list))
    return seed_lists, sets
        
#apply BFS on every key in seed_graph
#make a list of sets

def remove_duplicates(seed_lists, sets):
    res_sets = []
    res_idx = []
    for i in range(len(sets)):
        if sets[i] in res_sets:
            res_idx.append(i)
        if sets[i] not in res_sets:
            res_sets.append(sets[i])
    for idx in sorted(res_idx, reverse=True):
        del seed_lists[idx]
    return seed_lists, res_sets
#remove duplicates in both sets and seed_lists

def remove_subsets(seed_lists, res_sets):
    remove_idx = []
    k = 0
    j = 0
    while j < len(res_sets):
        k = j + 1
        while k < len(res_sets):
            if res_sets[j].issubset(res_sets[k]) == True:
                remove_idx.append(j)
            if res_sets[k].issubset(res_sets[j]) == True:
                remove_idx.append(k)
            k = k + 1
        j = j + 1
    #find the idx of the subsets needed to be removed
    res_remove_idx = []
    for idx in remove_idx:
        if idx not in res_remove_idx:
            res_remove_idx.append(idx)
    #remove duplicates from remove_idx
    res_remove_idx.sort(reverse=True)
    #sort res_remove_idx, values in reverse order
    for indx in res_remove_idx:
        del seed_lists[indx]
    #remove subsets from seed_lists
    return(seed_lists)

def seed_cluster(seed_graph):
    list_1, sets = seed_lists_sets(seed_graph)
    list_2, res_sets = remove_duplicates(list_1, sets)
    list_3 = remove_subsets(list_2, res_sets)
    return list_3

def find_all_paths(graph, start, end, path):
    path = path + [start]
    if start == end:
        return [path]
    if start not in graph:
        return []
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = find_all_paths(graph, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths

def remove_unavbl(seeds, mimo_seq):
    ini_seeds = []
    for seed in seeds:
        seed_seq = ''
        for s in seed:
            s_char = ''.join(list(filter(str.isalpha, s)))
            seed_seq = seed_seq + s_char
        if seed_seq in mimo_seq:
            ini_seeds.append(seed)
    return ini_seeds

def final_seeds(seed_lists, seed_graph, mimo_seq):
    seeds = []
    for seed_cluster in seed_lists:
        for i in range(len(seed_cluster)-1):
            for j in range(i+1, len(seed_cluster)):
                path = []
                paths = find_all_paths(seed_graph, seed_cluster[i], seed_cluster[j], path)
                for p in paths:
                    seeds.append(p)
    seeds = remove_unavbl(seeds, mimo_seq)
    seed_sets = []
    for s in seeds:
        seed_sets.append(set(s))
    seeds, seed_sets = remove_duplicates(seeds, seed_sets)
    seeds = remove_subsets(seeds, seed_sets)
    return seeds

def generate_seeds(mimo_seq, PD_cutoff, surface_dict):    
    matches = find_matches(mimo_seq, PD_cutoff, surface_dict)
    seed_graph = ini_seed_graph(matches, surface_dict)
    seed_graph = change_knv_name(matches, seed_graph, mimo_seq)
    seed_graph = sort_name(seed_graph, mimo_seq)
    seed_lists, sets = seed_lists_sets(seed_graph)
    seed_lists, res_sets = remove_duplicates(seed_lists, sets)
    seed_lists = remove_subsets(seed_lists, res_sets)
    seeds = final_seeds(seed_lists, seed_graph, mimo_seq)
    return seeds