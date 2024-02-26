def top_seeds(seed_s):
    filtered_seed_path = sorted(seed_s, key=len, reverse=True)[:30]
    return filtered_seed_path

#search from node0 in seed_graph
def cluster_BFS(node0, seed_graph):
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

def cluster_remove_duplicates(seed_lists, sets):
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

def cluster_remove_subsets(seed_lists, res_sets):
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

def cluster_lists_sets(seed_graph):
    sets = []
    seed_lists = []
    for key in seed_graph.keys():
        seed_lists.append(cluster_BFS(key, seed_graph))
    for seed_list in seed_lists:
        sets.append(set(seed_list))
    return seed_lists, sets

def sort_by_occurrence(alignments):
    final_alignments = []
    for z in alignments:
        for i in z:
            final_alignments.append(i)
    #final_alignments = sorted(final_alignments, key = final_alignments.count, reverse = True)
    set01 = set(final_alignments)
    dict01 = {item: final_alignments.count(item) for item in set01}
    sorted_x = sorted(dict01.items(), key=lambda a: a[1], reverse=True)
    return dict01, sorted_x

#could be added into a previous function (seeds finalizing)
def sort_seeds(all_seeds):
    all_seeds_list = []
    for g in range(len(all_seeds)):
        seed_g = [r[2:] for r in all_seeds[g]]
        seed_set = set(seed_g)
        if len(seed_set) == len(seed_g):
            all_seeds_list.append(all_seeds[g])
    return all_seeds_list

def initial_clusters(all_seeds_list):
    cluster_graph = {n:[] for n in range(len(all_seeds_list))}
    for i in range(len(all_seeds_list)-1):
        seed_i = [x[2:] for x in all_seeds_list[i]]
        set_i = set(seed_i)
        for j in range(i+1, len(all_seeds_list)):
            seed_j = [y[2:] for y in all_seeds_list[j]]
            set_j = set(seed_j)
            if len(list(set_i.intersection(set_j))) >= 2:
                cluster_graph[i].append(j)
    cluster_lists, cluster_sets = cluster_lists_sets(cluster_graph)
    cluster_lists, res_cluster_sets = cluster_remove_duplicates(cluster_lists, cluster_sets)
    cluster_lists = cluster_remove_subsets(cluster_lists, res_cluster_sets)
    ini_clusters = []
    for l in range(len(cluster_lists)):
        clusterl = []
        for e in cluster_lists[l]:
            for w in all_seeds_list[e]:
                w = w[2:]
                if w not in clusterl:
                    clusterl.append(w)
        ini_clusters.append(clusterl)
    return ini_clusters

def sort_filter_clusters(ini_clusters, all_seeds_list, dict01):
    clusters = []
    for i in range(len(ini_clusters)):
        sorted_cluster = [(icn,dict01[icn]) for icn in ini_clusters[i]]
        sorted_test = sorted(sorted_cluster, key=lambda a: a[1], reverse=True)
        cluster = []
        for a in sorted_test[0:2]:
            for t in all_seeds_list:
                if a[0] in [c[2:] for c in t]:
                    for residue in t:
                        residue = residue[2:]
                        if residue not in cluster:
                            cluster.append(residue)
        clusters.append(cluster)
    return clusters

