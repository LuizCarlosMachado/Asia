import sys
from collections import defaultdict
import pyslim
import tskit
from collections import Counter
import random
from random import randrange
import msprime
import math
from math import isclose
import numpy as np
from scipy.stats import fisher_exact
from math import comb

def calc_tstat(C1,C2,n,m,k1,k2):
    t1_max_likelihood = -(n - k1)*(math.log(C1/(n - k1))) - (n - k1)
    t2_max_likelihood = -(m - k2)*(math.log(C2/(m - k2))) - (m - k2)
    a = -(n + m - k1 - k2)
    b = math.log((C1 + C2)/(n + m - k1 - k2))
    c = (n + m - k1 - k2)
    t1_t2_max_likelihood = (a * b) - c
    T_statistic = (t1_max_likelihood + t2_max_likelihood) - t1_t2_max_likelihood
    return(T_statistic)

def subsample(ts,nsample):
    tree = ts.first()
    d = defaultdict(list)
    sample_ids = list(tree.samples())
    for s in sample_ids:
        ind = ts.node(s).individual
        d[ind].append(s)
    keys = random.sample(list(d),nsample)
    sub_samples = [d[k] for k in keys]
    flat_sub_samples = [item for sublist in sub_samples for item in sublist]
    ts = ts.simplify(samples=flat_sub_samples)
    return ts


def times2intervals(times_dict):
    dictlen = len(times_dict)
    cur_age = 0
    cur_node = 0
    coal_intervals = {}
    for i,(k,v) in enumerate(times_dict.items(),1):
        #print("i, k ,v ", i,k,v)
        if i == 1 and i != dictlen:
            cur_age = v
            cur_node = k
        elif i == dictlen and i == 1:
            #print("Last item")
            #coal_intervals[cur_node] = cur_age - v
            coal_intervals[k] = v
        if i != dictlen:
            #print("not First or last")
            coal_intervals[cur_node] = cur_age - v
            cur_age = v
            cur_node = k
        elif i == dictlen:
            #print("Last item")
            coal_intervals[cur_node] = cur_age - v
            coal_intervals[k] = v
        elif i == 1 and i != dictlen:
            #print("first")
            cur_age = v
            cur_node = k
    return coal_intervals


def log_pvalue(k, fk,N,fN):
    if fN < N and fk < k and fN > 0:
        px  = math.log(math.factorial(N-fN-1)) - math.log(math.factorial(k-fk-1)) - math.log(math.factorial(N-k+fk-fN))
        px += math.log(math.factorial(fN-1)) - math.log(math.factorial(fk-1)) - math.log(math.factorial(fN-fk))
        px -= math.log(math.factorial(N-1)) - math.log(math.factorial(k-1)) - math.log(math.factorial(N-k))
        logp = px
        x = fN - fk
        y = N - k
        c = N - 1
        while x < N-k:
            var  = fk + x
            px  += np.log( (y-x)/(x+1.0) * var / (c - var) )
            logp = np.log( 1.0 + np.exp( px - logp ) ) + logp;
            x += 1
        if logp > 0:
            logp = 0
        return(logp)
    else:
        print("Error with input quantities")

def children_of_node(node,tree):
    all_children = []
    node_list = []
    node_list.append(node)
    current_node = node
    while len(node_list) > 0:
        for n in node_list:
            children = tree.children(n)
            for c in children:
                node_list.append(c)
                all_children.append(c)
            node_list.remove(n)
    return all_children

def tips_of_node(node,tree):
    all_children = []
    node_list = []
    node_list.append(node)
    current_node = node
    while len(node_list) > 0:
        for n in node_list:
            children = tree.children(n)
            for c in children:
                node_list.append(c)
                all_children.append(c)
            node_list.remove(n)
    all_tips = []
    for c in all_children:
        if len(children_of_node(c,tree)) == 0:
           all_tips.append(c) 
    return all_tips

def ts_node_dict_anc(position,ts):
    # get genotypes at position
    for v in ts.variants():
        if v.site.position == position:
            sample_genotypes = list(v.genotypes)
            parent_node = v.site.mutations[0].node
    # get the tree at the selected position
    tree = ts.at(position)
    sel_ids = tips_of_node(parent_node,tree)
    node_age = tree.time(parent_node)
    #print("Pos: ", position, " Age: ", tree.time(parent_node))
    nodes = list(tree.nodes())
    leaves = list(tree.leaves())
    sel_leaves = list(tree.leaves(parent_node))
    anc_leaves = []
    for l in leaves:
        if l not in sel_leaves:
            anc_leaves.append(l)
    #print("len anc_leaves, anc_leaves:", len(anc_leaves),anc_leaves)
    #print("len sel_leaves, sel_leaves:", len(sel_leaves),sel_leaves)
   
    coal_ages = {}
    for n in nodes:
        age = tree.time(n)
        coal_ages[n] = age
    # sort by age
    coal_ages = dict(sorted(coal_ages.items(), key=lambda item: item[1],reverse=True))
    #print("Coal_ages len:",len(coal_ages),coal_ages)
    coal_intervals = times2intervals(coal_ages)
    #print("Coal intervals len: ", len(coal_intervals),coal_intervals)
    # sample IDs must be same as node IDs, which is default
    sample_ids = list(tree.samples())
    #print(sample_ids)
    #print(sample_genotypes)
    #sel_ids = []
    sel_dict = {}
    ind_dict = {}
 
    for s in sample_ids:
        if s in sel_ids:
            genotype = 1
        else:
            genotype = 0
        sel_dict[s] = genotype
        ind = ts.node(s).individual
        if ind in ind_dict:
            ind_dict[ind][s] = genotype
        else:
            ind_dict[ind] = {s:genotype}

    icounter = 0
    for i in ind_dict:
        if 1 in ind_dict[i].values():
            icounter += 1
    #Check if tip node based on len of children list of node 1
    # len(tree.children(1))
    roots = tree.roots
    k = 0
    # all coalescent intervals and ages
    coal_int = {}
    # selected coalescent intervals with all selected child nodes
    sel_int = []
    anc_int = []
    # loop over coal_intervals sorted by age
    for n in coal_intervals:
        #print("interval and age: ", n, tree.time(n))
        children = children_of_node(n,tree)
        tips = tips_of_node(n,tree)
        #print("n, children, tips", n,children,tips)
        if len(children) != 0:
            #age = tree.time(n)
            coal_int[n] = coal_intervals[n]
        # non-tip node with all mutated children
        if len(children) != 0 and set(tips).issubset(sel_ids):
            sel_int.append(n)
            #print("Selected node:",n)
            #if n in roots:
                #k += 1
            #    pass
            #else:
            #    parent = tree.parent(n)
            #    ptips = tips_of_node(parent,tree)
                #if not set(ptips).issubset(sel_ids):
                    #k += 1
                #    pass
        # ancestral node
        elif len(children) != 0 and not set(tips).issubset(sel_ids):
            #if n not in roots:
            #    print("Neutral node:",n, tree.time(n),tree.time(tree.parent(n)))
            #anc_int.append(n)
            # age of the node under which the mutation occured
            pre_mutation_node_age = tree.time(tree.parent(parent_node))
            # Coalescence is younger than derived allele
            #if tree.time(n) < node_age and tree.time(tree.parent(n)) < node_age:
            #    anc_int.append(n)
            #
            if tree.time(n) < pre_mutation_node_age:
                if len(children) != 0:   
                    anc_int.append(n)
                if tree.time(n) < node_age and tree.time(tree.parent(n)) > node_age:
                    k += 1
                    #print("Increment k2 for n, age_n, age_parent_n:", n, tree.time(n),tree.time(tree.parent(n)))
                #anc_int.append(n)
                #if tree.parent(n) in roots:
                #    k += 1
                # count two lineages for each ancestral coalescence contemporaneous with the derived mutation
                #else:
                #    k += 2
                #if tree.parent(n) not in anc_int:
                #    anc_int.append(tree.parent(n))
                #print("Neutral root node:",n)
        elif len(children) == 0 and not set([n]).issubset(sel_ids):
            if tree.time(n) < node_age and tree.time(tree.parent(n)) > node_age:
                k += 1
                #print("Increment k2 for n, age_n, age_parent_n:", n, tree.time(n),tree.time(tree.parent(n)))
            #print("Neutral node:",n, tree.time(n),tree.time(tree.parent(n)))
            #node_age

    anc_ages = {}
    sel_ages = {}
    for key,value in coal_ages.items():
        if key in anc_int:
            anc_ages[key] = value
        elif key in sel_int:
            sel_ages[key] = value
    anc_intervals = times2intervals(anc_ages)
    sel_intervals = times2intervals(sel_ages)
    #print("Converted anc ages to intervals:", anc_ages,anc_intervals)
    #print("anc_intervals:", len(anc_intervals), anc_intervals)    
    #print("coal_intervals:", len(coal_intervals),coal_intervals)
    #print("anc_int",len(anc_int), anc_int)
    #print("sel_int",len(sel_int), sel_int)
    #print("all_int including pre mutation", len(list(coal_intervals.keys())),list(coal_intervals.keys()))
    return [coal_int, sel_ids, k, anc_int, coal_ages,icounter,anc_leaves,sel_leaves,anc_intervals,sel_intervals,sel_int,parent_node]
 
def ts_node_dict_neut(position,ts):
    # get genotypes at position
    parent_node = None
    for v in ts.variants():
        if v.site.position == position:
            sample_genotypes = list(v.genotypes)
            parent_node = v.site.mutations[0].node
    # get the tree at the selected position
    tree = ts.at(position)
    nodes = list(tree.nodes())
    leaves = list(tree.leaves())
    sel_leaves = []
    sel_ids = []
    if parent_node:
        sel_ids = tips_of_node(parent_node,tree)
        node_age = tree.time(parent_node)
        sel_leaves = list(tree.leaves(parent_node))
    anc_leaves = []
    for l in leaves:
        if l not in sel_leaves:
            anc_leaves.append(l)
   
    coal_ages = {}
    for n in nodes:
        age = tree.time(n)
        coal_ages[n] = age
    # sort by age
    coal_ages = dict(sorted(coal_ages.items(), key=lambda item: item[1],reverse=True))
    #print("Coal_ages len:",len(coal_ages),coal_ages)
    coal_intervals = times2intervals(coal_ages)
    #print("Coal intervals len: ", len(coal_intervals),coal_intervals)
    # sample IDs must be same as node IDs, which is default
    sample_ids = list(tree.samples())
    #print(sample_ids)
    #print(sample_genotypes)
    #sel_ids = []
    sel_dict = {}
    ind_dict = {}
 
    #Check if tip node based on len of children list of node 1
    # len(tree.children(1))
    roots = tree.roots
    k = 0
    # all coalescent intervals and ages
    coal_int = {}
    # selected coalescent intervals with all selected child nodes
    sel_int = []
    anc_int = []
    # loop over coal_intervals sorted by age
    for n in coal_intervals:
        children = children_of_node(n,tree)
        tips = tips_of_node(n,tree)
        if len(children) != 0:
            coal_int[n] = coal_intervals[n]
        # non-tip node with all mutated children
        if len(sel_ids) == 0:
            anc_int.append(n)
        elif len(children) != 0 and set(tips).issubset(sel_ids):
            sel_int.append(n)
       # ancestral node
        elif len(children) != 0 and not set(tips).issubset(sel_ids):
            pre_mutation_node_age = tree.time(tree.parent(parent_node))
            if tree.time(n) < pre_mutation_node_age:
                if len(children) != 0:   
                    anc_int.append(n)
                if tree.time(n) < node_age and tree.time(tree.parent(n)) > node_age:
                    k += 1
        elif len(children) == 0 and not set([n]).issubset(sel_ids):
            if tree.time(n) < node_age and tree.time(tree.parent(n)) > node_age:
                k += 1

    anc_ages = {}
    sel_ages = {}
    for key,value in coal_ages.items():
        if key in anc_int:
            anc_ages[key] = value
        elif key in sel_int:
            sel_ages[key] = value
    anc_intervals = times2intervals(anc_ages)
    sel_intervals = times2intervals(sel_ages)
    #print("Converted anc ages to intervals:", anc_ages,anc_intervals)
    #print("Converted sel ages to intervals:", sel_ages,sel_intervals)
    #print("anc_intervals:", len(anc_intervals), anc_intervals)    
    #print("coal_intervals:", len(coal_intervals),coal_intervals)
    return [coal_int, coal_ages, anc_ages, sel_ages, anc_intervals, sel_intervals, sel_int, anc_int, anc_leaves, sel_leaves]
 
def ts_node_dict(position,ts):
    # get genotypes at position
    for v in ts.variants():
        if v.site.position == position:
            sample_genotypes = list(v.genotypes)
            parent_node = v.site.mutations[0].node
    # get the tree at the selected position
    tree = ts.at(position)
    sel_ids = tips_of_node(parent_node,tree)
    #print("sel_ids:",len(sel_ids),sel_ids)
    node_age = tree.time(parent_node)
    #print("Pos: ", position, " Age: ", tree.time(parent_node))
    nodes = list(tree.nodes())
    coal_ages = {}
    for n in nodes:
        age = tree.time(n)
        coal_ages[n] = age
    # sort by age
    coal_ages = dict(sorted(coal_ages.items(), key=lambda item: item[1],reverse=True))
    dictlen = len(coal_ages)
    cur_age = 0
    cur_node = 0
    coal_intervals = times2intervals(coal_ages)
    #print("coal_ages:",coal_ages)
#    for i,(k,v) in enumerate(coal_ages.items(),1):
#        if cur_age > 0 and i != dictlen:
#            #print("not First or last")
#            coal_intervals[cur_node] = cur_age - v
#            cur_age = v
#            cur_node = k
#        elif i == dictlen and i == 1:
#            #print("Last item")
#            #coal_intervals[cur_node] = cur_age - v
#            coal_intervals[k] = v
#        elif i == dictlen:
#            #print("Last item")
#            coal_intervals[cur_node] = cur_age - v
#            coal_intervals[k] = v
#        else:
#            #print("first")
#            cur_age = v
#            cur_node = k
    #print("coal_intervals:",coal_intervals)
    # sample IDs must be same as node IDs, which is default
    sample_ids = list(tree.samples())
    #print(sample_ids)
    #print(sample_genotypes)
    #sel_ids = []
    sel_dict = {}
    ind_dict = {}

 
    for s in sample_ids:
        if s in sel_ids:
            genotype = 1
        else:
            genotype = 0
        sel_dict[s] = genotype
        ind = ts.node(s).individual
        if ind in ind_dict:
            ind_dict[ind][s] = genotype
        else:
            ind_dict[ind] = {s:genotype}

    #for idx,g in enumerate(sample_genotypes):
    #    if g == 1:
    #        genotype = 1
    #        sel_ids.append(sample_ids[idx])
    #    else:
    #        genotype = 0
    #    sel_dict[sample_ids[idx]] = genotype
    #    ind = ts.node(sample_ids[idx]).individual
    #    if ind in ind_dict:
    #        ind_dict[ind][sample_ids[idx]] = genotype
    #    else:
    #        ind_dict[ind] = {sample_ids[idx]:genotype}
    
    #print(selnode)
    #sn_children = children_of_node(selnode,tree)
    #counter = 0
    #for c in sn_children:
    #    if c in sel_ids:
    #        counter += 1
    icounter = 0
    for i in ind_dict:
        if 1 in ind_dict[i].values():
            icounter += 1
    # A test suggests that the mutation node is not the parent of all of the variants, possibly due to recombination
    #print(counter, " of ",len(sel_ids)," selected genomes were found in the selected parent node")
    #print(icounter, " of ",len(ind_dict)," selected individuals were found")
    # Note that the TreeSeq includes each subgenome of a diploid as a separate sample
    # The number of samples with a mutation thus reflects the number of subgenoms that have it
    
    
    # When trees do not converge, or when mutated samples are scattered across the tree
    # we need to extract a list of subtrees based on the list of mutated genomes
    #for r in tree.roots:
    #    print(len(children_of_node(r,tree)))
    
    # Iterate over every non-tip node in one tree sequence
    # for each node, if all its children are mutant
    # add the node and its age to a dict
    # We also need to calculate k, the total number of independent lineages for the mutation
    # if all children of node are mutant but this is not true for parent node OR there is no parent node
    # then increment k by 1
    
    #Check if tip node based on len of children list of node 1
    # len(tree.children(1))
    roots = tree.roots
    k = 0
    # all coalescent intervals and ages
    coal_int = {}
    # selected coalescent intervals with all selected child nodes
    sel_int = []
    for n in coal_intervals:
        children = children_of_node(n,tree)
        tips = tips_of_node(n,tree)
        if len(children) != 0:
            #age = tree.time(n)
            coal_int[n] = coal_intervals[n]
        # non-tip node with all mutated children
        if len(children) != 0 and set(tips).issubset(sel_ids):
            sel_int.append(n)
            if n in roots:
                k += 1
            else:
                parent = tree.parent(n)
                ptips = tips_of_node(parent,tree)
                if not set(ptips).issubset(sel_ids):
                    k += 1
    #coal_int = dict(sorted(coal_int.items(), key=lambda item: item[1],reverse=True))
    #print("Coal_int:",coal_int)
    #print("Nodes: ", len(nodes))
    #print("Sel IDs" , len(sel_ids))
    #print(len(coal_int))
    #print(len(sel_int))
    #print(k)
    sel_ages = {}
    for key,value in coal_ages.items():
        if key in sel_int:
            sel_ages[key] = value
    sel_intervals = times2intervals(sel_ages)
    #print("sel_intervals:", sel_intervals)    

    return [coal_int, sel_ids, k, sel_int, coal_ages,icounter,sel_intervals]
   
   
def twosite_tstat(ts,sel_pos,ts2,neut_pos):
    #print("Analysing two sites!")
    neut_res = ts_node_dict_neut(neut_pos,ts2)
    anc_res = ts_node_dict_anc(sel_pos,ts)
    #print("coal_int, sel_ids, k, anc_int, node_age,icounter,anc_leaves,sel_leaves")
    #print("Derived tree:", anc_res)
    #[coal_int, coal_ages, anc_ages, sel_ages, anc_intervals, sel_intervals, sel_int, anc_int, anc_leaves, sel_leaves]
    #print("Neutral tree:", neut_res)
    t1_children = anc_res[10]
    tree1 = ts.at(sel_pos)
    # age of parent node for mutation
    t1_node_age = tree1.time(anc_res[11])

    #t2_children = neut_res[10]
    tree2 = ts2.at(neut_pos)
    #t2_node_age = neut_res[4]

    # numbers correspond to chromosomes NOT diploid individuals
    # derived extant tips    
    n = len(list(tree1.leaves()))
    # ancestral extant tips
    m = len(list(tree2.leaves()))
    #m = len(rand_neut_res[1])
    # subsampling can sometimes lead to loss of all lineages with the mutation
    if n == 0 or m == 0:
        return 

    #n = len(td1)
    #m = len(td2)
    
    # for soft sweeps, k1 will have to be calculated similarly to k2
    k1 = 1
    k2 = 1
    #k2 = rand_neut_res[2]
    
    #if k1 == 0:
    #    k1 = 1
    if k2 == 0:
        k2 = 1
    # Carry out Fisher's exact test
    twobytwo = np.array([[k1, k2], [n-k1, m-k2]])
    oddsr, fet_p = fisher_exact(twobytwo, alternative='two-sided')
    
    C1a = 0
    C1b = 0
    # subtree count vars
    der_intervals = anc_res[9]
    neut_intervals = neut_res[5]
    t1_count = 2
    t2_count = 2

    C1 = 0
    C2 = 0
    C1_full = 0
    C2_full = 0
   
    C_der = 0
    C_anc = 0
    C_anc_time = 0
    Csub_der = 0
    Csub_anc = 0
    Csub_anc_time = 0
    time_k = None
    C_all = 0
    for i, (key,value) in enumerate(anc_res[0].items(),1):
        C1_full += math.comb(i,2) * value
        if key in der_intervals:
            C_der += math.comb(i,2) * value
            Csub_der += math.comb(t1_count,2) * der_intervals[key]
            t1_count += 1
    #print("Selected mutation age:", t1_node_age)
    #print("Neutral coal ages: ", neut_res[1])
    #print("Neutral interval number and content: ",len(neut_intervals),neut_intervals)
    for i, (key,value) in enumerate(neut_res[0].items(),1):
        #print("i,key,value,sub_i:",i,key,value,t2_count)
        C2_full += math.comb(i,2) * value
        if key in neut_intervals:
            #print("Key in neut_intervals:", key, "with interval value: ", neut_intervals[key])
            C_anc += math.comb(i,2) * value
            Csub_anc += math.comb(t2_count,2) * neut_intervals[key]
            t2_count += 1

        if tree2.time(key) < t1_node_age and key in neut_intervals:
            #print("Key is in neut_intervals and younger than selected allele:", key, "with interval:",neut_intervals[key])
            if not time_k:
                # number of lineages present in neutral when selected var arose
                time_k = t2_count - 2
            C_anc_time += math.comb(i,2) * value
            Csub_anc_time += math.comb(t2_count,2) * neut_intervals[key]

    N_der_full = C1_full / (2*(n-k1))
    N_neut_full = C2_full / (2*(m-k2))
    N_der = C_der / (2*(n-k1))

    T_full = calc_tstat(C1_full,C2_full,n,m,k1,k2)
    try:
        T_part = calc_tstat(C_der,C_anc,n,m,k1,k2)
        T_part_time = calc_tstat(C_der,C_anc_time,n,m,k1,time_k)
        T_part_sub = calc_tstat(Csub_der,Csub_anc,n,m,k1,k2)
        T_part_time_sub = calc_tstat(Csub_der,Csub_anc_time,n,m,k1,time_k)

        N_neut = C_anc / (2*(m-k2))
        N_der_sub = Csub_der / (2*(n-k1))
        N_neut_sub = Csub_anc / (2*(m-k2))
        N_neut_time_sub = Csub_anc_time / (2*(m-time_k))

    except:
        print("Failed to calculate all T statistics - ensure both sites contain variants!")
        T_part = T_part_time = T_part_sub = T_part_time_sub = "NA"
        N_neut = N_der_sub = N_neut_sub = N_neut_time_sub = "NA"        
    #N = number of lineages in total at the final generation
    #fN = frequency of the derived allele at the final generation
    #k = number of lineages in total at generation N
    #fk = frequency of derived allele at generation N
    #logF = vector of precomputed log factorial values, e.g. logF[2] = log(2!)
    k = k1 + k2
    k3 = time_k
    fk = k1
    N = n + m
    fN = n
    #print("k,fk,N,fN:",k,fk,N,fN)
    relate_logp = log_pvalue(k,fk,N,fN)
    if time_k:
        relate_logp_post_sel = log_pvalue(k3,fk,N,fN)
    else:
        relate_logp_post_sel = "NA"
        time_k = "NA"
    #print("T_full,T_part,T_part_sub,T_part_time,T_part_time_sub:",T_full,T_part,T_part_sub,T_part_time,T_part_time_sub)
    outdict = {'position1': sel_pos, 
               'position2': neut_pos, 
               'OddsR': oddsr, 
               'fet_p': fet_p,
               'relate_logp': relate_logp,
               'relate_logp_post_sel': relate_logp_post_sel,
               'n': n,
               'k1': k1,
               'k2': k2,
               'm': m,
               'k3': time_k,
               'Tstat_full': T_full,
               'Tstat_derived': T_part,
               'Tstat_derived_sub': T_part_sub,
               'Tstat_derived_post_sel': T_part_time,
               'Tstat_derived_post_sel_sub': T_part_time_sub,
               'N_der_full':N_der_full,
               'N_neut_full':N_neut_full,
               'N_der':N_der,
               'N_neut':N_neut,
               'N_der_sub':N_der_sub,
               'N_neut_sub':N_neut_sub,
               'N_neut_sub_post_sel':N_neut_time_sub}

    #print(outdict)
    print(*list(outdict.keys()),sep="\t")
    print(*list(outdict.values()),sep="\t")

    #print("Coalescence intervals for derived subtree:",der_intervals)
    #print("Terminal DAF: ", n/(n+m))
    #print("Coalescence intervals for derived subtree list:",der_intervals.values())
    #print("Coalescence intervals for ancestral subtree:",anc_intervals)
    #print("Coalescence intervals for derived subtree:",der_intervals)
    #print("Coalescence ages for tree:",sorted(list(set((int(k) for k in t1_node_age.values())))))
 

def onesite_tstat(ts,sel_pos):
    
    anc_res = ts_node_dict_anc(sel_pos,ts)
    t1_children = anc_res[10]
    t2_children = anc_res[3]
 
    tree1 = ts.at(sel_pos)
    tree2 = ts.at(sel_pos)

    t1_node_age = anc_res[4]
    t2_node_age = anc_res[4]

    # numbers correspond to chromosomes NOT diploid individuals
    # derived extant tips    
    n = len(anc_res[7])
    # ancestral extant tips
    m = len(anc_res[6])
    #m = len(rand_neut_res[1])
    coal_ages = anc_res[0]

    der_nodes = anc_res[10]

    # Count coals post derived mutation
    der_coals = 0
    anc_coals = 0
    for coal_node,coal_age in coal_ages.items():
        if coal_node in der_nodes:
            der_coals += 1
        elif coal_age < tree1.time(anc_res[11]):
            anc_coals += 1
    #coals_der = len(anc_res[8])
    #coals_anc  = len(anc_res[9])

    # subsampling can sometimes lead to loss of all lineages with the mutation
    if n == 0 or m == 0:
        return 

    # for soft sweeps, k1 will have to be calculated similarly to k2
    k1 = 1
    k2 = anc_res[2]
    #k2 = rand_neut_res[2]
    
    #if k1 == 0:
    #    k1 = 1
    if k2 == 0:
        k2 = 1
    # Carry out Fisher's exact test
    twobytwo = np.array([[k1, k2], [n-k1, m-k2]])
    oddsr, fet_p = fisher_exact(twobytwo, alternative='two-sided')
    
    C1a = 0
    C1b = 0
    # subtree count vars
    der_intervals = anc_res[9]
    anc_intervals = anc_res[8]
    t1_count = 2
    t2_count = k2
    sub_C1a = 0
    sub_C1b = 0

    C2a = 0
    C2b = 0
    sub_C2a = 0
    sub_C2b = 0

    C1 = 0
    C2 = 0
   
    C_der = 0
    C_anc = 0
    Csub_der = 0
    Csub_anc = 0

    C_all = 0
    for i, (key,value) in enumerate(anc_res[0].items(),1):
        #if i != 1:
        C_all += math.comb(i,2) * value
        if key in der_intervals:
            #print(key, " is T_der child with weighting ",i," ",t1_count, " and interval ",value, "and C:", math.comb(i,2) * value)
            csum1a = (i*i) * value
            csum1b = i*value
            C1a += csum1a
            C1b += csum1b
            C_der += math.comb(i,2) * value
            Csub_der += math.comb(t1_count,2) * der_intervals[key] * (1-(n/(n+m)))
            #print("C old:",  (csum1a + csum1b) / 2)

            sub_value = der_intervals[key]
            sub_csum1a = (t1_count*t1_count) * sub_value
            sub_csum1b = t1_count * sub_value
            sub_C1a += sub_csum1a
            sub_C1b += sub_csum1b
            t1_count += 1

        elif key in anc_intervals:
            #print(key, " is T_anc child with weightings",i," ",t2_count, " and interval ",value, "and C:", math.comb(i,2) * value)
            csum2a = (i*i) * value
            csum2b = i * value
            C2a += csum2a
            C2b += csum2b
            C_anc += math.comb(i,2) * value
            Csub_anc += math.comb(t2_count,2) * anc_intervals[key] * ((n/(n+m)))
            #print("C old:",  (csum2a + csum2b) / 2)

            sub_value = anc_intervals[key]
            sub_csum2a = (t2_count*t2_count) * sub_value
            sub_csum2b = t2_count * sub_value
            sub_C2a += sub_csum2a
            sub_C2b += sub_csum2b
            t2_count += 1

    #C1 = (C1a + C1b) / 2
    C1 = C_der
    #sub_C1 = (sub_C1a - sub_C1b) / 2
    if C1 == 0:
        C1 = 1
    
    #C2 = (C2a + C2b) / 2
    C2 = C_anc
    #sub_C2 = (sub_C2a - sub_C2b) / 2
    if C2 == 0:
        C2 = 1
    N_all = C_all / (2*((n+m)-1))
    N_der = C1 / (2*(n-k1))
    N_anc = C2 / (2*(m-k2))

    #N_der = C1 / (2*(n-k1)) *(1-(n/(n+m)))
    #N_anc = C2 / (2*(m-k2)) *(n/(n+m))

    N_der_sub = Csub_der / (2*(n-k1))
    N_anc_sub = Csub_anc / (2*(m-k2))
    
    #N_der_sub = Csub_der / (2*(n-k1)) *(1-(n/(n+m)))
    #N_anc_sub = Csub_anc / (2*(m-k2)) *(n/(n+m))

    t1_max_likelihood = -(n - k1)*(math.log(C1/(n - k1))) - (n - k1)
    t2_max_likelihood = -(m - k2)*(math.log(C2/(m - k2))) - (m - k2)
    a = -(n + m - k1 - k2)
    b = math.log((C1 + C2)/(n + m - k1 - k2))
    c = (n + m - k1 - k2)
    t1_t2_max_likelihood = (a * b) - c

    t1_sub_max_likelihood = -(n - k1)*(math.log(Csub_der/(n - k1))) - (n - k1)
    t2_sub_max_likelihood = -(m - k2)*(math.log(Csub_anc/(m - k2))) - (m - k2)
    a = -(n + m - k1 - k2)
    b_sub = math.log((Csub_der + Csub_anc)/(n + m - k1 - k2))
    c = (n + m - k1 - k2)
    t1_t2_sub_max_likelihood = (a * b_sub) - c
 
    pi = (n - k1) / (n + m - k1 -k2)
    phi = C1 / (C1 + C2)
    #phi = C1*(1-(n/n+m)) / (C1*(1-(n/n+m)) + C2*((n/n+m)))
    
    
    
    
    
    
    
    #print("pi :",pi)
    #print("phi :", phi)
    # This simplification has an error so leaving the full terms for now
    #T_statistic_final = (n - k1) * math.log(pi/phi) + ((m - k2) * math.log((phi-pi)/(1-phi)))
    #Ta = (n - k1) * math.log(((C1*n)-(C1*k1)+(C2*n)-(C2*k1))/((C1*n)+(C1*m)-(C1*k1)-(C1*k2)))
    #Tb = (m - k2) * math.log(((C1*m)-(C1*k2)+(C2*m)-(C2*k2))/((C2*n)+(C2*m)-(C2*k1)-(C2*k2)))
    #T_statistic_final = Ta + Tb
    #print("T_statistic_final: ", T_statistic_final)
    
    #term1 = (n - k1) * math.log(((C1 + C2)*(n-k1))/(C1*(n + m - k1 - k2))) 
    #term2 = (m - k2) * math.log(((C1 + C2)*(m - k2))/(C2*(n + m - k1 - k2)))
    #T_statistic = term1 + term2
    #sub_term1 = (n - k1) * math.log(((sub_C1 + sub_C2)*(n-k1))/(sub_C1*(n + m - k1 - k2)))
    #sub_term2 = (m - k2) * math.log(((sub_C1 + sub_C2)*(m - k2))/(sub_C2*(n + m - k1 - k2)))
    #sub_T_statistic = sub_term1 + sub_term2
    #print("T_statistic: ", T_statistic, "sub_T_statistic: ", sub_T_statistic)
    
    T_statistic_basic = (t1_max_likelihood + t2_max_likelihood) - t1_t2_max_likelihood
    T_statistic_improved = T_statistic_basic -math.log(n) - math.log(comb(m,k2)) + math.log(comb(n+m, k1+k2))    
    sub_T_statistic = (t1_sub_max_likelihood + t2_sub_max_likelihood) - t1_t2_sub_max_likelihood
    sub_T_statistic_improved = sub_T_statistic -math.log(n) - math.log(comb(m,k2)) + math.log(comb(n+m, k1+k2))


    #Csub_der = (N_der_sub)*(2*(n-k1))
    #Csub_anc = (N_anc_sub)*(2*(m-k2))
    #phi_sub = Csub_der*(1-(n/n+m)) / (Csub_der*(1-(n/n+m)) + Csub_anc*((n/n+m)))
    #pi = (n - k1) / (n + m - k1 -k2)
    #tst_sub =  (n - k1) * math.log(pi / phi_sub) + (m - k2) * math.log((1 - pi) / (1 - phi_sub))
    #tst_sub_improved = tst_sub -math.log(n) - math.log(comb(m,k2)) + math.log(comb(n+m, k1+k2))


    #N = number of lineages in total at the final generation
    #fN = frequency of the derived allele at the final generation
    #k = number of lineages in total at generation N
    #fk = frequency of derived allele at generation N
    #logF = vector of precomputed log factorial values, e.g. logF[2] = log(2!)
    k = k1 + k2
    fk = k1
    N = n + m
    fN = n
    #print("k,fk,N,fN:",k,fk,N,fN)
    relate_logp = log_pvalue(k,fk,N,fN)
    outdict = {'tree':sys.argv[1], 'position1': sel_pos, 'position2': sel_pos, 'OddsR': oddsr, 'fet_p': fet_p, 'relate_logp':relate_logp, 'Tstatistic': T_statistic_basic, 'n': n, 'k1': k1, 'k2': k2, 'm': m,'sub_Tstatistic': sub_T_statistic,'sub_Tstatistic_improved': sub_T_statistic_improved,'Tstatistic_improved': T_statistic_improved, 'pi': pi, 'phi': phi, 'C1': C1, 'C2': C2, 'N_all':N_all,'N_der_sub':N_der_sub,'N_anc_sub':N_anc_sub, 'N_anc':N_anc,'N_der':N_der,'der_coals':der_coals,'anc_coals':anc_coals}
    #print(outdict)
    print(*list(outdict.keys()),sep="\t")
    print(*list(outdict.values()),sep="\t")


    #print("Coalescence intervals for derived subtree:",der_intervals)
    #print("Terminal DAF: ", n/(n+m))
    #print("Coalescence intervals for derived subtree list:",der_intervals.values())
    #print("Coalescence intervals for ancestral subtree:",anc_intervals)
    #print("Coalescence intervals for derived subtree:",der_intervals)
    #print("Coalescence ages for tree:",sorted(list(set((int(k) for k in t1_node_age.values())))))
    

    
ts = tskit.load(sys.argv[1])
sel_pos = int(sys.argv[2])
## Sampling disabled ##
#nsample = int(sys.argv[3])
#ts = subsample(ts,nsample)
if len(sys.argv) > 3:
    ts2 =  tskit.load(sys.argv[3])
    neut_pos = int(sys.argv[4])
    #ts2 = subsample(ts2,nsample)
    twosite_tstat(ts,sel_pos,ts2,neut_pos)
else:
    onesite_tstat(ts,sel_pos)





#print(sys.argv[1])
#print("Coalescence ages for tree:",sorted(list(set((int(k) for k in t1_node_age.values())))))
#tree = ts2.at(neut_pos)
#print(tree.draw_text())
