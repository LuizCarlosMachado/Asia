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


ts = tskit.load(sys.argv[1])

for v in ts.variants():
    var_tree = ts.at(v.site.position)
    parent_node = v.site.mutations[0].node
    root = var_tree.root
    tips = len(tips_of_node(parent_node,var_tree))
    leaves = len(list(var_tree.leaves()))
    freq = tips / leaves
    node_age = var_tree.time(parent_node)
    root_age = var_tree.time(root)
    delta_age = root_age - node_age
    if freq >0.20:
        print(*[sys.argv[1], int(v.site.position), freq, node_age, root_age, delta_age],sep="\t")
        #print(parent_node, sep="\t")














 
  
   
    
     
      
       
