import pandas as pd
from ete3 import Tree
import csv
import os
import argparse
import glob

name_conversions = pd.read_csv("SPECIES_NAME_CONVERSION.csv")
name_map = dict(zip(name_conversions['Species'], name_conversions['SpeciesName']))

def get_closest_relative(tree, leaf_name):
    node = tree.search_nodes(name=leaf_name)[0]
    parent = node.up
    if parent:
        relatives = [leaf.name for leaf in parent.get_leaves() if leaf.name != leaf_name]
        return set(relatives)
    return set()


gene_tree_file="/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/results/251021_by_POL/412/412.fa.pol.hmm.fa.algn_cleaned.fasta.treefile"
te=os.path.basename(os.path.dirname(gene_tree_file))
st = Tree("AMP_tree.tre", format=1, quoted_node_names=True)
gt = Tree(gene_tree_file, format=1, quoted_node_names=True)

seen_species = set()
nodes_to_remove = []
for leaf in gt:
    if "_" in leaf.name:
        leaf.name = leaf.name.split("_")[0]
    if leaf.name in name_map:
        leaf.name = name_map[leaf.name]
    if leaf.name in seen_species:
        nodes_to_remove.append(leaf)
    else:
        seen_species.add(leaf.name)
for node in nodes_to_remove:
    node.detach()

common = list(set(st.get_leaf_names()) & set(gt.get_leaf_names()))
st.prune(common)
gt.prune(common)

with open("HTT_"+te+".csv", mode='w', newline='') as f:
    writer = csv.writer(f)
    for sp in common:
        s_rel = get_closest_relative(st, sp)
        g_rel = get_closest_relative(gt, sp)
    
        if g_rel and not (g_rel & s_rel):
            normal = next(iter(s_rel)) if s_rel else "None"
            swapped = next(iter(g_rel))
            writer.writerow([sp, swapped, te])
            
        
