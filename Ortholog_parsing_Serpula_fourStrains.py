#!/usr/bin/env python

## Full script for parsing orthogroups into Arboretum format, using ETE3

from ete3 import PhyloTree

import re
import os

# read species tree
tree_file = open("/Users/jacky/Research/Projects/InHouse/ANALYSES/Gene_trees/species_tree.tree","r")
nw_spec = ""
for line in tree_file.readlines():
    nw_spec += line.strip()
tree_file.close()

# create new species tree
species_tree = PhyloTree(nw_spec)

# read gene trees
gene_tree_dir = "/Users/jacky/Research/Projects/InHouse/ANALYSES/Gene_trees/all_gene_trees"
tree_files = os.listdir(gene_tree_dir)

gene_trees = {}
for tfile in tree_files:
    if tfile.endswith(".tree"):
        tree_file = open(gene_tree_dir+"/"+tfile,"r")
        nw = ""
        p = re.compile('n\d*')
        for line in tree_file.readlines():
            nw += p.sub('',line.strip())
        tree_file.close()

        # create new gene tree
        gene_trees[tfile.split(".")[0]] = PhyloTree(nw)

# generate Arboretum input
fp_out = open("/Users/jacky/Research/Projects/InHouse/ANALYSES/Arboretum/OGid_members.txt","w")
for OG in gene_trees.keys():
    
    # reconcile with species tree
    gene_trees[OG].set_species_naming_function(lambda node: node.name[:4] )
    recon_tree, events = gene_trees[OG].reconcile(species_tree)

    # get speciation trees
    ntrees, ndups, sptrees = recon_tree.get_speciation_trees(autodetect_duplications=False)
    
    # write out OG info
    # lacE,lacJ,shas,shim
    dup_count = 1
    for spt in sptrees:
        leaves = {"lacE": "NONE",
                  "lacJ": "NONE",
                  "shas": "NONE",
                  "shim": "NONE"}
        for leaf in spt.get_leaves():
            if len(leaf.name) > 4:
                leaves[leaf.name[:4]] = leaf.name
        fp_out.write(OG+"_"+str(dup_count)+"\t"+leaves["lacE"]+","+leaves["lacJ"]+","+leaves["shas"]+","+leaves["shim"]+"\n")
        dup_count += 1
fp_out.close()
        