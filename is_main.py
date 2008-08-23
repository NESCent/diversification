#!/usr/bin/python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names, module and package names.
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

# Notes (******PLEASE READ*******)
# -----
#
# -1. In python anything within triple quotes """ is ignored. It is a
#    python convention to place a documentation string within triple
#    quotes at the beginning of each function.
#
# 0. IN PYTHON INDENTATION IS MEANINGFUL!!!
#    So a for loop to calculate the sum of first 10 natural numbers 
#    could be written as:
#    sum = 0
#    for j in range(0, 10):
#       sum = sum+j
#    print sum <---- the indention is the only indication that the for loop
#                    has ended at this point!!!!!!!!!!!!!!!!!!!
#    This may seem weird, but it's really easy to get used to.
#
# 1. ALL EQUATION, SECTION, PAGE NUMBERS ETC. REFER TO mossEqModel.pdf
#
# 2. In what follows, a "migration" means the change of state of the
#    phenotypic character. 
#
# 3. In what follows, a "lineage" is used to refer to a species in the given phylogenetic
#    tree (pardon the error, but note that I have not used this terminology
#    uniformly in mossEqModel.pdf. There I use lineages and species
#    interchangeably. - Ganesh)
#
# 4. Keep in mind that python list indices start from 0.
#
# 5. Note that the word 'state' is used in two different senses:
#       
#       (a) to denote the character state (0/1) of a lineage, for example if it
#           is boreal or neotropical. 
#
#       (b) as in when we speak of the 'state space' of a Markov chain. In
#           what follows whenever we speak of the state of the 'conditional
#           jump chain' this is what we mean.

# standard python imports
import os
import sys
import re
import levels

# Ganesh imports
import parse
#import tree_traverse
#import levels
#from is_classes import Node 
#from is_classes import Tree 

p_lambda = 0.1
p_b = 0.05
p_B = 1000
p_alpha = 0.1
p_mu = 0.11

sigma = (p_lambda, p_b, p_B, p_alpha, p_mu)

#tree_file = sys.argv[0]

trees_block = False
states_block = False
begin_trees = r'begin(\s)+trees:'
re_begin_trees = re.compile(begin_trees, re.IGNORECASE)  
begin_states = r'begin(\s)+states:'
re_begin_states = re.compile(begin_states, re.IGNORECASE)  
begin_taxa = r'begin(\s)+taxa:'
re_begin_taxa = re.compile(begin_taxa, re.IGNORECASE)  
end = r'end(;)?'
re_end = re.compile(end, re.IGNORECASE)

#both tables are dictionaries
taxon_table = {}
state_table = {}

taxa_block = False
states_block = False

TREE_FILE = open('onetree.nex', "r")

for line in TREE_FILE.readlines():
    # skip if line contains just whitespaces
    if len(line.split()) == 0:
        continue

    # if end of a block, set all block flags to false and continue
    x = re_end.match(line)
    if x != None:
        states_block = False
        taxa_block = False
        trees_block = False
        continue

    # get the name of the taxon, put it in the taxon table and move on
    # Note: a taxon block in a Nexus file may have slightly different
    # syntax. This is a place holder
    if taxa_block == True:
        number_name = line.split()
        number = number_name[0]
        name = number_name[1]
        taxon_table[number] = name
        continue

    # get the name of the area, put it in the area table and move on
    if states_block == True:
        if len(taxon_table) == 0:
            taxon_name_state = line.split()
            taxon_name = taxon_name_state[0]
            state = taxon_name_state[1]
            state_table[taxon_name] = state
        else:
            taxon_number_state = line.split()
            taxon_number = taxon_number_state[0]
            state = taxon_number_state[1]
            state_table[taxon_table[taxon_number]] = state
        continue

    if trees_block == True:
        # The name of a tree along with the equal to sign and any spaces around
        # the equal to sign. A tree name can be any sequence of alphanumeric
        # characters. Note that underscore is considered alphanumeric
        line = line.strip()
        tree_tuple = line.split('=')
        tree_name = tree_tuple[0].strip()
        tree_string = tree_tuple[1].strip()
        parsed_tree = parse.Read(tree_name, taxon_table, state_table, tree_string)
        parsed_tree.name = tree_name
        #tree_traverse.PostOrderTraverse(parsed_tree.root, tree_traverse.CalculateNodeAges)
        #tree_traverse.PrintLabel(parsed_tree.root)
        #tree_traverse.PostOrderTraverse(parsed_tree.root, tree_traverse.PrintSelfParentLabel)
        levels.DemarcateLevels(parsed_tree)
        #levels.PrintAllLevels(parsed_tree)
        #tree_traverse.TreeLeaves(parsed_tree)

        #transition_matrices = []
        #backpropose.MakeTransitionMatrices(parsed_tree.num_leaves, transition_matrices)
        #backpropose.Backpropose(parsed_tree, transition_matrices)
        #likelihood = is_likelihood.LikelihoodOfParameters(parsed_tree, state_table, sigma)
        #print likelihood
        continue
    
    #look for a trees block or an states block
    r = re_begin_trees.match(line) 
    s = re_begin_states.match(line)
    t = re_begin_taxa.match(line)
    if r != None:
        trees_block = True
        continue
    if s != None:
        states_block = True
        continue
    if t != None:
        taxa_block = True
   
TREE_FILE.close()
