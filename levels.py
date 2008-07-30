from is_classes import Node 
from is_classes import Tree 
from is_classes import Level 
import tree_traverse

def PrintLevel(level):
    print '+++++++++++++++++++'
    print 'Begin: ' + str(level.begin_time) + level.begin_node.label
    print 'End: ' + str(level.end_time) + level.end_node.label
    for l in level.lineages:
        tree_traverse.PrintLabel(l)

def PrintAllLevels(tree):
    for level in tree.levels:
        PrintLevel(level)

def PrintLevelHistory(level):
    for (time, state) in level.history:
        print state + ' at ' + '%4.2f' % time

def GetBifurcationAge(node):
    return node.bifurcation_age

#The demarcated levels are stored in the list levels[] in the tree
#instance. The levels are stored in ascending order according to their
#begin_time

def DemarcateLevels(tree):
    """ 
    Creates a level-by-level representation of the input tree.
        
    Input parameter
    ---------------

    tree    an instance of the class Tree (defined in is_classes.py)

    Output parameter  
    ---------------

    None.

    But the function populates the list tree.levels which is assumed to be 
    initially empty.

    If there are k leaves in the tree, the tree consists of levels k, k-1,
    ..., 2. Then, tree.levels[j] should be the level k-j of the tree, where 0 <= j
    <= n-2. 

    Each level should be an instance of the class Level (defined in
    is_classes.py)

    Remarks:
    -------

    1. If the input tree has less than two nodes, it is left unmodified.
    
    """

    #all the nodes of tree, including internal
    all_nodes = tree_traverse.TreePartialOrder(tree)

    #   the age here is the bifurcation_age. For leaves, bifurcation_age is 0,
    #   and for internal nodes, the bifurcation_age is the age of the
    #   more recent of its two end points.
    all_nodes_sorted_by_age = sorted(all_nodes, key=GetBifurcationAge)
    for node in all_nodes_sorted_by_age:
        print node.label+':'+str(node.bifurcation_age)
    
    n_nodes = len(all_nodes)
    print "N_NODES = ", n_nodes 

    #   curr_node <- earliest bifurcating node
    #   begin_time <- earliest bifurcating time
    curr_node = all_nodes_sorted_by_age[0]
    begin_node = curr_node
    begin_time = curr_node.bifurcation_age

    #   make the first (going from the tips of the tree to the root) 
    #   level. The first level consists of all terminal
    #   nodes, and is termintated by the first coalescence event
    level = Level()
    level.lineages.append(all_nodes_sorted_by_age[0])

    # The logic here is as follows: initialize begin_time to
    # earliest bifurcation_age. Then go through the rest of the sorted list
    # of nodes one by one. Whenever the previous bifurcation_age (stored
    # in begin_time) differs from the current node's bifurcation_age, it
    # denotes the end of an level. The level then gets its begin_time and
    # end_time set, and gets its list of lineages updated. The update
    # involves removing the two lineages that are involved in the
    # coalescence event and inserting the new node that's the result of the
    # coalescence event.
    # On the other hand, when the previous bifurcation age matches the
    # current node's bifurcation age, it could be we are still going
    # through the leaves (or) as is more unlikely, there could be two
    # coalescences at the same time. WE IGNORE THIS POSSIBILITY. 
    # In that case we just append the
    # current node to the list of this level's lineages [NOTE LINEAGES AND
    # NODES ARE BEING USED INTERCHANGEABLY]. Except if this the last node,
    # then this is the last level; This last level's duration happens to be
    # zero time units. We update this last levels begin, end times and its
    # lineage list and bail out.
   
    for i in range(1, n_nodes):
        curr_node = all_nodes_sorted_by_age[i]
        curr_time  = curr_node.bifurcation_age
        if (begin_time == curr_time):
            level.lineages.append(curr_node)
            if (i == n_nodes-1):
                level.begin_time = begin_time
                level.begin_node = begin_node
                level.end_time = curr_time
                level.end_node = curr_node
                tree.levels.append(level)
                assert len(level.lineages) == 2
            continue

        level.begin_time = begin_time
        level.end_time = curr_time
        level.begin_node = begin_node
        level.end_node = curr_node
        tree.levels.append(level)
        if (i == n_nodes-1):
            assert len(level.lineages) == 2
            continue

        begin_time = curr_time
        begin_node = curr_node

        new_level = Level()
        new_level.lineages = list(level.lineages)  #this is a copy operation
        for child in curr_node.children:
            new_level.lineages.remove(child)
        new_level.lineages.append(curr_node)
        level = new_level
