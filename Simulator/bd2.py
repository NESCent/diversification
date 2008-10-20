import math
import random

# The nodes map is a map from a unique Node index to the array representing that node:
# [Time of Event, Left Child ,Right Child, Parent node, region, migration_time].
# migration_time is -1 unless a migration has occured. Lineages and Nodes are the same
# (the lineage corresponds to the node at which it was created)
# Left child is the index of the parent node
# Right child is the index of the new node and therefore is the same as the key
# Parent node is the index of the parent node
# Therefore, Left Child is the same as the parent, and Right child is the same as
# the key. However, these are included incase an extinction has occured
# then the extinct node/lineage will be -1.
# Also, migrations will only occur once. Therefore, the birth region of a lineage can be
# inferred from the current region and whether or not time = -1.
# Assume there is at first one lineage, lineage 0.


def Migration(y,time,nodes, direction):
    """
    Input Paramters
    ---------------
    y           the index of the node migrating
    time        the time of migration
    nodes       the node map from node index to
    [ Time of Event (0), Left Child (1),Right Child (2), Parent node (3), region (4), migration_time (5)]
    direction   if 1, the migration is occuring boreal to tropical, if 0, the migration is occuring tropical to boreal

    Return Value(s)
    --------------
    updated events
    """
    if direction==nodes[y][4]:
        print "error, migrating to the state already in"
            
    nodes[y][4]=direction;
    nodes[y][5]=time;
    
    return(nodes)

def Birth(y, time, nodes):
    """
    Input Parameters
    ----------------
    y           the index of the node giving birth
    time        the time of birth
    nodes       the node map from node index to
                [ Time of Event (0), Left Child (1),Right Child (2), Parent node (3), region (4), migration_time (5)]
    Return Value(s)
    --------------
    updated nodes
    """
    #len(nodes) will be the index of the new node, ie. the right child of y, the node giving birth
    #the region of the new node will have to be the same as the region of the parent node
    new_index=len(nodes)
    new_node=[time,y,new_index,y,nodes[y][4], -1]
    nodes.update({new_index:new_node})
    return nodes

def Death(y,nodes):
    """
    Input Parameters
    ----------------
    y           the index of the node giving dying
    time        the time of death
    nodes       the node map from node index to
                [ Time of Event (0), Left Child (1),Right Child (2), Parent node (3), region (4), migration_time (5)]
    Return Value(s)
    --------------
    updated nodes
    """
    return nodes
def CreatePopulation(birth_0, death_0, birth_1, death_1, migration_1, migration_0, dt, intitial_0, initial_1, runs):
    """
    Starting with initial_0 unconnected lineages in the boreal region, and initial_1 unconnected lineages in
    the tropical region, this function simulates the birth-death-migration process until all the lineages left are related
    to eachother
    Input Parameters
    ----------------
    birth_0         the per lineage birth rate in state 0
    death_0         the per lineage death rate in state 0
    birth_1         the per lineage birth rate in state 1
    death_1         the per lineage death rate in state 1
    migration_1     the effective migration rate from 0 to 1
    migration_0     the effective migration rate from 1 to 0
    initial_0       the size of the initial population in state 0
    initial_1       the size of the initial population in state 1
    dt              the time step to be used
    runs            the number of times to run the burn in simulation
        
    Return Values
    -------------
    nodes       the node map
                [ Time of Event (0), Left Child (1),Right Child (2), Parent node (3), region (4), migration_time (5)]
                of the final population
    final_0   the number of lineages in the boreal region at the tips of the tree (after all remaining lineages are related to eachother)
    final_1   the number of lineages in the tropical region at the tips of the tree (after all remaining are related to eachother)
    """
    return [nodes, final_0, final_1]

def Simulate(birth_0, death_0, birth_1, death_1, migration_1, migration_0, dt, intitial_0, initial_1, runs):
    """
    Input Parameters
    ----------------
    birth_0         the per lineage birth rate in state 0
    death_0         the per lineage death rate in state 0
    birth_1         the per lineage birth rate in state 1
    death_1         the per lineage death rate in state 1
    migration_1     the effective migration rate from 0 to 1
    migration_0     the effective migration rate from 1 to 0
    initial_0       the size of the initial population in state 0
    initial_1       the size of the initial population in state 1
    dt              the time step to be used
    runs            the number of times to run the burn in simulation
        
    Return Values
    -------------
    nodes       the node map
                [ Time of Event (0), Left Child (1),Right Child (2), Parent node (3), region (4), migration_time (5)]
                of the final population
    """
    #start with initial_0 in the boreal, initial_1 in the tropical, and run the birth-death, character change (only one allowed
    #per lineage) simulator until the remaining lineages all have the same parent. Then use each of these lineages as initial_0,
    #initial_1 (ie no longer connected) until the remaining lineages all have the same parent.
    #Do this the desired number of times. The last run will be nodes array to be returned.
    #This last run will produce a tree with a certain age (from root to tips). The output tree cannot be shorter in
    #time-span, but it could be longer. As of now, I am not implementing this feature. In order for this simulation
    #to work, initial_0 and initial_1 have to be large, since death will be greater than birth, and therefore,
    #after each run, the total number of lineages will decrease (most likely). The extent of total lineage decline is
    #not known, so it will be best to at first have the birth rate only slightly smaller than the death rate,
    #and have run equal 2 or 3.
    for i in range(0,runs):
        [nodes,initial_0,initial_1]=CreatePopulation(birth_0, death_0, birth_1, death_1, migration_1, migration_0, dt, intitial_0, initial_1, runs)
    return nodes