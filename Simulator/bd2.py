import math
import random

# The nodes map is a map from a unique Node index to the array representing that node:
# [Time of Birth, Parent node, region, migration_times, Time of Death,related_to].
# migration_times is [] unless a migration has occured. Time of Death is -1 unless the lineage
# has died. Lineages and Nodes are the same. (the lineage corresponds to the node at which it was created)
# Region is the current region of the node.
# Migrations can occur multiple times, and one lineage can
# migrate from one trait to the other multiple times (ie 0->1, 0->1 are allowed).
# The migration_times will be a list with each value in the list being
# a list with the first component being the migration time, and the second component the
# state to which the lineage migrated.
# Assume there is at first one lineage, lineage 0.

def Migration(y,time,nodes, direction):
    """
    Input Paramters
    ---------------
    y           the index of the node migrating
    time        the time of migration
    nodes       the node map from node index to
    [Time of Birth (0), Parent node (1), region (2), migration_times (3), Time of Death (4), related_to (5)]
    direction   if 1, the migration is occuring boreal to tropical, if 0, the migration is occuring tropical to boreal
        
    Return Value(s)
    --------------
    updated nodes
    """
    nodes[y][2]=direction;
    times=nodes[y][3];
    times.append([time,direction]);
    return(nodes)

def Birth(y, time, nodes):
    """
    Input Parameters
    ----------------
    y           the index of the node giving birth
    time        the time of birth
    nodes       the node map from node index to
                [Time of Birth (0), Parent node (1), region (2), migration_times (3), Time of Death (4), related_to (5)]
    
    Return Value(s)
    --------------
    updated nodes
    the new index
    """
    #len(nodes) will be the index of the new node, ie. the right child of y, the node giving birth
    #the region of the new node will have to be the same as the region of the parent node
    related_to=nodes[y][5];
    new_index=len(nodes)
    new_related_to=related_to.copy()
    new_related_to.add(y)
    new_node=[time,y,nodes[y][2],[], -1, new_related_to]
    nodes.update({new_index:new_node})
    return (nodes, new_index)

def Death(y,time,nodes):
    """
    Input Parameters
    ----------------
    y           the index of the node giving dieing
    time        the time of death
    nodes       the node map from node index to
               [Time of Birth (0), Parent node (1), region (2), migration_times (3), Time of Death (4), related_to (5)]
    Return Value(s)
    --------------
    updated nodes
    """
    nodes[y][4]=time;
        
    return nodes

def CreatePopulation(birth_0, death_0, birth_1, death_1, migration_0, migration_1, initial_0, initial_1):
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
    migration_1     the per lineage migration rate from 0 to 1
    migration_0     the per lineage migration rate from 1 to 0
    initial_0       the size of the initial population in state 0
    initial_1       the size of the initial population in state 1
        
    Return Values
    -------------
    nodes       the node map
                [Time of Birth (0), Parent node (1), region (2), migration_times (3), Time of Death (4), related_to (5)]
                of the final population
    final_0   the number of lineages in the boreal region at the tips of the tree (after all remaining lineages are related to eachother)
    final_1   the number of lineages in the tropical region at the tips of the tree (after all remaining are related to eachother)
    common_rel  the common ancestors of the remaining lineages
    """
    all_related=False;
    time=0;
    nodes={};
    in_0=set(range(initial_0));
    in_1=range(initial_1);
    for i in range(initial_1):
        in_1[i]=in_1[i]+initial_0;
    in_1=set(in_1);
        
    for i in range(initial_0):
        nodes.update({i : [0,-1,0,[],-1,set([i])]});
    for i in range(initial_1) :
        nodes.update({i+initial_0 : [0,-1,1,[],-1,set([i+initial_0])]});

      
    j=0;
    time=0;
    while(all_related==False):
        j=j+1;
        #print(nodes);
        total_0=len(in_0);
        total_1=len(in_1);
        total_events=birth_0+birth_1+death_0+death_1+migration_0+migration_1;
        time_till_event=random.gammavariate(1,total_events);
        time=time+time_till_event;
        #time=time+dt;
        #prob_X is at least one X event in the time dt (1-Poisson probability no events), where X is birth_0, death_0
        #birth_1, death_1, migration_0 or migration_1. Therefore, dt needs to be small enough
        #such that the probability of two events is much smaller than the probability of one event in dt.
        #prob_birth_0=1-math.exp(-birth_0*total_0*dt);
        #prob_death_0=1-math.exp(-death_0*total_0*dt);
        #prob_birth_1=1-math.exp(-birth_1*total_1*dt);
        #prob_death_1=1-math.exp(-death_1*total_1*dt);
        #prob_migration_0=1-math.exp(-migration_0*total_0*dt);
        #prob_migration_1=1-math.exp(-migration_1*total_1*dt);
        #Poisson distribution, probability of no events in time dt
        #prob_of_no_event=math.exp((-birth_0*dt+-death_0*dt+-migration_0*dt)*total_0+(-birth_1*dt+-death_1*dt+-migration_1*dt)*total_1);
        in_0_copy=in_0.copy()
        in_1_copy=in_1.copy()
        all_lineages=in_0_copy.union(in_1_copy);
        all_lineages_copy=all_lineages.copy();
        print("all lineages "+str(in_0)+str(in_1));
        if(len(all_lineages)==0):
            print "all died"
            return -1;
        random_list=range(len(in_0));
        for i in range(len(in_0)):
            random_list[i]=in_0_copy.pop();
        random.shuffle(random_list);
        if len(random_list)==0:
            random_index_0=-1;
        else:            
            random_index_0=random_list[0];
        
        random_list=range(len(in_1));
        for i in range(len(in_1)):
            random_list[i]=in_1_copy.pop();
        random.shuffle(random_list);
        if len(random_list)==0:
            random_index_1=-1;
        else:            
            random_index_1=random_list[0];

        random_list=range(len(all_lineages ))     
        for i in range(len(all_lineages)):
            random_list[i]=all_lineages_copy.pop();
        random.shuffle(random_list);
        random_index_all=random_list[0];
        #row_for_choosing 0:birth_0, 1:death_0, 2:birth_1, 3:death_1, 4:migration_0, 5:migration_1    
        row_for_choosing=range(6);
        #special case
        if(len(in_0)==0):
            row_for_choosing[0]=0;
            row_for_choosing[1]=0;
            row_for_choosing[2]=0;
        #normal
        else:
            row_for_choosing[0]=birth_0;
            row_for_choosing[1]=row_for_choosing[0]+death_0;
            row_for_choosing[2]=row_for_choosing[1]+migration_0;
        #special case
        if(len(in_1)==0):
            row_for_choosing[3]=0+row_for_choosing[0];
            row_for_choosing[4]=0+row_for_choosing[0];
            row_for_choosing[5]=0+row_for_choosing[0];
        #normal
        else:
            row_for_choosing[3]=row_for_choosing[2]+birth_1;
            row_for_choosing[4]=row_for_choosing[3]+death_1;
            row_for_choosing[5]=row_for_choosing[4]+migration_1;

        #normalize
        sum=0;
        for i in range(len(row_for_choosing)):
            sum=sum+row_for_choosing[i];
        #print("not normalized "+str(row_for_choosing));
        print("normalizing constant "+str(row_for_choosing[len(row_for_choosing)-1]));
        for i in range(len(row_for_choosing)):
            row_for_choosing[i]=row_for_choosing[i]/row_for_choosing[len(row_for_choosing)-1];
        #print(row_for_choosing[len(row_for_choosing)-1]);
        #print(row_for_choosing);       
        event_number=GetEvent(row_for_choosing);
        print "event number "+str(event_number);    
        
        if(event_number==0):#birth 0
            #print "birth 0";
            index_to_birth=random_index_0;
            (nodes, new_index)=Birth(index_to_birth,time,nodes);
            in_0.update([new_index]);
        if(event_number==1):#death 0
            #print "death 0";
            index_to_die=random_index_0;
            nodes=Death(index_to_die,time,nodes);
            in_0.remove(index_to_die);
        if(event_number==3):#birth 1
            #print "birth 1";
            index_to_birth=random_index_1;
            (nodes, new_index)=Birth(index_to_birth,time,nodes);
            in_1.update([new_index]);
        if(event_number==4):#death 1
            #print "death 1";
            index_to_die=random_index_1;
            nodes=Death(index_to_die,time,nodes);
            in_1.remove(index_to_die);
        if(event_number==2):#migration 0
            #print "migration 0";
            index_to_migrate=random_index_all;
            nodes=Migration(index_to_migrate,time,nodes,0)
            if len(in_1.intersection(set([index_to_migrate])))!=0:
                in_1.remove(index_to_migrate);
                in_0.update([index_to_migrate]);
        if(event_number==5):#migration 1
            #print "migration 1";
            index_to_migrate=random_index_all;
            nodes=Migration(index_to_migrate,time,nodes,1)
            if len(in_0.intersection(set([index_to_migrate])))!=0:
                in_0.remove(index_to_migrate);
                in_1.update([index_to_migrate]);

        #print("nodes");
        #print(nodes);
        #print("in 0")
        #print(in_0);
        #print("in 1");
        #print(in_1);
        all_lineages=in_0.union(in_1);
        length=len(all_lineages);
        all_related=True;
        common_rel=nodes[all_lineages.pop()][5]
        #print("common");
        #print(length-1);
        while(len(all_lineages)>0):
            common_rel=nodes[all_lineages.pop()][5].intersection(common_rel);
            #print(common_rel);
            if(len(common_rel)==0):
                all_related=False;
        #print("j "+str(j));
        if j==10000:
            all_related=True;
            print("forced quit");
        
    #return [nodes, final_0, final_1]
    return [nodes, in_0, in_1, len(in_0), len(in_1), common_rel]    

def GetEvent(row_for_choosing):
    """
    actual probabilities
    [.1, .2, .4, .3]
    
    row for choosing
    [.1, .3, .7, 1]

    """
    n=random.uniform(0,1);
    print(str(n)+"random")
    for i in range(len(row_for_choosing)):
        if(i==0):
            if(0<=n<row_for_choosing[i]):
                return (i)
        if(i==(len(row_for_choosing)-1)):
            if(row_for_choosing[i]<=n<=1):
                return (i)
        if(row_for_choosing[i-1]<=n<row_for_choosing[i]):
            return (i)
    return -1;            

def Simulate(birth_0, death_0, birth_1, death_1, migration_1, migration_0, intitial_0, initial_1, runs):
    """
    Input Parameters
    ----------------
    birth_0         the per lineage birth rate in state 0
    death_0         the per lineage death rate in state 0
    birth_1         the per lineage birth rate in state 1
    death_1         the per lineage death rate in state 1
    migration_1     the per lineage migration rate from 0 to 1
    migration_0     the per lineage migration rate from 1 to 0
    initial_0       the size of the initial population in state 0
    initial_1       the size of the initial population in state 1
    runs            the number of times to run the simulation
        
    Return Values
    -------------
    nodes       the node map
                [Time of Birth (0), Parent node (1), region (2), migration_times (3), Time of Death (4)]
                of the final population
    common_rel  the common ancestors of the remaining lineages
    in_0   the number of lineages in the boreal region at the tips of the tree (after all remaining lineages are related to eachother)
    in_1   the number of lineages in the tropical region at the tips of the tree (after all remaining are related to eachother)
    """
    #start with initial_0 in the boreal, initial_1 in the tropical, and run the birth-death, character change simulator
    #until the remaining lineages all have the same parent. Then use each of these lineages as initial_0,
    #initial_1 (ie no longer connected) until the remaining lineages all have the same parent.
    #Do this the desired number of times. The last run will be nodes array to be returned.
    #This last run will produce a tree with a certain age (from root to tips). The extent of total lineage decline is
    #not known, so it will be best to at first have the birth rate only slightly smaller than the death rate,
    #and have run equal to one.
    for i in range(0,runs):
        [nodes,in_0, in_1, initial_0,initial_1, common_rel]=CreatePopulation(birth_0, death_0, birth_1, death_1, migration_1, migration_0, intitial_0, initial_1)
        #common_rel=CreatePopulation(birth_0, death_0, birth_1, death_1, migration_1, migration_0, dt, intitial_0, initial_1)
    return [nodes, in_0, in_1, common_rel]

def GetCounts(nodes_real):
    counts=range(len(nodes_real));
    for i in range(len(counts)):
     	counts[i]=0;
     	
    nodes=nodes_real.copy();
    for i in range(len(nodes)):
        size=len(nodes[i][5]);
        if(nodes[i][4]==-1):#only count a node if it is not dead
            #print("not dead "+str(i));
            #print("size of ancestors "+str(size));
            for j in range(size):
                index=nodes[i][5].pop()
                #print("parent node "+str(index));
                counts[index]=counts[index]+1;  
    return counts

def PrintNodes(nodes, common_rel, in_0, in_1, total_begin):
    """
    Input:
    ------
    nodes       the node map
                [Time of Birth (0), Parent node (1), region (2), migration_times (3), Time of Death (4)]
    common_rel
    total_begin the total number of beginning lineages

    Output:
    ------
    nodes_map:
                parent: [children], a map of each each parent node to an array of all its immediate children
    """
    #common=common_rel.copy();
    #common_rel=common.pop();
    nodes_map={};
    full_final_map_0={};
    full_final_map_1={};
    #iterate through all nodes and assign them to a parent
    for child in range(len(nodes)):
        parent_node=nodes[child][1];
        if(nodes_map.has_key(parent_node)):
            children_array=nodes_map[parent_node];
            children_array.append(child);
        else:
            nodes_map.update({parent_node:[child]});
    print("real test"+str(nodes_map.has_key(25)));

##    for node in in_0:
##        if(not nodes_map.has_key(node)):#only search for lineages at leaves
##            final_map_0={};
##            (final_map_0, found)=Dfs(node, nodes_map, final_map_0,common_rel,'FALSE');
##            full_final_map_0[node]=final_map_0;
##    #for node in in_1:
##    for node in in_1:
##        if(not nodes_map.has_key(node)):#only search for lineages at leaves
##            final_map_1={};
##            (final_map_1, found)=Dfs(node,nodes_map, final_map_1,common_rel,'FALSE');
##            full_final_map_1[node]=final_map_1;
    node=common_rel.intersection(set(range(total_begin))).pop();
    print node
    final_map={};
    (final_map, found)=Dfs(node, nodes_map, final_map, node, nodes, 'FALSE');
    full_final_map=final_map

    return [nodes_map, full_final_map]            

def Dfs(begin_node, nodes_map, final_map, current_node, nodes, found):
    print "current node";
    print current_node
    print "node";
    print begin_node
    print("test "+str(nodes_map.has_key(25)));
    if(nodes_map.has_key(current_node)):
        children=nodes_map[current_node];
        for child in children:
            (final_map,found)=Dfs(begin_node, nodes_map, final_map, child, nodes, found)
            if(found=='TRUE'):
                print("update")
                if(final_map.has_key(current_node)):
                   final_map[current_node].append(child);
                   print(str(current_node)+":"+str(final_map[current_node]));
                else:
                    final_map.update({current_node:[child]});
                    print(str(current_node)+":"+str([child]));
##            print "child"
##            print child
##            if(child==end_node):
##                final_map.update({current_node:child});
##                print "found"
##                print "updated"
##                print(str(current_node)+":"+str(child));
##                found='TRUE';
##            else:
##                print "else"
##                print("current node "+str(current_node));
##                print(str(child)+" "+str(end_node));
##                (final_map, found)=Dfs(end_node,nodes_map, final_map, child, found);
##                if(found=='TRUE'):
##                    print "found true within else"
##                    final_map.update({current_node:child});
##                    print "updated"
##                    print(str(current_node)+":"+str(child));
    else:#got to a leaf
        if(current_node==25):
            print("25" +str(nodes_map.has_key(current_node)));
        print("got to leaf with current node "+str(current_node));
        if(nodes[current_node][4]==-1):#leaf node is not dead
            found='TRUE'
            print("found")
##        else:
##            found='FALSE'
##            print("not found")
        
    return (final_map, found)               
            

    

    