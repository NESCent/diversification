#! /usr/bin/env python

# This code creates a population simulation for 2 geographical regions: boreal and neotropical. 
# Author: NIHSHANKA DEBROY
# Initially, I generate a starting population in each region. The BOREAL region can undergo 
# TURNOVER events or MIGRATION events. In a turnover event, a random speciation event is followed
# by a random extinction event. In a migration event, an organism moves to the neotropical region,
# leaving a copy in the boreal region. The population size of the boreal region remains constant in 
# this manner. The neotropical region can have SPECIATIONS or EXTINCTIONS, in addition to the migration
# events from the boreal region.  So the population size in the neotropical region can change.
# The OUTPUT of the simulation is the SET OF ORGANISMS in each population region and the EVENT HISTORY.

# importing  packages
import random
import math

# Function to find the MAXIMUM of two integers, a and b, provided as input
def Maximum(a,b):
	if a>b:
		return a
	else:
		return b
# end of function MAXIMUM

# Function for MIGRATION event: Takes as INPUT array holding boreal population (c_pop), index for the
# species undergoing migration (y), time of event (f_times), array holding
# event history (I_event) for boreal region, array holding neotropical population (c_pop1), array holding
# event history (I_event1) for neotropical region.
# left_lineage=right_lineage=-2 is a flag for migration
# RETURNS updated population and event history arrays for boreal and neotropical regions.
#def Migration(c_pop,y,f_times,I_event,c_pop1,I_event1):
def Migration(population,y,time,event):
	    #event # and node # are the same
	    # For a description of what the entries in the I_event array represent, look in "MAIN CODE"
		# Adding the event to the boreal region as [event #,time of event,-1,-1,-1].
	    #I_event.append([len(I_event)+1,f_times,-2,-2,migrate_parent])
	    # Add a copy of the migrating organism to the neotropical region
	    event[y][5]=1;
	    event[y][6]=time;
	    #c_pop1.append([c_pop[y-1][0],-1,0])
	    population[y][0]=1
	    # In the event history of the neotropical region, add an event as [event #, time of event, -1,-1,-1].
        # return the updated population and event history arrays for both regions
	    #return (c_pop, I_event,c_pop1,I_event1)
        return(population,event)
# end of function Migration

# Function for BIRTH event: Takes as INPUT array holding current population (c_pop), index for the
# species undergoing birth (y), name of species born (i_count+1), time of event (f_times), array holding
# event history (I_event), ch_pop is equal to "bor" if the birth is happening in the boreal region and "neo"
# if it is happening in the neotropical region.
#  RETURNS updated population array (c_pop) and event history array (I_event)
def Birth(c_pop, y, i_count_species, i_count_events, f_times,I_event,ch_pop):
	# i_temp is the number of events so far - 1 
	i_temp = len(I_event)-1
	# If there is at least one species already in the population c_pop,
	# add new species to c_pop.
	if len(c_pop) > 0:
		print c_pop[y-1][0],'undergoes birth.'
		if ch_pop == "bor":
			c_pop.append([i_count_species + 1,-1,0])
		else:
			c_pop.append([i_count_species + 1,-1,1]) 
		# Otherwise, put the first element into c_pop
		print 'A birth took place' 
	else:
		if ch_pop == "bor":
			c_pop = [[i_count_species+1,-1,0]]
		else:
		    c_pop = [[i_count_species+1,-1,1]] 
		    #print 'Inter-event time:',f_times,'seconds. Species:', c_pop, '.'
		    # Update I_event to have event time, left and right children
		    # if first event (that is, time of latest event is 0), 
		    if I_event[0][1] == 0:
		    	# set event time to be f_times
		    	I_event[0][1]=f_times;
                # set left child of event to be first species in c_pop
                I_event[0][2]=c_pop[0][0]
                # if c_pop has at least 2 species, set right child of event to be second species in c_pop
		if len(c_pop) > 1:
			I_event[0][3]=c_pop[1][0]
	else:
		# otherwise, append the event to the event history:
		# i_iter  is the number of events so far - 1
		i_iter=len(I_event)-1
		i_check=0
		i_found=-1
		if (y-1) >=0: 
			# Find the most recent event that has the species undergoing birth as its left/right child
			# Set i_found to point to that event - this is the index of the parent node in the tree for
			# the latest event due to the birth. If there is no parent node, i_found stays as -1.
			while i_iter >=0  and i_check==0:
				if I_event[i_iter][2]==c_pop[y-1][0] or I_event[i_iter][3]==c_pop[y-1][0]:
					i_found = I_event[i_iter][0]
                    i_check=1

				else:
					i_iter=i_iter-1
		#print 'Event history',I_event[i_count-1][1]+f_times,'c_pop',c_pop[y-1][0]
		
		# Append event: ([event #, event time, left child = y-1 (species undergoing birth), right child = i_count+1, parent node (i_found)])
		I_event.append([i_count_events+1,f_times,c_pop[y-1][0],i_count+1,i_found])

        #print 'Event history:',I_event,'\n'
	    #return the updated population and event history arrays
	   	return (c_pop, I_event)
# end of function BIRTH

# Function for DEATH event: Takes as INPUT array holding current population (c_pop), index for the
# species undergoing death (y), i_count (event # = i_count+1), inter-event waiting time (f_times), array holding
# event history (I_event). RETURNS updated population array (c_pop) and event history array (I_event)
def Death(c_pop,y,i_count,f_times,I_event,ch_pop): 
	   # i_iter is the number of events so far - 1
           i_iter = len(I_event)-1
	   # i_check keeps track of whether species dying out is found while traversing tree of events
           i_check = 0

	   # delete the most recent edge corresponding to the dying species
           while i_iter >=0  and i_check==0:
		# if there is at least one species in the population
                if len(c_pop)>0:
			# if left child of an event points to species dying out
			if I_event[i_iter][2]==c_pop[y-1][0]:
				# set the left child to -1
				I_event[i_iter][2]=-1
                        	i_check=1
			else:
				# if right child of an event points to species dying out
				if I_event[i_iter][3]==c_pop[y-1][0]:
					# set the right child to -1
					I_event[i_iter][3]=-1
					i_check=1
				else:
					i_iter=i_iter-1
                else:
                	i_check=1;

           # if the population only has one organism and i_count is 0
           if len(c_pop) == 1 and i_count==0:
		# set the event time = f_times
           	I_event[0][1]=f_times
	   else:		
		# otherwise, append the event as ([event #,event time, -1, -1, -1])
           	I_event.append([i_count+1,f_times,-1,-1,-1])
	
	   if len(c_pop) > 0:
           	print c_pop[y-1][0],'undergoes death. \n'#Inter-event time:', f_times,'seconds.'
                # removing the organism from the population
	   	c_pop.remove(c_pop[y-1])
	   else:
		# This happens in the rare case that there is a death, but the population is already 0. 
		# So the death does not result in any change in this case.
		print 'No species remaining. Population unchanged.'
	   #print 'Species:',c_pop,'.'	   
	   #print 'Event history:',I_event,'\n'
	   # return the updated population and event history arrays
           return (c_pop,I_event)       
# end of function DEATH

# Function for POPULATION CREATION: Takes as INPUT the current population array (c_pop), 
# desired population size (i_pop_size), birth rate(f_brate), death rate(f_drate), event history (I_event).
# RETURNS updated population and event history arrays
def Create_Population(c_pop,i_pop_size,f_brate,f_drate,I_event,ch_pop):

	# Pr(T<t) = 1 - exp(- ((f_brate+f_drate)*(i_count+1)) * t), random.expovariate returns the waiting times
	for i_count in range(0,i_pop_size):
		# f_times is the inter-event waiting time
		f_times=random.expovariate((f_brate*(i_count+1))+(f_drate*(i_count+1)))
                if len(I_event) > 0:
			f_maxtime = I_event[len(I_event)-1][1]
		else:
			f_maxtime=0;

		f_tot=f_times+f_maxtime
      		# For each event, consider an individual at random: 
        	x= random.random()
		# y is the index of the organism undergoing birth
        	y= int(math.ceil(x*len(c_pop)))
                # Note: There is a very small, but non-zero, probability that an individual could undergo birth more than once
        	# simulating births and deaths randomly 
        	f_ran = random.random()
        	f_temp = (f_brate*(i_count+1))/(((f_brate+f_drate)*(i_count+1)))

		# BIRTH
        	if f_ran <= f_temp:
			if ch_pop == "bor":
				c_pop, I_event = Birth(c_pop,y,i_count,f_tot,I_event,"bor")
			else:
				c_pop, I_event = Birth(c_pop,y,i_count,f_tot,I_event,"neo") 
                        #print 'Birth'
        	# DEATH
        	else:
			if ch_pop == "bor":
				c_pop, I_event = Death(c_pop,y,i_count,f_tot,I_event,"bor")
			else:
				c_pop, I_event = Death(c_pop,y,i_count,f_tot,I_event,"neo")
			#print 'Death'
		#print 'Species:',c_pop,'Number of species:', len(c_pop)
		#print c_pop,"\n", I_event
	# return the updated population and event history arrays
	return (c_pop,I_event)
# end of function CREATE_POPULATION

# function to print out the ancestry of a species in the population
# Input: Population array (c_pop), event history (I_event), index of species of interest (i_index)
def Print_Ancestry(c_pop, I_event,i_index):
	i_stop=0
	i_term=0
	while((i_stop==0) and (i_term<10)):
		i_term=i_term+1
		if (i_index==0):
			print "Species:",c_pop[i_index][0]," No ancestor in current population."
			i_stop=1
		else:
			#while(i_index!=0):
			print "Species:",c_pop[i_index][0]
			i_search= len(I_event)-1
			# loop to find event # in history that has fourth element (right child) == species of interest
			while((I_event[i_search][3] != c_pop[i_index][0]) and (i_search>=0)):
				i_search=i_search-1
			#print i_search
        		#print c_pop[i_index][0],"was born from", I_event[i_search][2],"at time",I_event[i_search][1]
			if i_search==-1:
				print c_pop[i_index][0],"was not born from some species - it was there at the start"
				i_term=11
			else:
				# if left child is not -1
				if I_event[i_search][2]!=-1:
					# if left child is not equal to right child
					if I_event[i_search][2]!=I_event[i_search][3]:
        					print c_pop[i_index][0],"was born from", I_event[i_search][2],"at time",I_event[i_search][1]
						i_test=0
						i_found=0
						while((i_test<len(c_pop))and(i_found==0)):
							# search for left edge (species) in population array
							if(c_pop[i_test][0]==I_event[i_search][2]):
								i_found=1
								i_index=i_test
							else:
								i_test=i_test+1
						#print i_test
						if(i_found==0):
							i_term=11
					else:
        					print c_pop[i_index][0],"was born at time",I_event[i_search][1]
						print '- no ancestor, since population was empty before it was born'
				else:
					print c_pop[i_index][0],"has no living parent"
					i_term=11
        				#print c_pop[i_index][0],"was born from", I_event[I_event[i_search][4]-1][3],"at time",I_event[i_search][1]
		
	return 0

# end of function PRINT_ANCESTRY

# ***************************************          MAIN CODE   ************************************
# ***************************************          MAIN CODE   ************************************

# Inputting the speciation, extinction, migration and turnover rates
f_brate=raw_input("Please enter the speciation rate:")
f_brate=float(f_brate)
f_drate=raw_input("Please enter the extinction rate:")
f_drate=float(f_drate)
f_mrate=raw_input("Please enter the migration rate:");
f_mrate=float(f_mrate)
f_turnover=raw_input("Please enter the turnover rate:");
f_turnover=float(f_turnover)

# Inputting the desired STARTING population sizes (BEFORE the simulation begins)
i_bor_size=raw_input("Please enter desired number of events to create STARTING population in boreal region:")
i_bor_size=int(i_bor_size)
i_neo_size=raw_input("Please enter desired number of events to create STARTING population in neotropical region:")
i_neo_size=int(i_neo_size)

# Inputting the number of events to be carried out in the simulation
i_sim_size=raw_input("Please enter desired number of events for simulation:")
i_sim_size=int(i_sim_size)

# VARIABLES and INITIALIZATION

# Population of individuals
# The -1 flags indicate that the organisms have not undergone migration yet
# Boreal (denoted by third argument - 0)
c_bor = [[0,-1,0]]
# Neotropical (denoted by third argument - 1)
c_neo = [[0,-1,1]]

# Event histories - first no. is event #, second no. is event time, third no. is left-child, fourth no. is right-child.
# fifth no. is parent (all pointers initialized to -1)
# For boreal region
I_event_bor = [[1,0,-1,-1,-1]]
# For neotropical region
I_event_neo = [[1,0,-1,-1,-1]]

# ALGORITHM: 
# Step 1: Generate the initial boreal and neotropical populations, using the Create_Population function.
# Step 2: As many times as we want to simulate, generate inter-arrival waiting times 
# from an exponential distribution. 
# Step 3: For each waiting time(event), a turnover, migration, birth or death occurs (probabilities computed below).
#  Turnovers and migrations occur in the boreal region. Turnovers involve a speciation followed by an extinction.
#  In migration, the species leaves a copy in the boreal region and moves to the neotropical region.
#  Births and deaths occur in the neotropical region. Note: For a stationary distribution to be reached, the death
#  rate should exceed the birth rate. 
# Building the tree -  The vertices/nodes in the tree are events, and the edges correspond to species. Each time a
# birth takes place, two edges come out of the corresponding event node - the organism undergoing birth/speciation and 
# the child. The edge lengths correspond to the waiting times, which can be obtained from the event history record. 
# Each time a death takes place, the most recent edge corresponding to the extinct organism gets removed. In this manner, 
# the edges from the leaf nodes correspond to existing organisms.

print '\t\t\t\t\t\t\t\t**********************************************\n'
print '\t\t\t\t\t\t\t\t**********************************************\n'

print '\nThe existing species and the event history for each population region are printed out. Each event has five pieces of information'
print ' [event #, event time, left edge, right edge, parent node(event)]. If there is no left/right edge, that corresponding entry'
print '  is set to -1. So, for example, [2,.5,2,-1,1] would refer to event 2 (that has event 1 as its parent node), with species 2'
print ' as left edge and no right edge, occurring at time .5 seconds. The event time displayed is the event time in the cumulative sense,'
print ' not the inter-event time. When a death occurs, an event is added to the history as [event #, event time, -1,-1,-1]. \n\n The simulation'
print 'starts out with species [0] in existence. The third entry for each species indicates whether it is from the boreal or neotropical regions.'
print ' So [x,y,0] means that species x is in the boreal region, [x,y,1] means it is in the neotropical region. The second entry for each species'
print ' indicates whether it has migrated yet. So [x,-1,0] means that species x has not migrated yet.'

print '\n\nBirth rate:',f_brate,'Death rate:',f_drate,'Migration rate:',f_mrate,'Turnover rate:',f_turnover
print '\n\nGenerating starting population in boreal region...'

# Create the initial population and event history in boreal region
(c_bor,I_event_bor) = Create_Population(c_bor,i_bor_size,f_brate,f_drate,I_event_bor,"bor")
print '\n\nGenerating starting population in neotropical region...'
# Do the same for the neotropical region
(c_neo,I_event_neo) = Create_Population(c_neo,i_neo_size,f_brate,f_drate,I_event_neo,"neo")
print '\nPrinting out the starting populations and event histories for the simulations...\n'
print 'Boreal:\n','\nPopulation:',c_bor,'\nEvent history:',I_event_bor,'\n\nNeotropical:\n','\nPopulation:',c_neo,'\nEvent history:',I_event_neo,'\n'

# The events take place as follows. Let P = f_turnover+((i_count+1)*(f_brate+f_drate))+f_mrate. With probability,
# 1) f_turnover/P, turnover event
# 2) f_mrate/P, migration event
# 3) i_count*f_brate/P, birth event
# 4) i_count*f_drate/P, exinction event
print 'SIMULATION...\n'

for i_count in range(0,i_sim_size):    
    # inter-event waiting time
    P = f_turnover +((i_count+1)*(f_brate+f_drate)) + f_mrate
    f_times=random.expovariate(P)
    f_maxtime=Maximum(I_event_bor[len(I_event_bor)-1][1],I_event_neo[len(I_event_neo)-1][1])
    f_tot=f_times+f_maxtime
	#print f_maxtime
        #print f_times
    f_rand = random.random()
	# if random is b/w 0 and f_turnover/P, let's say it's a turnover event
    if f_rand <= f_turnover/P:
		# Pick a species randomly, it undergoes birth.
        x= random.random()
		# y is the index of the species undergoing birth
        y= int(math.ceil(x*len(c_bor)))
        print '\nTurnover event\n'
        #if len(c_bor)>0:
        #	print c_bor[y-1],' undergoes birth\n'
        # call Birth function
        (c_bor, I_event_bor) = Birth(c_bor,y,len(I_event_bor),f_tot,I_event_bor,"bor")
		#print c_bor, I_event_bor
		# Pick a species randomly, it undergoes death 
        x= random.random()
        y= int(math.ceil(x*len(c_bor)))
		# y is the index of the species undergoing death
		#if len(c_bor)>0:
		#	print c_bor[y-1],' undergoes death\n'
		# call Death function
	 	#print c_bor, y, i_count, I_event_bor
        (c_bor, I_event_bor) = Death(c_bor,y,len(I_event_bor),f_tot,I_event_bor,"bor")
		#print c_bor, I_event_bor
	# if random is b/w f_turnover/P and (f_turnover+f_mrate)/P, let's say it's a migration event
    elif (f_rand > f_turnover/P) and (f_rand <= (f_turnover+f_mrate)/P):
		# Pick a species randomly, copy/migrate it to the neotropical region
		x= random.random()
		# y is the index of the species undergoing migration
		y= int(math.ceil(x*len(c_bor)))
		if len(c_bor)>0 and c_bor[y-1][1]==-1:
			print '\nMigration event:\nSpecies ', c_bor[y-1][0],'moves over from the boreal region\n'
			# call Migration function
			(c_bor, I_event_bor,c_neo, I_event_neo) = Migration(c_bor,y,f_tot,I_event_bor, c_neo, I_event_neo)
			#if len(c_bor)>0:	
			#	print c_bor[y-1],' (from boreal to neotropical, leaving a copy in boreal)\n'
			c_bor[y-1][1]=1
			#print c_bor, I_event_bor,'\n', c_neo, I_event_neo		
		else:
			print 'This species has already migrated or there are no species to migrate. No change in populations.\n'
	# if random is b/w (f_turnover+f_mrate)/P and (f_turnover+((i_count+1)*f_brate)+f_mrate)/P, 
	# it's a birth event in the neotropical region
    elif f_rand > (f_turnover+f_mrate)/P and f_rand <= (f_turnover+((i_count+1)*f_brate)+f_mrate)/P:
		# Pick a species randomly, it undergoes birth
        x= random.random()
		# y is the index of the species undergoing birth
        y= int(math.ceil(x*len(c_neo)))
        if len(c_neo)>0:
        	print '\nSpeciation event: '# undergone by',c_neo[y-1],'\n'
		# call Birth function
		(c_neo, I_event_neo) = Birth(c_neo,y,len(I_event_neo),f_tot,I_event_neo,"neo")
		#print c_neo, I_event_neo
	# else, it's a death event in the neotropical region 
	else:
		# Pick a species randomly, it undergoes extinction
		x= random.random()
		# y is the index of the species undergoing extinction
        y= int(math.ceil(x*len(c_neo)))
        if len(c_neo) > 0:
			print '\nExtinction event: '# undergone by',c_neo[y-1],'\n'
		# call Death function
        (c_neo, I_event_neo) = Death(c_neo,y,len(I_event_neo),f_tot,I_event_neo,"neo")
		#print c_neo, I_event_neo

print "\nFinal population in boreal region: ",c_bor,"\n\nFinal population in neotropical region:",c_neo
c_choice=raw_input("\nDo you wish to print out the lineage of a species (y/n)?")
#print c_choice
while c_choice == 'y':
	c_region=raw_input("From boreal or neotropical region (b/n)?")
	#print c_region
	if c_region =='b':
		print "You chose the BOREAL region."
		i_index=raw_input("Enter index of species of interest in boreal pop. array (indices start from 0):")
		i_index=int(i_index)
		if ((len(c_bor)-1)>=i_index):
			print "You picked species",c_bor[i_index][0]
			# call function
			Print_Ancestry(c_bor, I_event_bor,i_index)
		else:
			print "Wrong input - index exceeds size of array"
	else:
		print "You chose the NEOTROPICAL region"
		i_index=raw_input("Enter index of species of interest in neotropical pop. array (indices start from 0):")
		i_index=int(i_index)
		if ((len(c_neo)-1)>=i_index):
			print "You picked species",c_neo[i_index][0]
			# call function
			# if it was born in the neotropical region...
			if (c_neo[i_index][2]==1):
				Print_Ancestry(c_neo, I_event_neo,i_index)
			else:
			# it migrated from the boreal region
				print "This migrated from the boreal region. Moving to Boreal region..:"
				i_find=0
				i_trav=0
				# look for the migrated species in boreal region
				while((i_find==0)and(i_trav < len(c_bor))):
					if c_neo[i_index][0]==c_bor[i_trav][0]:
						i_find=1
					else:
						i_trav=i_trav+1
				if(i_find==0):
					print "Migrated species has become extinct in boreal region."
				else:
					Print_Ancestry(c_bor, I_event_bor,i_trav)
		else:
			print "Wrong input - index exceeds size of array"
		
	c_choice=raw_input("\nDo you wish to print out the lineage of a species (y/n)?")