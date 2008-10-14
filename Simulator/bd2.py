import math
import random

# The event map is a map from a unique Node to the array
# [Time of Event, Left Child ,Right Child, Parent node, region, migration_time]
# migration_time is -1 unless a migration has occured. Lineages and Nodes are the same
# (the lineage corresponds to the node at which it was created)
# Therefore, Left Child is the same as the parent, and Right child is the same as
# the key. However, these are included incase an extinction has occured
# then the extinct node/lineage will be -1.
# Also, migrations will only occur once. Therefore, the birth region of a lineage can be
# inferred from the current region and whether or not time = -1.
# Assume there is at first one lineage, lineage 0.


def Migration(y,time,events, direction=1):
    """
    Input Paramters
    ---------------
    y           the node migrating
    time        the time of migration
    events      the event array
                [ Time of Event (0), Left Child (1),Right Child (2), Parent node (3), region (4), migration_time (5)]
    direction   if 1, the migration is occuring boreal to tropical, if 0, the migration is occuring tropical to boreal

    Return Value(s)
    --------------
    updated events
    """
    if direction==events[y][4]:
        print "error, migrating to the state already in"
            
    events[y][4]=direction;
    events[y][5]=time;
    
    return(events)