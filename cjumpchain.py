#!/usr/bin/python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names, module and package names.
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

#   levels=tree.levels
#    current_level=levels[level_number]
#    event_history=current_level.event_history
#    lineage_names=event_history[0].keys()b
#   j=level_num
#    allcomb=[]
#    next_coalesce=CombinationChoose2(lineage_names,allcomb,j)

from is_classes import Level
from is_classes import Tree
from is_classes import Node
from levels import DemarcateLevels
from numpy import *           #sum, zeros, eye, array, allclose
from scipy.misc import * #comb, factorial 
from scipy.linalg import solve
import random
import sys
import math

#These strings are supposed to be constants. 
SMALLER_S_BT = "smaller_s_bt"
SMALLER_S_TT = "smaller_s_tt"
SMALLER_S_BB = "smaller_s_bb"
KAPPA = "kappa"
M_1 = "m_1"
M_2 = "m_2"
SMALLER_S_B_ARROW_T = "smaller_s_b_arrow_t"
S_BT = "s_bt"
S_BB = "s_bb"
S_TT = "s_tt"
S_B_ARROW_T_GLOBAL = "s_b_arrow_t_global"
 
#sys.path=[sys.path,'/Network/Servers/orrorin.nescent.org/Volumes/xraid/home/alexandrabalaban/Documents/Summer2008']
#sys.path=[sys.path,'C:\Documents and Settings\Perry Zheng\workspace\Summer2008\src']    
#sys.path.append('C:\Users\Lonnie\Desktop\Summer2008_2\Summer2008\src')
#sys.path.append('C:\Users\Lonnie\Desktop\Summer2008_2\Summer2008\src')
#sys.path.append('C:\Documents and Settings\Owner\Desktop\Summer2008\src')
def Equal(a, b):
    """
    For testing purposes only, equal if float a and float b are close 
    enough they're identically equal
    """
    return abs(a-b) < a*1e-7
    
def CombinationChoose2(list,allcomb,j):
    """
    Returns all the j choose 2 combinations of a list of j
    items.

    Input parameters
    ---------------
    list        the list from which to create combinations
    allcomb     an empty list which the function will return filled
                with the combinations
    j           the length of list

    Return value
    ------------
    a list containing all possible j choose 2 combinations
    
    Details
    -------
        Input: list=[1,2,3], allperm=[], j=3
        Output: [[1, 3], [2, 3], [1, 2]]
    """
    if j==1:
        return allcomb
    j=j-1     
    for i in range(j):
        allcomb.append([list[i],list[j]])
        CombinationChoose2(list,allcomb,j)

    return -1        
    
def LogFactorialOfNegative(z,terms=20):
    """
    Returns an approximation of the log(factorial) any real number that is not a non-positive integer
    Input paramters
    ---------------
    z       any number that is not a non-positive integer
    terms   the number of terms for which to carry out the approximation

    Return value
    ------------
    An approximation of the factorial.  Uses the Euler infinite product
    definition of the Gamma function.  Returns the log(factorial(z))

    Details
    -------
    If z+1 is negative: Make z positive for the calculations. Then if floor(z) is even, the approxiimation
    would have been negative if carried out to a number of terms greater than floor(z),
    so return 1/approx. Otherwise, if floor(z) is odd, the approximation would have
    been positive if carried out to a number of terms greater than floor(z), so return approx.
    I.E. log((-a)^even_num)=even_num*log(a); log((-a)^odd_num)=1/(odd_num*log(a))
   """

    z=float(z)
    z=z+1
    if(z>0):
        fact=math.log(1/(z),10)
        for n in range(terms):
            n=n+1
            n=float(n)
            fact=fact+math.log(float(pow((1+1/n),z))/float(1+z/n),10)
        return fact
    if(z==0):
        return 1
    if(z<0):
        z=-z
        #above line * for future reference
        fact=math.log(1/(z),10)
        for n in range(terms):
            n=n+1
            n=float(n)
            fact=fact+math.log(float(pow((1+1/n),z))/float(1+z/n),10)
        if(math.fmod(math.floor(z),2)==0):
            return 1/fact
            #because there is one more negative term in *
        else:
            return fact
            #because * makes the total number of negative #s multiplied even

def CombOfNegative(n,k,terms=20):
    """
    An approximation of the Combination function using FactorialOfNegative

    Input
    -----
    n   any real number that is not a non-negative integer and such that n-k is not
        a non-negative integer
    k   any real number that is not a non-negative integer and such that n-k is not
        a non-negative integer
        
    
    Output
    ------
    n choose k

        
    """
    log_n_fact=LogFactorialOfNegative(n,terms)
    log_k_fact=LogFactorialOfNegative(k,terms)
    log_n_fact_over_k_fact=log_n_fact-log_k_fact
    log_n_minus_k_fact=LogFactorialOfNegative(n-k,terms)
    #print(str(pow(10,log_n_fact_over_k_fact-log_n_minus_k_fact))+" comb of Neg");
    return  pow(10,log_n_fact_over_k_fact-log_n_minus_k_fact)

def GetEvents():
    """
    Returns sets of various types of events (or equivalently, state
    transitions) that could happen in the conditional jump chain described
    in Section 5.3 of mossEqModel.pdf. 

    No Input Parameters

    Return value
    ------------
    a 3-tuple consisting of 3 sets: (bigbig_K, bigbig_M, bigbig_H).
    bigbig_K is the set of various types of speciation events.
    bigbig_M is the set of various types of migration events.
    bigbig_H is the union of the above two sets.
    See Page 19, segment entitled "Events as operators on states". 
    Also see 'Details' below.

    Details
    -------

    Possible transitions (i.e., events) of the condition jump chain. 
    
    See Page 19, segment entitled "Events as operators on states". 
    The events can be one of the following:
   
    1. a speciation event involving the two next-coalescing lineages
    
    2a, 2b, 2c. a speciation that does not involve one or both of the
                next-coalescing lineages. These can be three types,
                as noted in Section 3. (Note that in page 19, segment
                titled 'Events as operators on states, 2a, 2b and 2c
                are all lumped into one category, which is not quite
                accurate. My apologies - Ganesh).
   
    3a, 3b. a migration (i.e., state transition) of one of the two
            next-coalescing lineages.
   
    4. Any other migration event.

    1, 2a, 2b, 2c together will make bigbig_K. 
    3a, 3b and 4 will make bigbig_M 
    bigbig_H is the union of the above two sets.
    """
    #sys.path=[sys.path,'/Network/Servers/orrorin.nescent.org/Volumes/xraid/home/alexandrabalaban/Documents/Summer2008']
    #sys.path=[sys.path,'C:\Documents and Settings\Perry Zheng\workspace\Summer2008\src']    
    # 'kappa' is the speciation
    # involving the two next-coalescing species. The other three follow
    # notation from Section 3, bullet point entitled 'Events', page 5.
    # Also see Figure 1, page 6.
    # set() is a python command for creating a set from a list.
    
    #Notice that I erased 's_tt', 's_bb', 's_bt' b/c
    #kappa+smaller_s_tt = s_tt
    #kappa+smaller_s_bt = s_bt
    #kappa+smaller_s_bb = s_bb  
    bigbig_K = set([KAPPA, SMALLER_S_TT, SMALLER_S_BT, SMALLER_S_BB])

    # M_1 is the migration of next_coalescing_lineage1: i.e., lineage whose
    # state is state_of_next_coalescing_lineages_1;
    # M_2 is the migration of next_coalescing_lineage1: i.e., lineage whose
    # state is state_of_next_coalescing_lineages_2;
    # 's_b_arrow_t' is any other migration
    #Again, SMALLER_S_B_ARROW_T = m_1+m_2
    bigbig_M = set([M_1, M_2, SMALLER_S_B_ARROW_T])

    # bigbig_H = bigbig_K union bigbig_M. 
    bigbig_H = bigbig_K.union(bigbig_M)

    return((bigbig_K, bigbig_M, bigbig_H))

def MakeTupleList(level_number):
    """
    Creates a list of all the possible 4-tuples states for a given level
    and then creates a one-to-one correspondence between state numbers and
    the all possible 4-tuple states
 

    Input parameters
    ----------------  
    1. level number    The level for which the dictionary is computed

    
    Return value (s)
    ----------------    
    A list of 4-element tuple (q_t, r_t, x_1, x_2)
    
    """
   
    tuples=range(4*(level_number-1))
    j=0
    for i in range(level_number+1):
        if((i!=0)&(i!=1)):
            tuples[i-2]=(i,level_number-i,0,0)
            #print(i-2)

    for i in range(level_number+1):
        if((i!=0)&(i!=level_number)):
            j=(i-1)+1*(level_number-1)
            #print(j)
            tuples[j]=(i,level_number-i,0,1)
        
    for i in range(level_number+1):
        if((i!=level_number-1)&(i!=level_number)):
            j=i+2*(level_number-1)
            #print(j)
            tuples[j]=(i,level_number-i,1,1)

    for i in range(level_number+1):
        if((i!=0)&(i!=level_number)):
            j=(i-1)+3*(level_number-1)
            #print(j)
            tuples[j]=(i,level_number-i,1,0)
        
    return(tuples)

def MakeState2IndexDictionary(level_number):
    """
    Returns a dictionary (or map) between the state (key) and the index (value)
    for a given level
    The state is the four-element tuple (q_t, r_t, x_1, x_2), represented by (). 
    Since tuples are immutable objects, they can be used as keys. 
    
    The index is the numerical index used in MakeTupleList.
     
    Input parameters
    ----------------

    1.level number    The level for which the state to index dictionary
                        is computed

    Details
    --------
    So when this method is called: 
    x = MakeState2IndexDictionary(2)
    print x      #Prints {'0,2,0,1': 3, '0,2,0,0': 0, '2,0,1,0': 11, 
                '2,0,1,1': 8, '1,1,0,1': 4, '1,1,0,0': 1, '0,2,1,0': 9, 
                '0,2,1,1': 6, '1,1,1,0': 10, '1,1,1,1': 7, '2,0,0,1': 5, 
                '2,0,0,0': 2}

    y = x['1,1,1,0'] 
    print y      #Prints 10
    
    Return value (s)
    ---------------

    A dictionary matching tuples representing states of the conditional
    jump chain to indices in the transition matrix
    """
    tuple_list=MakeTupleList(level_number)
    state_to_index_dictionary={}
    for i in range(len(tuple_list)):
        state_to_index_dictionary[tuple_list[i]]=i        
        
    return(state_to_index_dictionary)

def ReverseMap(level_number):
    """
    Returns a dictionary (or map) between the index (key) and the state (value)
    for a given level
    The index is the numerical index used in the list MakeTupleList.
    The state is the four-element tuple (q_t, r_t, x_1, x_2)
     
    Makes and returns the index to state dictionary for a given level

    Input parameters
    ----------------

    1.level number    The level for which the state to index dictionary
                        is computed

    
    Return value (s)
    ---------------

    A dictionary matching indices in the transition matrix to tuples representing states of
    the conditional jump chain
    """
    tuple_list=MakeTupleList(level_number)
    index_to_state_dictionary={}
    for i in range(len(tuple_list)):    
        index_to_state_dictionary[i] = tuple_list[i]
        
    return(index_to_state_dictionary)

def GetCondJumpChainStateSpace(level_number):
    """
    Input parameters
    -----
    level_number   The level for which the conditional jump chain space is calculated

    Output
    ------
    A list of all the possible states, where the state's index in the array corresponds to the state's index
    to the four-element tuple (called "state", i.e. (q_t,r_t,x_1,x_2) 
    in state_to_index_in_transistion_matrix dictionary. 
    
    Details
    -------
    Note in here the "state" means the four-element tuple (q_t,r_t,x_1,x_2) where
    q_t and r_t are the species in boreal and neotropical regions, respectively. And 
    x_1 and x_2 are the "transition states" of the two next-coalescing species. 

    Thus "transition states" here refer to boreal (state 0) and neotropical (state 1)
    regions. Going FORWARD in time, only transition from state 0 to state 1 is 
    possible. Going BACKWARD in time, only transition from state 1 to state 0 is 
    possible.
    
    From now on we will call such state the "four-element tuple" to avoid confusion
    with "transition states." 
    """
    
    state_space=range(4*(level_number-1))
    index_to_state_map=ReverseMap(level_number)
    for i in range(len(state_space)):
        state_space[i]=index_to_state_map[i]
    
    return(state_space)

def IsValidState(z,level):
    """
    Tests whether the four-element tuple is a state in the jump chain. 
    
    Input Parameters
    ----------------
    z        A state in the jump chain in the current_level, 
             or the four-element tuple, (q_t, r_t, x_1, x_2). See ApplyEvent()
             for more details.
    level    The current level
             
    Details
    --------
    Some sample states that are NOT valid include: 
    (5,0,1,0),  (4,1,1,1), (0,2,0,1), (5,1,1,1)
    
    Note we assume x_1 and x_2 are correctly from (0,1). This "validity"
    refers more to states that not possible due to configuration, not
    due to misspelling or typos. For instance, (5,0,1,0) is an invalid state
    b/c next-coalescing lineage1's state is 1 and there is no state 1 overall.
    
    Return Value
    ------------
    True if the state is a valid state, false otherwise. 
    """
    q_t = z[0]
    r_t = z[1]
    x_1 = z[2]
    x_2 = z[3]
    if (x_1==1 and x_2==1 and r_t<2):
        return False
    if (x_1==0 and x_2==0 and q_t<2):
        return False
    if (x_1==1 and r_t<1):
        return False
    if (x_2==1 and r_t<1):
        return False
    if (x_1==0 and q_t<1):
        return False
    if (x_2==0 and q_t<1):
        return False
    if (x_1>1 or x_1<0 or x_2>1 or x_2<0):
        return False
    if(q_t<0 or r_t<0):
        return False
    if(q_t+r_t>level):
        return False
    return True;

def ApplyEvent(z,w,level):
    """
    Returns the state that results when event w is applied to state z, i.e., 
    w(z) in eqn (24). The result is also a four-element tuple. 
    
    Input Parameters
    ----------------
    level    The current level
    z        A state in the jump chain in the current_level, 
             or the four-element tuple, (q_t, r_t, x_1, x_2) where
             q_t is the number of current lineages in boreal region (state 0)
             r_t is the number of current lineages in neotropics region (state 1)
             x_1 is the state transition (0 or 1) of the first of the 
                     next-coalescing lineages pair
             x_2 is the state transition (0 or 1) of the second of the 
                     next-coalescenct lineages pair                     
    w        w is an event that acts on the state in the jump chain. As discussed in 
             page 20 eqn (24), we could view an event as an operator on the 
             state. w could be: 
             Speciation Events (collectively known as bigbig_K): 
             KAPPA     The coalescence (convergence going backwards in time, 
                         or speciation going forward in time) of next-coalescing 
                         lineages
             SMALLER_S_TT      The coalescence of lineages in (t,t,), (b,b) and (b,t) which aren't the next coalescing lineages.
             SMALLER_S_BB      Note1: (b,t) is really a pseudo speciation event.
             SMALLER_S_BT      Note2: these events cannot occur in the current model
             
             Migration Events (collectively known as bigbig_M): 
             M_1           migration of next-coalescing lineage1, the first lineage
                             in the next-coalescing lineage pair
             M_2           migration of next-coalescing lineage2
             "s_b_arrow_t"   any other migration
             
             We denote bigbig_H = bigbig_K UNION bigbig_M
             So w could be any of the events in bigbig_H.             
            
    Details
    --------
    Note: This is NOT P(F | w(z))!!!. This is simply the state w(z) in the 
    expression. We calculate probability based on what we know about w(z). 
    
    When the method returns (-1, -1, -1, -1) it means it w(z) is not possible -
    it is not possible to go from z to w(z). 
    
    Return values
    -------------
    The state that results when an event happens to state z. Also a four-element
    tuple, just like z.
    
    Details
    -------
    IMPORTANT:SMALLER_S_BT, SMALLER_S_BB, and SMALLER_S_TT can never occur in the backwards probability calculation.
    Therefore, return (-1,-1,-1,-1) as a signifier of an error. Also, return (-1,-1,-1,-1)
    if the input current state, z, is not valid, or if w, the transition, is not a valid event.
    """
    q_t = z[0]
    r_t = z[1]
    x_1 = z[2]
    x_2 = z[3]
    if (IsValidState(z,level)==False):
        return (-1,-1,-1,-1)
    
    if (w==SMALLER_S_BT or w==SMALLER_S_BB or w==SMALLER_S_TT):
        #Since these events while going backwards have probability zero in the current model, return (-1,-1,-1,-1)
        return (-1,-1,-1,-1)
    
    if (w==KAPPA):
        if (x_1==1 and x_2==1):
            new_state= (q_t, r_t-1, -1, -1)
            #(q_t, r_t, x_1, x_2) --> (q_t, r_t-1, UNDEFINED, UNDEFINED)
            if(q_t>-1 and r_t-1 >-1):
                return new_state
            else:
                return (-1,-1,-1,-1)
    
        if ((x_1==1 and x_2==0) or (x_1==0 and x_2==1)):
            new_state=(q_t,r_t-1, -1, -1)
            if(q_t>-1 and r_t-1 >-1):
                return new_state
            else:
                return (-1,-1,-1,-1)

        if (x_1==0 and x_2==0):
            new_state= (q_t-1, r_t, -1, -1)
            if(q_t-1 >-1 and r_t >-1):
                return new_state
            else:
                return (-1,-1,-1,-1)
    
    if(w==M_1):
        if(x_1==1):
        #(q_t,r_t,x_1,x_2) --> (q_t,r_t,0,x_2)
            tuple = (q_t+1, r_t-1, 0, x_2)
            if (IsValidState(tuple,level)==False):
                return (-1,-1,-1,-1)
            else:
                return tuple
        else:
            #x_1 must have state 1
            return (-1,-1,-1,-1)
        
    if(w==M_2):
        if(x_2==1):
            tuple = (q_t+1,r_t-1,x_1,0)
            if (IsValidState(tuple,level)==False):
                return (-1,-1,-1,-1)
            else:
                return tuple
        else:
            #x_2 must have state 1
            return (-1,-1,-1,-1)            
            
    if(w==SMALLER_S_B_ARROW_T):
        #IMPORTANT!!!!
        #Here I take SMALLER_S_B_ARROW_T to mean NOT M_1 NOR M_2, that is: 
        #M_1, M_2 and "smaller_s_b_arrow_t" form an disjoint set!
        #Only true for conditional probability. For uncondtional ones,
        #M_1 and M_2 ARE s_b_arrow_t events. But for 
        #conditional cases, "s_b_arrow_t" really means those "s_b_arrow_t"
        #events EXCLUDING M_1 and M_2. 
        
        tuple = (q_t+1,r_t-1,x_1,x_2)
        if (IsValidState(tuple,level)==False):
            return (-1,-1,-1,-1)
        return tuple
       
    return (-1,-1,-1,-1)
    # w was not a valid event

def Make2DimensionalArray(dimension):
    """
    Make a 2-dimensional square array of dimensions with all
    slots initialized to 0.0. Number of rows should equal to 
    the number of columns in the array. 
    
    Input Parameters
    ----------------
    dimension        An integer, indicating the "width" of the square array
    
    Details
    --------
    This is used to initialize the transition_matrix used in 
    MakeTransitionMAtrixForLevel method. 
    
    An two-dimension array is really a list of a list. 
    Return Value
    ----------------
    A two-dimensional square array with all slots initialized to 0.0. 
    """
    #Old implementation, now obsolete but still want to keep the refernece:
    #Taken from http://mail.python.org/pipermail/python-list/2000-December/063456.html    
    #matrix= [([0.0] * dimension)[:] for i in range(dimension)]  
    
    #New implementation, using zeros method in Numpy
    matrix = zeros( (dimension,dimension) )  
    #Explanation: [([2] * 5)[:] for i in range(3)] creates an the matrix 
    #[[2, 2, 2, 2, 2], [2, 2, 2, 2, 2], [2, 2, 2, 2, 2]]
    #that is, a 3-element list, each of whose element is a list [2]*5. 
    #Note that this creation ensures that each [2,2,2,2,2] is independent of 
    #each other.
    return matrix

def GetLinearEquations(state_space_at_current_level,sigma, level):
    """
    Output a constant vector (made up of kappa_ij 's), and 
    a dimension*dimension matrix each row of which 
    composes of coefficients of the unknown (which in this case is
    P(F|z)
    
    Input Parameters
    -----------------
    state_space_at_current_level    
                                A list containing
                                all the states (4-tuples)
                                at the current level (e.g. level 5), 
                                really generated by the
                                GetCondJumpChainState(level) method
                                See comments in MakeTransitionMatrixForLevel
                                for more details.
    level                       the current level
    

    Output value
    -------------
    Return (A,b) where A is a level_number*level_number matrix, and b is a constant vector
    (i,j) of A represents a transition from state with index i to state with index j
    
    Details
    --------
    Get the necessary ingredients about the systems
    of equations to plug into the LinearEquationSolver
    
    The systems of equations are:
    
    for all z in Z_k, 
    P(F|z) = sum_{w in M} P(w|z)*P(F|w(z)) + P(k_ij | z)
    
    where P(F|z) and P(F|w(z) are unknowns (could be different or same), 
    and P(k_ij | z) and P(w|z) are knowns
    
    Now, by our CONSTRUCTION (more like assumption), 
    (M_1, M_2, "s_b_arrow_t") are DISJOINT SETS in 
    conditional probability cases, we have, after rearranging terms:
     
    P(k_ij | z) = P(F|z) - sum_{w in M} P(w|z) * P(F|w(z)) for all z in Z_k
                = P(F|z) - P(M_1|z)*P(F|M_1(z)) - P(M_2|z)*P(F|M_2(z))
                    - P("s_b_arrow_t"|z)*P(F|"s_b_arrow_t"(z))
    
    There are 4(k-1) such z's, because w(z) WILL be inside
    M by our construction, and if it's not, we simply ignore it in our 
    computation (b/c we would get (-1,-1,-1,-1) if it's not
    in the current state or an invalid state from ApplyEvent()). And 
    when m(z) is undefined, we would treat P(w|z) as 0 
    
    So we can create a 4(k-1) by 4(k-1) matrix called A, 
    filled with coefficients of various P(F|z)'s in each row. 
    
    How do we fill out matrix A? 
    
    For instance, suppose there are 5 equations of the following format
    x_0 = 4 + 3x_0 + 11x_3 + 8x_4         
    x_1 = 8 + 2x_1 + 10x_3
    x_2 = 17
    x_3 = 25 + 3x_1 + 12x_3 + 22x_2
    x_4 = 0 + 3x_3
    
    Rearranging, we get:
    4 = x_0 - 3x_0 - 11x_3 - 8x_4 
    8 = x_1 - 2x_1 - 10x_3
    17 = x_2
    25 = x_3 - 3x_1 - 12x_3- 22x_2
    0 = x_4 - 3x_3
    
    We want to eventually turn it into the form that: 
    [4,8,17,25,0].transpose() = [[-1, 0, 0, -11, -8],  * [x_0,x_1,x_2,x_3,x_4].transpose
                                 [0, -1, 0, -10, 0],
                                 [0, 0, 1, 0, 0],
                                 [0, -3, -22, -11, 0],
                                 [0, 0, 0, -3, 1]]
                                 
    We call the LHS b, and the 5x5 matrix A, we would output (A, b)                             
    """ 
    length = len(state_space_at_current_level)
    constant_vector = [UnconditionalTransitionProbability(KAPPA,z,sigma) for z in state_space_at_current_level]
    #when level_number=2, 
    #we have state_space_at_current_level = [(2, 0, 0, 0), (1, 1, 0, 1), (0, 2, 1, 1), (1, 1, 1, 0)]
    #and constant_vector = [1.0, 0.00049975012493753122, 0.058823846204939162, 0.00049975012493753122]
    
    #First create a level by level matrix with diagonal
    #initialized to 1's and rest 0's Fill in the diagonals with 1's
    #We really want the coefficients to be ints. 
    matrix = eye(length,length, dtype=float)
    
    #Fill in the coefficient P(M_1|z) into cell matrix[index of z][index of w(z)]

    #First, create a Probability of w given Z vector
    prob_m1_vector = [-UnconditionalTransitionProbability(M_1, z, sigma) for z in state_space_at_current_level]
    prob_m2_vector = [-UnconditionalTransitionProbability(M_2, z, sigma) for z in state_space_at_current_level]
    prob_s_b_arrow_t_vector = [-UnconditionalTransitionProbability(SMALLER_S_B_ARROW_T, z, sigma) for z in state_space_at_current_level]
    #print "prob_m1_vector is: " + str(prob_m1_vector)
    #print "prob_m2_vector is: " + str(prob_m2_vector)
    #print "prob_s_b_arrow_t_vector is: " + str(prob_s_b_arrow_t_vector)
    
    #For instance, when level=2, prob_s_b_arrow_t is
    #[0.0, 0.99950024987506236, 0.94117615379506092, 0.99950024987506236]
    #Remember prob_w_given_z_vector does not depend on (x_1,x_2)
    #This way, prob_w_given_z[3] refers to the probability of
    #w_given_z_vector for the third element in state_space_at_current_level, 
    #or really state_space_at_current_level[2]
    
    m1_vector = [ApplyEvent(z,M_1,level) for z in state_space_at_current_level]
    m2_vector = [ApplyEvent(z,M_2,level) for z in state_space_at_current_level]
    s_b_arrow_t_vector = [ApplyEvent(z,SMALLER_S_B_ARROW_T,level) for z in state_space_at_current_level]
    #print m1_vector
    #print m2_vector
    #print s_b_arrow_t_vector
    #For instance, the following code
    #state_space_at_current_level = GetCondJumpChainStateSpace(2)
    #print state_space_at_current_level
    #s_b_arrow_t_vector = [ApplyEvent(z,"s_b_arrow_t") for z in state_space_at_current_level]
    #generates:
    #state_space_at_current_level = [(2, 0, 0, 0), (1, 1, 0, 1), (0, 2, 1, 1), (1, 1, 1, 0)]
    #s_b_arrow_t_vector = [(-1, -1, -1, -1), (2, 0, 0, 1), (1, 1, -1, -1), (2, 0, 1, 0)]
    #So s_b_arrow_t_vector is an "ordered" list of the w(z)'s. 
    
    level_number = length/4 + 1
    state_to_index_map = MakeState2IndexDictionary(level_number)
    #print "state_to_index_map[(1,2,0,1)] is: " + str(state_to_index_map[x])
    assert(length == len(state_space_at_current_level) == len(prob_m1_vector) == len(m1_vector))
    
    for row in range(length):
        #for each row in the matrix, only fill in 3 numbers
        #So this implementation is faster than going through the
        #entire matrix
        if (m1_vector[row]!=(-1,-1,-1,-1)):
            #Otherwise, keep it 0
            index_of_m1 = state_to_index_map[m1_vector[row]]
            matrix[row,index_of_m1] = prob_m1_vector[row]
        
        if (m2_vector[row]!=(-1,-1,-1,-1)):
            #Otherwise, keep it 0
            index_of_m2 = state_to_index_map[m2_vector[row]]
            matrix[row,index_of_m2] = prob_m2_vector[row]

        if (s_b_arrow_t_vector[row]!=(-1,-1,-1,-1)):
            #Otherwise, keep it 0
            index_of_s_b_arrow_t = state_to_index_map[s_b_arrow_t_vector[row]]                                      
            matrix[row,index_of_s_b_arrow_t] = prob_s_b_arrow_t_vector[row]
    
    return (matrix, constant_vector)

def LinearEquationSolver(A, b):
    """
    Solve system of linear equations given A*x = b, outputs x in vector
    form. 
    
    Input parameters
    ----------------
    A        A square matrix 
    b        A column vector
    
    Output Value
    ------------
    Solution to the system of equations A*x = b. Outputs x which is a 
    vector. 
    """
    return solve(A,b)

def CalculatePi(j_species,sigma):
    """
    Calculates p_j, or the probability that the equilibrium number of species equals j given
    sigma
    
    Input Parameters
    ----------------
    @param j_species    j species; a nonnegative integer
    @param sigma    set of all parameters for the model
                    For the model in mossEqModel.pdf, 
                    sigma = (lambda, b, B, alpha, mu), where 
                    
                    lambda is the turnover rate in state 0 (the boreal region), 
                    b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                    mu = the per lineage extinction rate in state 1 (the neotropical region).

    Return Value
    ----------
    p_j, or the probability that the equilibrium number of species equals j given
    
    Details
    --------    
    Only 3 components of sigma are needed:
    alpha     sigma[3] = per species birth rate in neotropical region; positive constant
    b         sigma[1] = migration rate
    mu        sigma[4] = death rate of species in neotropical region when number of species is i; positive constant

    if log_pi=-1, then pi=0 (since the log of 0 is undefined).
    """
    log_pi=CalculateLogPi(j_species,sigma)
    alpha = sigma[3]
    mu = sigma[4]
    if (log_pi == -1): 
        return 0;      
    return(pow(alpha/mu,log_pi))
    
def CalculateLogPi(j_species, sigma):
    """
    Calculates p_j, or the probability that the equilibrium number of species equals j given
    sigma
    
    Input Parameters
    ----------------
    @param j_species        
                    j species; a nonnegative integer
    @param sigma    set of all parameters for the model
                    For the model in mossEqModel.pdf, 
                    sigma = (lambda, b, B, alpha, mu), where 
                    
                    lambda is the turnover rate in state 0 (the boreal region), 
                    b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                    mu = the per lineage extinction rate in state 1 (the neotropical region).
    
    Return value
    -------------
    The log base alpha/mu of p_j, where p_j is defined as
    p_j = C(b/alpha + j - 1,j)*(alpha/mu)^j*(1-alpha/mu)^(-b/alpha)
    if alpha=mu, -1, a flag is returned (since log of zero, and log base 1 is undefined)

    Details
    --------    
    Only 3 components of sigma are needed:
    alpha     sigma[3] = per species birth rate in neotropical region; positive constant
    b         sigma[1] = migration rate
    mu        sigma[4] = death rate of species in neotropical region when number of species is i; positive constant
    

    if alpha=mu, pi will be always zero, and CalculateLogPi returns -1, a flag that alpha=mu
    
    We want to calculate pi_j = C(b/a+j-1,j)*(a/u)^j*(1-a/u)^(-b/a)
    so we take log_a/u of pi_j, which prevents (a/u)^j from becoming
    0 when j is large (i.e. when j=1500). 
    
    Thus x = log_{a/u} pi_j = log_{a/u} C(b/a+j-1,j) + j + log_{a/u} (1-a/u)^(-b/a)
    
    This functions returns x. 
    
    CalculatePi returns pi_j = (a/u)^x, by simply calling
    this helper class.  
    """
    alpha = sigma[3]
    b = sigma[1]
    mu = sigma[4]
    coefficient=1;
    if(alpha>mu):

        #print("Error, alpha is greater than mu")        
        return -2
    if (alpha== mu):
        return -1
    if(b<alpha):
        combination= float(CombOfNegative(b/alpha+j_species-1,j_species))
    else:
        combination = float(comb(b/alpha+j_species-1,j_species))

    #print(combination);
    #print(alpha/mu);
    try:
        coefficient = math.log(combination,alpha/mu)
    except OverflowError:
        print("beta "+str(b));
        print("j_species "+str(j_species));
        print("alpha "+str(alpha));
        print("mu "+str(mu));
        print("combination "+str(combination)+" base "+str(alpha/mu));
        
    ratio = float(j_species)
    ratio2 = -b/alpha*math.log(pow(float(1-alpha/mu),1),alpha/mu)
    return coefficient+ratio+ratio2

#I got this function from R documentation
def LogSum(logX,logY,base):
    """
    Input Parameters
    ----------------
    logY   the log of x with base=base, i.e., log_base X
    logY   the log of y with base=base, i.e., log_base Y
    base    the base of the logs of x and y

    Return Value
    ------------
    log_base (X+Y)

    Details
    -------
    Note: log_x and log_y must share a common base. 
    For instance, if logX = 3 = log_2 8
                     logY = 4 = log_2 16
    then LogSum(3,4,2) = log_2 (24) = 4.585
    
    This equality and function were obtained from the R function logSum in the package TileHMM 
    """
    temp = 1 + pow(base,logY-logX)
    return logX + math.log(temp,base)

def CalculateLogPiStar(n_t, current_state_for_uncond_probs, sigma, num_sum=20):
    """
    Calculates p(n_t | q_t, r_t) = pi_star, or eqn (10) in mossEqModel.pdf
    
    Input parameters
    -----------------
    @param n_t      total number of neotropical species in the POPULATION
                    at time t
    @param current_state_for_uncond_probs    
                    the current state for unconditional probability
                    current_state_for_uncond_probs[0] = q_t : number of ancestral species in boreal in sample history at time t
                    current_state_for_uncond_probs[1] = r_t : number of ancestral species in neotropics in sameple history at time t              
                    Note: only [1]=r_t is needed.
    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum
    @param sigma    set of all parameters for the model
                    sigma = (lambda, b, B, alpha, mu)
                    See CalculateLogPi for more details. 
                    
    Details
    --------
    Only 3 components of sigma are used:
    alpha    sigma[3] = the speciation rate in state 1 (the neotropical region)
    b        sigma[1] = the migration rate
    mu       sigma[4] = the extinction rate in state 1 (the neotropical region)
    
    See CalculateLogPi for theoretical details. 
    
    Here we are calculating 
    LogPiStar = log_{a/u} pi_j - log_{a/u} (sum of a bunch of pi's)
              = log_{a/u} pi_j / (sum of a bunch of pi's)
    
    Return value
    -------------
    The log of piStar, where piStar is defined as:
    p(n_t | q_t, r_t)= pi_n_t / sum_{k>=r_t}(pi_k); 
    Note that piStar is called the stead-state frequency of n_t, normalized to condition
    on the fact that at least r_t neotropical species are known to exist at time t. 
    """
    r_t = current_state_for_uncond_probs[1]
    alpha = sigma[3]
    mu = sigma[4]
    log_numerator = CalculateLogPi(n_t,sigma)
    log_denominator = CalculateLogPi(r_t,sigma)
    for i in range(r_t+1, r_t+num_sum):
        log_denominator=LogSum(log_denominator,CalculateLogPi(i,sigma),alpha/mu)
    
    return log_numerator-log_denominator

def CalculatePiStar(n_t, current_state_for_uncond_probs, sigma, num_sum=20):
    """
    Calculates p(n_t | q_t, r_t) = pi_star, or eqn (10) in mossEqModel.pdf
    
    Input parameters
    -----------------
    @param n_t      total number of neotropical species in the POPULATION
                    at time t
    @param current_state_for_uncond_probs    
                    the current state for unconditional probability
                    current_state_for_uncond_probs[0] = q_t
                    current_state_for_uncond_probs[1] = r_t              
                    Note: only [1]=r_t is needed.
    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum
    @param sigma    set of all parameters for the model
                    sigma = (lambda, b, B, alpha, mu)
                    See CalculatePi for more details. 
                    
    Details
    --------
    Only 3 components of sigma are used:
    alpha    sigma[3] = the speciation rate in state 1 (the neotropical region)
    b        sigma[1] = the migration rate
    mu       sigma[4] = the extinction rate in state 1 (the neotropical region)
    
    This returns the piStar value using the helper method 
    CalculateLogPiStar
    
    Return value
    -------------
    p(n_t | q_t, r_t)= pi_n_t / sum_{k>=r_t}(pi_k); a float
    This is called the stead-state frequency of n_t, normalized to condition
    on the fact that at least r_t neotropical species are known to exist at time t. 
    """  
    alpha = sigma[3]
    mu = sigma[4]
    logPiStar = CalculateLogPiStar(n_t, current_state_for_uncond_probs, sigma, num_sum=20)
    return pow(alpha/mu,logPiStar)

    """
    The old way of calculating PiStar
    Note that CalculatePiStar is the "new", logarithmic way of testing
    
    For testing purposes only. 
    
    Details
    ------------
    current_state_for_uncond_probs = [1500,1500]
    num_sum=20
    #sigma = (lambda, b, B, alpha, mu)
    sigma=(.1, .1, 1000, .05, .1)
    print CalculateLogPi(1500,sigma)
    numerator = CalculatePi(1500,sigma)
    print numerator==0    //True! That's why we use logs instead of 
                            adding up all the Pi's, b/c numerator
                            and denominators are all 0s. 
    """
    
    r_t = current_state_for_uncond_probs[1]
    
    numerator = CalculatePi(n_t,sigma)
    denominator=0.0
    
    for i in range(r_t, r_t+num_sum): 
        denominator=denominator+CalculatePi(i,sigma)
    
    return numerator/denominator     
    
def UncondProbSTTGlobal(current_state_for_uncond_probs, sigma, num_sum=20):
    """
    The s_tt global, that includes s_tt and kappa. 
    Calculates p(S_TT|q_t,r_t), probability of speciation event within 
    the neotropics that is captured in the history of the sample, 
    conditioned on the fact that at time t the numbers of ancestral 
    species in boreal and neotrpoics regions are q_t and r_t, respectively. (actually, only r_t matters)
    See eqn (13) in mossEqnModel.pdf
    
    Input parameters
    ----------------
    @param current_state_for_uncond_probs    
                    the current state for unconditional probability
                    current_state_for_uncond_probs[0] = q_t 
                    current_state_for_uncond_probs[1] = r_t              
                    Only [1] = r_t is needed for this method. 
    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum
    @param sigma    sigma = (lambda, b, B, alpha, mu), where  
                    lambda is the turnover rate in state 0 (the boreal region), 
                    b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                    mu = the per lineage extinction rate in state 1 (the neotropical region).
 
    Details
    ---------
    Only 3 components of sigma are used (some are used ONLY to call CalculatePiStar and are
    marked with *): 
    alpha    sigma[3] = the speciation rate in state 1 (the neotropical region)
    b        *sigma[1] = the migration rate
    mu       *sigma[4] = the extinction rate in state 1 (the neotropical region)
       
    Return value
    ------------ 
    p(S_TT|q_t,r_t) = 2*alpha*C(r_t,2)*\sum_{k>=r_t}(PiStar_{k-1}/k)
    """
    r_t=current_state_for_uncond_probs[1]
    alpha=sigma[3]
    #Trivia: Note that comb changes the number from type "float" to "numpy.float64"
    #because comb is a method of numpy class
    coefficient = 2*alpha*comb(r_t,2)
    sum=0
    for k in range(r_t,r_t+num_sum):
        if (k!=0):
            #In case k=0, we don't want division by 0
            sum=sum+ float(CalculatePiStar(k-1,current_state_for_uncond_probs,sigma,num_sum))/k
    return coefficient*sum

def UncondProbSBBGlobal(current_state_for_uncond_probs, sigma, num_sum=20):
    """
    The s_bb global that incluces regular s_bb and kappa
    Calculates p(S_BB | q_t), the probability of speciation event in the boreal region 
    that is captured in the history of the sample. This is simple 
    since boreal region maintains a constant fixed number of species, 
    conditioned on the fact that at time t, the number of ancestral species 
    in boreal region is q_t. 
    See eqn (17) in mossEqnModel.pdf
                
    Input parameters
    ----------------  
    @param current_state_for_uncond_probs    
                     [0] = q_t and [1] = r_t, see UncondProbSTT for details
    @param num_sum   r_t+num_sum is the upper limit in summation series, this is NOT needed. 
                     I created this just so the method parameters are identical
    @param sigma     set of all parameters for the model
                     sigma = (lambda, b, B, alpha, mu)
                     See UncondProbSTT for more details
                     
    Details
    --------
    Only 2 components of sigma are used: 
    lam      sigma[0] = rate of turnover in boreal region
    B        sigma[2] = the total number of species in boreal community (or state 0); an assumption

    Return value
    -------------
    p(S_BB | q_t) = 2*lam*C(q_t,2)/B^2
    """
    
    q_t=current_state_for_uncond_probs[0]
    lam=sigma[0]
    B=sigma[2]
    return 2*lam*comb(q_t,2)/(B*B)

def UncondProbSBTGlobal(current_state_for_uncond_probs, sigma, num_sum=20):
    """
    Calculates p(S_{BT} | n_t,q_t,r_t), the probability of migration 
    of species from the boreal to the neotropics region when the duplicate 
    in the boreal region also appears in the sample history, aka a 
    "pseudo-speciation" event in the history of the sample, conditioned 
    on n_t, q_t, r_t. 
    See eqn (22) in mossEqnModel.pdf
    
    Input Parameters
    ----------------
    @param current_state_for_uncond_probs    
                     [0] = q_t and [1] = r_t, see UncondProbSTT for details
                     Both q_t and r_t are needed for this method. 
                     Recall:
                     q_t = number of ancestral species in boreal region in the sample
                           history at time t
                     r_t = number of ancestral species in the neotropics
                           region in the sample history at time t  
    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum
    @param sigma     set of all parameters for the model
                     sigma = (lambda, b, B, alpha, mu)
                     See uncondProbSTT for more details
                
    Details
    --------
    Only 4 components of sigma are needed (for instance, some are used only to call CalculatePi; they're marked with *. Those
    marked with ** are used for both CalculatePi and UncondProbSBT): 
    B        sigma[2] = the total number of species in boreal community (or state 0); an assumption
    alpha    *sigma[3] = per species birth rate of in neotropical region when the number of species is i; positive constant
    b        **sigma[1] = effective migration rate (used in both methods)
    mu       *sigma[4] = per species death rate in neotropical region when number of species is i; positive constant
     
    Return Value 
    -----------_
    Returns p{S_{BT} | n_t,q_t,r_t} = b*r_t*q_t/B * \sum_{k>=r_t}(pi_{k-1}/k)
    """    
    q_t=current_state_for_uncond_probs[0] 
    r_t = current_state_for_uncond_probs[1]
    B=sigma[2]
    b=sigma[1]
    
    coef = b*r_t*q_t/float(B)
    sum=0
    for k in range(r_t,r_t+num_sum):
        if (k!=0):
            #if k==0, then we don't want to divide by 0.
            sum=sum+CalculatePi(k-1,sigma)/k
            
    return coef*sum

def UncondProbSBarrowTGlobal(current_state_for_uncond_probs, sigma, num_sum=20):    
    """
    Calculates p(S_{B->T}), the probability of migration of species 
    from the boreal to the neotropics region when the duplicate in the
    boreal region does NOT appear in the sample.
    See eqn (21) in mossEqModel.pdf
    
    Input Parameters
    ----------------
    @param current_state_for_uncond_probs    
                     [0] = q_t and [1] = r_t, see UncondProbSTT for details
                     Both q_t and r_t are needed for this method. 
                     Recall:
                     q_t = number of ancestral species in boreal region in the sample
                           history at time t
                     r_t = number of ancestral species in the neotropics
                           region in the sample history at time t  
    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum
    @param sigma     set of all parameters for the model
                     sigma = (lambda, b, B, alpha, mu)
                     See UncondProbSTT for more details
                
    Details
    --------
    Only 4 components of sigma are needed (for instance, some are used only to call CalculatePi; they're marked with *): 
    B        sigma[2] = the total number of species in boreal community (or state 0); an assumption
    alpha    *sigma[3] = per species birth rate of in neotropical region when the number of species is i; positive constant
    b        **sigma[1] = effective migration rate (used in both methods)
    mu       *sigma[4] = per species death rate of in neotropical region when number of species is i; positive constant
     
    Return Value 
    -----------_
    Returns p(S_{B->T}) = b*r_t*(1-q_t/B)*\sum_{k>=r_t}(pi_{k-1}/k) 
    """
    q_t = current_state_for_uncond_probs[0]
    r_t = current_state_for_uncond_probs[1]
    B = sigma[2]
    b = sigma[1]
    coef = b*r_t*(1-float(q_t/B))
    sum=0
    for k in range(r_t,r_t+num_sum):
        if (k!=0):
            #In case k==0, we don't want divide by 0
            sum=sum+CalculatePi(k-1,sigma)/k
    return coef*sum 

def UncondProbSTT(current_state_for_uncond_probs, sigma, num_sum=20):
    """
    s_tt that is exclusive of kappa!
    """
    tt_global = UncondProbSTTGlobal(current_state_for_uncond_probs, sigma, num_sum=20)
    kappa = UncondProbKappa("s_tt", current_state_for_uncond_probs, sigma, num_sum)
    return tt_global - kappa

def UncondProbSBB(current_state_for_uncond_probs, sigma, num_sum=20):  
    """
    s_bb that is exclusive of kappa!
    """
    bb_global = UncondProbSBBGlobal(current_state_for_uncond_probs, sigma, num_sum=20)
    kappa = UncondProbKappa("s_bb", current_state_for_uncond_probs, sigma, num_sum)
    return bb_global - kappa

def UncondProbSBT(current_state_for_uncond_probs, sigma, num_sum=20):
    """
    s_bt that is exclusive of kappa!
    """
    bt_global = UncondProbSBTGlobal(current_state_for_uncond_probs, sigma, num_sum=20)
    kappa = UncondProbKappa(S_BT, current_state_for_uncond_probs, sigma, num_sum)
    return bt_global - kappa
        
def UncondProbSBarrowT(current_state_for_uncond_probs, sigma, num_sum=20):
    """
    Calculates p(S_{B->T}), the probability one of the not-next-coalescing-lineage's
    migration from the boreal to the neotropics region. The lineage's duplicate in the
    boreal region does NOT appear in the sample.
    See eqn (21) in mossEqModel.pdf
    
    Input Parameters
    ----------------
    @param current_state_for_uncond_probs    
                     [0] = q_t and [1] = r_t, see UncondProbSTT for details
                     Both q_t and r_t are needed for this method. 
                     Recall:
                     q_t = number of ancestral species in boreal region in the sample
                           history at time t
                     r_t = number of ancestral species in the neotropics
                           region in the sample history at time t  
    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum
    @param sigma     set of all parameters for the model
                     sigma = (lambda, b, B, alpha, mu)
                     See uncondProbSTT for more details
                
    Details
    --------
    Only 4 components of sigma are needed (for instance, some are used only to call CalculatePi; they're marked with *): 
    B        sigma[2] = the total number of species in boreal community (or state 0); an assumption
    alpha    *sigma[3] = per species birth rate of in neotropical region when the number of species is i; positive constant
    b        **sigma[1] = effective migration rate (used in both methods)
    mu       *sigma[4] = per species death rate in neotropical region when number of species is i; positive constant

    The only possible direction of migration when going backwards in time is 1->0

    p(S_{B->T}) is multiplied by (r_t - 2) / r_t or r_t - 1) / r_t to reflect the ratio of lineages that could
    undergo migration from 1->0 that are not the next coalescing lineages to all lineages that could undergo migration
    1->0
     
    Return Value 
    -----------_
    Returns p(S_{B->T})*various coefficients
    """


    s_b_arrow_t = UncondProbSBarrowTGlobal(current_state_for_uncond_probs, sigma, num_sum=20)
    x_1=current_state_for_uncond_probs[2]
    x_2=current_state_for_uncond_probs[3]
    r_t=current_state_for_uncond_probs[1]
    if(r_t==0):
        return s_b_arrow_t
    if(x_1==1 and x_2==1):
        return s_b_arrow_t * (r_t - 2) / r_t
    if(x_1==0 and x_2==0):
        return s_b_arrow_t
    if((x_1==1 and x_2==0) or (x_1==0 and x_2==1)):
        return s_b_arrow_t * (r_t - 1) / r_t
    
def UncondProbM1(current_state_for_uncond_probs, sigma, num_sum=20):
    """
    Calculates p(S_{B->T}), the probability one of the next-coalescing-lineage's
    migration from the boreal to the neotropics region. The lineage's duplicate in the
    boreal region does NOT appear in the sample.
    See eqn (21) in mossEqModel.pdf
    
    Input Parameters
    ----------------
   
    @param current_state_of_cond_jump_chain    Really just a 4-tuple (Not a list of tuples - PZ)  
    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum
    @param sigma     set of all parameters for the model
                     sigma = (lambda, b, B, alpha, mu)
                     See uncondProbSTT for more details
                
    Details
    --------
    Only 4 components of sigma are needed (for instance, some are used only to call CalculatePi; they're marked with *): 
    B        sigma[2] = the total number of species in boreal community (or state 0); an assumption
    alpha    *sigma[3] = per species birth rate of in neotropical region when the number of species is i; positive constant
    b        **sigma[1] = effective migration rate (used in both methods)
    mu       *sigma[4] = per species death rate in neotropical region when number of species is i; positive constant

    The only possible direction of migration when going backwards in time is 1->0

    p(S_{B->T}) is multiplied by 1/r_t to reflect the ratio of next coalescing lineages marked x_1 in the tuple
    able to undergo migration to the total number of lineages in state 1
     
    Return Value 
    -----------_
    Returns p(S_{B->T})*various coefficients
    """
    s_b_arrow_t = UncondProbSBarrowTGlobal(current_state_for_uncond_probs, sigma, num_sum=20)
    x_1= current_state_for_uncond_probs[2]
    x_2= current_state_for_uncond_probs[3]
    r_t= current_state_for_uncond_probs[1]
    if (r_t==0):
        return 0
    if(x_1==1):
        return s_b_arrow_t / r_t
    return 0

def UncondProbM2(current_state_for_uncond_probs, sigma,num_sum=20):
    """
    Calculates p(S_{B->T}), the probability one of the next-coalescing-lineage's
    migration from the boreal to the neotropics region. The lineage's duplicate in the
    boreal region does NOT appear in the sample.
    See eqn (21) in mossEqModel.pdf
    
    Input Parameters
    ----------------
    @param current_state_for_uncond_probs    Just a 4-tuple (Not a list of tuples - PZ) 
    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum
    @param sigma     set of all parameters for the model
                     sigma = (lambda, b, B, alpha, mu)
                     See uncondProbSTT for more details
                
    Details
    --------
    Only 4 components of sigma are needed (for instance, some are used only to call CalculatePi; they're marked with *): 
    B        sigma[2] = the total number of species in boreal community (or state 0); an assumption
    alpha    *sigma[3] = per species birth rate of in neotropical region when the number of species is i; positive constant
    b        **sigma[1] = effective migration rate (used in both methods)
    mu       *sigma[4] = per species death rate in neotropical region when number of species is i; positive constant

    The only possible direction of migration when going backwards in time is 1->0

    p(S_{B->T}) is multiplied by 1/r_t to reflect the ratio of next coalescing lineages marked x_2 in the tuple
    able to undergo migration to the total number of lineages in state 1
     
    Return Value 
    -----------_
    Returns p(S_{B->T})*various coefficients
    """
    s_b_arrow_t = UncondProbSBarrowTGlobal(current_state_for_uncond_probs,sigma, num_sum=20)    
    x_1= current_state_for_uncond_probs[2]
    x_2= current_state_for_uncond_probs[3]
    r_t= current_state_for_uncond_probs[1]
    if(r_t==0):
        return 0
    if(x_2==1):
        return s_b_arrow_t / r_t
    return 0

def ReinterpretUncondEvent(event, z):
    """
    To calculate unconditional probabilities (not based on
    each coealescing lineages' present states), M_1 and M_2 are 
    treated as "s_b_arrow_t", and "kappa" is treated as 
    either "s_bb", "s_tt", or S_BT.
    
    This method "reinterprets" the event for unconditional probabilities. 
        
    Input parameters
    ----------------
    event        An event operator, UnconditionalTransitionProbability and
                 MaketransitionMatrixForLevel for more details. 
                 An event could be either one of the following:
                 M_1, M_2, "s_b_arrow_t", KAPPA, "s_bb", "s_tt",
                 or S_BT
                 
    z            a 4-tuple (q_t, r_t, x_1, x_2)
    
    Details
    -------
    In unconditional probability calculations, 
    
    event 'kappa' is treated as (a) 's_bb' if the state of the
                                     next-coalescing lineages is (0, 0)
                                (b) 's_tt' if the state of the two
                                     next-coalescing lineages is (1, 1).
                                (c) S_BT is the state of the
                                     next-coalescing lineages is (1, 0) or
                                     (0, 1)
    For instance, if event kappa operates on z=(5,2,1,0), then
    this must be a S_BT speciation event. 
                                     
    Return value
    ------------
    A string (which is a constant in this program)
    as "s_bb", "s_tt", S_BT, or "s_b_arrow_t", M_1,or M_2
    """
    if (event=="kappa"):
        (x_1,x_2) = z[2:4]
        map = {(0,0):"s_bb", (1,1):"s_tt", (1,0):S_BT, (0,1):S_BT}
        return map[(x_1,x_2)]
    
    return event;

def UncondProbKappa(event, current_state_for_uncond_probs, sigma, num_sum=20):
    """
    Calculates unconditional probability of Kappa, which is
    a very, very small subset of s_tt, s_bb, or s_bt. 
    """
    q_t = current_state_for_uncond_probs[0]
    r_t = current_state_for_uncond_probs[1]
    x_1 = current_state_for_uncond_probs[2]
    x_2 = current_state_for_uncond_probs[3]
    if (event==S_BT):
        if((x_1==0 and x_2==1) or (x_1==1 and x_2==0)):
            temp = UncondProbSBTGlobal(current_state_for_uncond_probs, sigma, num_sum)
            if (q_t!=0 and r_t!=0):
                return temp / (q_t*r_t)
    if (event=="s_bb"):
        if(x_1==0 and x_2==0):
            if(comb(q_t,2)!=0):
                temp = UncondProbSBBGlobal(current_state_for_uncond_probs, sigma, num_sum)
                return temp / comb(q_t, 2)
    if (event=="s_tt"):
        if(x_1==1 and x_2==1):
            if(comb(r_t,2)!=0):
                temp = UncondProbSTTGlobal(current_state_for_uncond_probs, sigma, num_sum)
                return temp / comb(r_t,2)
    return 0

def UnconditionalTransitionProbability(event, current_state_of_cond_jump_chain, sigma):
    """
    Calculates the unconditioned probability (i.e., not conditioned on the
    tree) of an event at the current state, based on Equation 23 (page 18).

    Input parameters
    ----------------
    
    current_state_of_cond_jump_chain    is a 4-tuple (n_lineages_in_state_zero, n_lineages_in_state_one, state_of_next_coalescing_lineages_1, state_of_next_coalescing_lineages_2)
                                        
                                        equivalent to (q, r, x_1, x_2) in
                                        paragraph 3, page 19.

                                        The next-coalescing
                                        lineages are the two that are involved in the
                                        next (going back in time) speciation event,
                                        based on the input phylogenetic tree. See
                                        paragraph entitled 'State spaces of conditional
                                        jump chains' on page 19.

                                        BUT NOTE THAT ONLY THE FIRST TWO
                                        ELEMENTS, n_lineages_in_state_zero
                                        and n_lineages_in_state_one, ARE
                                        RELEVANT TO THIS FUNCTION, SINCE
                                        ONLY UNCONDITIONAL PROBABILITIES
                                        ARE BEING CALCULATED HERE.


    event       is one of the events in bigbig_H (see function GetEvents,
                and also segment entitled "Events as operators", page 19.

                basically, an event could be either:

                'kappa'     denoting a speciation involving the two
                            next-coalescing lineages. 
                
                's_bb'      speciation events not involving both the
                S_BT      next-coalescing lineages (see Section 3, bullet
                's_tt'      point entitled "Events", and also Figure 1 on page 6

                M_1       migration events involving the first and the
                M_2       second of the next-coalescing lineages
                            respectively.
                
                's_b_arrow_t'   any other migration event. The model in
                                mossEqModel allows only migrations of a
                                lineage from character state 0 to character

    @param sigma    sigma =(lambda, b, B, alpha, mu), where  
                    lambda is the turnover rate in state 0 (the boreal region), 
                    b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                    mu = the per lineage extinction rate in state 1 (the neotropical region).


    
    Return value
    ------------
    
    The probability of the given event happening when the chain is in state
    current_state_of_cond_jump_chain, computed based on Equation 23 on page
    18.

    But note that the probability calculated is Pr(event | n_lineages_in_state_zero, n_lineages_in_state_one)
    i.e., the probability is not conditioned on the state of the next-coalescing lineages.


    Details
    -------

    Since the function calculates unconditioned probabilities, 
    
    event 'kappa' is treated as (a) 's_tt' if the state of the two
                                     next-coalescing lineages is (1, 1).
                                (b) 's_bb' if the state of the
                                     next-coalescing lineages is (0, 0)
                                (c) S_BT is the state of the
                                     next-coalescing lineages is (1, 0) or
                                     (0, 1)
    
    events M_1, M_2, and 's_b_arrow_t' are not treated as disjoint
    """

    # current state that's relevent for the computation of unconditional
    # probabilities. Assuming n_lineages_in_state_zero and
    # n_lineages_in_state_one are the first two elements in
    # current_state_of_cond_jump_chain
    ##I fixed a typo
    event = ReinterpretUncondEvent(event, current_state_of_cond_jump_chain)
 
    # In python the set() statement makes a set from a list.
    set_of_events = set(['s_bb', 's_tt', S_BT, M_1,M_2,SMALLER_S_B_ARROW_T])

    # UncondProbSTT implements Equation 13.
    # UncondProbSBB implements Equation 17.
    # UncondProbSBarrowT implements Equation 21.
    # UncondProbSBT implements Equation 22.
    #I create two maps, the numerator map maps event to its "instantaneous",
    #rates, such as KAPPA, M_1, and SMALLER_S_B_ARROW_T. The denominator
    #map maps s_tt to s_tt, s_bb to s_bb, etc. But
    #since M_1+M_2+SMALLER_S_B_ARROW_T = s_b_arrow_t_global,  
    #and also m_1 and m_2 exist in set_of_events, I keep them
    #in the denominator map. 
    numerator_event_to_function = {'s_tt': UncondProbKappa, 's_bb': UncondProbKappa, SMALLER_S_B_ARROW_T: UncondProbSBarrowT, S_BT: UncondProbKappa, M_1:UncondProbM1, M_2:UncondProbM2}
    
    denominator_event_to_function = {'s_tt': UncondProbSTTGlobal, 's_bb': UncondProbSBBGlobal, SMALLER_S_B_ARROW_T: UncondProbSBarrowT, S_BT: UncondProbSBTGlobal, M_1:UncondProbM1, M_2:UncondProbM2}
    # e loops over events in set(['s_tt', 's_bb',  S_BT, M_1, M_2,SMALLER_S_B_ARROW_T])
    # we will also need the sum of all rates, for calculating the
    # probabilities from rates (basically, the sum is the denominator is
    # Eq. 23).
    numerator = 0.0
    sum_of_rates = 0.0
    for e in set_of_events:
        # Call UncondProbSTT, UncondProbSBB, UncondProbSBarrowT or
        # UncondProbSBT depending on what the event is.
        #We only calculate numerator once, when event==e
        #We then divide the methods in the map into if/else based
        #on their two different parameters
        if (event==e):
            if(e==S_BT or e=="s_bb" or e=="s_tt"):
                numerator = numerator_event_to_function[e](event, current_state_of_cond_jump_chain, sigma, num_sum=20)
            else:
                #Must be M_1, M_2, or s_b_arrow_t
                numerator = numerator_event_to_function[e](current_state_of_cond_jump_chain, sigma, num_sum=20)
        #All methods in the denominator map have the same parameters
        sum_of_rates += denominator_event_to_function[e](current_state_of_cond_jump_chain, sigma, num_sum=20)
    #print "sum of rates is: " + str(sum_of_rates)    
    # This is basically Eqn. 23, with the numerator and denominator.
    probability_of_event = numerator/sum_of_rates
    return(probability_of_event)

def MakeTransitionMatrixForLevel(level_number, sigma):
    """
    Makes and returns the transition probability matrix for a given level 

    Section 5 of mossEqModel.pdf is devoted to this.

    Input parameters
    ----------------
    
    1. level number    The level for which the transition matrix is
                       computed.
    


    2. sigma  is the set of model all parameters. 

              for the model in mossEqModel.pdf, sigma = (lambda, b, B, alpha, mu), where lambda is the
              turnover rate, b is the migration rate, B is the total number of
              species in state 0 (the boreal region, according to
              mossEqModel.pdf), and alpha and mu are the speciation and
              extinction rates, respectively, in state 1 (i.e., the
              neotropical region).
    

    Return value(s)
    ---------------

    a 3-tuple: (transition_matrix, state_to_index_in_transition_matrix, index_in_transition_matrix_to_state)
    
    transition_matrix                       4(k-1)+1 x 4(k-1)+1 matrix of floats, where 
                                            k = level number

    state_to_index_in_transition_matrix     is a map from the states of the
                                            jump chain (but only those that are relevant 
                                            at the given le
                                            vel) to
                                            integers, i.e., basically a numbering of the states. 
                                            This is so that the transition matrix can be 
                                            indexed with integers.

    index_in_transition_matrix_to_state     is the inverse map of  
                                            state_to_index_in_transition_matrix
    


    Details
    -------

    When the conditional jump chain is in level k (i.e., the total
    number of lineages = k), the number of states at that level is 4(k-1) (see
    paragraph 3, page 19). Besides transitioning within level k, i.e.,
    among its 4(k-1) states, the chain can also effect precisely one
    transition to a state in level k-1, representing the speciation involving the two
    next_coalescing_lineages whose character states (whether 0/1) are
    represented in state_of_cond_jump_chain. We add this state as the
    4(k-1)+1-th state in the transition matrix.
    
    The precise identity of this state in level k-1 does not matter. 
    All that we need is to be able calculate the probability of the 
    above-mentioned speciation event  from each of the 4(k-1) states in level k.

    """

    # k is just an alias for level_number
    k = level_number

    # first create a map between the states and their indices in the
    # transition matrix. In python, this can be implemented as a
    # dictionary. All the states in a level can be figured out from the
    # level number.
    state_to_index_in_transition_matrix = MakeState2IndexDictionary(level_number)
    # also make the reverse map.
    ##I changed input of the next function from state_to_index_in_transition_matrix to level_number
    index_in_transition_matrix_to_state = ReverseMap(level_number)

    # The set of possible types of events.
    # bigbig_K is the set of various types of speciation events.
    # bigbig_M is the set of various types of migration events.
    # bigbig_H is the union of the above two sets.
    # See function GetEvents() for more details.
    # Also see Page 19, segment entitled "Events as operators on states". 
    (bigbig_K, bigbig_M, bigbig_H) = GetEvents()

    # ========================================================================
    # Begin Task: Solving the system of equations represented by Equation 25 #
    # ========================================================================

    # Equation 25 is a system of 4(k+1) equations in 4(k+1) unknowns. 
    # The unknowns are Pr(F|z) for each state z in the level-k state space,
    # where F is the event that no speciation events occur prior to the
    # one involving the two known next-coalescing lineages. 
    #
    # 
    # For ease of description, let the states in Z_k be states z_0
    # through z_{4(k+1)-1}. Now, if we look at Equation 25, each unknown z_i
    # occurs with a coefficient of 1 in exactly one equation. In all
    # other equations, its coefficient is either 0 or negative. We will
    # call the one equation where the coefficient of z_i is positive
    # "the equation corresponding to z_i" or simply "equation z_i." In
    # Equation 25, what comes after "for all z in Z_k" is the equation
    # corresponding to z.
    #               _0 through z_{4(k+1)-1} and for each
    # state we generate its corresponding equation. We will then
    # collect all the equations and stick them in an equations-solver
    # and get out the answers.
    
    # first get the entire state space of the conditional jump chain at
    # level k. The state space depends only on k, the total number of lineages, and
    # can be calculated from state_to_index_in_transition_matrix
    #
    # Note: state_space_at_the_current_level must be such that:
    #   1. length(state_space_at_the_current_level) = 4(k+1)
    #   2. if state_space_at_the_current_level[j] = z, then state_to_index_in_transition_matrix[z] = j
    #      Basically, the states should be listed in
    #      state_space_at_the_current_level according to their indices in
    #      state_to_index_in_transition_matrix. 
    #   The asserts following this statement verifies this.
    #   If the condition inside the assert statement fails, the program
    #   will exit at this point.
    state_space_at_the_current_level = GetCondJumpChainStateSpace(level_number)
    assert(len(state_space_at_the_current_level) == 4*(k-1))
    # In python, range(0, n) = [0, 1, 2, ..., n-1]
    for j in range(0, 4*(k-1)):
        assert(state_to_index_in_transition_matrix[state_space_at_the_current_level[j]] == j)
    
    # I am using the python map syntax. 
    # For example, if I say squares = [j**2 for j in [0, 1, 2, 3]], then
    # squares = [0, 1, 4, 9]. Similarly, system_of_equations will be a
    # list of equations one for each unknown z_i. 
    # GetLinearEquation(z, state_of_cond_jump_chain)
    # could (I am saying 'could' because this is a suggested representation, and the actual
    # implementation can change - Ganesh) return an
    # equation in the form of a vector of length 4(k+1)+1, where the
    # i-th element of the vector is the coefficient of the unknown z_i,
    # and the last element is the constant in the equation. 
    matrix,b = GetLinearEquations(state_space_at_the_current_level, sigma, level_number)
    solution = LinearEquationSolver(matrix, b)
    # ======================================================================
    # End Task: Solving the system of equations represented by Equation 25 #
    # ======================================================================


    # =================================================================
    # Begin Task: Calculating the transition matrix using Equation 24 #
    # =================================================================

    # make a 2-dimensional array of dimensions 4(k+1)+1 x 4(k+1)+1
    # Make2DimensionalArray should initialize the array with 0.0
    transition_matrix = Make2DimensionalArray(4*(k-1)+1)
        
    # This is to make sure that the array is initialized properly.
    # In python, range(0, n) = [0, 1, 2, ..., n-1]
    for j in range(0, 4*(k-1)+1):
        # assert that each row sums to 0.0. If the assert fails, the
        # program will exit at this point
        assert(sum(transition_matrix[j]) == 0.0)

    # The following loop fills out the transition matrix state by state. 
    # Each state corresponds to a row in the matrix, and the value of each index in the
    # row corresponds to the transition probability of the row's state to the
    # column's state.
    # Thus, each iteration of the following for loop handles one row of the
    # transition matrix.
    # Remember that state_space_at_the_current_level contains the 4(k-1) states at
    # the level k. It does not include the one state at level k-1 the chain
    # might transition to. But we need not fill the row corresponding to
    # that state (the one in level k-1) because no transitions are possible *from* that state,
    # and appropriately the probabilities in that row are already initialized to 0.0
    # 
    # Variables z and w have the same meaning that they have in Equation
    # 24. z loops over all states in Z_k, and w loops over all events in
    # bigbig_H
    for z in state_space_at_the_current_level:
        row = state_to_index_in_transition_matrix[z]
        # from z, the chain can only transition to those states x such that
        # x = w(z) for some event w. So we don't have to visit all the
        # columns and fill in transition_matrix[row][column]; we visit only
        # those states that are one event away from z. The other transition
        # probabilities are already correctly set to 0.0. 
        for w in bigbig_H:
            # ApplyEvent will take the current state (z, in Equation 24) and an
            # event (w, in Eq. 24) and return the state w(z). That is, the state
            # that the chain transitions when event w happens when it (the chain)
            # is in state z.
            w_of_z = ApplyEvent(z, w,level_number)
            #If column == (-1,-1,-1,-1), which means the event
            #cannot be applied to z, then P(w|z)=0, so we don't need to 
            #fill out this entry anyhow. 
            if (w_of_z != (-1,-1,-1,-1)):
                #Now, if column != (,,-1,-1),that is, 
                #when the event is "kappa", we place it in the last column
                if(w_of_z[2]==-1 and w_of_z[3]==-1):
                    column = 4*(k-1)
                else: 
                    column = state_to_index_in_transition_matrix[w_of_z]
   
            # In the numerator of the right side of Eq. 24, 
            #   1. Pr(F | w(z)) is in terms of the program's variables,
            #      solution[column], since
            #           a. solution[j] = Pr(F | z) where z = state_space_at_the_current_level[j] 
            #           b. state_to_index_in_transition_matrix[w_of_z] = column, 
            #              and this implies w_of_z = state_space_at_the_current_level[column]
            #           
            #   2. Pr(w | z) is to be calculated using the function
            #      UnconditionalTransitionProbability(z, w, sigma) which will implement
            #      Equation 23 (it also has to implement all equations in
            #      Section 5.1 in order to do so)
            #
            #   The denominator of Eq. 24 is a normalization constant. We
            #   can normalize the whole row after filling the entries out.
            #   Normalization means: suppose we have a vector [1, 2, 3, 4],
            #   after normalization it will become [1/(1+2+3+4), 2/10, 3/10, 4/10]
            #   After normalization the elements in the vector should add
            #   up to 1.
            #   transition_matrix[row][column] = numerator in eqn(24)
            #    If event is kappa, place it in the last column
                if (column == 4*(k-1)):
                    transition_matrix[row,column] = UnconditionalTransitionProbability("kappa", z, sigma)
                else: 
                    transition_matrix[row,column] = solution[column] * UnconditionalTransitionProbability(w, z, sigma)
        #Testing if sum of rows = P(F|z), really normalizing against P(F|z)
        testing = transition_matrix[row].copy()
        testing2 = transition_matrix[row].copy()
        #Equal is defined in this class to test two numbers are so close
        #they're considered "equal" - use this is b/c == is hard
        #to apply to float numbers. 
        #assert(Equal(testing.sum(), solution[row]))
        for j in range(len(testing)):
            testing[j] = testing[j]/solution[row]
        Normalize(testing2)  
#        assert(allclose(testing,testing2))
        
        Normalize(transition_matrix[row])
    # =================================================================
    # End Task: Calculating the transition matrix using Equation 24 #
    # =================================================================
    #Note that last row of transition_matrix is always 0, b/c we don't
    #even touch it! 
    return((transition_matrix, state_to_index_in_transition_matrix, index_in_transition_matrix_to_state))

def Normalize(row):
    """
    Input parameter
    ----------------
    row        an array of class ndarray, such as ndarray(2, 3, 5, 6)
    
    Details
    -------
    If the input is (2,3,5,6), return (2/(2+3+5+6), 3/(16), 5/16, 6/16)
    
    Output value
    ------------
    Return a normalized version of the array. 
    """
    #Scipy function calculating the sum of the entries of the array
    sum = row.sum()  
    if sum==0:
        #in case it's a row of zeros, don't want to divide by 0.
        return row
    for i in range(len(row)):
        row[i] = float(row[i])/sum
    return row

def MakeTransitionMatricesbyLevels(G, delta, sigma):
    """
    Makes transition matrices, one for each level for the input tree, for
    the conditional jump chain described in Section 5.3.

    Input parameters
    ----------------
    
    G, delta, sigma         see function LikelihoodOfParameters

    Return value
    ------------

    A 3-tuple:
    (transition_matrices, state_to_index_in_transition_matrices, index_in_transition_matrices_to_state)
   
    transition_matrices[k] should contain the transition matrix for
    level n-k, where n is the number of leaves in G, and 0 <= k <= n-1
    (A level is the portion of G between two
    speciation events.) 

    state_to_index_in_transition_matrices[k] is a dictionary for level
    n-k, which maps each state (of the jump chain) in level k to its index
    in that level's transition matrix.
   
    index_in_transition_matrices_to_state[k] is the reverse map of
    state_to_index_in_transition_matrices[k]
    """

    # level numbers range from k, ...., 2, 1
    n_leaves = G.num_leaves
    
    transition_matrices = []
    state_to_index_in_transition_matrices = []
    index_in_transition_matrices_to_state = []

    # range(n_leaves, 1, -1) is [n_leaves, n_leaves-1, ..., 2]
    for j in range(n_leaves, 1, -1):

        (x, y, z) = MakeTransitionMatrixForLevel(j, sigma)
        transition_matrices.append(x)
        state_to_index_in_transition_matrices.append(y)
        index_in_transition_matrices_to_state.append(z)

    return((transition_matrices, state_to_index_in_transition_matrices, index_in_transition_matrices_to_state))        
        
def GetInitialStatesofLineages(delta):
    """
    returns (#lineages in state 0, #lineages in state 1) based on what's in
    delta.

    Input parameters
    ----------------

    delta  is a map from lineages to 0 or 1. That is, delta
           assigns the state 0 or 1 to lineages. The state could be
           thought of as a phenotypic character or a geographic location.
           In the manuscript mossEqModel.pdf, state 0 is assumed to be the
           boreal region and state 1 the neotropical region.

    Return value(s)
    ---------------
    
    a tuple: (#lineages in state 0, #lineages in state 1). 
    """
    
    # function to be implemented. Basically go through delta and calculate
    # the number of lineages in state_0 and the number of lineages in state
    # 1.
    n_lineages_in_state_zero = 0
    n_lineages_in_state_one = 0
    for i in delta:
        if delta[i] == 0:
            n_lineages_in_state_zero = n_lineages_in_state_zero + 1 
        if delta[i] == 1:
            n_lineages_in_state_one = n_lineages_in_state_one + 1        
    return(n_lineages_in_state_zero, n_lineages_in_state_one)


def GetInitialStateofCondJumpChain(G, delta):
    """
    return the initial state of the conditional jump chain, by looking at
    delta and the topology of tree G; delta assigns character states (0/1) to
    the tips of G.

    Input parameters
    ----------------
    G       is a rooted phylogenetic tree with known time-order of
            branching events. G is an instance of class Tree defined in
            is_classes.py
    
    delta  is a map from the tips of G to 0 or 1. That is, delta
           assigns the state 0 or 1 to each tip of G. The state could be
           thought of as a phenotypic character or a geographic location.
           In the manuscript mossEqModel.pdf, state 0 is assumed to be the
           boreal region and state 1 the neotropical region.

    Return value(s)
    ---------------

    a 4-tuple: (n_lineages_in_state_zero, n_lineages_in_state_one, state_of_next_coalescing_lineages_1, state_of_next_coalescing_lineages_2)
               For the meaning of the variables in the tuple, see 'Details'
               below.

    Details
    -------

    There are 4 state variables, described below.
    They are equivalent, respectively, to q, r, x_1 and x_2 in page
    19, paragraph 3, line 4.
    
    IMPORTANT: A LINEAGE IS A SPECIES IN THE SAMPLE OR AN ANCESTOR OF
               A SPECIES IN THE SAMPLE.
   
    n_lineages_in_state_zero              the number of lineages in state 0
    n_lineages_in_state_one               the number of lineages in state 1
    
                                          |the states of the two lineages 
    state_of_next_coalescing_lineages_1   |involved in the next (going
    state_of_next_coalescing_lineages_2   |backwards in time) speciation
                                          |event.
    
    The unconditional jump chain (see Section 5.2) uses the first two;
    and the conditional jump chain uses all the four state variables.
    For more details, see Section 5.2 and the paragraph titled 
    "State spaces of conditional jump chains", page 19.

    The first two state variables can be inferred from delta.
    The next two state variables can be inferred from G

    """

    # Function to be implemented.
    
    # As mentioned in "Details", get the first two states from delta.
    (n_lineages_in_state_zero, n_lineages_in_state_one) = GetInitialStatesofLineages(delta)
    

    # go through G and get state_of_next_coalescing_lineages_1 and
    # state_of_next_coalescing_lineages_2
    #
    # Basically, G would have all its levels demarcated. G.levels[0]
    # contains the initial level, i.e., the level with all the tips. 
    # G.levels[0].end_node is the parent of the first two coalescing
    # lineages.
    [next_coalescing_lineages_1, next_coalescing_lineages_2] = G.levels[0].end_node.children

    #In this new implementation, next_coalescing_lineage_1 is of 
    #class Node, so to get the index, we need to do next_coalescing_lineage_1.index
    state_of_next_coalescing_lineages_1 = delta[next_coalescing_lineages_1.index]
    state_of_next_coalescing_lineages_2 = delta[next_coalescing_lineages_2.index]

    return((n_lineages_in_state_zero, n_lineages_in_state_one, state_of_next_coalescing_lineages_1, state_of_next_coalescing_lineages_2))
        
def PickNextStateofChain(row_of_transition_matrix):
    """
    Pick the index of next state from a row of transition matrix, along with its
    probability of transition (which is really row_of_transition_matrix[x])
    
    For instance, if the row_of_transition_matrix is [0,0,.33,0,.54,.12]
    then a possible output is (4,.54). Another possible output is (2,.33)
    """
#    Got the idea from http://snippets.dzone.com/posts/show/732
#    Randomly choose next state based on their weights. For instance, 
#    array = [0,0,.33,0,.54,.12] would have weights 0<=x<=.33 corrspond to 
#    array[2], .33<=x<.88 correspond to array[4], and .88<=x<1 corresponds 
#    to array[5]

    #don't have to worry about row_of_transition_matrix=[.25,.25,.5,0], row_for_choosing=[.25,.5,1,1]
    #If n randomly equals 1,
    #and the the third "if statement" is 1<=1<=1, you might think the code falsly chooses the 3rd index, 4th spot in
    #the row_of_transition_matrix, which has probability 0. But the code would have gotten to the second "if statement"
    #first, .5<=1<1, correctly choosing the 2nd index, 3rd spot in row_of_transition_matrix. Or say,
    #row_of_transition_matrix=[.25,.25,0,.5], row_for_choosing=[.25,.5,.5,1] and n=.5. The code would return
    #(1, row_of_transition_matrix[1]) before it got to i=2, and even if it did get to i=2, it is not true that .5<=.5<.5. 

    n = random.uniform(0,1)
   # for i in range(len(row_of_transition_matrix)): 
   #     if (n < row_of_transition_matrix[i]):
   #         return (i, row_of_transition_matrix[i])
   #     n = n - row_of_transition_matrix[i]
   # return n
   
    spot=0
    row_for_choosing=range(len(row_of_transition_matrix))
    for i in range(len(row_of_transition_matrix)):
        spot+=row_of_transition_matrix[i]
        row_for_choosing[i]=spot

    for i in range(len(row_for_choosing)):
        if(i==0):
            if(0<=n<row_for_choosing[i]):
                return (i, row_of_transition_matrix[i])
        if(i==(len(row_for_choosing)-1)):
            if(row_for_choosing[i]<=n<=1):
                return (i, row_of_transition_matrix[i])
        if(row_for_choosing[i-1]<=n<row_for_choosing[i]):
                return (i, row_of_transition_matrix[i])
            
    #Tally: Perry, the code I wrote already works the way you explained below. 
    #I changed this back to my implementation b/c I was
    #getting errors with the other implementation. I am adding this comment
    #in the log. From now on, please put comments in the log.
    
    #Perry: Tally, random(0,1) actually generates a number 
    #from 0<= x < 1. See: http://docs.python.org/lib/module-random.html
    #Furthermore, if row=[.25,.25,0,.5], and n=.5, then we 
    #should return (3,.5) b/c
    #0<=x<.25 go into slot 0, .25<=x<.5 go into slot 1, and .5<=x<1 go into 
    #slot 3. This way it's truly "uniformly" distributed.           

        
def WhetherInTheSameLevel(state1, state2):
    """
    Test if state1 and state 2 are in the same level, where level is 
    defined as the portion between two bifurcation events
    """
    if state1[0]+state1[1] != state2[0]+state2[1]:
        return False
    return True

def MigrationType(z1, z2):
    """
    Identify migration type from z1 to z2. Return whether it's
    M_1, M_2, or SMALLER_S_B_ARROW_T
    
    Precondition: Assumes z1 and z2 are in same level. 
    """
    ret = ""
    if (z1[2]==z2[2] and z1[3]==z2[3]):
        ret = SMALLER_S_B_ARROW_T
    if (z1[2]==z2[2] and z1[3]!=z2[3]):
        ret = M_2
    if (z1[2]!=z2[2] and z1[3]==z2[3]):
        ret = M_1
    return ret

def ChooseLineageandUpdateDelta(G, level_number, migration_type, current_delta, all_delta_earlier):
    """
    Choose a lineage to migrate based on migration_type and return the migrated lineages and updated deltas
    
    Input Values
    ------------
    G               of type Tree
    level_number    the current level
    migration_type  either s_b_arrow_t or m1 or m2
    current_delta   a map of all the existing lineages at the level by index to their current states
    all_delta_earlier       a map of all existing lineages at all levels by index to their current states

    
    Update all_delta based on migration_type, which could 
    be M_1, M_2, SMALLER_S_B_ARROW_T. Also return which lineage
    to migrate. 
    
    Need to know also the next-coalescing lineages. 
    Say level is at 5, and next-coalescing lineages are in species 3 and 4
    So say delta=[1,1,0,1,1] and migration_type is M_2, then 
    this method should change delta to [1,1,0,1,0], and return 4, since 
    4th lineage is the lineage to migrate.

    Return Values
    -------------
     migrating_lineage  the lineage choosen to migrate
     current_delta      updated
     all_delta_earlier  updated

    """
    #This returns the "indices" of the next_colalescing_lineages
    #remember "level_number" is really calculated backwards in G.levels
    #where G[0] is the most recent level, and G[n-1] is the first level
    next_coalescing_lineages = G.levels[G.num_leaves - level_number].end_node.children
    lineage1=-1
    lineage2=-1
    if (migration_type ==M_1):
        lineage1 = next_coalescing_lineages[0].index
        #migration always goes from 1 --> 0 backwards in time
        assert(current_delta[lineage1]==1)
        current_delta[lineage1] = 0
        all_delta_earlier[lineage1] = 0
        return (lineage1,current_delta,all_delta_earlier)
    if (migration_type==M_2):
        lineage2 = next_coalescing_lineages[1].index
        assert(current_delta[lineage2]==1)
        current_delta[lineage2] = 0
        all_delta_earlier[lineage2] = 0
        return (lineage2,current_delta,all_delta_earlier)
    if (migration_type==SMALLER_S_B_ARROW_T):
        #need to randomly pick lineages whose deltas==1 to migrate, EXCLUDING
        #those involving next_coalescing_lineages. 
        #So first delete next_coalescing_lineages from consideration and those lineages
        #whose state is already 0o
        random_list = []
        for i in current_delta.keys():
            if (i != next_coalescing_lineages[0].index and i != next_coalescing_lineages[1].index and current_delta[i]!=0):
                #print "this is next_coalescing lineage" + str(i)
                random_list.append(i)
        #so in our example, newlist becomes [0,1] b/c besides
        #the next_coalescing_lineages 3 and 4, lineages 0 and 1 would
        #fall into SMALLER_S_B_ARROW_T
        #We now "randomly shuffle" it and select the first element, that's
        #where we want to migrate to
        #print "before shuffling: " + str(random_list)
        random.shuffle(random_list)
        #print "after shuffling: " + str(random_list)
        #update current_delta first
        lineage_to_migrate = random_list[0]
        current_delta[lineage_to_migrate] = 0
        all_delta_earlier[lineage_to_migrate] = 0
        return (random_list[0],current_delta,all_delta_earlier)

def MovetoNextLevel(current_state, current_delta, all_delta_earlier, all_delta_later, current_level, next_level):
    """
    (a) update current_delta,
    (b) return the initial state for the next level's conditional jump chain
        recall this has to be a kappa event.
    (c) initialize the event history for the next level with the updated delta
    (d) update all_delta_earlier
    (e) update all_delta_later
    
    Input parameters
    -----------------
    current_delta           a map of lineage indices to character states, these are the character states
                            at the "beginning" i.e. later in time part of the level
    all_delta_earlier       a map of all existing lineages at all levels by index to their current states
    all_delta_later         a map of all existing lineages at all levels by index to their most recent states
    current_level           is the of class Level
    next_level              is also of class Level

    Return Value
    ------------
    initial_state_in_the_next_level     initial state for the next (earlier) level's jumpchain
    current_delta                       updated
    all_delta_earlier                           updated

    Details
    -------
    
    Say next-coalescing lineages are 3 and 4, which coalesce to form lineage 6. 
    Say delta = {1:1, 2:1, 3:0, 4:0, 5:0}, and current_level=5, next_level=4. 
    The updated delta should be {1:1, 2:1, 6:0, 5:0} Note: this is a map so order
    does not matter.  all_delta_earlier would be {1:1, 2:1, 6:0, 5:0, 3:0, 4:0}
    
    Now when kappa happens to (1,4,1,1), it becomes 
    (1,3,-1,-1), the last two elements are -1 b/c we don't know what are
    "its" next-coelescing lineages. But if we have the tree, we would know: 
    it's just G.levels[num_leaves - (level_number - 1)].end_node.children
    """
    q_t = current_state[0]
    r_t = current_state[1]
    x_1 = current_state[2]
    x_2 = current_state[3]
    current_x_1_index = current_level.end_node.children[0].index
    current_x_2_index = current_level.end_node.children[1].index
    if(len(next_level.lineages)==1):
        next_x_1_index = -1
        next_x_2_index = -1
    else:        
        next_x_1_index = next_level.end_node.children[0].index
        next_x_2_index = next_level.end_node.children[1].index
    next__level_speciating_lineage_index=next_level.begin_node.index
    assert(next_level.begin_node.index==current_level.end_node.index)#just a check
    new_index=next__level_speciating_lineage_index
    
    #print current_x_1_index
    #print current_x_2_index
    #print next_x_1_index
    #print next_x_2_index
    #print next__level_speciating_lineage_index
    #print new_index
    
    current_delta.pop(current_x_1_index)
    current_delta.pop(current_x_2_index)
    
    
    if (x_1==0 and x_2==0): #sbb
        current_delta.update({new_index:0})
        all_delta_earlier.update({new_index:0})
        all_delta_later.update({new_index:0})

        assert(q_t>=1)#this should be a given by the structure of the transition matrix
        if (q_t >= 1):
            if(len(next_level.lineages)==1):
                return ((q_t-1, r_t, -1, -1),current_delta,all_delta_earlier,all_delta_later);
            else:
                return ((q_t-1, r_t, current_delta[next_x_1_index], current_delta[next_x_2_index]),current_delta,all_delta_earlier,all_delta_later)
    if (x_1==1 and x_2==1):#stt
        current_delta.update({new_index:1})
        all_delta_earlier.update({new_index:1})
        all_delta_later.update({new_index:1})

        assert(r_t>=1)#this should be a given by the structure of the transition matrix
        if (r_t >= 1):
            if(len(next_level.lineages)==1):
                return ((q_t, r_t-1, -1, -1),current_delta,all_delta_earlier,all_delta_later);
            else:
                return ((q_t, r_t-1, current_delta[next_x_1_index], current_delta[next_x_2_index]),current_delta,all_delta_earlier,all_delta_later)
    if ((x_1==1 and x_2==0) or(x_1==0 and x_2==1)):#sbt
        current_delta.update({new_index:0})
        all_delta_earlier.update({new_index:0})
        all_delta_later.update({new_index:0})

        #(1,0) --> 0 so r_t decreases by 1
        assert(r_t>=1)
        if (r_t >=1):#this should be a given by the structure of the transition matrix
            if(len(next_level.lineages)==1):
                return ((q_t, r_t-1, -1,-1),current_delta,all_delta_earlier,all_delta_later);
            else:
                return ((q_t, r_t-1, current_delta[next_x_1_index], current_delta[next_x_2_index]),current_delta,all_delta_earlier,all_delta_later)

def Test(G):
    sigma=(0.10000000000000001, 0.040000000000000012, 50, 0.090000000000000003, 0.10000000000000001)
    delta={0:0,1:1,2:1,3:1}
    alpha=sigma[3]
    mu=sigma[4]
    if(alpha>mu):
        return('error, alpha is greater than mu')
    (transition_matrices, state_to_index_in_transition_matrices, index_in_transition_matrices_to_state)=MakeTransitionMatricesbyLevels(G,delta,sigma)
    
    return SampleFromIS(G,delta,sigma,(transition_matrices, state_to_index_in_transition_matrices, index_in_transition_matrices_to_state))

def SampleFromIS(G, delta, sigma, (transition_matrices, state_to_index_in_transition_matrices, index_in_transition_matrices_to_state)):
    """
    IMPORTANT: DO NOT ALLOW MIGRATION IN ROOT LINEAGE
    This function assigns migrations (i.e., character state changes) to the
    branches of G, with migrations being drawn from a probability
    distribition such that the density of one assignment of migrations M is
    Pr(M | G, delta, sigma).

    migrations will be assigned to branches with known time-order, but 
    without specifying the actual times of migrations. "all_delta_earlier" will be created.
    "all_delta_earlier" is dictionary mapping the indices of all lineages in all levels to their character
    states right before the lineages were created in a branching event. Also, "all_delta_later" will be created.
    "all_delta_later" is dictionary mapping the indices of all lineages in all levels to their most recent character
    states.
    

    Input parameters
    ----------------

    G, delta, sigma         see function LikelihoodOfParameters


    The fourth parameter is a 3-tuple. Each of its elements is described below.

    transition_matrices     The migrations are assigned by running a
                            discrete-time Markov chain (the conditional
                            jump chain) once. 
                            
                            transition_matrices[k] should contain the transition matrix for
                            level n-k, where n is the number of leaves in G, and 0 <= k <= n-1
                            (A level is the portion of G between two
                            speciation events.) 

    state_to_index_in_transition_matrices           is a list of
                                                    dictionaries, one for each level.
                                                    
                                                    state_to_index_in_transition_matrices[k] 
                                                    is a dictionary for level
                                                    n-k, which maps each state (of the jump chain) 
                                                    in level k to its index
                                                    in that level's transition matrix.
   
    index_in_transition_matrices_to_state           is a list of
                                                    dictionaries, one for each level.
                                                    
                                                    is the reverse map of
                                                    state_to_index_in_transition_matrices[k]
    
    Return value(s)
    ---------------
    
    a real number density_A, where A is a a time-order event history A (for use in the right side of
                        Equation 8)
    all_delta_earlier   a dictionary mapping the indices of all lineages in all levels to their character
                        states right before the lineages coalesce

    all_delta_later     a dictionary mapping the indices of all lineages in all levels to their most recent
                        character states

    -  The time-order event history A is generated
       using a probabilistic process, and hence along with A
       the probability of A is also returned, which we denote density_A.

    -   essentially, A is the input phylogenetic tree G 
        augmented so that migrations (i.e., character state transitions) 
        are assigned to each branch of the phylogenetic tree. The order of
        occurrence (i.e., the time order) of the migrations will be known, but not
        the actual times. See Figure 3.

        * all classes are defined in is_classes.py*

        G is an instance of class Tree, and has a data attribute
        levels[], which contains a level-by-level representation of the
        tree (G.levels[k] contains level n-k, where n is the total
        number of leaves in G, and 0 <= k <= n-1). 
        
        Each level is an instance of class Level, and has a data
        attribute event_history which is left unpopulated when the
        level instance is created. This (the current) function
        populates the event history.

        Suppose level = 4, and that there are 4 lineages 1, 2, 3, & 4.
        And let the character state assignments be {1:0, 2:1, 3:1,
        4:1}. And let following events happen before the chain goes to
        level 3:

        a. lineage 2 goes from character state 1 to character state 0
        b. lineage 3 goes from character state 1 to character state 0
        c. lineages 2 & 3 coalesce to form a parent lineage 5; and this
           ends level 4.

        The event history for level 4 will be a list: [{1:0, 2:1, 3:1,
        4:1}, 2, 3]. The first element is the character state at the
        beginning. The other elements are the lineages that undergo
        migration. Given that lineages only migrate between 2 states,
        it is sufficient to list the lineages that migrated (their
        states can be inferred).

        For the next level, (i.e., level 3) the initial character state
        assignment will be: {1:0, 4:1, 5:0} (since 2 & 3 were in state
        0, their parent 5 also will be in state 0).

    Details: 
    -------

    if the event_history is     [{0: 1, 1: 1, 2: 1, 3: 1}]
                                [{2: 1, 3: 1, 4: 1}, 4]
                                [{4: 0, 5: 1}, 5]
                                
    all_delta_earlier is        {0: 1, 1: 1, 2: 1, 3: 1, 4: 0, 5: 0}
    
    all_delta_later is          {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1}

    
    The input phylogenetic tree can be viewed as a tuple (tau, BRL(tau))
    (see also page 7, bullet point entitled "Sample history")
    where tau is the tree topology with *with known time order of branching
    events*, BRL(tau) are the actual time durations assigned to the
    branches of tau. 

    The function follows the following basic scheme in assigning 
    migrations to the branches of G as follows:

    G = (tau, BRL(tau)) ===> tau ===> (tau, MIG(tau)) ===> (tau, MIG(tau), BRL(tau))
                                  /           |        1
                                 /            |        1---------------> here we are just sticking BRL(tau) back.
                                /             V        
                               /        the migrations assigned to
                              /            branches of tau  
                             |
                             |
                             V
                 this is a stochastic process (a discrete-time Markov
                 chain, to be precise);
                 MIG(tau) is sampled from a probability distribution such
                 that its density under the distribution is Pr(MIG(tau) |
                 tau, delta, sigma). Sections 5.1-5.3 are devoted to this.
                 
                 THE DISCRETE-TIME MARKOV CHAIN INVOLVED IN THIS IS THE
                 CONDITIONAL JUMP CHAIN OF SECTION 5.3

    # A detailed description of the conditional jump  chain(PLEASE READ, HIGHLY IMPORTANT)
    # ------------------------------------------------------------------------------------
    #  
    # IMPORTANT: A LINEAGE IS A SPECIES IN THE SAMPLE OR AN ANCESTOR OF
    #            A SPECIES IN THE SAMPLE.
    #
    # 1. The conditional jump chain has the following properties.
    #        ==========================
    #   
    #       a. it is characterized by a 4-tuple: (set of states, 
    #                                             designated initial state,
    #                                             designated set of terminal states,
    #                                             matrix of probabilities of transitions between states
    #                                             )
    #
    #       b.  each state is a 4-tuple:
    #
    #                n_lineages_in_state_zero              the number of lineages in state 0
    #                n_lineages_in_state_one               the number of lineages in state 1
    #    
    #                                                      |the states of the two lineages 
    #                state_of_next_coalescing_lineages_1   |involved in the next (going
    #                state_of_next_coalescing_lineages_2   |backwards in time) speciation
    #                                                      |event.
    #
    #                The next_coalescing_lineages are known since the tree
    #                is known. Thus, this jump chain is constructed
    #                conditional on the tree being known (and hence its
    #                name).
    #
    #           The 4 variables are equivalent, respectively, to q, r, x_1 and x_2 in page
    #           19, paragraph 3, line 4.
    #
    #           When there there are k lineages (i.e.,
    #           n_lineages_in_state_zero + n_lineages_in_state_one = k), the
    #           chain is said to be in level-k. 
    #                                  ======= 
    #   
    #       c.  its state transitions represent either migration (i.e., change
    #           of state of lineages), or speciation events. Migrations leave
    #           the number of lineages unchanged, while the speciation events
    #           reduce the number by 1. Thus migrations are "within-level"
    #                                                        =============
    #           transitions and speciations are "between-level" transitions. 
    #                                           ===============
    #
    #       d.  the initial state of the chain reflects the situation at the
    #           tips of the input phylogenetic tree and delta: n_lineages_in_state_zero
    #           and n_lineages_in_state_one depend on delta; the next-coalescing
    #           lineages are the two lineages involved in the most recent
    #           speciation event in the tree.
    #
    #       e.  the terminal state(s) of the chain are such that
    #           n_lineages_in_state_zero + n_lineages_in_state_one = 1.
    #      
    #       f.  the chain is designed (i.e., the state transition matrix is
    #           set up) in such a way that the within-level
    #           transitions (i.e., speciation) are exactly those that occur
    #           in the input phylogenetic trees, and occur in the same
    #           order.
    #
    #
    # 2. The transition matrix  and the operation of the conditional jump chain.
    #        =================
    #
    #    The conditional jump chain starts from the initial state, and repeatedly
    #    transitions to other states until a terminal state (i.e., one
    #    where the total number of lineages is 1) is reached. As noted
    #    above in (b), each state transition represents an event - either
    #    migration or speciation - in the history of the sample. 
    #
    #    Generally, if a discrete-Time Markov chain has N states, the
    #    transition matrix is an N x N matrix T such that T[i, j] is the
    #    probability that the next state is j, given that the current state
    #    is i. 
    #    
    #    In our case, we don't define the entire transition matrix all at
    #    once, since it is convenient to do so level-by-level. At each
    #    level, we define a matrix for within-level transitions, and also
    #    calculate the probability of the one between-level transition
    #    that's allowed (namely, that of the speciation involving the two
    #    next-coalescing lineages), and a within-level or between-level
    #    transition is made based on the probabilities. 
    #
    #    The whole operation of the Markov chain can be characterized as
    #    follows:
    #       
    #       level <- n, where n is the number of tips of the input
    #                phylogenetic tree
    #       
    #       UNTIL level 1 is reached, do:
    #
    #           while remaining within the level:
    #               choose either a within-level transition or between-level
    #               transition based on their probabilities.
    #           
    #           level <- level - 1
    """
    alpha=sigma[3]
    mu=sigma[4]
    if(alpha>mu):
        #print(alpha);
        #print(mu);
        return("!!error, alpha is greater than mu")
    # In what follows, uncond is short for unconditional (i.e., not
    # conditioned on the tree), and cond is short for conditional (i.e., 
    # conditioned on the tree). And  prob, of course, is short for probability.
    state_of_cond_jump_chain = GetInitialStateofCondJumpChain(G, delta)

    n_leaves = G.num_leaves
    # at first we are at the n-th level, where n is the number of leaves in # G.
    current_level_number = n_leaves

    # We will generate a time-order event history backwards.
    # This history is a list of states the chain passes through, 
    # The probability_of_history is initially set to 1.
    # the 1.0 will ensure that probability_of_history is a float.
    # Otherwise, python will assume it is an int.
    probability_of_history = 1.0
    
    # with each event, the assignment states to
    # lineages changes. current_delta keeps track of the current assignment
    # of states to lineages.

    # VERY IMPORTANT: delta is a dictionary. And 
    # "current_delta = delta.copy()"  is very different from saying
    # current_delta = delta (this holds for any *mutable* (that is,
    # changeable) object like lists or tuples.). The former, called a "deep
    # copy" creates a copy of delta and maps the name current_delta. 
    # Thus after this assignment,
    # current_delta and delta are pointing to two different objects.
    # Whereas after saying "current_delta = delta", both the names point to
    # the same object - the one which delta was pointing to before the
    # assignment.
    # I am doing a deep copy since I don't want to mess with the original
    # delta.
    current_delta = delta.copy()
    
    all_delta_earlier=delta.copy()
    all_delta_later=delta.copy()

    while not current_level_number == 1:
        #print("level "+str(current_level_number))
        # Note: the following, and all assignments in fact, are copies by
        # reference, since all an assignment does in python 
        # is to map a name to an object.
        current_level = G.levels[n_leaves - current_level_number]
        current_level.event_history.append(current_delta.copy())                  
        # when the conditional jump chain is in level k (i.e., the total
        # number of lineages = k), the number of states at that level 
        # is 4(k-1) (see paragraph 3, page 19). Further, there is one
        # possible transition to level k-1 (which is the speciation
        # involving the next-coalesceing lineages). 
        # Thus, # transition_matrix_for_the_level will be a 4(k-1)+1 x 4(k-1)+1 # matrix.  
        #
        # NOTE: it is not known in advance which state in level k-1 the
        # chain will transition to from level k. But that knowledge is not
        # necessary. From each of the 4(k+1) states in level k, it is
        # possible to calculate the probability of transitioning to level
        # k-1. In fact, that's precisely the reason why we include the
        # character state of the next coalescing lineages in description of
        # the jump chain states. 

        # The left side variable is the state-to-index dictionary for the
        # current level.
        state_to_index_in_transition_matrix = state_to_index_in_transition_matrices[n_leaves - current_level_number]
        # also get the reverse map.
        # The left side variable is the index-to-state dictionary for the
        # current level.
        index_in_transition_matrix_to_state = index_in_transition_matrices_to_state[n_leaves - current_level_number]

        index_of_current_state = state_to_index_in_transition_matrix[state_of_cond_jump_chain]
        transition_matrix_for_the_level = transition_matrices[n_leaves - current_level_number]
    
        # transition_matrix_for_the_level[index_of_current_state] is a row
        # of probabilities of transitions from the current state. 
        # the probabilities in the row must sum up to 1. We can verify that
        # by placing an assert. If the condition fails, the program exits
        # with a message. In python, sum(some_list) returns the sum of the
        # elements of  some_list, if list consists of real numbers.
        assert(1.0-1e5<(sum(transition_matrix_for_the_level[index_of_current_state]))< 1.0+1e-5)

        whether_in_the_same_level = "True"
        while whether_in_the_same_level == "True":
            # pick a next state to transition to such that 
            # Pr(index of next state = j | index of current state = index_of_current_state) = transition_matrix_for_the_level[index_of_current_state][j]
            (index_of_next_state, probability_of_transition) = PickNextStateofChain(transition_matrix_for_the_level[index_of_current_state])
            #print("prob 1 "+str(index_of_next_state)+" "+str(probability_of_transition))
            #print("current_state "+str(state_of_cond_jump_chain))
            if(index_of_next_state==len(index_in_transition_matrix_to_state)):
                whether_in_the_same_level = "False"
            # check if current state and the proposed next state are in the same level
            #whether_in_the_same_level = WhetherInTheSameLevel(state_of_cond_jump_chain, next_state_of_cond_jump_chain)
            if whether_in_the_same_level == "True":
                # migration event is happening. So pick a lineage to
                # migrate. For our specific model we need to pick one of
                # the 1 (neotropical lineages) to migrate to 0 (remember
                # that we are going back in time, so backward migrations
                # are from 1->0. The function ChooseLineageandUpdateDelta
                # updates current_delta and also returns a lineage that
                # migrates. The function should choose one lineage
                # uniformly at random among
                # all lineages that *could* migrate and return it.
                next_state_of_cond_jump_chain = index_in_transition_matrix_to_state[index_of_next_state]
                #print("next state "+str(next_state_of_cond_jump_chain))
                migration_type = MigrationType(state_of_cond_jump_chain, next_state_of_cond_jump_chain)
                (migrating_lineage, current_delta, all_delta_earlier) = ChooseLineageandUpdateDelta(G, current_level_number, migration_type, current_delta, all_delta_earlier)
                current_level.event_history.append(migrating_lineage)
                #print("appended migrating lineage "+str(migrating_lineage))
                
                # updating the state of the chain.
                state_of_cond_jump_chain = next_state_of_cond_jump_chain
                index_of_current_state = index_of_next_state
            else:
                next_level_number = current_level_number-1
                if(next_level_number!=0):
                    next_level = G.levels[n_leaves - next_level_number]
                    # The chain wants to move to the next level (to level k-1 if
                    # the current level is k). 
                    # Basically, the chain wants to effect a speciation
                    # involving the next_coalescing_lineages, and if we know
                    # the character state (0/1) of the next_coalescing_lineages
                    # in the current state of the cond. jump chain, we can use
                    # that to figure out the state of the parent lineage. 
                    # Look at Figure 1 (page
                    # 6). The the two daughter lineages are in states (0, 0),
                    # the parent will be in in state 0. Similaly, (1, 1) -> 1,
                    # and (0, 1) -> 0. 
                    # 
                    # Note 1: The next coalescing lineages are the children of from
                    #         current_level.end_node
                    #
                    # Note 2. The character state of the next coalescing
                    #         lineages can be got from current_delta
                    #
                    # The function MovetoNextLevel will apply these rules and 
                    # (a) update current_delta, 
                    # (b) return the initial state for the next level's conditional jump chain 
                    # (c) initialize the event history for the next level with
                    #     the updated delta
                    (initial_state_in_the_next_level,current_delta,all_delta_earlier, all_delta_later) = MovetoNextLevel(state_of_cond_jump_chain, current_delta, all_delta_earlier, all_delta_later, current_level, next_level)
                    #print((initial_state_in_the_next_level,current_delta,all_delta_earlier, all_delta_later));
                    #print(current_level_number);
                    #print(next_level_number);
                    # and *now* update the state of the chain
                    state_of_cond_jump_chain = initial_state_in_the_next_level
                    # update current level number
                current_level_number=next_level_number
    
            probability_of_history = probability_of_history * probability_of_transition
    current_level = G.levels[n_leaves - 1]
    current_level.event_history.append(current_delta)   

    return(probability_of_history,all_delta_earlier, all_delta_later)            

def PrepareTree():
    """
                  |
                  n6
               /     \
              /       \
             /         n5
            /          /\
           n4         /  \
          /  \       /    \
        n0   n1     n2    n3

    4 levels: from top: 1,2,3,4

    num_leaves=4

    so tree.levels[j] is level 4-j, j=(0,1,2,3)
    """
    #G.levels[0].end_node is the parent of the first two coalescing
    # lineages.
    n0 = Node()
    n0.children = []
    n0.all_leaves = []
    n0.partial_order = []
    n0.num_leaves = 0
    n0.label = "node0"
    n0.index = 0
    
    n1 = Node()
    n1.children = []
    n1.all_leaves = []
    n1.partial_order = []
    n1.num_leaves = 0
    n1.label = "node1"
    n1.index = 1
    
    n2 = Node()
    n2.children = []
    n2.all_leaves = []
    n2.partial_order = []
    n2.num_leaves = 0
    n2.label = "node2"
    n2.index = 2
    
    n3 = Node()
    n3.children = []
    n3.all_leaves = []
    n3.partial_order = []
    n3.num_leaves = 0
    n3.label = "node3"
    n3.index = 3

    n4 = Node()
    n4.children = [n0,n1]
    n4.all_leaves = [n0,n1]
    n4.partial_order = [n0,n1,n4]
    n4.num_leaves = 2
    n4.label = "node4"
    n4.index = 4    

    n5 = Node()
    n5.children = [n2,n3]
    n5.all_leaves = [n2,n3]
    n5.partial_order = [n2,n3,n5]
    n5.num_leaves = 2
    n5.label = "node5"
    n5.index = 5

    n6 = Node()
    n6.children = [n4,n5]
    n6.all_leaves = [n0,n1,n2,n3]
    n6.partial_order = [n0,n1,n4,n2,n3,n5,n6]
    n6.num_leaves = 4
    n6.label = "node6"
    n6.index = 6     

    n0.parent = n4
    n1.parent = n4
    n2.parent = n5
    n3.parent = n5
    n4.parent = n6
    n5.parent = n6
    n6.parent = "None"

    l1 = Level()
    l1.begin_time=-5
    l1.end_time=-4
    l1.lineages = [n6]
    l1.begin_node = n6
    l1.end_node = -1
    l1.event_history = []
    
    l2 = Level()
    l2.begin_time=-4
    l2.end_time=-3
    l2.lineages = [n4,n5]
    l2.begin_node = n5
    l2.end_node = n6
    l2.event_history = []
    
    l3 = Level()
    l3.begin_time=-3
    l3.end_time=-2
    l3.lineages = [n4,n2,n3]
    l3.begin_node = n4
    l3.end_node = n5
    l3.event_history = []
    
    l4 = Level()             #the most recent level, aka the tips
    l4.begin_time=-2
    l4.end_time=0
    l4.lineages = [n0,n1,n2,n3]
    l4.begin_node = n0
    l4.end_node = n4
    l4.event_history = []
    
    levels = [l4,l3,l2,l1] #so levels[0] is the most recent level, aka the tips
    
    x = Tree()
    x.all_leaves = [n0,n1,n2,n3]
    x.partial_order = [n0,n1,n2,n3,n4,n5,n6]
    x.height = 2
    x.num_leaves = 4
    x.root = n6
    x.name = "MyTree"
    x.levels = levels
        
    return x
    
if __name__ == "__main__":
    G = PrepareTree()
    print G.all_leaves
    level_number = 4
    migration_type = M_2
    current_delta = [1,1,0,1]
    all_delta_earlier = [1,1,0,1]
    lineage = ChooseLineageandUpdateDelta(G, level_number, migration_type, current_delta, all_delta_earlier)
    #we expect current_delta to be [1,0,0,1]
    #and lineage = 2
    #print lineage
    
    current_state = (5, 2, 1, 1)
    current_delta = [1,1,0,1]
    all_delta_earlier = {0:1, 1:1, 2:0, 3:1}
    all_delta_later = {0:1, 1:1, 2:0, 3:1}
    current_level = G.levels[0]
    next_level = G.levels[1]
    x = 3
    print "nice" , x ,"nice"