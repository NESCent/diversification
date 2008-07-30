#!/usr/bin/python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names, module and package names.
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

#   levels=tree.levels
#    current_level=levels[level_number]
#    event_history=current_level.event_history
#    lineage_names=event_history[0].keys()
#   j=level_num
#    allcomb=[]
#    next_coalesce=CombinationChoose2(lineage_names,allcomb,j)

from is_classes import Level
from is_classes import Tree
from is_classes import Node
from levels import DemarcateLevels
from scipy.misc import comb
from scipy.misc import factorial
import sys
import math

#sys.path=[sys.path,'/Network/Servers/orrorin.nescent.org/Volumes/xraid/home/alexandrabalaban/Documents/Summer2008']
#sys.path=[sys.path,'C:\Documents and Settings\Perry Zheng\workspace\Summer2008\src']    
#sys.path.append('C:\Users\Lonnie\Desktop\Summer2008_2\Summer2008\src')
#sys.path.append('C:\Users\Lonnie\Desktop\Summer2008_2\Summer2008\src')
sys.path.append('C:\Documents and Settings\Owner\Desktop\Summer2008\src')
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

def LogFactorialOfNegative(z,terms=100):
    """
    Returns an approximation of the factorial of a negative number less than -1.

    Input paramters
    ---------------
    n       the negative number less than one
    terms   the number of terms for which to carry out the approximation

    Return value
    ------------
    An approximation of the factorial. Uses the Euler infinite product
    definition of the Gamma function
    """
    z=float(z)
    fact=0
    for n in range(terms):
        n=n+1
        n=float(n)
        fact=fact+math.log(float(pow((1+1/n),z))/float(1+z/n))
    return(fact)

def CombOfNegative(n,k):
    """
    An approximation of the Combination function using FactorialOfNegative\

    Input
    -----
    n,k     Any real number that is not a negative integer

    Output
    ------
    n choose k

        
    """
    log_n=LogFactorialOfNegative(n)
    log_k=LogFactorialOfNegative(k)
    n_fact_over_k_fact=pow(10,log_n-log_k)
    n_minus_k_fact=pow(10,LogFactorialOfNegative(n-k))
    return n_fact_over_k_fact/n_minus_k_fact
    
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
    
       
    bigbig_K = set(['kappa', 's_tt', 's_bb', 's_bt'])

    # m_1 is the migration of next_coalescing_lineage1: i.e., lineage whose
    # state is state_of_next_coalescing_lineages_1;
    # m_2 is the migration of next_coalescing_lineage1: i.e., lineage whose
    # state is state_of_next_coalescing_lineages_2;
    # 's_b_arrow_t' is any other migration
    bigbig_M = set(['m_1', 'm_2', 's_b_arrow_t'])

    # bigbig_H = bigbig_K union bigbig_M. 
    bigbig_H = bigbig_K.union(bigbig_M)

    return((bigbig_K, bigbig_M, bigbig_H))

def GetCondJumpChainStateSpace(level_number):
    """
    Input parameters
    -----
    level_number   The level for which the conditional jump chain space is calculated

    Output
    ------
    A list of all the possible states, where the state's index in the array corresponds to the state's index
    in the state_to_index_in_transistion_matrix dictionary
    """
    state_space=range(4*(level_number+1))
    state_map=ReverseMap(level_number)
    for i in range(len(state_space)):
        state_space[i]=state_map[i]
    
    return(state_space)

def GetLinearEquation(z, state_of_cond_jump_chain):
    """
    """
    q_t=z[0]
    r_t=z[1]
    x1=z[2]
    x2=z[3]
    speciation=''
    if ((x1==1)&(x2==0))|((x1==0)&(x2==1)):
        speciation='s_bt'

    if(x1==1):
        speciation='s_tt'
    if(x1==0):
        speciation='s_bb'

    return([])

def LinearEquationSolver(system_of_equations):
    """
    """
    return([])

def Make2DimensionalArray(dimension):
    """
    """
    return([[]])

def ApplyEvent(z, w):
    """
    """
    return(())

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
    
    Details
    --------    
    Only 3 components of sigma are needed:
    alpha     sigma[3] = per species birth rate in neotropical region; positive constant
    b         sigma[1] = migration rate
    mu        sigma[4] = death rate of species in neotropical region when number of species is i; positive constant
    
    Important!! b/a+j-1 must be nonegative, so should be j

    if alpha=mu, pi will be always zero, and CalculateLogPi returns -1, a flag that alpha=mu
    
    Returns value
    -------------
    p_j = C(b/alpha + j - 1,j)*(alpha/mu)^j*(1-alpha/mu)^(-b/alpha); a float
    if alpha=mu, -1, a flag is returned
    """
    alpha = sigma[3]
    b = sigma[1]
    mu = sigma[4]
    if (alpha==mu):
        return -1
    elif (b/alpha<1):
        coefficient = math.log(float(CombOfNegative(b/alpha+j_species-1,j_species)),alpha/mu)
        ratio = float(j_species)
        ratio2 = math.log(pow(float(1-alpha/mu),-b/alpha),alpha/mu)
        return coefficient+ratio+ratio2
    else:
        coefficient = math.log(float(comb(b/alpha+j_species-1,j_species,exact=0)),alpha/mu)
        ratio = float(j_species)
        ratio2 = math.log(pow(float(1-alpha/mu),-b/alpha),alpha/mu)
        return coefficient+ratio+ratio2
    
def CalculatePi(j_species, sigma):
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
                    
                    lambda is the turnover rate, 
                    b = the effective migration rate, 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the speciation rate in state 1 (the neotropical region), and 
                    mu = the extinction rate in state 1 (the neotropical region).
    
    Details
    --------
    Only 3 components of sigma are needed:
    alpha     sigma[3] = per species birth rate in neotropical region; positive constant
    b         sigma[1] = migration rate
    mu        sigma[4] = death rate of species in neotropical region when number of species is i; positive constant
    
    Important!! b/a+j-1 must be nonegative, so should be j

    
    Returns value
    -------------
    pi_j = C(b/alpha + j - 1,j)*(alpha/mu)^j*(1-alpha/mu)^(-b/alpha); a float
    
    """
    alpha = sigma[3]
    b = sigma[1]
    mu = sigma[4]
    coefficient = comb(b/alpha+j_species-1,j_species,exact=0)
    ratio = pow(float(alpha/mu),j_species)
    ratio2 = pow(float(1-alpha/mu),-b/alpha)
    return coefficient*ratio*ratio2

def CalculatePiStar(n_t, current_state_for_uncond_probs, upper, sigma):
    """
    Calculates p(n_t | q_t, r_t) = pi_star, or eqn (10) in mossEqModel.pdf
    
    Input parameters
    -----------------
    @param n_t      total number of neotropical species in the POPULATION
                    at time t
    @param current_state_for_uncond_probs    
                    the current state for unconditional probability
                    current_state_for_uncond_probs[0] = q_t - number of ancestral species in boreal in sample history at time t
                    current_state_for_uncond_probs[1] = r_t - number of ancestral species in neotropics in sameple history at time t              
                    Note: only [1]=r_t is needed.
    @param upper    an upper limit in the summation
    @param sigma    set of all parameters for the model
                    sigma = (lambda, b, B, alpha, mu)
                    See CalculatePi for more details. 
                    
    Details
    --------
    Only 3 components of sigma are used:
    alpha    sigma[3] = the speciation rate in state 1 (the neotropical region)
    b        sigma[1] = the migration rate
    mu       sigma[4] = the extinction rate in state 1 (the neotropical region)
    
    Return value
    -------------
    p(n_t | q_t, r_t)= pi_n_t / \sum_{k>=r_t}(pi_k); a float
    This is called the stead-state frequency of n_t, normalized to condition
    on the fact that at least r_t neotropical species are known to exist at time t. 
    """
    r_t = current_state_for_uncond_probs[1]
    
    #numerator = CalculatePi(n_t,sigma)
    #denominator=0.0
    
    #for i in range(r_t, upper): 
    #    denominator=denominator+CalculatePi(i,sigma)
    #return numerator/denominator 
    alpha = sigma[3]
    mu = sigma[4]
    log_numerator = CalculateLogPi(n_t,sigma)
    log_denominator=CalculateLogPi(r_t,sigma)
    if(log_numerator==-1):
    #the flag that both the numerator and denominator will be zero
        return 0.0
    else:
        for i in range(r_t,upper):
            log_denominator=LogSum(CalculateLogPi(i,sigma),log_denominator,alpha/mu)

    temp=log_numerator-log_denominator
        
    return pow(alpha/mu,temp)

#I got this function from R documentation
def LogSum(log_x,log_y,base):
    """
    Input Parameters
    ----------------
    log_x   the log of x with base=base 
    log_y   the log of y with base=base
    base    the base of the logs of x and y

    Return Value
    ------------
    log(x+y) with base=base

    Details
    -------
    This equality and function were obtained from the R function logSum in the package TileHMM 
    """
    temp=1+pow(base,log_y-log_x)
    return log_x+math.log(temp,base)

def UncondProbSTT(current_state_for_uncond_probs, upper, sigma):
    """
    Calculates p(S_TT|q_t,r_t), probability of speciation event within 
    the neotropics that is captured in the history of the sample, 
    conditioned on the fact that at time t the numbers of ancestral 
    species in boreal and neotrpoics regions are q_t and r_t, respectively. (actually, only r_t matters)
    See eqn (13) in mossEqnModel.pdf
    
    Input parameters
    ----------------
    @param current_state_for_uncond_probs    
                    the current state for unconditional probability
                    current_state_for_uncond_probs[0] = q_t - number of ancestral species in boreal in sample history at time t
                    current_state_for_uncond_probs[1] = r_t - number of ancestral species in neotropics in sameple history at time t              
                    Only [1] = r_t is needed for this method. 
    @param upper    an upper limit in the summation series
    @param sigma    set of all parameters for the model
                    for the model in mossEqModel.pdf, sigma = (lambda, b, B, alpha, mu), where 
                    lambda is the turnover rate, 
                    b = the migration rate, 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the speciation rate in state 1 (the neotropical region), and 
                    mu = the extinction rates in state 1 (the neotropical region).
 
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
    coefficient = 2*alpha*comb(r_t,2)
    sum=0
    for k in range(r_t,upper):
        sum=sum+ float(CalculatePiStar(k-1,current_state_for_uncond_probs,upper,sigma))/k
    return coefficient*sum

def UncondProbSBB(current_state_for_uncond_probs, upper, sigma):
    """
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
    @param upper     upper limit in summation series, this is NOT needed. 
                     I created this just so the method parameters are identical
    @param sigma     set of all parameters for the model
                     sigma = (lambda, b, B, alpha, mu)
                     See uncondProbSTT for more details
                     
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

def UncondProbSBarrowT(current_state_for_uncond_probs, upper, sigma):    
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
    @param upper    an upper limit in the summation series
    @param sigma     set of all parameters for the model
                     sigma = (lambda, b, B, alpha, mu)
                     See uncondProbSTT for more details
                
    Details
    --------
    Only 4 components of sigma are needed (for instance, some are used only to call CalculatePi; they're marked with *): 
    B        sigma[2] = the total number of species in boreal community (or state 0); an assumption
    alpha    *sigma[3] = birth rate of species in neotropical region when the number of species is i; positive constant
    b        **sigma[1] = migration rate (used in both methods)
    mu       *sigma[4] = death rate of species in neotropical region when number of species is i; positive constant
     
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
    for k in range(r_t,upper):
        sum=sum+CalculatePi(k-1,sigma)/k
    return coef*sum 

def UncondProbSBT(current_state_for_uncond_probs, upper, sigma):
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
    @param upper    an upper limit in the summation series
    @param sigma     set of all parameters for the model
                     sigma = (lambda, b, B, alpha, mu)
                     See uncondProbSTT for more details
                
    Details
    --------
    Only 4 components of sigma are needed (for instance, some are used only to call CalculatePi; they're marked with *. Those
    marked with ** are used for both CalculatePi and UncondProbSBT): 
    B        sigma[2] = the total number of species in boreal community (or state 0); an assumption
    alpha    *sigma[3] = birth rate of species in neotropical region when the number of species is i; positive constant
    b        **sigma[1] = migration rate (used in both methods)
    mu       *sigma[4] = death rate of species in neotropical region when number of species is i; positive constant
     
    Return Value 
    -----------_
    Returns p{S_{BT} | n_t,q_t,r_t} = b*r_t*q_t/B * \sum_{k>=r_t}(pi_{k-1}/k)
    """    
    q_t=current_state_for_uncond_probs[0]
    r_t=current_state_for_uncond_probs[1]
    B=sigma[2]
    b=sigma[1]
    
    coef = b*r_t*q_t/float(B)
    sum=0
    for k in range(r_t,r_t+ upper):
        sum=sum+CalculatePi(k-1,sigma)/k
    return coef*sum  

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
                's_bt'      next-coalescing lineages (see Section 3, bullet
                's_tt'      point entitled "Events", and also Figure 1 on page 6

                'm_1'       migration events involving the first and the
                'm_2'       second of the next-coalescing lineages
                            respectively.
                
                's_b_arrow_t'   any other migration event. The model in
                                mossEqModel allows only migrations of a
                                lineage from character state 0 to character

    sigma     is the set of model all parameters. 

              for the model in mossEqModel.pdf, sigma = (lambda, b, B, alpha, mu), where lambda is the
              turnover rate, b is the migration rate, B is the total number of
              species in state 0 (the boreal region, according to
              mossEqModel.pdf), and alpha and mu are the speciation and
              extinction rates, respectively, in state 1 (i.e., the
              neotropical region).

    
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
                                (c) 's_bt' is the state of the
                                     next-coalescing lineages is (1, 0) or
                                     (0, 1)
    
    events 'm_1' and 'm_2' are both treated as 's_b_arrow_t'.
    """

    # current state that's relevent for the computation of unconditional
    # probabilities. Assuming n_lineages_in_state_zero and
    # n_lineages_in_state_one are the first two elements in
    # current_state_of_cond_jump_chain
    ##I fixed a typo
    current_state_for_uncond_probs = current_state_of_cond_jump_chain[0:2]


    # In python the set() statement makes a set from a list.
    set_of_events = set(['s_tt', 's_bb', 's_b_arrow_t', 's_bt'])

    # This is a map from events to *instantaneous rates*
    # I am suggesting a python dictionary implementation.
    # This is just an initialization. Eventually 
    # these are to be calculated using Equations 13, 17, 21, 22 in Section 5.1
    dictionary_event_to_instantaneous_rate = {'s_tt': 1.0, 's_bb': 1.0, 's_b_arrow_t': 1.0, 's_bt': 1.0}

    # dictionary_event_to_function is a map from events to the functions that will calculate the
    # instantaneous rate for the event.
    # python allows these kind of dictionaries. But any way, I am using
    # features that are really specific to python because it is easy to
    # express the intuition in python syntax. But it doesn't mean we have
    # to stick to this implementation.
    # UncondProbSTT implements Equation 13.
    # UncondProbSBB implements Equation 17.
    # UncondProbSBarrowT implements Equation 21.
    # UncondProbSBT implements Equation 22.
    dictionary_event_to_function = {'s_tt': UncondProbSTT, 's_bb': UncondProbSBB, 's_b_arrow_t': UncondProbSBarrowT, 's_bt': UncondProbSBT}

    # e loops over events in set(['s_tt', 's_bb', 's_b_arrow_t', 's_bt'])
    # we will also need the sum of all rates, for calculating the
    # probabilities from rates (basically, the sum is the denominator is
    # Eq. 23).
    sum_of_rates = 0.0
    for e in set_of_events:
        # Call UncondProbSTT, UncondProbSBB, UncondProbSBarrowT or
        # UncondProbSBT depending on what the event is.
        uncond_rate = dictionary_event_to_function[e](current_state_for_uncond_probs, sigma)
        dictionary_event_to_instantaneous_rate[e] = uncond_rate
        sum_of_rates += uncond_rate

    # This is basically Eqn. 23, with the numerator and denominator.
    probability_of_event = dictionary_event_to_instantaneous_rate[event]/sum_of_rates
    return(probability_of_event)

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

    A dictionary matching strings representing states of the conditional
    jump chain to indices in the transition matrix

    Details
    -------

    Here, strings represent states which correspond to 4-tuples.
    For example, state=4-tuple=[0,3,0,1] is represented by '[0,3,0,1]'
    """
    tuple_list=MakeTupleList(level_number)
    tuple_string_dictionary={}
    for i in range(len(tuple_list)):
        tuple_string_dictionary[tuple_list[i]]=i        
        
    return(tuple_string_dictionary)

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
    A array matching indices to 4-tuple states
    
    """
   
    tuples=range(4*(level_number+1))
    j=0
    for i in range(level_number+1):
        tuples[i]=(i,level_number-i,0,0) #A tuple is represented by (), not [], as before

    for i in range(level_number+1):
        j=i+1*(level_number+1)
        tuples[j]=(i,level_number-i,0,1)
        
    for i in range(level_number+1):
        j=i+2*(level_number+1) 
        tuples[j]=(i,level_number-i,1,1)

    for i in range(level_number+1):
        j=i+3*(level_number+1)
        tuples[j]=(i,level_number-i,1,0)
        
    return(tuples)
        
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

    A dictionary matching indices in the transition matrix to strings representing states of
    the conditional jump chain
    ***in this case indices could be matched to lists representing states of the conditional
    jump chain instead of strings representing states of the conditional jump chain

    Details
    -------

    Here, strings represent states which correspond to 4-tuples.
    For example, state=4-tuple=[0,3,0,1] is represented by '0,3,0,1'
    """
    tuple_list=MakeTupleList(level_number)
    tuple_string_dictionary={}
    for i in range(len(tuple_list)):     
        tuple_string_dictionary[i] = tuple_list[i]
        
    return(tuple_string_dictionary)

def PickNextStateofChain(row_of_transition_matrix):
    """
    """
    return((0, 0.0))

def WhetherInTheSameLevel(state1, state2):
    """
    """
    return(False)

def ChooseLineageandUpdateDelta(current_delta):
    """
    """
    return(None)

def MovetoNextLevel(current_delta, current_level, next_level):
    """
    """
    return(())

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
    
    transition_matrix                       4(k+1)+1 x 4(k+1)+1 matrix of floats, where 
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
    number of lineages = k), the number of states at that level is 4(k+1) (see
    paragraph 3, page 19). Besides transitioning within level k, i.e.,
    among its 4(k+1) states, the chain can also effect precisely one
    transition to a state in level k-1, representing the speciation involving the two
    next_coalescing_lineages whose character states (whether 0/1) are
    represented in state_of_cond_jump_chain. We add this state as the
    4(k+1)+1-th state in the transition matrix.
    
    The precise identity of this state in level k-1 does not matter. 
    All that we need is to be able calculate the probability of the 
    above-mentioned speciation event  from each of the 4(k+1) states in level k. 
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
    assert(len(state_space_at_the_current_level) == 4*(k+1))
    # In python, range(0, n) = [0, 1, 2, ..., n-1]
    for j in range(0, 4*(k+1)):
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
    system_of_equations = [GetLinearEquation(z, state_of_cond_jump_chain) for z in state_space_at_the_current_level]
        
    # May be we can use an equation solver from scipy? 
    # In any case, solution *must* be a vector of length 4(k+1), such that
    # solution[j] = Pr(F | z) where z = state_space_at_the_current_level[j]
    solution = LinearEquationSolver(system_of_equations)
    # ======================================================================
    # End Task: Solving the system of equations represented by Equation 25 #
    # ======================================================================


    # =================================================================
    # Begin Task: Calculating the transition matrix using Equation 24 #
    # =================================================================

    # make a 2-dimensional array of dimensions 4(k+1)+1 x 4(k+1)+1
    # Make2DimensionalArray should initialize the array with 0.0
    transition_matrix = Make2DimensionalArray(4*(k+1)+1)
        
    # This is to make sure that the array is initialized properly.
    # In python, range(0, n) = [0, 1, 2, ..., n-1]
    for j in range(0, 4*(k+1)+1):
        # assert that each row sums to 0.0. If the assert fails, the
        # program will exit at this point
        assert(sum(transition_matrix[j]) == 0.0)

    # The following loop fills out the transition matrix state by state. 
    # Each state corresponds to a row in the matrix. 
    # Thus, each iteration of the following for loop handles one row of the
    # transition matrix.
    #
    # Remember that state_space_at_the_current_level contains the 4(k+1) states at
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
        # probabilities are already correctly set to 0.0
        for w in bigbig_H:
            # ApplyEvent will take the current state (z, in Equation 24) and an
            # event (w, in Eq. 24) and return the state w(z). That is, the state
            # that the chain transitions when event w happens when it (the chain)
            # is in state z.
            w_of_z = ApplyEvent(z, w)
            column = state_to_index_in_transition_matrix[w_of_z]
            # In the numerator of the right side of Eq. 24, 
            #   1. Pr(F | w(z)) is, in terms of the program's variables,
            #      solution[column], since
            #           a. solution[j] = Pr(F | z) where z = state_space_at_the_current_level[j] 
            #           b. state_space_at_the_current_level[w_of_z] = column, 
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
            transition_matrix[row][column] = solution[column] * UnconditionalTransitionProbability(w, z, sigma)
            
        Normalize(transition_matrix[row])
    # =================================================================
    # End Task: Calculating the transition matrix using Equation 24 #
    # =================================================================
    return((transition_matrix, state_to_index_in_transition_matrix, index_in_transition_matrix_to_state))
        
        
                     

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

    state_of_next_coalescing_lineages_1 = delta[next_coalescing_lineages_1]
    state_of_next_coalescing_lineages_2 = delta[next_coalescing_lineages_2]

    return((n_lineages_in_state_zero, n_lineages_in_state_one, state_of_next_coalescing_lineages_1, state_of_next_coalescing_lineages_2))

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
    n_leaves = G.num_leaves()
    
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
        

def SampleFromIS(G, delta, sigma, (transition_matrices, state_to_index_in_transition_matrices, index_in_transition_matrices_to_state)):
    """
    This function assigns migrations (i.e., character state changes) to the
    branches of G, with migrations being drawn from a probability
    distribition such that the density of one assignment of migrations M is
    Pr(M | G, delta, sigma).

    migrations will be assigned to branches with known time-order, but 
    without specifying the actual times of migrations.
    

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
    
    a real number density_A, where A is a 
    a time-order event history A (for use in the right side of
                                     Equation 8)

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

        Suppose level = 4, and that there are 4 lineaages 1, 2, 3, & 4.
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
        assignment will be: {1:0, 4:1, 5:1} (since 2 & 3 were in state
        1, their parent 5 also will be in state 1).

    

    Details: 
    -------

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

    while not current_level_number == 1:
        # Note: the following, and all assignments in fact, are copies by
        # reference, since all an assignment does in python 
        # is to map a name to an object.
        current_level = G.levels[n_leaves - current_level_number]    
        current_level.event_history.append[curr_delta]

        # when the conditional jump chain is in level k (i.e., the total
        # number of lineages = k), the number of states at that level 
        # is 4(k+1) (see paragraph 3, page 19). Further, there is one
        # possible transition to level k-1 (which is the speciation
        # involving the next-coalesceing lineages). 
        # Thus, # transition_matrix_for_the_level will be a 4(k+1)+1 x 4(k+1)+1 # matrix.  
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
        assert(sum(transition_matrix_for_the_level[index_of_current_state]) == 1.0)

        whether_in_the_same_level = True
        while whether_in_the_same_level == True:
            # pick a next state to transition to such that 
            # Pr(index of next state = j | index of current state = index_of_current_state) = transition_matrix_for_the_level[index_of_current_state][j]
            (index_of_next_state, probability_of_transition) = PickNextStateofChain(transition_matrix_for_the_level[index_of_current_state])
            next_state_of_cond_jump_chain = index_in_transition_matrix_to_state[index_of_next_state]

            # check if current state and the proposed next state are in the same level
            whether_in_the_same_level = WhetherInTheSameLevel(state_of_cond_jump_chain, next_state_of_cond_jump_chain)
            if whether_in_the_same_level == True:
                # migration event is happening. So pick a lineage to
                # migrate. For our specific model we need to pick one of
                # the 1 (neotropical lineages) to migrate to 0 (remember
                # that we are going back in time, so backward migrations
                # are from 1->0. The function ChooseLineageandUpdateDelta
                # updates current_delta and also returns a lineage that
                # migrates. The function should choose one lineage
                # uniformly at random among
                # all lineages that *could* migrate and return it.
                migrating_lineage = ChooseLineageandUpdateDelta(current_delta)
                current_level.event_history.append(migrating_lineage)
                
                # updating the state of the chain.
                state_of_cond_jump_chain = next_state_of_cond_jump_chain
                index_of_current_state = index_of_next_state
            else:
                next_level_number - current_level_number-1
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
                initial_state_in_the_next_level = MovetoNextLevel(current_delta, current_level, next_level)

                # and *now* update the state of the chain
                state_of_cond_jump_chain = initial_state_in_the_next_level
                # update current level number
                current_level_number - current_level_number-1

            probability_of_history = probability_of_history * probability_of_transition

    return(probability_of_history)            

if __name__ == "__main__":
    print MakeState2IndexDictionary(2)
    print MakeTupleList(2)