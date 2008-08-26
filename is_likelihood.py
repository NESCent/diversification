  #!/usr/bin/python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names, module and package names.
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

# standard python imports
import os
import sys
import re
from scipy.misc import comb
import cjumpchain
import math


def GetNumberOfLevels(A):
    return len(A.levels)

def GetLevelLength(A,level):
    n=A.num_leaves
    current_level=A.levels[n-level]
    branch_length=current_level.begin_time-current_level.end_time
    return branch_length

def GetNumberOfMigrationEvents(A,level):
    n=A.num_leaves
    current_level=A.levels[n-level]
    current_event_history=current_level.event_history
    return len(event_history)-1

def GetQ_bbAndQ_bt(A,level):
    n=A.num_leaves
    current_level=A.levels[n-level]
    current_event_history=current_level.event_history
    current_lineages=current_level.lineages
    q_bb=0
    q_bt=0
    for current_lineage in level_lineages:
        current_lineage_index=current_lineage.index
        #see if the lineages is a boreal lineage and not a migrating one
        if(current_event_history[0][current_lineage_index]==0 and (current_lineage_index not in current_event_history)):
            children=current_lineage.children
            child_1=children[0]
            child_2=children[1]
            child_1_index=child_1.index
            child_2_index=child_2.index
            A_levels=A.levels
            child_1_level_index=-1
            child_2_level_index=-1
            found="FALSE"
            #find the levels at which the two children lineages/nodes first appear
            for i in range(1,GetNumberOfLevels(A)+1):
                current_event_history=A_levels[i].event_history
                if(current_event_history[0].has_key(child_1_index)):
                    child_1_level_index=i
                    found="TRUE"
            if(found=="TRUE"):
                break

            for i in range(1,GetNumberOfLevels(A)+1):
                current_event_history=A_levels[i].event_history
                if(current_event_history[0].has_key(child_2_index)):
                    child_2_level_index=i
                    found="TRUE"
            if(found=="TRUE"):
                break            
                    
            child_1_level=A.levels[child_1_level_index]
            child_2_level=A.levels[child_2_level_index]
            child_1_state=child_1_level.event_history[0][child_1_index]
            child_2_state=child_2_level.event_history[0][child_2_index]
            if(child_1_state==1 or child_2_state==1):
                q_bt+=1

            else:
                q_bb+=1

        #end if statement making sure lineages are both boreal and non-migrating
    #end for-loop looping through all lineages in the specified level
    return(q_bt, q_bb)                


def ForwardProbMIG(A,sigma):
    """
    Calculates the probability of all SBArrowT events which is equivalent to
    P(MIG(theta)|sigma)
    Input Parameters
    ----------------
    @param A        the tuple (theta,BRL(theta),MIG(theta))
    @param level    the level in the tree
    """
    # level is 1 here for no particular reason, just need n's value which is independent of level
    n=GetNumberOfLevels(A)
    for level in range(2,n+1):
        l=GetLevelLength(A,level)
        k=GetNumberOfMigrationEvents(A,level)
        (q_bb,q_bt)=GetQ_bbAndQ_bt(A,level)
        mult=mult*ProbKMigrationInL(level,l,k,q_bb,q_bt,sigma)

    return mult        

def GetHBar(i,A,sigma):
    """
    gives the sum of the rates of all possible speciation events in the level i

    Input Parameters
    ----------------
    i       the level
    A       the tree G plus character state assignments, of class Tree
    sigma   see ForwardRateSTT
    Return Value
    -----------
    Sum of GetH(X,i,A,sigma) for X in {'s_tt,s_bb,s_bt'}
    """
    return GetH('s_tt',i,A,sigma)+GetH('s_bb',i,A,sigma)+GetH('s_bt',i,A,sigma)

def GetYi(A,i):
    """
    returns the speciation event that ends level i

    Input
    -----
    @param i    the level
    @param A    The tree G plus character assignments, of class Tree

    Return value
    ------------
    the speciation event that ends level i in Tree A, either s_tt, s_bb, or s_bt
    """
    n=A.num_leaves
    level_i_event_history=A.levels[n-i].eventhistory
    next_level_event_history=A.levels[n-(i+1)].eventhistory
    set_lineage_indices_level_i=set(level_i_event_history[0].keys())
    set_lineage_indices_next_level=set(next_level_event_history[0].keys())
    set_new_lineages_indices=set_lineage_indices_next_level.difference(set_lineage_indices_level_i)
    set_bifurcated_lineage_index=set_lineage_indices_level_i.difference(set_lineage_indices_next_level)
    new_lineages_indices=[set_new_lineages.pop(),set_new_lineages.pop()]
    bifurcated_lineage_index=[set_bifurcated_lineage_index.pop()]

    character_bifurcated=level_i_event_history[0][bifurcated_lineage_index]
    character_new_1=next_level_event_history[0][new_lineages_indices[0]]
    character_new_2=next_level_event_history[0][new_lineages_indices[1]]

    if(character_bifurcated==1):
        return 's_tt'
    if(character_new_1==1 or character_new_2==1):
        return 's_bt'
    else:
        return 's_bb'
    return ()    
    
    return ()
def GetH(X,i,A,sigma):
    """
    Gives the rate of event X in the ith level
    
    Input Parameters
    ----------------
    @param X        either 's_tt', 's_bb' or 's_bt'
    @param i        the level
    @param A        the tree G plus character state assignments, of class Tree
    @param sigma    see ForwardRateSTT

    Return Value
    ------------
    the rate of event X in the ith level, where X is either STT, SBB or SBT
    """
    (q_bb,q_bt)=GetQ_bbAndQ_bt(A,i)
    q_b_arrow_t=len(A.levels[i].event_history)-1
    r_t=i-q_b_arrow_t-q_bb-q_bt
    if(X=='s_tt'):
        return ForwardRateSTT((r_t,q_bb,q_b_arrow_t,q_bt),sigma)
    if(X=='s_bb'):
        return ForwardRateSBB((r_t,q_bb,q_b_arrow_t,q_bt),sigma)
    if(X=='s_bt'):
        return ForwardRateSBT((r_t,q_bb,q_b_arrow_t,q_bt),sigma)
    return ()

def ForwardProbThetaAndBRLGivenMIG(A,sigma):
    """
    """
    n=GetNumberOfLevels(A)
    for i in range(1,n):
        l=GetLevelLength(A,i)
        hbar_i=GetHBar(i,A,sigma)
        Yi=GetYi(A,i)
        mult=mult*hbar_i*math.exp(-hbar_i*l_i)*GetH(Y_i,i,A,sigma)/hbar_i

    l_n=GetLevelLength(A,n)
    hbar_n=GetHBar(n,A,sigma)
    mult=mult*hbar_n*math.exp(-hbar_n*l_n)
    return mult

def ForwardProbA(A,sigma):
    """
    """
    return ForwardProbThetaAndBRLGivenMIG(A,sigma)*ForwardProbMIG(A,sigma)
        
        
    
def ProbKMigrationInL(level,l,k,q_bb,q_bt,sigma):
    """
    Returns the forward probability of the k migrations that occur in the level. Defined in equation 29.
    Input Parameters
    ---------------
    @param level                        the level number
    @param l                            the duration of the level
    
    @param k                            the number of migrations in the level
    
    @param q_bb                         the number of lineages in the level whose next event is a SBB event
    @param q_bt                         the number of lineages in the level whose next event is a SBT event
                                        
    @param sigma                        sigma =(lambda, b, B, alpha, mu), see ForwarRateSTT for more details
    """
    n=A.num_leaves
    coef=1
    q_b_arrow_t=len(A.levels[n-i].event_history)-1
    r_t=level-q_b_arrow_t-q_bb-q_bt
    for i in range(1,k+1):
        coef=coef*Phi(i,(r_t,q_bb,q_b_arrow_t,q_bt),sigma)
    sum=0
    for j in range(1,k+1):
        sum=sum+PhiJK(j,k,(r_t,q_bb,q_b_arrow_t,q_bt),sigma)*math.exp(-Phi(j,(r_t,q_bb,q_b_arrow_t,q_bt),sigma)*l)                                                                   
    return coef*sum

def PhiJK(j,k,(r_t,q_bb,q_b_arrow_t,q_bt),sigma):
    """
    PhiJK is defined in Equation 29
    This is a helper function to ProbKMigrationInL

    Input Parameters
    ----------------
    @param j                    the input to Phi will be j-1
    @param k                    The number of migrations within the level
    @param (r_t,q_bb,q_b_arrow_t,q_bt)  the (r_t,q_bb,q_b_arrow_t,q_bt) tuple before the first migration within the level
                                        See ForwardRateSTT for more details                      
    @param sigma                sigma =(lambda, b, B, alpha, mu), where  
                                lambda is the turnover rate in state 0 (the boreal region), 
                                b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                                B = the total number of species in state 0 (the boreal region)
                                alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                                mu = the per lineage extinction rate in state 1 (the neotropical region).
    """
    product=1
    if(j==1):
        product=1
    else:
        for i in range(1,j):
            product=product*(Phi(i,(r_t,q_bb,q_b_arrow_t,q_bt),sigma)-Phi(j,(r_t,q_bb,q_b_arrow_t,q_bt),sigma))
  
    if(j!=k):
        for i in range(j+1,k+1):
            product=product*(Phi(i,(r_t,q_bb,q_b_arrow_t,q_bt),sigma)-Phi(j,(r_t,q_bb,q_b_arrow_t,q_bt),sigma))
    
    return pow(product,-1)        

def Phi(i,(r_t,q_bb,q_b_arrow_t,q_bt),sigma):
    """
    Phi(i) is the forward rate of SBarrowT right before the i-th migration

    Input Parameters
    ----------------
    @i                                      represents the ith migration
    @param (r_t,q_bb,q_b_arrow_t,q_bt)    the  (r_t,q_bb,q_b_arrow_t,q_bt) tuple before the first migration within the level
                                          See ForwardRateSTT for more details
                            
    @param sigma    sigma =(lambda, b, B, alpha, mu), See ForwardRateSTT for more details
    
    Return Value
    ------------
    ForwardRateSBArrowT(current_state,sigma), the forward rate of SBArrowT right before the ith migration within
    the level
    """
    return ForwardRateSBArrowT((r_t+i-1,q_bb,q_b_arrow_t-i+1,q_bt),sigma)


def ForwardRateSTT(current_state,sigma,num_sum=20):
    """
    Input Parameters
    ----------------

    @param current_state   (r_t,q_bb,q_b_arrow_t,q_bt)
                            Where
                            r_t is the number of ancestral tropical lineages at the current level and time within the level,
                            q_bb is the number of boreal lineages at the current level and time within the level
                            that will undergo a SBB speciation,
                            q_b_arrow_t is the number of boreal lineages at the current level and time within the level
                            that will undergo a SBArrowT migration
                            q_bt is the number of boreal lineages at the current level and time within the level that
                            will undergo a SBT speciation
    
    @param sigma    sigma =(lambda, b, B, alpha, mu), where  
                    lambda is the turnover rate in state 0 (the boreal region), 
                    b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                    mu = the per lineage extinction rate in state 1 (the neotropical region).

    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum

    Return value(s)
    --------------
    The forward probability of the event STT as described by equation 26    
          
    """
    r_t=current_state[0]
    alpha=sigma[3]
    coefficient = 2*alpha*comb(r_t,2)
    sum=0
    current_state_for_uncond_probs=(current_state[1]+current_state[2]+current_state[3],current_state[0])
    for k in range(r_t,r_t+num_sum):
        if (k!=0):
            #In case k=0, we don't want division by 0
            sum=sum+ float(cjumpchain.CalculatePiStar(k-1,current_state_for_uncond_probs,sigma,num_sum))/k
    return coefficient*sum

def ForwardRateSBB(current_state,sigma,num_sum=20):
    """
    Input Parameters
    ----------------

    @param current_state   (r_t,q_bb,q_b_arrow_t,q_bt) See ForwardRateSTT for more information
    
    @param sigma    sigma =(lambda, b, B, alpha, mu), where  
                    lambda is the turnover rate in state 0 (the boreal region), 
                    b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                    mu = the per lineage extinction rate in state 1 (the neotropical region).

    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum

    Return value(s)
    --------------
    The forward probability of the event SBB as described by equation 26    
          
    """
    q_bb=current_state[1]
    lam=sigma[0]
    B=sigma[2]
    return 2*lam*comb(q_bb,2)/(B*B)

def ForwardRateSBT(current_state,sigma,num_sum=20):
    """
    Input Parameters
    ----------------

    @param current_state   (r_t,q_bb,q_b_arrow_t,q_bt) See ForwardRateSTT for more information
    
    @param sigma    sigma =(lambda, b, B, alpha, mu), where  
                    lambda is the turnover rate in state 0 (the boreal region), 
                    b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                    mu = the per lineage extinction rate in state 1 (the neotropical region).

    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum

    Return value(s)
    --------------
    The forward probability of the event SBT as described by equation 26 and 27
          
    """
    q_bt=current_state[3] 
    B=sigma[2]
    b=sigma[1]

    return b*q_bt/float(B)

def ForwardRateSBArrowT(current_state,sigma,num_sum=20):
    """
    Input Parameters
    ----------------

    @param current_state   (r_t,q_bb,q_b_arrow_t,q_bt) See ForwardRateSTT for more information
    
    @param sigma    sigma =(lambda, b, B, alpha, mu), where  
                    lambda is the turnover rate in state 0 (the boreal region), 
                    b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                    B = the total number of species in state 0 (the boreal region)
                    alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                    mu = the per lineage extinction rate in state 1 (the neotropical region).

    @param num_sum  the number of summations to be done, the upper limit of summation is r_t+num_sum

    Return value(s)
    --------------
    The forward probability of the event SBarrowT as described by equation 26 and 27    
          
    """
    q_b_arrow_t=current_state[2] 
    B=sigma[2]
    b=sigma[1]

    return b*q_b_arrow_t/float(B)



def LikelihoodOfParameters(G, delta, sigma):
    """
    This function computes the likelihood of the set of parameters sigma,
    i.e., Pr(G, delta | sigma), according to Equation 8.

    Input parameters
    ----------------

    G      is a rooted phylogenetic tree with branch lengths.
    
    delta  is a map from the tips of G to 0 or 1. That is, delta
           assigns the state 0 or 1 to each tip of G. The state could be
           thought of as a phenotypic character or a geographic location.
           In the manuscript mossEqModel.pdf, state 0 is assumed to be the
           boreal region and state 1 the neotropical region.
    
    sigma  is the set of model all parameters. 

           for the model in mossEqModel.pdf, sigma = (lambda, b, B, alpha, mu), where lambda is the
           turnover rate, b is the migration rate, B is the total number of
           species in state 0 (the boreal region, according to
           mossEqModel.pdf), and alpha and mu are the speciation and
           extinction rates, respectively, in state 1 (i.e., the
           neotropical region).

    Return value(s)
    ---------------
    Pr(G, delta | sigma), computed according to Equation 8.
    """


    # The likelihood is calculated by evaluating the right side of Equation 8.

    # 1. First, get the transition matrices of of the conditional jump chain. 
    #   
    #       a. The conditional jump chain is a discrete-time Markov chain, 
    #          set up so that each run or "realization" of the chain would 
    #          give one A_i in the right side of equation 8.
    # 
    #       b. The transition matrices are computed one for each level in G. A level
    #          is the portion of G between two speciation events. The matrix for a
    #          level depends only on the number of lineages in the level and
    #          sigma (and so all the matrices can be computed given only the number of 
    #          leaves in G and sigma). 
    #
    #       c. transition_matrices[k] should contain the transition matrix for
    #          level n-k, where n is the number of leaves in G, and 0 <= k <= n-2
    #
    #       d. state_to_index_in_transition_matrices[k] is a dictionary for level
    #           n-k, which maps each state (of the jump chain) in level k to its index
    #           in that level's transition matrix.
    #
    #       e.  index_in_transition_matrices_to_state[k] is the reverse map of
    #           state_to_index_in_transition_matrices[k]
    
    (transition_matrices, state_to_index_in_transition_matrices, index_in_transition_matrices_to_state) = cjumpchain.MakeTransitionMatricesbyLevels(G, delta, sigma)

    # very_large_number is analogous to k in the right side of Equation 8.
    very_large_number = 1000000
    sum_of_importance_sampling_weights = 0

    # This loop basically evaluates the rhs of Equation 8. 
    # every iteration of the loop computes an "importance_sampling_weight",
    # and what we want is the average importance_sampling_weight as k goes
    # to infinity in Equation 8.
    for index in range(0, very_large_number):
        # sample 'A' from the importance-sampling distribution IS(.) by
        # running the conditional jump chain once. IS(A) is the density of 'A'
        # under the importance sampling distribution: i.e., the probability
        # that one run of the chain results in 'A'.
        # 
        # For example, IS could be the
        # distribution Upsilon in page 12, paragraph 2, line 3.
        density_of_A = cjumpchain.SampleFromIS(G, delta, sigma,
                transition_matrices, state_to_index_in_transition_matrices,index_in_transition_matrices_to_state) 
                    
        # SampleFromIS augments G with migration events. So A can now be
        # used in place of G.
        # calculate forward probability of A given sigma, Pr(A | sigma) 
        # according to the procedure described in Section 5.4.
        forward_probability_of_A_given_sigma = ForwardProbA(A,sigma)

        # It is guaranteed that density_of_A > 0
        importance_sampling_weight = forward_probability_of_A_given_sigma/density_of_A
        sum_of_importance_sampling_weights += importance_sampling_weight

    prob_g_and_delta_given_sigma = sum_of_importance_sampling_weights/very_large_number
    return(prob_g_and_delta_given_sigma)

