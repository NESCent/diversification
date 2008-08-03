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


def GetInfoForForwardProbMIG(A,level):
    """
    A helper function to get inputs for ProbMIG
    
    Input Parameters
    ----------------
    @param A        the tuple (theta,BRL(theta),MIG(theta))
    @param level    the level in the tree

    Return Values
    -------------
    k                           the number of migration (SBArrowT) events that will occur within the level
    states_within_level         a k+1 length array of all the (r_t,q_bb,q_b_arrow_t,q_bt) states within the level,
                                where states_within_level[i] is the state right after the ith migration
    l                           the time span of level "level"
    n                           the number of levels in A

    Details
    ------
    There are k+1 states within a level since there are k migration events within a level, and at each migration
    the values (r_t,q_bb,q_b_arrow_t,q_bt) change.
    states_within_level[i]= the tuple (r_t,q_bb,q_b_arrow_t,q_bt) representing the state after the i-th migration
    where i ranges from 0 to k.
    """
    ##right now I am returning one value in order to test out other functions
    ##here the level is 6+5+4+3=18
    k=2
    l=100
    states_within_level=range(3)
    states_within_level[0]=(0,5,4,3)
    states_within_level[1]=(1,5,3,3)
    states_within_level[2]=(2,5,2,3)
    return (l,k,states_within_level,1)

def ForwardProbMIG(A,sigma):
    """
    Calculates the probability of all SBArrowT events which is equivalent to
    P(MIG(theta)|sigma)
    Input Parameters
    ----------------
    @param A        the tuple (theta,BRL(theta),MIG(theta))
    @param level    the level in the tree
    """
    # level is 1 here for no particular reason, just need n
    (l,k,states_within_level,n)= GetInfoForForwardProbMIG(A,1)
    mult=1
    for level in range(1,n+1):
        (l,k,states_within_level,n)= GetInfoForForwardProbMIG(A,level)
        mult=mult*ProbKMigrationInL(l,k,states_within_level,sigma)

    return mult        
    

def ProbKMigrationInL(l,k,states_within_level,sigma):
    """
    Returns the forward probability of the k migrations that occur in the level lev. Defined in equation 29.
    Input Parameters
    ---------------
    @param l                            the duration of the level
    
    @param k                            the number of migrations in the level
    
    @param states_within_level          a k+1 length array of all the (r_t,q_bb,q_b_arrow_t,q_bt) states within the level,
                                        where states_within_level[i] is the state right after the ith migration
                                        
    @param sigma                        sigma =(lambda, b, B, alpha, mu), see ForwarRateSTT for more details
    """
    coef=1
    for i in range(1,k+1):
        #states_within_level[i-1] is the state right before the ith migration within the level
        coef=coef*Phi(states_within_level[i-1],sigma)

    sum=0
    #if the input to PhiJK is (j,k,...) the input to Phi, a term in PhiJK, is j-1
    for j in range(1,k+1):
        sum=sum+PhiJK(j,k,states_within_level,sigma)*math.exp(-Phi(states_within_level[j-1],sigma)*l)        
                                                             
    return coef*sum

def PhiJK(j,k,states_within_level,sigma):
    """
    PhiJK is defined in Equation 29
    This is a helper function to ProbKMigrationInL

    Input Parameters
    ----------------
    @param j                    the input to Phi will be j-1
    @param k                    The number of migrations within the level
    @param states_within_level  a k+1 length array of all the (r_t,q_bb,q_b_arrow_t,q_bt) states within the level,
                                where states_within_level[i] is the state right after the ith migration                       
    @param sigma                sigma =(lambda, b, B, alpha, mu), where  
                                lambda is the turnover rate in state 0 (the boreal region), 
                                b = the effective migration rate from state 0 to 1 (boreal to tropical), 
                                B = the total number of species in state 0 (the boreal region)
                                alpha = the per lineage speciation rate in state 1 (the neotropical region), and 
                                mu = the per lineage extinction rate in state 1 (the neotropical region).
    Details
    -------
    #states_within_level[i-1] is the state right before the ith migration within the level
    """
    product=1
    #print("j is" + str(j))
    if(j==1):
        product=1
        #print("j equals 1")

    else:
        for i in range(1,j):
            current_state_i=states_within_level[i-1]
            current_state_j=states_within_level[j-1]
            product=product*(Phi(current_state_i,sigma)-Phi(current_state_j,sigma))
            #print("in first")
            #print("range "+str(1)+" "+str(k))
            #print("i is " +str(i))
            #print(product)

    #if(j==k):
        #print("j equals k")    
    if(j!=k):
        for i in range(j+1,k+1):
            current_state_i=states_within_level[i-1]
            current_state_j=states_within_level[j-1]
            product=product*(Phi(current_state_i,sigma)-Phi(current_state_j,sigma))
            #print("in second")
            #print("range "+str(j+1)+" "+str(k))
            #print("i is " +str(i))
            #print(product)
    
    return pow(product,-1)        

def Phi(current_state,sigma):
    """
    Phi(i) is the forward rate of SBarrowT right before the i-th migration

    Input Parameters
    ----------------
    @param current_state    the current (r_t,q_bb,q_b_arrow_t,q_bt) tuple before the ith migration within the level
                            See ForwardRateSTT for more details
                            
    @param sigma    sigma =(lambda, b, B, alpha, mu), See ForwardRateSTT for more details
    
    Return Value
    ------------
    ForwardRateSBArrowT(current_state,sigma), the forward rate of SBArrowT right before the ith migration within
    the level
    """
    return ForwardRateSBArrowT(current_state,sigma)


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
        forward_probability_of_A_given_sigma = ForwardProbMIG(A, sigma)

        # It is guaranteed that density_of_A > 0
        importance_sampling_weight = forward_probability_of_A_given_sigma/density_of_A
        sum_of_importance_sampling_weights += importance_sampling_weight

    prob_g_and_delta_given_sigma = sum_of_importance_sampling_weights/very_large_number
    return(prob_g_and_delta_given_sigma)

