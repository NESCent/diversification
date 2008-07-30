#!/usr/bin/python

# Using CapWords for function, method, and class names
# Using underscored_names for variable names, module and package names.
# Using ALL_CAPS_WITH_UNDERSCORES for file handles

# standard python imports
import os
import sys
import re


def ForwardProbability(A, sigma):
    """
    """
    return(1.0)

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
                (transition_matrices, state_to_index_in_transition_matrices,index_in_transition_matrices_to_state) 
                    
        # SampleFromIS augments G with migration events. So G can now be
        # used in place of A.
        # calculate forward probability of A given sigma, Pr(A | sigma) 
        # according to the procedure described in Section 5.4.
        forward_probability_of_A_given_sigma = ForwardProbability(G, sigma)

        # It is guaranteed that density_of_A > 0
        importance_sampling_weight = forward_probability_of_A_given_sigma/density_of_A
        sum_of_importance_sampling_weights += importance_sampling_weight

    prob_g_and_delta_given_sigma = sum_of_importance_sampling_weights/very_large_number
    return(prob_g_and_delta_given_sigma)

