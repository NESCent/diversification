import pymc as py
import numpy as nu

def MetropolisHastings(G, delta, sigma):
    #First sample (sigma*, G*) from Q(. | sigma, G)
    S = MCMC(G, db='ram')
    S.sample(1000, length=None, verbose=0)
    GStar = uniform(S)
    #proposal_sig = P(G*, delta | sigma) * W(sigma*)
    #              ------------------------------
    #               P(G, delta| sigma) * W(sigma) 
    ProbGStarGivenSigma = LikelihoodOfParameters(GStar, delta, sigma)
    ProbGGivenSigma = LikelihoodOfParameters(G, delta, sigma)  

    alpha = sigma[3]
    beta = sigma[1]
    WeightSigmaStar = CalculateGamma(Gstar, alpha, beta)
    WeightSigma = CalculateGamma(G, alpha, beta)
    proposal_sig = ProbGStarGivenSigma * WeightSigmaStar / (ProbGGivenSigma * WeightSigma)
    
    M = Metropolis(stochastic, scale=1, sig=proposal_sig, dist=None, verbose=0)
    
def CalculateGamma(x, alpha, beta):
    #CalculateGamma(x; alpha, beta) = beta^alpha * e^(-beta*x) * x^(alpha-1)
    #                               --------------------------------------
    #                                        Big-Gamma(alpha)
    return beta**(alpha) * exp(-beta*x) * x**(alpha-1) / nu.gamma(alpha)