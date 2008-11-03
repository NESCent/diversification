import pymc as pymc
import numpy as numpy
import random as random

def MetropolisHastings(G, delta, sigma):
    #First sample (sigma*, G*) from Q(. | sigma, G)
    #an Independence Chain Metropolis-Hastings since the proposal distribution does not depend on the
    #current proposed value of sigma
    
    x = random()
    lambdaStar = x.gammavariate(.3,.5)
    betaStar = x.gammavariate(.3,.5)
    BStar = x.gammavariate(.3.5)
    alphaStar = x.gammavariate(.3,.5)
    muStar = x.gammavariate(.3,.5)
    
    sigmaStar= [randomLambda, randomBeta, randomB, randomalpha, randomMu]
    
    probBetaStar=CalculateGamma(lambdaStar,.3,.5)
    probLambdaStar=CalculateGamma(betaStar,.3,.5)
    probBStar=CalculateGamma(BStar,.3,.5)
    probAlphaStar=CalculateGamma(alphaStar,.3,.5)
    probMuStar=CalculateGamma(muStar,.3,.5)

    #wStarOfSigma is the probability of sigma star under the prior distribution, W    
    wStarOfSigma = probLambda*probBeta*probB*probAlpha*probMu

    lam = sigma[0]
    beta = sigma[1]
    B = sigma[2]
    alpha = sigma[3]
    mu = sigma[4]
    probBeta=CalculateGamma(lam,.3,.5)
    probLambda=CalculateGamma(beta,.3,.5)
    probB=CalculateGamma(B,.3,.5)
    probAlpha=CalculateGamma(alpha,.3,.5)
    probMu=CalculateGamma(mu,.3,.5)
    wOfSigma = probBeta*probLambda*probB*probAlpha*probMu    
    
    S = MCMC(G, db='ram')
    S.sample(1000, length=None, verbose=0)
    ran = uniform(1)
    #proposal_sig = P(G*, delta | sigma) * W(sigma*)
    #              ------------------------------
    #               P(G, delta| sigma) * W(sigma) 
    ProbGStarGivenSigma = LikelihoodOfParameters(GStar, delta, sigma)
    ProbGGivenSigma = LikelihoodOfParameters(G, delta, sigma)  

    alpha = sigma[3]
    beta = sigma[1]
    #Right now I am just using GStar as a placeholder for SigmaStar
    #WeightSigmaStar = CalculateGamma(Gstar, alpha, beta)
    #WeightSigma = CalculateGamma(G, alpha, beta)
    hastingsRatio= ProbGStarGivenSigma *wStarOfSigma/ (ProbGGivenSigma * wOfSigma)
    
    #M = Metropolis(stochastic, scale=1, sig=proposal_sig, dist=None, verbose=0)
    
def CalculateGamma(x, alpha, beta):
    #CalculateGamma(x; alpha, beta) = beta^alpha * e^(-beta*x) * x^(alpha-1)
    #                               --------------------------------------
    #                                        Big-Gamma(alpha)
    return beta**(alpha) * exp(-beta*x) * x**(alpha-1) / nu.gamma(alpha)