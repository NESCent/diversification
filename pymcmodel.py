#import pymc as pymc
import numpy as numpy
import random
#from random import gammavariate
import math
from scipy.special import *

def MetropolisHastings(G, delta, sigma):
    #First sample (sigma*, G*) from Q(. | sigma, G)
    #an Independence Chain Metropolis-Hastings since the proposal distribution does not depend on the
    #current proposed value of sigma
    sigmas=range(1,2);
    random.randint(0,1);
    random.gammavariate(.3,.5);
    for i in range(1,2):    
        lambdaStar = random.gammavariate(.3,.5)
        betaStar = random.gammavariate(.3,.5)
        BStar = random.gammavariate(.3,.5)
        alphaStar = random.gammavariate(.3,.5)
        muStar = random.gammavariate(.3,.5)
        
        sigmaStar= [lambdaStar, betaStar, BStar, alphaStar, muStar]
        
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
        

        ProbGStarGivenSigma = is_likelihood.LikelihoodOfParameters(GStar, delta, sigma)
        ProbGGivenSigma = is_likelihood.LikelihoodOfParameters(G, delta, sigma)  
        hastingsRatio= ProbGStarGivenSigma *wStarOfSigma/ (ProbGGivenSigma * wOfSigma)
        u=random.uniform(0,1);
        if(u<hastingsRatio):
            sigma=sigmaStar;
        sigmas[i]=sigma;            
    
    
def CalculateGamma(x, alpha, beta):
    #CalculateGamma(x; alpha, beta) = beta^alpha * e^(-beta*x) * x^(alpha-1)
    #                               --------------------------------------
    #                                        Big-Gamma(alpha)
    return beta**(alpha) * math.exp(-beta*x) * x**(alpha-1) / scipy.special.gamma(alpha)