#import pymc as pymc
#import numpy as numpy
#from random import gammavariate
import math
import scipy.special 
import random
import is_likelihood
import cjumpchain

#def MetropolisHastings(G, delta, sigma):
def MetropolisHastings(delta, runs, a):
    #First sample (sigma*, G*) from Q(. | sigma, G)
    #an Independence Chain Metropolis-Hastings since the proposal distribution does not depend on the
    #current proposed value of sigma
    #sigma = (lambda,b,B,alpha,mu)
    #mean of a gamma distribution is alpha*beta, sd is sqrt(alpha)*beta, so make alpha one for now
    sigma = (0.1, 0.2, 100, 0.10, 0.11)
    lam = sigma[0]
    beta = sigma[1]
    B = sigma[2]
    alpha = sigma[3]
    mu = sigma[4]
    G=cjumpchain.PrepareTree();
    sigmas=range(runs);
    for i in range(runs):
        sigmas[i]=range(len(sigma));
    for i in range(runs):
        lambdaStar = random.gammavariate(1,lam)
        betaStar = random.gammavariate(1,beta)
        #next is to create smaller variance while keeping the same mean of B
        BStar = random.gammavariate(float(10000000),float(B)/float(10000000))
        less_than=False;
        while(not(less_than)):
            alphaStar = random.gammavariate(1,alpha)
            muStar = random.gammavariate(1,mu)
            if(alphaStar<muStar):
                less_than=True;
        #print("alphaStar "+str(alphaStar));
        #print("muStar "+str(muStar));
        #print("betaStar "+str(betaStar));
        
        sigmaStar= [lambdaStar, betaStar, BStar, alphaStar, muStar]
        #print("BStar "+str(BStar));
        
        probBetaStar=CalculateGamma(betaStar,a,beta/a)
        #print("probBetaStar "+str(probBetaStar));
        probLambdaStar=CalculateGamma(lambdaStar,a,lam/a)
        #print("probLambdaStar "+str(probLambdaStar));
        #probBStar=CalculateGamma(BStar,1,B)
        #print("probBStar "+str(probBStar));
        probAlphaStar=CalculateGamma(alphaStar,a,alpha/a)
        #print("probAlphaStar "+str(probAlphaStar));
        probMuStar=CalculateGamma(muStar,a,mu/a)
        #print("probMuStar "+str(probMuStar));

        #wStarOfSigma is the probability of sigma star under the prior distribution, W    
        #wStarOfSigma = probLambdaStar*probBetaStar*probBStar*probAlphaStar*probMuStar
        wStarOfSigma = probLambdaStar*probBetaStar*probAlphaStar*probMuStar
        
        probLambda=CalculateGamma(lam,a,lam/a)
        #print("probLambda "+str(probLambda));
        probBeta=CalculateGamma(beta,a,beta/a)
        #print("probBeta "+str(probBeta));
        #probB=CalculateGamma(B,1,B)
        #print("probB "+str(probB));
        probAlpha=CalculateGamma(alpha,a,alpha/a)
        #print("probAlpha "+str(probAlpha));
        probMu=CalculateGamma(mu,a,mu/a)
        #print("probMu "+str(probMu));
        #wOfSigma = probBeta*probLambda*probB*probAlpha*probMu
        wOfSigma = probBeta*probLambda*probAlpha*probMu

        ratio_of_B=math.exp(-B*(BStar-B));
        #print("ratio of B "+str(ratio_of_B));
        
        ProbGGivenSigmaStar = is_likelihood.LikelihoodOfParameters(G, delta, sigmaStar)
        print("ProbGGivenSigmaStar "+str(ProbGGivenSigmaStar));
        print("wStarOfSigma "+str(wStarOfSigma));
        ProbGGivenSigma = is_likelihood.LikelihoodOfParameters(G, delta, sigma)
        print("ProbGGivenSigma "+str(ProbGGivenSigma));
        print("wOfSigma "+str(wOfSigma));
        hastingsRatio= ProbGGivenSigmaStar *wStarOfSigma/ (ProbGGivenSigma * wOfSigma)*ratio_of_B
        print("hasting Ratio "+str(hastingsRatio));
        u=random.uniform(0,1);
        if(u<hastingsRatio):
            sigma=sigmaStar;
        #print(i)
        print(sigma);
        for j in range(len(sigma)):
            sigmas[i][j]=sigma[j];
            
    return(sigmas)        
    
    
def CalculateGamma(x, alpha, beta):
    #CalculateGamma(x; alpha, beta) = beta^alpha * e^(-beta*x) * x^(alpha-1)
    #                               --------------------------------------
    #                                        Big-Gamma(alpha)
    return beta**(alpha) * math.exp(-beta*x) * x**(alpha-1) / scipy.special.gamma(alpha)