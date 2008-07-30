from scipy.misc import comb
import time

def calcPi(j,a,b,u):
    """
    Calculates pi_j. 
    
    @param a - birth rate of species when number of species 
                is i; positive constant
    @param b - a fixed increase in # of species
    @param j - state j; a nonnegative integer
    @param u - death rate of species when number of species
                is i; positive constant, a<u
    Note: b/a+j-1 must be nonnegative, so is j
    -----------
    Returns pi_j = C(b/a+j-1,j)*(a/u)^j*(1-a/u)^(-b/a)
    """
    coefficient = comb(b/a+j-1,j,exact=0)
    ratio = pow(float(a/u),j)
    ratio2 = pow(float(1-a/u),-b/a)
    return coefficient*ratio*ratio2

def calcPi_Star(n_t,r_t,upper,a,b,u):
    """
    Calculates p(n_t | q_t, r_t) = pi_star. 
    
    @param n_t - total number of neotropical species in the population
                at time t
    @param r_t - number of ancestral species in neotropics 
                in the sample history at time t
    @param upper - an upper limit in the summation
    @param a - birth rate
    @param b - fixed birth constant
    @param u - death rate
    
    -----------
    Returns p(n_t | q_t, r_t)= pi_n_t / \sum_{k>=r_t}(pi_k). 
    """
    numerator = calcPi(n_t,a,b,u)
    denominator=0
    for i in range(r_t, upper): 
        denominator=denominator+calcPi(i,a,b,u)
    return numerator/denominator    

#def CalculatePiStar2(n_t,r_t,upper,alpha,b,mu):
    """
    Calculates p(n_t | q_t, r_t) = pi_star. 

    
    Input parameters
    -------------------    
    @param n_t - total number of neotropical species in the population
                at time t
    @param r_t - number of ancestral species in neotropics 
                in the sample history at time t
    @param upper - an upper limit in the summation
    @param alpha - birth rate
    @param b - fixed birth constant
    @param mu - death rate
    
    Return value
    ------------
    p(n_t | q_t, r_t)= pi_n_t / \sum_{k>=r_t}(pi_k). 

    numerator = CalculatePi(n_t,alpha,b,mu)
    denominator=0.0
    for i in range(r_t, upper): 
        denominator=denominator+CalculatePi(i,alpha,b,mu)
    return numerator/denominator 
    """
    
def calcS_TT(r_t,upper,a,b,u):
    """
    Calculates p(S_TT|q_t,r_t), probability of speciation event within 
    the neotropics that is captured in the history of the sample, 
    conditioned on the fact that at time t the numbers of ancestral 
    species in boreal and neotrpics regions are q_t and r_t, respectively
    
    @param r_t - number of ancestral species in 
            neotropics in sample history at time t
    @param upper - an upper limit in the summation series
    @param a - birth rate
    @param b - fixed birth constant
    @param u - death rate 
    
    -----------
    Returns p(S_TT|q_t,r_t) = 2*a*C(r_t,2)*\sum_{k>=r_t}(piStar_{k-1}/k)
    """
    coefficient = 2*a*comb(r_t,2)
    sum=0
    for k in range(r_t,upper):
        sum=sum+calcPi_Star(k-1,r_t,upper,a,b,u)/k
    return coefficient*sum

def calcS_BB(lam,q_t,B):
    """
    Calculates p(S_BB | q_t), the probability of speciation event
                in the boreal region that is captured in the history
                of the sample. This is simple since boreal region
                maintains a constant fixed number of species, 
                conditioned on the fact that at time t, the number
                of ancestral species in boreal region is q_t. 
                
    @param lam - lambda - rate of turnover in boreal region
    @param q_t - number of ancestral species in boreal region at time t
    @param B - an assumption, a fixed number of species in the boreal
                community
    ------------
    Returns p(S_BB | q_t) = 2*lam*C(q_t,2)/B^2
    """
    return 2*lam*comb(q_t,2)/(B*B)


def calcS_B_T(q_t, B, r_t, upper, a, b, u):    
    """
    Calculates p(S_{B->T}), the probability of migration of species 
    from the boreal to the neotropics region when the duplicate in the
    boreal region does NOT appear in the sample.

    @param q_t - number of ancestral species in the boreal region
                in the sample history at time t
    @param B - fixed number of species in the boreal community
    @param r_t - number of ancestral species in 
                neotropics in sample history at time t
    @param upper - an upper limit in the summation series
    @param a - birth rate
    @param b - fixed birth constant
    @param u - death rate 
    -----------_
    Returns p(S_{B->T}) = b*r_t*(1-q_t/B)*\sum_{k>=r_t}(pi_{k-1}/k) 
    """
    coef = b*r_t*(1-float(q_t/B))
    sum=0
    for k in range(r_t,upper):
        sum=sum+calcPi(k-1,a,b,u)/k
    return coef*sum    

def calcS_BT(q_t,B,r_t,upper,a,b,u):  
    """
    Calculates p(S_{BT} | n_t,q_t,r_t), the probability of migration 
    of species from the boreal to the neotropics region when the duplicate 
    in the boreal region also appears in the sample history, aka a 
    "pseudo-speciation" event in the history of the sample, conditioned 
    on n_t, q_t, r_t. 
    
    @param q_t - number of ancestral species in the boreal region
                in the sample history at time t
    @param B - fixed number of species in the boreal community
    @param r_t - number of ancestral species in 
                neotropics in sample history at time t
    @param upper - an upper limit in the summation series
    @param a - birth rate
    @param b - fixed birth constant
    @param u - death rate 
    -----------
    Returns p{S_{BT} | n_t,q_t,r_t} = b*r_t*q_t/B * \sum_{k>=r_t}(pi_{k-1}/k)
    """      
    coef = b*r_t*q_t/float(B)
    sum=0
    for k in range(r_t,upper):
        sum=sum+calcPi(k-1,a,b,u)/k
    return coef*sum  


if __name__ == "__main__":
    print time.time()
    print calcS_TT(1000,1020,.05,.1,.1)
    print time.time()