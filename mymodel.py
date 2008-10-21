import pymc
import numpy

#Some data
n = 5*numpy.ones(4, dtype=int)
x = numpy.array([-.86, -.3, -.05, .73])

#Priors on unknown parameters
alpha = pymc.Normal('alpha', mu=0, tau=.01)
beta = pymc.Normal('beta', mu=0, tau=.01)

#Arbitrary deterministic function of parameters
def theta(a=alpha, b=beta, d=dose):
    """theta = logit^{-1}(a+b)"""
    return pymc.invlogit(a+b*d)

#Binomial likelihood for data
d = Binomial('d', n=n, p=theta, value=numpy.array([0., 1., 3., 5.]), isdata=True)
