import unittest
import cjumpchain
from scipy.misc import comb
from numpy import matrix
from scipy.linalg import solve
from numpy import allclose, array, zeros

class TestCjumpchain(unittest.TestCase):
    def testMakeTupleList(self):
        x = cjumpchain.MakeTupleList(2)
        self.assertEqual(x,[(2, 0, 0, 0), (1, 1, 0, 1), (0, 2, 1, 1), (1, 1, 1, 0)])
        
    def testMakeState2IndexDictionary(self):
        x = cjumpchain.MakeState2IndexDictionary(2)
        self.assertEqual(x, {(1, 1, 0, 1): 1, (1, 1, 1, 0): 3, (2, 0, 0, 0): 0, (0, 2, 1, 1): 2})
    
    def testReverseMap(self):
        x = cjumpchain.ReverseMap(2)
        self.assertEqual(x, {0: (2, 0, 0, 0), 1: (1, 1, 0, 1), 2: (0, 2, 1, 1), 3: (1, 1, 1, 0)})   
        #^Interesting in this case the the map is "ordered," but it's usually not.
        
    def testGetCondJumpChainStateSpace(self):
        x = cjumpchain.GetCondJumpChainStateSpace(2)
        self.assertEqual(x, [(2, 0, 0, 0), (1, 1, 0, 1), (0, 2, 1, 1), (1, 1, 1, 0)])
        
        #Testing Ganesh's comments in MakeTransitionMatrixForLevel method
        state_space_at_the_current_level = cjumpchain.GetCondJumpChainStateSpace(3)
        self.assertEqual(len(state_space_at_the_current_level),4*(3-1))
        
        state_to_index_in_transition_matrix=cjumpchain.MakeState2IndexDictionary(3)

        for j in range(0,4*(3-1)):
            self.assertEqual(j, state_to_index_in_transition_matrix[state_space_at_the_current_level[j]])
    
    def testIsValidState(self): 
        self.assertEqual(False,cjumpchain.IsValidState((5,0,1,0)))                           
        self.assertEqual(False,cjumpchain.IsValidState((4,1,1,1)))
        self.assertEqual(False,cjumpchain.IsValidState((0,2,0,1)))
        self.assertEqual(False,cjumpchain.IsValidState((5,1,1,1)))
        self.assertEqual(False,cjumpchain.IsValidState((5,2,3,0))) 
        self.assertEqual(True,cjumpchain.IsValidState((5,2,1,1)))
        self.assertEqual(True,cjumpchain.IsValidState((6,10,0,1)))
        
    def testApplyEvent(self):
        x = cjumpchain.ApplyEvent((3,2,1,1), "kappa")
        self.assertEqual(x, (3,1,-1,-1))
        
        x = cjumpchain.ApplyEvent((1,2,1,0),"kappa")
        self.assertEqual(x,(1,1,-1,-1))
        
        x = cjumpchain.ApplyEvent((1,0,0,1), "kappa")
        self.assertEqual(x,(-1,-1,-1,-1))
        
        x = cjumpchain.ApplyEvent((5,2,1,1), "m_1")
        self.assertEqual(x, (6,1,0,1))
        
        x = cjumpchain.ApplyEvent((0,6,0,1), "m_2")
        self.assertEqual(x,(-1,-1,-1,-1))
        
        x = cjumpchain.ApplyEvent((2,8,1,1), "s_b_arrow_t")
        self.assertEqual(x,(3,7,1,1))
        
        x = cjumpchain.ApplyEvent((0,1,1,1), "m_1")
        self.assertEqual(x, (-1,-1,-1,-1)) #invalid input
    
        x = cjumpchain.ApplyEvent((2,0,0,0), "s_b_arrow_t")
        #cannot operate "s_b_arrow_t" on (2,0,0,0)
        self.assertEqual(x, (-1,-1,-1,-1))
        
    def testMake2DimensionalArray(self):
        x = cjumpchain.Make2DimensionalArray(5)
        #allclose is a method of testing if two vectors/arrays are equal
        self.assertTrue(allclose(x,[[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0]]))
        
        x[0][3] = 17
        #Test for correct referencing
        #allclose is a method of testing if two vectors/arrays are equal
        self.assertTrue(allclose(x,[[0.0, 0.0, 0.0, 17, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0]]))
        
        #Check if each slot is a float
        self.assertTrue(isinstance(x[0][2],float))
        
        #Test the comment under "Calculating the transition matrix using Equation 24"
        transition_matrix = cjumpchain.Make2DimensionalArray(4*(5-1)+1)
        for j in range(0,4*(5-1)+1):
            self.assertTrue(sum(transition_matrix[j]) == 0.0)
        
    def testCalculateLogPi(self):
        current_state_for_uncond_probs = [1000,1000]
        num_sum=20
        #sigma = (lambda,b,B,alpha,mu)
        sigma = (0.1, 0.2, 2000, 0.10, 0.11)                
        x = cjumpchain.CalculateLogPi(700,sigma)
        self.assertAlmostEqual(x,580.932986844)
        
        x = cjumpchain.CalculateLogPi(999,sigma)
        self.assertAlmostEqual(x,876.20571036)
        
        x = cjumpchain.CalculateLogPi(1001, sigma)
        self.assertAlmostEqual(x, 878.184747198)
        
        sigma = [.01, 8, 1000, 4, 2]
        x = cjumpchain.CalculateLogPi(10, sigma)
        self.assertAlmostEqual(x,-1.0) 
        #since alpha>mu, LogPi should give an error message
        
        sigma = [.01, .1, 1000, .05, .1] #Used in Tally's Mathematica calculation
        x = cjumpchain.CalculateLogPi(10, sigma)
        self.assertAlmostEqual(x, 4.54056838136)
        
        sigma=(.1,.1,1000,.05,.1)
        x = cjumpchain.CalculateLogPi(1500, sigma)
        self.assertAlmostEqual(x, 1487.44829174)

    def testCalculatePi(self):
        current_state_for_uncond_probs = [1000,1000]
        num_sum=20
        #sigma = (lambda, b, B, alpha, mu)
        sigma = (0.1, 0.2, 2000, 0.10, 0.11)                
        x = cjumpchain.CalculatePi(700, sigma)
        self.assertAlmostEqual(x, 8.98718698348e-025)  
        
        x = cjumpchain.CalculatePi(999, sigma)
        self.assertAlmostEqual(x, 5.38881032071e-037)
        
        x = cjumpchain.CalculatePi(1001, sigma)
        self.assertAlmostEqual(x, 4.46246937301e-037)
        
        j = 10
        #sigma = [lambda, b, B, alpha, mu]
        #Now we redefine sigma to be something else, this is allowable in Python
        sigma = [.01, 8, 1000, 4, 2]
        x = cjumpchain.CalculatePi(j, sigma)
        self.assertAlmostEqual(x,0.0) #if alpha>mu, pi ought to be undefined
        
        sigma = [.01, .1, 1000, .05, .1] #Used in Tally's Mathematica calculation
        x = cjumpchain.CalculatePi(j, sigma)
        self.assertAlmostEqual(x, 0.04296875)
        
        sigma=(.1,.1,1000,.05,.1)
        x = cjumpchain.CalculatePi(1500, sigma)
        self.assertAlmostEqual(x, 0)

    def testCalculatePiStar(self):
        #test1
        current_state_for_uncond_probs = [1000,1000]
        num_sum=20
        #sigma = (lambda, b, B, alpha, mu)
        sigma = (0.1, 0.2, 2000, 0.10, 0.11)
        x = cjumpchain.CalculatePiStar(1001, current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(x, 0.0965433633826)
        #test2
        current_state_for_uncond_probs = [1000,1000]
        num_sum=20
        #sigma = (lambda, b, B, alpha, mu)
        sigma=(0.1, .1, 1000, .05, .1)
        x = cjumpchain.CalculatePiStar(1001, current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(x, 0.250000243178)

    def testLogFactorialOfNegative(self):
        x=pow(10,cjumpchain.LogFactorialOfNegative(4))
        y=4*3*2*1
        self.assertAlmostEqual(x,y,1)

    def testCombOfNegative(self):
        x=cjumpchain.CombOfNegative(5,2)
        y=comb(5,2)
        self.assertAlmostEqual(x,y,1)
        
    def testUncondProbSTTGlobal(self):
        current_state_for_uncond_probs = [1000, 1000]
        num_sum=20
        #sigma is defined as = (lambda, b, B, alpha, mu)
        sigma = (.01, .1, 1000, .05, .1)
        x = cjumpchain.UncondProbSTTGlobal(current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(x, 99.7006007003,6)
        #99.70 is the value calculated by Tallly using Mathematica
        
        sigma = (0.1, 0.2, 2000, 0.10, 0.11) #redefine sigma
        x = cjumpchain.UncondProbSTTGlobal(current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(x, 109.0710861,6)
        
    def testUncondProbSBBGlobal(self):
        current_state_for_uncond_probs = [1000,1000]
        num_sum = 20
        sigma = [.01, .1, 1000, .05, .1]
        x = cjumpchain.UncondProbSBBGlobal(current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(x, 0.00999,6)
        
        sigma = (0.1, 0.2, 2000, 0.10, 0.11) #redefine sigma
        x = cjumpchain.UncondProbSBBGlobal(current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(x, 0.024975,6)
        
    def testUncondProbSBarrowTGlobal(self):
        current_state_for_uncond_probs = [1000,1000,0,1]
        num_sum=20
        sigma = [.01, .1, 2000, .05, .1]
        x = cjumpchain.UncondProbSBarrowTGlobal(current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(x, 1.00931556086e-036,6)
    
    def testUncondProbSBarrowTGlobal(self):
        current_state_for_uncond_probs = [1,1,0,1]
        num_sum=20
        sigma = [.01, .1, 2000, .05, .1]
        s_b_arrow_t_global = cjumpchain.UncondProbSBarrowTGlobal(current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(s_b_arrow_t_global, 0.799999237061)        
        
        #Testing sum of m_1+m_2+s_b_arrow_t = s_b_arrow_t_global
        m_1 = cjumpchain.UncondProbM1(current_state_for_uncond_probs, sigma, num_sum)
        m_2 = cjumpchain.UncondProbM2(current_state_for_uncond_probs, sigma, num_sum)
        s_b_arrow_t = cjumpchain.UncondProbSBarrowT(current_state_for_uncond_probs, sigma, num_sum)
        self.assertEqual(m_1, 0)
        self.assertAlmostEqual(m_2, 0.799999237061)
        self.assertEqual(s_b_arrow_t, 0)
        self.assertEqual(m_1+m_2+s_b_arrow_t, s_b_arrow_t_global)
        
    def testUncondProbST(self):
        current_state_for_uncond_probs = [1000,1000,0,1]
        num_sum = 20
        sigma = [.01, .1, 2000, .05, .1]
        x = cjumpchain.UncondProbSBTGlobal(current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(x, 5.04657780429e-037,6)            
    
    def testLogSum(self):
        x = cjumpchain.LogSum(3,4,2)
        #Suppose to get log_2 (24) = 4.585
        self.assertAlmostEqual(x,4.58496250072) 
        
        x = cjumpchain.LogSum(4,4,2)
        self.assertAlmostEqual(x,5)
        
    def testReiterpretUncondEvent(self):
        x = cjumpchain.ReinterpretUncondEvent("kappa", (5,2,1,0))
        self.assertEqual(x, "s_bt")
        
        x = cjumpchain.ReinterpretUncondEvent("s_b_arrow_t", (5,2,1,0))
        self.assertEqual(x, "s_b_arrow_t")
        
        x = cjumpchain.ReinterpretUncondEvent("s_b_arrow_t", (2,0,0,0))
        self.assertEqual(x, "s_b_arrow_t")
    
    def testUncondProbKappa(self):
        current_state_for_uncond_probs = [1,1,0,1]
        num_sum=20
        sigma = [.01, .1, 2000, .05, .1]
        x = cjumpchain.UncondProbKappa("s_bt", current_state_for_uncond_probs, sigma, num_sum=20)
        self.assertAlmostEqual(x, 0.00039999961853)
        
        x = cjumpchain.UncondProbKappa("s_bb", current_state_for_uncond_probs, sigma, num_sum=20)
        self.assertEqual(x,0)
        
        x = cjumpchain.UncondProbKappa("s_t", current_state_for_uncond_probs, sigma, num_sum=20)
        self.assertEqual(x,0)

    def testUnconditionalTransitionProbability(self):
        current_state_of_cond_jump_chain = [1000,1000,1,0]
        sigma = [.01, .1, 2000, .05, .1]
        num_sum=20
        x = cjumpchain.UnconditionalTransitionProbability("kappa", current_state_of_cond_jump_chain, sigma)
        self.assertAlmostEqual(x,7.48833482867e-301)
        #Verifying each of the instantaneous
        #s_bt, s_tt, s_b_arrow_t, s_bb rates
        
        current_state_for_uncond_probs = current_state_of_cond_jump_chain[0:2]
        s_bb = cjumpchain.UncondProbSBBGlobal(current_state_for_uncond_probs, sigma, num_sum)
        self.assertAlmostEqual(s_bb, 0.0024975)
        
        s_tt = cjumpchain.UncondProbSTTGlobal(current_state_of_cond_jump_chain, sigma, num_sum)
        self.assertAlmostEqual(s_tt, 99.7006007003) #Tally's mathematica calculation
        
        s_b_arrow_t = cjumpchain.UncondProbSBarrowTGlobal(current_state_of_cond_jump_chain, sigma, num_sum)
        self.assertAlmostEqual(s_b_arrow_t, 1.49322036556e-298)
        
        s_bt = cjumpchain.UncondProbSBTGlobal(current_state_of_cond_jump_chain, sigma, num_sum)
        self.assertAlmostEqual(s_bt, 7.46610182779e-299)
        
        sum_of_rates = s_bb + s_tt + s_b_arrow_t + s_bt
        self.assertAlmostEqual(sum_of_rates, 99.7030982003)
        
        #Tests below are ENCOURAGING evidence that this method is "really"
        #correct
        current_state_of_cond_jump_chain = [2,0,0,0]
        sigma = [.01, .0001, 2000, .05, .1]
        x = cjumpchain.UnconditionalTransitionProbability("kappa", current_state_of_cond_jump_chain, sigma)
        #Since there are only two species, the "only" possible
        #event is "kappa", which in this case is "s_bb"
        #Furthermore, s_bt=s_tt=s_b_arrow_t=0.0, and s_bb=5e-009, 
        #sum_of_rates=5e-009
        self.assertEqual(x, 1.0)

    def testUncondProbSTT(self):
        current_state_for_uncond_probs = [6,1,0,1]
        num_sum=20
        #sigma is defined as = (lambda, b, B, alpha, mu)
        sigma = (.01, .1, 1000, .05, .1)
        tt = cjumpchain.UncondProbSTTGlobal(current_state_for_uncond_probs, sigma, num_sum)
        bt = cjumpchain.UncondProbSBTGlobal(current_state_for_uncond_probs, sigma, num_sum)
        bb =  cjumpchain.UncondProbSBBGlobal(current_state_for_uncond_probs, sigma, num_sum)
    
        smalltt = cjumpchain.UncondProbSTT(current_state_for_uncond_probs, sigma, num_sum)
        smallbt = cjumpchain.UncondProbSBT(current_state_for_uncond_probs, sigma, num_sum)
        smallbb = cjumpchain.UncondProbSBB(current_state_for_uncond_probs, sigma, num_sum)
        kappatt = cjumpchain.UncondProbKappa("s_tt", current_state_for_uncond_probs, sigma, num_sum)
        kappabt = cjumpchain.UncondProbKappa("s_bt", current_state_for_uncond_probs, sigma, num_sum)
        kappabb = cjumpchain.UncondProbKappa("s_bb", current_state_for_uncond_probs, sigma, num_sum)
        
        #Testing if tt_global is really the sum of smalltt+kappa
        self.assertEqual(tt,smalltt+kappatt)
        self.assertAlmostEqual(bt,smallbt+kappabt)
        self.assertEqual(bb,smallbb+kappabb)
        
    def testLinearEquationSolver(self):
        A=matrix([[1,1,1],[4,4,3],[7,8,5]]) # 3 lines 3 rows
        b = matrix([1,2,1]).transpose()     # 3 lines 1 rows.
        x = cjumpchain.LinearEquationSolver(A, b)
        #all close is a way of testing whether two vectors are equal
        self.assertTrue(allclose(x,matrix([1.,-2.,2.]).transpose()))
    
    def testNormalize(self):
        row = array([1,2,3,4], dtype=float)
        x = cjumpchain.Normalize(row)
        #allclose is a method of testing if two vectors/arrays are equal
        self.assertTrue(allclose(x, [1/10.,2/10.,3/10.,4/10.]))
        
        #Test if it can normalize a zero vector
        row = zeros(4)
        x = cjumpchain.Normalize(row)
        self.assertTrue(allclose(x,[0,0,0,0]))
                                             
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCjumpchain)
    unittest.TextTestRunner(verbosity=2).run(suite)
