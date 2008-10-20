import unittest
import bd2

class Testbd2(unittest.TestCase):
    #the node array
    #[ Time of Event (0), Left Child (1),Right Child (2), Parent node (3), region (4), migration_time (5)]
    def testMigration(self):
        nodes={0:[0,0,1,0,0,-1],1:[10,0,1,0,0,-1]}
        nodes=bd2.Migration(1,20,nodes)
        self.assertEqual(nodes,{0:[0,0,1,0,0,-1],1:[10,0,1,0,1,20]})

    def testBirth(self):
        nodes={0:[0,0,1,0,0,-1],1:[10,0,1,0,0,-1]}
        nodes=bd2.Birth(1,20,nodes)
        self.assertEqual(nodes,{0:[0,0,1,0,0,-1],1:[10,0,1,0,0,-1],2:[20,1,2,1,0,-1]})

    def testDeath(self):
        nodes={0:[0,0,1,0,0,-1],1:[10,0,1,0,0,-1]}
        nodes=bd2.Death(1,nodes)
        self.assertEqual(nodes,{0:[0,0,1,0,0,-1],1:[10,-1,-1,0,0,-1]})
                                            
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Testbd2)
    unittest.TextTestRunner(verbosity=2).run(suite)
