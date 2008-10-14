import unittest
import bd2

class Testbd2(unittest.TestCase):
    def testMigration(self):
        events={0:[0,0,1,0,0,-1],1:[10,0,1,0,0,-1]}
        events=bd2.Migration(1,20,events)
        print events[1]
        self.assertEqual(events,{0:[0,0,1,0,0,-1],1:[10,0,1,0,1,20]})
        
                                                 
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Testbd2)
    unittest.TextTestRunner(verbosity=2).run(suite)
