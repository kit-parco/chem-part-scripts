import unittest
import chempartlib2 as chempartlib
import random
import math

from networkit import *

class TestChempartlib(unittest.TestCase):
	
	def test_dpPartition(self):
		runs = 100
		for run in range(runs):
			G = generators.ErdosRenyiGenerator(random.randint(10,100), random.random(), False).generate()
			if run == 0:
				G = readGraph('ubiquitin_complete.graph', Format.METIS)

			n = G.numberOfNodes()
			k = random.randint(2,int(n/2))
			epsilon = random.random()

			part1 = chempartlib.dpPartition(G, k, epsilon, [])
			sizes = part1.subsetSizes()
			self.assertTrue(max(sizes) <= math.ceil(n/k)*(1+epsilon))
			self.assertEqual(len(sizes), k)

			part2 = chempartlib.dpPartition(G, k, epsilon, [], True)
			sizes = part2.subsetSizes()
			self.assertTrue(max(sizes) <= math.ceil(n/k)*(1+epsilon))
			self.assertEqual(len(sizes), k)
			self.assertTrue(min(sizes) >= math.floor(n / k)*(1-epsilon))

	def test_dpPartitionInputCheck(self):
		G = generators.ErdosRenyiGenerator(100, 0.1, False).generate()
		with self.assertRaises(AssertionError):
			chempartlib.dpPartition(G, 1, 0)

		with self.assertRaises(AssertionError):
			chempartlib.dpPartition(G, 101, 0)

		with self.assertRaises(AssertionError):
			chempartlib.dpPartition(G, 10, -1)



if __name__ == '__main__':
    unittest.main()



   #dpPartition(G, k, imbalance, isCharged=[], useLowerBounds=False)