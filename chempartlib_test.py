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
			self.assertEqual(part1.numberOfSubsets(), k)
			self.assertTrue(partitioning.computeImbalance(part1, G) <= 1+epsilon)# if this fails after the previous size check worked, your version of networkit is probably outdated
			self.assertTrue(chempartlib.partitionValid(G, part1, math.ceil(n/k)*(1+epsilon)))

			naive = chempartlib.naivePartition(G, k)
			if naive.numberOfSubsets() == k:
				self.assertTrue(partitioning.computeEdgeCut(part1, G) <= partitioning.computeEdgeCut(naive, G))

			part2 = chempartlib.dpPartition(G, k, epsilon, [], True)
			sizes = part2.subsetSizes()
			self.assertTrue(max(sizes) <= math.ceil(n/k)*(1+epsilon))
			self.assertEqual(len(sizes), k)
			self.assertTrue(min(sizes) >= math.floor(n / k)*(1-epsilon))
			self.assertTrue(partitioning.computeImbalance(part2, G) <= 1+epsilon)
			self.assertTrue(chempartlib.partitionValid(G, part2, math.ceil(n/k)*(1+epsilon)))

			naiveSizes = naive.subsetSizes()
			if len(naiveSizes) == k and min(naiveSizes) >= math.floor(n / k)*(1-epsilon):
				self.assertTrue(partitioning.computeEdgeCut(part2, G) <= partitioning.computeEdgeCut(naive, G))

			self.assertTrue(partitioning.computeEdgeCut(part1, G) <= partitioning.computeEdgeCut(part2, G))

	def test_dpPartitionCharged(self):
		runs = 100
		for run in range(runs):
			G = generators.ErdosRenyiGenerator(random.randint(10,100), random.random(), False).generate()
			n = G.numberOfNodes()
			k = random.randint(2,int(n/2))
			minGapSize = random.randint(0,3)

			epsilon = random.random()

			chargedNodes = random.sample(G.nodes(), int(k*0.8))
			isCharged = [v in chargedNodes for v in G.nodes()]

			try:
				part = chempartlib.dpPartition(G, k, epsilon, isCharged)
				for a in chargedNodes:
					for b in chargedNodes:
						if part[a] == part[b]:
							self.assertEqual(a,b)
				self.assertTrue(chempartlib.partitionValid(G, part, math.ceil(n/k)*(1+epsilon), isCharged, minGapSize))
			except ValueError as e:
				pass

	def test_dpPartitionNegativeWeights(self):
		runs = 100
		absSum = 0
		absSumNaive = 0
		for run in range(runs):
			G = generators.ErdosRenyiGenerator(random.randint(10,100), random.random(), False).generate()
			if run == 0:
				G = readGraph('ubiquitin_complete.graph', Format.METIS)

			n = G.numberOfNodes()
			k = random.randint(2,int(n/2))
			epsilon = random.random()

			GWeighted = Graph(n, True)
			for edge in G.edges():
				randomWeight = random.random()*6-3
				GWeighted.addEdge(edge[0], edge[1], randomWeight)

			G = GWeighted

			part1 = chempartlib.dpPartition(G, k, epsilon, [])
			sizes = part1.subsetSizes()
			self.assertTrue(max(sizes) <= math.ceil(n/k)*(1+epsilon))
			self.assertEqual(len(sizes), k)
			self.assertEqual(part1.numberOfSubsets(), k)
			self.assertTrue(partitioning.computeImbalance(part1, G) <= 1+epsilon)# if this fails after the previous size check worked, your version of networkit is probably outdated
			self.assertTrue(chempartlib.partitionValid(G, part1, math.ceil(n/k)*(1+epsilon)))

			naive = chempartlib.naivePartition(G, k)
			absSum += abs(partitioning.computeEdgeCut(part1, G))
			absSumNaive += abs(partitioning.computeEdgeCut(naive, G))
			#print("Run",run,": Naive", partitioning.computeEdgeCut(naive, G), ", DP:", partitioning.computeEdgeCut(part1, G))
			#if naive.numberOfSubsets() == k:
			#	self.assertTrue(abs(partitioning.computeEdgeCut(part1, G)) <= abs(partitioning.computeEdgeCut(naive, G)))

			part2 = chempartlib.dpPartition(G, k, epsilon, [], True)
			sizes = part2.subsetSizes()
			self.assertTrue(max(sizes) <= math.ceil(n/k)*(1+epsilon))
			self.assertEqual(len(sizes), k)
			self.assertTrue(min(sizes) >= math.floor(n / k)*(1-epsilon))
			self.assertTrue(partitioning.computeImbalance(part2, G) <= 1+epsilon)
			self.assertTrue(chempartlib.partitionValid(G, part2, math.ceil(n/k)*(1+epsilon)))

			naiveSizes = naive.subsetSizes()
			#if len(naiveSizes) == k and min(naiveSizes) >= math.floor(n / k)*(1-epsilon):
			#	self.assertTrue(abs(partitioning.computeEdgeCut(part2, G)) <= abs(partitioning.computeEdgeCut(naive, G)))
		print("AbsSum:",absSum, ", AbsSumNaive:",absSumNaive)


	def test_dpPartitionInputCheck(self):
		G = generators.ErdosRenyiGenerator(100, 0.1, False).generate()
		with self.assertRaises(AssertionError):
			chempartlib.dpPartition(G, 1, 0)

		with self.assertRaises(AssertionError):
			chempartlib.dpPartition(G, 101, 0)

		with self.assertRaises(AssertionError):
			chempartlib.dpPartition(G, 10, -1)

	def test_naivePartition(self):
		runs = 100
		for run in range(runs):
			G = generators.ErdosRenyiGenerator(random.randint(10,100), random.random(), False).generate()
			n = G.numberOfNodes()
			k = random.randint(2,int(n/2))

			part = chempartlib.naivePartition(G, k)
			sizes = part.subsetSizes()
			self.assertTrue(max(sizes) <= math.ceil(n/k))
			self.assertTrue(len(sizes) <= k)
			# in edge cases, the naive partition contains less than k fragments

	def test_greedyPartition(self):
		runs = 100
		for run in range(runs):
			G = generators.ErdosRenyiGenerator(random.randint(10,100), random.random(), False).generate()
			numComp = components.ParallelConnectedComponents(G).run().numberOfComponents()
			if numComp > 1:
				continue
				
			n = G.numberOfNodes()
			k = random.randint(2,int(n/2))
			epsilon = random.random()

			chargedNodes = random.sample(G.nodes(), int(k*0.8))
			isCharged = [v in chargedNodes for v in G.nodes()]
			part = chempartlib.greedyPartition(G, k, epsilon, isCharged)
			self.assertTrue(chempartlib.partitionValid(G, part, math.ceil(n/k)*(1+epsilon), isCharged))
			self.assertTrue(part.numberOfSubsets() >= k)

	def test_mlPartition(self):
		runs = 100
		for run in range(runs):
			G = generators.ErdosRenyiGenerator(random.randint(10,100), random.random(), False).generate()
			numComp = components.ParallelConnectedComponents(G).run().numberOfComponents()

			n = G.numberOfNodes()
			k = random.randint(2,int(n/2))
			minGapSize = random.randint(0,3)
			epsilon = random.random()

			chargedNodes = random.sample(G.nodes(), int(k*0.8))
			isCharged = [v in chargedNodes for v in G.nodes()]
			if numComp > 1:
				with self.assertRaises(RuntimeError):
					part = chempartlib.mlPartition(G, k, epsilon, isCharged)
				continue
			else:
				part = chempartlib.mlPartition(G, k, epsilon, isCharged, False, minGapSize)
			self.assertTrue(chempartlib.partitionValid(G, part, math.ceil(n/k)*(1+epsilon), isCharged, minGapSize))
			self.assertEqual(part.numberOfSubsets(), k)

	def test_fmPartition(self):
		runs = 100
		for run in range(runs):
			G = generators.ErdosRenyiGenerator(random.randint(10,100), random.random(), False).generate()
			numComp = components.ParallelConnectedComponents(G).run().numberOfComponents()

			n = G.numberOfNodes()
			k = random.randint(2,int(n/2))
			minGapSize = random.randint(0,3)
			epsilon = random.random()

			chargedNodes = random.sample(G.nodes(), int(k*0.8))
			isCharged = [v in chargedNodes for v in G.nodes()]
			if numComp > 1:
				with self.assertRaises(RuntimeError):
					part = chempartlib.mlPartition(G, k, epsilon, isCharged)
				continue
			else:
				try:
					dynamic = chempartlib.dpPartition(G, k, epsilon, isCharged)
				except ValueError:
					# if no continuous partition is possible, no good starting point for FM exists
					continue
				part = chempartlib.fmPartition(G, k, epsilon, isCharged, minGapSize)
			self.assertTrue(chempartlib.partitionValid(G, part, math.ceil(n/k)*(1+epsilon), isCharged, minGapSize))
			self.assertEqual(part.numberOfSubsets(), k)

	def test_repairPartition(self):
		runs = 100
		for run in range(runs):
			# generating random graph
			G = generators.ErdosRenyiGenerator(random.randint(10,100), random.random(), False).generate()
			n = G.numberOfNodes()
			k = random.randint(2,int(n/2))
			epsilon = random.random()

			#generating random partition, possibly invalid
			badPartition = partitioning.Partition(n)
			badPartition.setUpperBound(k)

			nodesLeft = set(G.nodes())
			
			for f in range(k):
				fragmentSize = random.randint(1, len(nodesLeft) - (k-f) + 1) if f < k - 1 else len(nodesLeft)
				fragment = random.sample(nodesLeft, fragmentSize)
				for i in fragment:
					badPartition.moveToSubset(f, i)
					nodesLeft.remove(i)

			badPartition.compact()
			self.assertTrue(badPartition.numberOfSubsets() == k)
			
			repaired = chempartlib.repairPartition(G, badPartition, epsilon)
			self.assertTrue(chempartlib.partitionValid(G, repaired, math.ceil(n/k)*(1+epsilon)))

if __name__ == '__main__':
    unittest.main()