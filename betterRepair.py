from heapq import *
import math

def repairPartition(G, partition, imbalance = 0.2, isCharged = []):
	n = G.numberOfNodes()

	partition.compact()
	fragmentSet = set(partition)
	k = len(fragmentSet)

	fragmentSizes = [0 for f in fragmentSet]
	fragmentCharges = [[] for f in fragmentSet]
	edgeCuts = [[0 for f in fragmentSet] for v in G.nodes()]

	gapsFound = False

	for v in G.nodes():
		fragmentSizes[partition[v]] += 1
		if isCharged[v]:
			fragmentCharges[partition[v]].append(v)
		if G.hasNode(v+2) and partition[v+2] == partition[v] and G.hasNode(v+1) and partition[v+1] != partition[v]:
			gapsFound = True

		for u in G.neighbors(v):
			edgeCuts[v][partition[u]] += G.weight(v, u)

	maxBlockSize = int(math.ceil(n / k)*(1+imbalance))
	if max(fragmentSizes) <= maxBlockSize and max([len(group) for group in fragmentCharges]) <= 1 and not gapsFound:
		return partition

	maxGain = [- float('inf') for v in G.nodes()]
	maxTarget = [-1 for v in G.nodes()]
	heap = []

	def allowed(v, target):
		allowed = True
		if isCharged[v] and ((len(fragmentCharges[target]) > 0 and v not in fragmentCharges[target]) or len(fragmentCharges[target]) > 1):
			allowed = False
		if fragmentSizes[target] > maxBlockSize:
			allowed = False
		if fragmentSizes[target] == maxBlockSize and partition[v] != target:
			allowed = False
		if G.hasNode(v+2) and partition[v+2] == target and partition[v+1] != target:
			allowed = False
		if v > 0 and G.hasNode(v-1) and G.hasNode(v+1) and partition[v-1] == partition[v+1] and target != partition[v-1]:
			allowed = False
		return allowed


	for v in G.nodes():
		for neighbor in G.neighbors(v):
			target = partition[neighbor]# actually, I only need to iterate over the partitions. This is more for an easier asymptotic running time

			if allowed(v, target) and edgeCuts[v][target] - edgeCuts[v][partition[v]] > maxGain[v]:
				maxGain[v] = edgeCuts[v][target] - edgeCuts[v][partition[v]]
				maxTarget[v] = target

		heappush(heap, (maxGain[v], v))

	visited = [False for v in range(n)]
	assert(len(heap) == n)
	i = 0

	while len(heap) > 0:
		assert(len(heap) +  i == n)
		(key, v) = heappop(heap)
		i += 1
		fragment = partition[v]
		visited[v] = True

		def gapAt(v):
			return v >= 0 and G.hasNode(v) and G.hasNode(v+2) and partition[v] == partition[v+2] and partition[v] != partition[v+1]

		# if fragment of v is alright, skip node
		if fragmentSizes[fragment] <= maxBlockSize and (not isCharged[v] or len(fragmentCharges[fragment]) <= 1) and not gapAt(v-2) and not gapAt(v-1) and not gapAt(v):
			continue

		if key == -float('inf'):
			# new partition necessary
			if partition.upperBound() <= max(fragmentSet)+1:
				partition.setUpperBound(max(fragmentSet)+2)
				fragmentSizes.append(0)
				fragmentCharges.append([])
				for u in G.nodes():
					edgeCuts[u].append(0)
			maxTarget[v] = max(fragmentSet)+1		

		assert(maxTarget[v] >= 0)
		assert(maxTarget[v] < partition.upperBound())

		# otherwise, move v to best allowed fragment and update data structures
		fragmentSizes[partition[v]] -= 1
		fragmentSizes[maxTarget[v]] += 1

		if isCharged[v]:
			fragmentCharges[partition[v]].remove(v)
			fragmentCharges[maxTarget[v]].append(v)
	
		for neighbor in G.neighbors(v):
			edgeCuts[neighbor][partition[v]] -= G.weight(neighbor, v)
			edgeCuts[neighbor][maxTarget[v]] += G.weight(neighbor, v)

		partition[v] = maxTarget[v]

		# update max gains and queue positions of neighbors
		for neighbor in G.neighbors(v):
			if visited[neighbor]:
				continue

			oldKey = maxGain[neighbor]
			for target in fragmentSet:
				if allowed(neighbor, target) and edgeCuts[neighbor][target] - edgeCuts[neighbor][partition[neighbor]] > maxGain[neighbor]:
					maxGain[neighbor] = edgeCuts[neighbor][target] - edgeCuts[neighbor][partition[neighbor]]
					maxTarget[neighbor] = target

			if maxGain[neighbor] != oldKey:
				heap.remove((oldKey, neighbor))
				heapify(heap)
				heappush(heap, (maxGain[neighbor], neighbor))

	return partition