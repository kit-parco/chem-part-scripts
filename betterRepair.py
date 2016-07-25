from heapq import *
import math

def repairPartition(G, partition, imbalance = 0.2, isCharged = []):
	n = G.numberOfNodes()
	z = G.upperNodeIdBound()
	if len(isCharged) > 0 and len(isCharged) != z:
		raise ValueError("If charges are given, charge array must have the same size as graph")

	partition.compact()
	fragmentSet = set(partition)
	k = len(fragmentSet)
	maxBlockSize = int(math.ceil(n / k)*(1+imbalance))

	fragmentSizes = [0 for f in fragmentSet]
	fragmentCharges = [[] for f in fragmentSet]
	edgeCuts = [[0 for f in fragmentSet] for v in G.nodes()]

	gapsFound = False

	def gapAt(v, target):
		if not G.hasNode(v):
			return False

		# check whether v is in the middle of a gap
		if v >= 1 and G.hasNode(v-1) and G.hasNode(v+1) and partition[v-1] == partition[v+1] and partition[v-1] != target:
			return True

		#check whether v is directly left of a gap
		if G.hasNode(v+1) and G.hasNode(v+2) and target == partition[v+2] and partition[v+1] != target:
			return True

		#check whether v is directly right of a gap
		if v >= 2 and G.hasNode(v-2) and G.hasNode(v-1) and partition[v-2] == target and partition[v-1] != target:
			return True

		return False

	def sizeAllowed(v, target):
		return fragmentSizes[target] < maxBlockSize or (fragmentSizes[target] == maxBlockSize and partition[v] == target)

	def chargeAllowed(v, target):
		numCharged = len(fragmentCharges[target])
		return len(isCharged) == 0 or not isCharged[v] or numCharged == 0 or numCharged == 1 or fragmentCharges[target] == {v}

	def allowed(v, target):
		return chargeAllowed(v, target) and sizeAllowed(v, target) and not gapAt(v, target)

	def createNewFragment():
		if partition.upperBound() <= max(fragmentSet)+1:
			partition.setUpperBound(max(fragmentSet)+2)
			fragmentSizes.append(0)
			fragmentCharges.append([])
			for u in G.nodes():
				edgeCuts[u].append(0)
		newfrag = max(fragmentSet)+1
		fragmentSet.add(newfrag)
		return newfrag

	# check if already valid and prepare data structures
	for v in G.nodes():
		fragmentSizes[partition[v]] += 1
		if len(isCharged) > 0 and isCharged[v]:
			fragmentCharges[partition[v]].append(v)
		if gapAt(v, partition[v]):
			gapsFound = True

		for u in G.neighbors(v):
			edgeCuts[v][partition[u]] += G.weight(v, u)

	# if partition is already valid, return it unchanged
	if max(fragmentSizes) <= maxBlockSize and max([len(group) for group in fragmentCharges]) <= 1 and not gapsFound:
		return partition

	del(fragmentSizes)

	#first handle charged nodes
	for fragment in fragmentSet:
		while len(fragmentCharges[fragment]) > 1:
			# charged node must be moved. We don't care about the size or gap constraints here, these can be handled later.
			bestMovementCandidate = fragmentCharges[fragment][0]
			bestTargetFragment = -1
			bestGain = -float("inf")

			for chargedNode in fragmentCharges[fragment]:
				for target in fragmentSet:
					gain = edgeCuts[chargedNode][target] - edgeCuts[chargedNode][fragment]
					if chargeAllowed(chargedNode, target) and gain > bestGain:
						bestGain = gain
						bestTargetFragment = target
						bestMovementCandidate = chargedNode

			if bestTargetFragment == -1:
				raise ValueError("Input partition contains multiple charges per fragment and one of them cannot be moved.")
	
			fragmentCharges[fragment].remove(bestMovementCandidate)
			fragmentCharges[bestTargetFragment].append(bestMovementCandidate)
		
			for neighbor in G.neighbors(bestMovementCandidate):
				edgeCuts[neighbor][fragment] -= G.weight(neighbor, bestMovementCandidate)
				edgeCuts[neighbor][bestTargetFragment] += G.weight(neighbor, bestMovementCandidate)

			partition[v] = bestTargetFragment

	#then handle gaps
	for v in G.nodes():
		if v > 0 and G.hasNode(v-1) and G.hasNode(v+1) and partition[v-1] == partition[v+1] and partition[v] != partition[v+1]:
			#we have a gap here. If possible, switch move v to surrounding block, not considering size constraints
			if len(isCharged) == 0 or not isCharged[v] or len(fragmentCharges[partition[v+1]]) == 0:
				#move to surrounding block.
				partition[v] = partition[v+1]
				
			else:
				#one of the neighbours must be uncharged, since they belong to the same partition
				assert(not isCharged[v-1] or not isCharged[v+1])
				if not isCharged[v+1]:
					partition[v+1] = partition[v]
				else:
					#switch blocks of with right neighbor
					ownFragment = partition[v]
					partition[v] = partition[v+1]
					partition[v+1] = ownFragment

	#rebuild indices of fragment sizes and charges
	fragmentSizes = [0 for f in fragmentSet]
	fragmentCharges = [[] for f in fragmentSet]
	edgeCuts = [[0 for f in fragmentSet] for v in G.nodes()]

	for v in G.nodes():
		fragmentSizes[partition[v]] += 1
		if len(isCharged) > 0 and isCharged[v]:
			fragmentCharges[partition[v]].append(v)

		for u in G.neighbors(v):
			edgeCuts[v][partition[u]] += G.weight(v, u)

		#charges should be still valid
		assert(chargeAllowed(v,partition[v]))
		#no gaps should be left
		assert(not gapAt(v,partition[v]))

	assert(sum(fragmentSizes) == G.numberOfNodes())
	assert(max([len(chargeList) for chargeList in fragmentCharges]) <= 1)

	#now, build heap of all other nodes and handle size constraints
	maxGain = [- float('inf') for v in G.nodes()]
	maxTarget = [-1 for v in G.nodes()]
	heap = []

	for v in G.nodes():
		for target in fragmentSet:
			if allowed(v, target) and edgeCuts[v][target] - edgeCuts[v][partition[v]] > maxGain[v]:
				maxGain[v] = edgeCuts[v][target] - edgeCuts[v][partition[v]]
				maxTarget[v] = target

		heappush(heap, (-maxGain[v], v))

	visited = [False for v in range(n)]
	assert(len(heap) == n)
	i = 0
	heapify(heap)

	while len(heap) > 0:
		assert(len(heap) +  i == n)
		(key, v) = heappop(heap)
		key *= -1
		#print("i:",i,",key:",key,",node:", v)
		i += 1
		fragment = partition[v]
		visited[v] = True

		# if fragment of v is alright, skip node
		if fragmentSizes[fragment] <= maxBlockSize and (len(isCharged) == 0 or not isCharged[v] or len(fragmentCharges[fragment]) <= 1) and not gapAt(v, partition[v]):
			continue

		if key == -float('inf'):
			#recompute if still the case
			for target in fragmentSet:
				if allowed(v, target) and edgeCuts[v][target] - edgeCuts[v][partition[v]] > maxGain[v]:
					maxGain[v] = edgeCuts[v][target] - edgeCuts[v][partition[v]]
					maxTarget[v] = target
			if maxGain[v] == -float('inf'):
				#now we have a problem. 
				raise RuntimeError("k:"+str(k)+",maxBlockSize:"+str(maxBlockSize)+",v:"+str(v)+", partition"+str(partition))

			## new partition necessary
			#maxTarget[v] = createNewFragment()

		assert(maxTarget[v] >= 0)
		assert(maxTarget[v] < partition.upperBound())
		if not allowed(v, maxTarget[v]):
			errorString = "Node "+str(v)+" cannot be moved to block "+str(maxTarget[v])+" of size "+str(fragmentSizes[maxTarget[v]])
			#print("Node ", v, " cannot be moved to block", maxTarget[v], " of size ", fragmentSizes[maxTarget[v]])
			if not chargeAllowed(v, maxTarget[v]):
				errorString += "\nNode"+str(v)+"is charged and block"+str(maxTarget[v])+"already contains"+str(len(fragmentCharges[maxTarget[v]]))+"charged nodes"
			if not sizeAllowed(v, maxTarget[v]):
				errorString += "\nThe maximum block size is"+str(maxBlockSize)
			if gapAt(v, maxTarget[v]):
				errorString+="\nA gap would result."
			raise RuntimeError(errorString)

		# move v to best allowed fragment and update data structures
		fragmentSizes[partition[v]] -= 1
		fragmentSizes[maxTarget[v]] += 1

		if len(isCharged) > 0 and isCharged[v]:
			fragmentCharges[partition[v]].remove(v)
			fragmentCharges[maxTarget[v]].append(v)
	
		for neighbor in G.neighbors(v):
			edgeCuts[neighbor][partition[v]] -= G.weight(neighbor, v)
			edgeCuts[neighbor][maxTarget[v]] += G.weight(neighbor, v)

		partition[v] = maxTarget[v]

		# update max gains and queue positions of other nodes
		for node in G.nodes():
			if visited[node]:
				continue

			oldKey = maxGain[node]
			maxGain[node] = - float('inf')# reset, since the old target might not be valid any more
			for target in fragmentSet:
				if allowed(node, target) and edgeCuts[node][target] - edgeCuts[node][partition[node]] > maxGain[node]:
					maxGain[node] = edgeCuts[node][target] - edgeCuts[node][partition[node]]
					maxTarget[node] = target

			if maxGain[node] != oldKey:
				heap.remove((-oldKey, node))
				heapify(heap)
				heappush(heap, (-maxGain[node], node))

	assert(max(fragmentSizes) <= maxBlockSize)
	assert(max([len(chargeList) for chargeList in fragmentCharges]) <= 1)
	return partition