
# coding: utf-8

# In[1]:

# %load chempartlib.py
from networkit import *
import math, sys, subprocess
import time


# In[2]:

def getCutWeight(G, part, v, block):
    n = G.numberOfNodes()
    z = G.upperNodeIdBound()
    assert(len(part) == n)
    assert(v < n)
    assert(block in part)
    
    return sum([G.weight(v, u) for u in G.nodes() if G.hasEdge(v,u) and part[u] == block])


# In[3]:

def continuousPartition(G, X, tolerance, k=0, isCharged=[]):
    """
    Partition G into subsets of size X +/- tolerance and with consecutive node ids.
    Charged nodes are not grouped into the same subset.
    """

    # validate input
    n = G.numberOfNodes()
    if len(isCharged) > 0:
        assert(len(isCharged)==n)
        if k > 0:
            assert(sum(isCharged) <= k)
    else:
        isCharged = [False for i in range(n)]
    assert(X > 0)
    assert(tolerance >= 0)
    if (n % X > 0 or (k>0 and k*X != n)) and tolerance == 0:
        if (k > 0):
            print("Creating", k, "partitions of size", X, "from a graph with", n, " nodes and tolerance 0 is impossible.")
        else:
            print("Creating partitions of size", X, "from a graph with", n, " nodes and tolerance 0 is impossible.")
        print("Set tolerance to 1.")
        tolerance = 1
        
    if k > 0:
        assert(n <= (X + tolerance)*k)
        assert(n >= (X - tolerance)*k)
        
    maxNumberOfPartitions = int(n / max(X - tolerance, 1)) if (k == 0) else k
    
    def w(i,l):
        """
        Weight of cutting after node i.
        """
        assert(i >= 0)
        assert(i < n)
        weightSum = 0
        for k in range(l+1,i+1):
            for neighbor in G.neighbors(k):
                if neighbor <= l or neighbor > i:
                    weightSum += G.weight(k, neighbor)
        return weightSum
    
    # allocate cut and predecessor table
    table = [[float("inf") for j in range(maxNumberOfPartitions)] for i in range(n)]
    pred = [[-1 for j in range(maxNumberOfPartitions)] for i in range(n)]
    
    # fill table 
    for i in range(n):
        for j in range(maxNumberOfPartitions):
            if (j == 0):
                if abs(i+1-X) <= tolerance:
                    table[i][j] = w(i,-1)
            elif (i >= (X-tolerance)):
                windowStart = max(i-(X+tolerance),0)
                
                # make sure that no two charged nodes are in the same partition
                chargeEncountered = False
                for l in reversed(range(windowStart, i+1)):
                    assert(l >= windowStart)
                    if isCharged[l]:
                        if chargeEncountered:
                            windowStart = l
                            break
                        else:
                            chargeEncountered = True
                
                predList = [table[l][j-1] + w(i,l) for l in range(windowStart, min(i,i-(X-tolerance)+1))]
                if (len(predList) > 0):
                    minPred = min(predList)
                    table[i][j] = minPred
                    pred[i][j] = predList.index(minPred) + windowStart
                                      
    # get result from table
    resultlist = [table[n-1][j] for j in range(maxNumberOfPartitions)]
    if len(resultlist) == 0:
        raise ValueError("Combination of parameters allows no partition!")
    
    # if k is given, select best solution with k subsets. If not, select best overall solution
    if (k == 0):
        bestCutValue = min(resultlist)
        k = resultlist.index(bestCutValue) + 1
    else:
        bestCutValue = table[n-1][k-1]
        
    if (bestCutValue == float("inf")):
        raise ValueError("Combination of parameters allows no partition!")
    result = partitioning.Partition(n)
    result.setUpperBound(k)
    
    # search best path backwards
    j = k-1
    i = n-1
    c = bestCutValue
    
    while (j > 0):
        nextI = pred[i][j]
        assert(nextI >= 0)
        # assign partitions to nodes
        for l in range(nextI+1, i+1):
            result[l] = j
        j -= 1
        c -=w(i,nextI)
        i = nextI
        
    # assign partitions to first nodes not covered by previous loop
    for l in range(0, nextI+1):
        result[l] = 0
        
    # check results:
    for i in range(n):
        assert(result[i] >= 0)
        assert(result[i] < k)
        
    for size in result.subsetSizes():
        if (abs(size-X) > tolerance):
            print("For n=", n, ", k=", k, ", X=", X, ", tolerance=", tolerance, ", ", size, " is wrong.")
        assert(abs(size-X) <= tolerance)
    
    return result


# In[4]:

def spiralLayout(G, k, rowheight = 10, colwidth = 10):
    """
    Return two lists, of x and y coordinates for a spiral layout of G.

    k nodes are put in one row, keywords rowheight and colwidth determine spacing
    """
    n = G.numberOfNodes()
    z = G.upperNodeIdBound()
    x = [0 for i in range(z)]
    y = [0 for i in range(z)]
    for i in range(z):
        if G.hasNode(i):
            if int(i / k) % 2 > 0:
                x[i] = colwidth*(k-(i % k)-1)
            else:
                x[i] = colwidth*(i % k)
            
            y[i] = rowheight*int(i / k)
            
            # adapt coordinates for rounded bends
            
            ydelta = int(rowheight / 4)
            xdelta = colwidth*(1-math.cos(math.pi/3))
            rightwards = int(i / k) % 2 == 0
    
            if i % k == k-1:
                y[i] += ydelta
                x[i] = x[i] - xdelta if rightwards else x[i] + xdelta
            if i > 0 and i % k == 0:
                y[i] -= ydelta
                x[i] = x[i] - xdelta if not rightwards else x[i] + xdelta
        
    for i in range(z):
        x[i] += 1# gephi ignores coordinates with value 0
        y[i] += 1
    return x, y


# In[5]:

def naivePartition(G, X):
    """
    Chop a new fragment off G every X nodes
    """
    n = G.numberOfNodes()
    naivePart = partitioning.Partition(n)
    naivePart.allToSingletons()
    for i in range(n):
        naivePart.moveToSubset(int(i/X), i)
    return naivePart


# In[6]:

def mlPartition(G, k, imbalance, isCharged=[]):
    """
    Use a multi-level approach with Fiduccia-Matheyses to partition G.

    Subsets have size at most (1+imbalance)*ceil(n/k)
    """
    n = G.numberOfNodes()
    if len(isCharged) > 0:
        assert(len(isCharged)==n)
        if k > 0:
            assert(sum(isCharged) <= k)
    else:
        isCharged = [False for i in range(n)]
    
    listOfChargedNodes = [i for i in range(n) if isCharged[i]]
    mlp = partitioning.MultiLevelPartitioner(G, k, imbalance, False, listOfChargedNodes, True)
    mlp.run()
    return mlp.getPartition()


# In[7]:

def greedyPartition(G, sizelimit, isCharged=[]):
    """
    Starting with singleton clusters, greedily merge the heaviest edge as long as smaller than sizelimit.
    """
    n = G.numberOfNodes()
    if len(isCharged) > 0:
        assert(len(isCharged)==n)
    else:
        isCharged = [False for i in range(n)]
    n = G.numberOfNodes()
    part = partitioning.Partition(n)
    part.allToSingletons()
    chargedPartitions = set([part.subsetOf(i) for i in range(n) if isCharged[i]])
    
    def getWeight(edge):
        return G.weight(edge[0], edge[1])
    
    sortedEdges = sorted(G.edges(), key=getWeight)
    
    # merge heaviest edge, as long as allowed
    while len(sortedEdges) > 0:
        allowed = True
        heaviestEdge = sortedEdges.pop()
        firstPart = part.subsetOf(heaviestEdge[0])
        secondPart = part.subsetOf(heaviestEdge[1])
        if firstPart in chargedPartitions and secondPart in chargedPartitions:
            allowed = False
        sizeMap = part.subsetSizeMap()
        if sizeMap[firstPart] + sizeMap[secondPart] > sizelimit:
            allowed = False
        partSet = {firstPart, secondPart}
        for i in range(n-2):
            if part[i] in partSet and part[i+2] in partSet and not part[i+1] in partSet:
                allowed = False #otherwise, would create single embedded node
        if allowed:
            part.mergeSubsets(firstPart, secondPart)
            if firstPart in chargedPartitions or secondPart in chargedPartitions:
                chargedPartitions.add(part.subsetOf(heaviestEdge[0]))
    
    part.compact()
    return part


# In[8]:

def exportToGephi(G, xcoords, ycoords, part):
    """
    Export graph to Gephi, along with coordinates and partition
    """
    client = gephi.streaming.GephiStreamingClient()
    client.clearGraph()
    client.exportGraph(G)
    client.exportNodeValues(G, part, "partition")
    client.exportNodeValues(G, xcoords, 'x')
    client.exportNodeValues(G, [-elem for elem in ycoords], 'y')


# In[9]:

def moveAllowed(G, partition, v, targetBlock, maxBlockSize = 0, isCharged = []):
    """
    Returns False/True, depending on whether it is allowed to move node v to block targetBlock without violating the charge and gap constraints
    """
    
    z = G.upperNodeIdBound()
    n = G.numberOfNodes()
    assert(n <= z)
    assert(len(partition) == z)
    assert(G.hasNode(v))
    assert(len(isCharged) == 0 or len(isCharged) == z)
    
    if maxBlockSize == 0:
        maxBlockSize = n
        
    targetBlockSize = sum([partition[u] == targetBlock for u in G.nodes()])
    
    if len(isCharged) == 0:
        isCharged = [False for i in range(z)]
    
    # move would be forbidden if node and target partition are already charged
    if isCharged[v] and any([isCharged[u] for u in G.nodes() if partition[u] == targetBlock and u != v]):
        return False
    
    # move would also be forbidden if gap in main chain of size 1 is created
    if v-2 >= 0 and G.hasNode(v-2) and partition[v-2] == targetBlock and G.hasNode(v-1) and partition[v-1] != targetBlock:
        return False
    
    if G.hasNode(v+2) and partition[v+2] == targetBlock and G.hasNode(v+1) and partition[v+1] != targetBlock:
        return False
    
    # move is forbidden if targetBlock is already at maximum capacity
    if targetBlockSize > maxBlockSize or (targetBlockSize == maxBlockSize and partition[v] != targetBlock):
        return False
    
    # move is allowed!
    return True
    


# In[10]:

def getCutWeights(G, part, v):
    """
    Returns a dictionary with block numbers as keys and cut weights as values
    """
    n = G.numberOfNodes()
    assert(len(part) == n)
    assert(v < n)
    
    sums = {}
    
    for block in set(part):
        sums[block] = 0
        
    for u in G.neighbors(v):
        sums[part[v]] += G.weight(u,v)
        
    return sums


# In[11]:

def repairPartition(G, partition, maxBlockSize = 0, isCharged = []):
    """
    Repairs errors in partition, for example multiple charged nodes in the same partition or gaps of size 1
    """
    z = G.upperNodeIdBound()
    assert(len(partition) == z)
    assert(len(isCharged) == 0 or len(isCharged) == z)
    if len(isCharged) == 0:
        isCharged = [False for i in range(z)]
    
    fragmentSet = {partition[v] for v in range(len(partition))}
    
    def findBestMove(u):
        allowedMoves = [f for f in fragmentSet if moveAllowed(G, partition, u, f, maxBlockSize, isCharged)]
        weights = getCutWeights(G, partition, u)
        if len(allowedMoves) == 0:
            print("No allowed target partition found for node", u, ", creating new partition.")
            newFragmentID= max(fragmentSet)+1
            fragmentSet.add(newFragmentID)
            return newFragmentID

        assert(len(allowedMoves) > 0)
        bestMove = allowedMoves[0]
        for move in allowedMoves:
            if weights[move] > weights[bestMove]:
                #we want to minimize the total weight of cut edges.
                #Thus we select as target partition the one with the highest edge weight to the current node
                bestMove = move
        return bestMove
    
    for u in G.nodes():
        #can the node stay where it is?
        if not moveAllowed(G, partition, u, partition[u], maxBlockSize, isCharged):
            # if not, find the best move for it
            partition[u] = findBestMove(u)
    return partition


# In[12]:

def kaHiPWrapper(G, k, imbalance = 3, pathToKaHiP = '/home/moritzl/Gadgets/KaHIP/deploy/kaffpa', multiple=False):
    """
    Calls KaHiP, an external partitioner.
    """
    
    tempFileName = 'tempForKaHiP.graph'
    outputFileName = 'tmppartition'+str(k)
    n = G.numberOfNodes()
    
    maxWeight = max([G.weight(u,v) for (u,v) in G.edges()])
    
    """
    KaHiP only accepts integer weights, thus we scale and round them.
    Weights must be under 1 million, otherwise the METIS graph writer switches to scientific notation,
    which confuses KaHiP
    """
    scalingFactor = int((10**6-1)/maxWeight)
    
    
    #copy and scale graph
    Gscaled = G.copyNodes()
    for (u,v) in G.edges():
        Gscaled.addEdge(u,v,int(G.weight(u,v)*scalingFactor))
    
    # write out temporary file
    writeGraph(Gscaled, tempFileName, Format.METIS)
    
    # call KaHIP
    callList = [pathToKaHiP, '--k='+str(k), '--imbalance='+str(imbalance), '--preconfiguration=strong']
    if multiple:
        callList.append('--time_limit=1')
    callList.append(tempFileName)
    subprocess.call(callList)
    
    # read in partition
    part = community.PartitionReader().read(outputFileName)
    
    # remove temporary files
    subprocess.call(['rm', tempFileName])
    subprocess.call(['rm', outputFileName])
    
    return part


# In[13]:

def partitionValid(G, partition, maxBlockSize = 0, isCharged = []):
    z = G.upperNodeIdBound()
    n = G.numberOfNodes()
    if len(partition) != z:
        return False
    
    if len(isCharged) != 0 and len(isCharged) != z:
        return False
    
    if len(isCharged) == 0:
        isCharged = [False for i in range(z)]
        
    if maxBlockSize == 0:
        maxBlockSize = n
        
    chargedFragments = set()
    
    fragmentSizes = {}
    
    for v in range(G.numberOfNodes()):
        if not G.hasNode(v):
            print("Node ", v, " not in graph.")
            return False
            
        # partition invalid if two charged nodes in same fragment
        if isCharged[v]:
            if partition[v] in chargedFragments:
                print("Node", v, " is charged, but fragment", partition[v], "already has a charged node.")
                return False
            else:
                chargedFragments.add(partition[v])
    
        # partition also invalid if gaps of size 1 exist
        if G.hasNode(v+2) and partition[v+2] == partition[v] and G.hasNode(v+1) and partition[v+1] != partition[v]:
            print("Nodes", v, "and", v+2, "are in fragment", partition[v], "but", v+1, "is in fragment", partition[v+1])
            return False
    
        # partition invalid if fragment is larger than allowed
        if not partition[v] in fragmentSizes:
            fragmentSizes[partition[v]] = 1
        else:
            fragmentSizes[partition[v]] += 1
        if fragmentSizes[partition[v]] > maxBlockSize:
            print("Fragment", partition[v], "contains", fragmentSizes[partition[v]], "nodes, more than", maxBlockSize)
            return False
    
    # no reason to complain found, partition is valid
    return True


# In[19]:

def comparePartitionQuality(G, k, imbalance, chargedNodes = set()):
    n = G.numberOfNodes()
    
    isCharged = [v in chargedNodes for v in range(G.numberOfNodes())]
    sizelimit = int(math.ceil(n / k)*(1+imbalance))
    result = {}
    
    ml = mlPartition(G, k, imbalance, isCharged)
    print("MultiLevel:", partitioning.computeEdgeCut(ml, G))
    if not partitionValid(G, ml, sizelimit, isCharged):
        ml = repairPartition(G, ml, sizelimit, isCharged)
        print("Repaired Multilevel:", partitioning.computeEdgeCut(ml, G))
        assert(partitionValid(G, ml, sizelimit, isCharged))

    print()
    result['ml'] = partitioning.computeEdgeCut(ml, G)
    
    greedy = greedyPartition(G, sizelimit, isCharged)
    print("Greedy:", partitioning.computeEdgeCut(greedy, G))
    if not partitionValid(G, greedy, sizelimit, isCharged):
        greedy = repairPartition(G, greedy, sizelimit, isCharged)
        print("Repaired Greedy:", partitioning.computeEdgeCut(greedy, G))
        assert(partitionValid(G, greedy, sizelimit, isCharged))
        
    print()
    result['greedy'] = partitioning.computeEdgeCut(greedy, G)
    
    X = int(n / k)
    tolerance = int(math.ceil(n / k)*(1+imbalance)) - X
    
    try:
        cont = continuousPartition(G, X, tolerance, k, isCharged)
        print("Dynamic Programming:", partitioning.computeEdgeCut(cont, G))
        if not partitionValid(G, cont, sizelimit, isCharged):
            cont = repairPartition(G, cont, sizelimit, isCharged)
            print("Repaired Dynamic:", partitioning.computeEdgeCut(cont, G))
            assert(partitionValid(G, cont, sizelimit, isCharged))
        result['cont'] = partitioning.computeEdgeCut(cont, G)
        
    except ValueError as e:
        print(e)

    print()
    
    
    result['bestOfThree'] = min([result[key] for key in result])
    
    ka = kaHiPWrapper(G, k, imbalance*100)
    print("Raw KaHIP:", partitioning.computeEdgeCut(ka, G))
    if not partitionValid(G, ka, sizelimit, isCharged):
        ka = repairPartition(G, ka, sizelimit, isCharged)
        print("Repaired KaHiP:", partitioning.computeEdgeCut(ka, G))
        assert(partitionValid(G, ka, sizelimit, isCharged))
    print()
    result['ka'] = partitioning.computeEdgeCut(ka, G)

    
    naive = naivePartition(G, X)
    print("Naive:", partitioning.computeEdgeCut(naive, G))
    if not partitionValid(G, naive, sizelimit, isCharged):
        naive = repairPartition(G, naive, sizelimit, isCharged)        
        print("Repaired Naive:", partitioning.computeEdgeCut(naive, G))
        assert(partitionValid(G, naive, sizelimit, isCharged))
    result['naive'] = partitioning.computeEdgeCut(naive, G)
    result['bestOfFive'] = min([result[key] for key in result])

    
    return result


# In[20]:

def main(numberOfChargedNodes = 5):
    pathPrefix = "/home/moritzl/NetworKit/chemfork/NetworKit-chemfork/input/"
    Gnames = ["ubiquitin", "bacteriorhodopsin-10-2.5", "bacteriorhodopsin-10-5", "bacteriorhodopsin-10-10", "bubble"]
    scores = []
    
    for Gname in Gnames:
        G = readGraph(pathPrefix + Gname + ".graph", Format.METIS)

        chargedNodes = set()
        for v in range(numberOfChargedNodes):
            chargedNodes.add(G.randomNode())
            
        print("Graph:", Gname)
        print("chargedNodes =", chargedNodes)
        qualities = comparePartitionQuality(G, 10, 0.1, chargedNodes)
        scores.append(qualities)
        print('------------------------------------------------------------------')
    
    results = {}
    for score in scores:
        for key in score:
            if not key in results:
                results[key] = []
            results[key].append(score[key])
            
    for method in results:
        print(method, ':', str(sum(results[method]) / len(results[method])))


# In[21]:

if __name__ == "__main__":
    # execute only if run as a script
    main(0)


# In[ ]:



