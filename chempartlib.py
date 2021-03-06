from networkit import *
import math, sys

def continuousPartition(G, X, tolerance, k=0, isCharged={}):
    """
    Partition G into subsets of size X +/- tolerance and with consecutive node ids. Charged nodes are not grouped into the same subset.
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
    
    def w(i):
        """
        Weight of cutting after node i.
        TODO: consider other weights beside the single cut edge
        """
        assert(i >= 0)
        assert(i < n)
        if (i == n-1):
            return 0
        else:
            return G.weight(i,i+1)
    
    # allocate cut and predecessor table
    table = [[float("inf") for j in range(maxNumberOfPartitions)] for i in range(n)]
    pred = [[-1 for j in range(maxNumberOfPartitions)] for i in range(n)]
    
    # fill table 
    for i in range(n):
        for j in range(maxNumberOfPartitions):
            if (j == 0):
                if abs(i+1-X) <= tolerance:
                    table[i][j] = w(i)
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
                
                predList = [table[l][j-1] for l in range(windowStart, min(i,i-(X-tolerance)+1))]
                if (len(predList) > 0):
                    minPred = min(predList)
                    table[i][j] = minPred + w(i)
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
        c -=w(i)
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


def mlPartition(G, k, imbalance, isCharged={}):
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


def greedyPartition(G, sizelimit, isCharged={}):
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
        for i in range(n-1):
            if part[i] in partSet and part[i+2] in partSet and not part[i+1] in partSet:
                allowed = False #otherwise, would create single embedded node
        if allowed:
            part.mergeSubsets(firstPart, secondPart)
            if firstPart in chargedPartitions or secondPart in chargedPartitions:
                chargedPartitions.add(part.subsetOf(heaviestEdge[0]))
    
    part.compact()
    return part

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
