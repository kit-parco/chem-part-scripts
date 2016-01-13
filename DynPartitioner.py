from networkit import *
import math, sys

def continuousPartition(G, X, tolerance, k=0, isCharged={}):
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
                    
    #for i in range(n):
    #    row = table[i]
    #    print(i, ': [',end="")
    #    for elem in row:
    #        if math.isinf(elem):
    #            print("inf, ", end="")
    #        else:
    #            print(int(elem*100)/100, ", ", end="")
    #    print("]")
                    
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


filename = sys.argv[1]
X = int(sys.argv[2])
tolerance = int(sys.argv[3])
k = int(sys.argv[4]) if len(sys.argv) >= 5 else 0
chargedNodes = []
for i in range(5, len(sys.argv)):
    chargedNodes.append(int(sys.argv[i]))
    
G = readGraph(filename, Format.METIS)
n = G.numberOfNodes()
isCharged = [False for i in range(n)]
for c in chargedNodes:
    assert(c < n)
    isCharged[c] = True

part = continuousPartition(G, X, tolerance, k, isCharged)
for a in chargedNodes:
	for b in chargedNodes:
		if a != b:
			assert(part[a] != part[b])

for i in range(n):
	print(part[i])
#community.PartitionWriter().write(part, filename+'.part')

