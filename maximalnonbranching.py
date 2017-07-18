def outdegree(g,node):
    """Returns the number of outgoing edges from a node"""
    try:
        return len(g[node])
    except:
        return 0
        
def indegree(g,node):
    """Returns the number of incoming edges for a node"""
    incoming = [k for k,v in g.iteritems() if node in v]
    return len(incoming)
    
def oneInOneOutNode(g, node):
        return indegree(g,node) == outdegree(g,node) == 1
        
def outgoingEdges(g, node):
    outgoing_edges = []
    for neighbour in g[node]:
        outgoing_edges.append((node, neighbour))
    return outgoing_edges
    
def findIsolatedCycles(g, nodes):
    paths = []
    cycle = [nodes[0]]
    while len(nodes) > 0:
        source = cycle[-1]
        cycle.append(g[source].pop())
        if len(g[source]) == 0: del g[source]
        nodes.remove(source)
        if cycle[0] == cycle[-1]:
            paths.append(cycle)
            if len(nodes) > 0:
                cycle = [nodes[0]]
    return paths    
                                
    
def maximalNonBranchingPaths(g):
    """ Input: The adjacency list of a graph whose nodes are integers.
     Output: The collection of all maximal nonbranching paths in this graph."""
    paths = []
    nodes = g.keys()
    print 'all_nodes', nodes
    for node in g:
        if oneInOneOutNode(g, node) is False:
            if outdegree(g, node) > 0:
                nodes.remove(node)
                for (v,w) in outgoingEdges(g, node):
                    nonbranchingPath = [v,w]
                    while oneInOneOutNode(g, w):
                        nonbranchingPath.extend(g[w])
                        nodes.remove(w)
                        w = g[w]
                    paths.append(nonbranchingPath)
    
    print 'paths', paths
    print 'one-in-one-out-nodes', nodes
    if len(nodes) > 0:
        paths.extend(findIsolatedCycles(g, nodes))
    return paths          
    
def read_input(filename):
    f = open(filename, 'r')
    g = dict()
    for line in f:
        line = line.strip()
        line = line.split(' -> ')
        g[int(line[0])] = [int(x) for x in line[1].split(',')]
    return g
        
g = {1:[2], 2:[3], 3:[4], 4:[21], 6:[7], 7:[6], 10:[11], 11:[12], 12:[10], 21:[23], 23:[25], 25:[8,9]}
maximalNonBranchingPaths(g)
           