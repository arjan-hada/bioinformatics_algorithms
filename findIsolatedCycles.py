def findIsolatedCycles(g, nodes):
    """Find isolated cycles given a list of one-in-one-out nodes"""
    paths = []
    for node in nodes:
        isolatedCycle = []
        def visit(n):
            while len(g[n]) > 0:
                nodes.remove(n)
                dst = g[n].pop()
                visit(dst)
            isolatedCycle.append(n)
        visit(node)
        isolatedCycle = isolatedCycle[::-1] # reverse and then take all but last node
        paths.append(isolatedCycle)
    return paths  

#g = {1:[2], 2:[3], 3:[4,5], 6:[7], 7:[8], 8:[6]}
g = {1:[2], 2:[3], 3:[4,5], 6:[7], 7:[6]}
#g = {1:[2], 2:[3], 3:[4,5], 6:[6]}
nodes = [6,7]
print findIsolatedCycles(g, nodes)