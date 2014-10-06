# Jeremy Wang
# Jan 9, 2012
#
# Minimum Vertex Cover, with vertex weights
# -- this does an exhaustive recursive search --
#
# node:      a vertex name (v) to start at
# weights:   a dictionary of {node:weight, node:weight ...}
# E:         a list of edges [(v, v), (v, v) ...]
# rootseen:  the vertices we have already seen and do not need to revisit
#
# This is tricky-ass code. Do not take it for granted.
#

import time # for testing

# computes MVC for the connected subgraph starting at node
# aggressive BRANCH AND BOUND
# returns (cost, seen, removed, done) where done can be a value 0 (impossible), 1 (possible), or 2 (achieved/done)
def MVC(node, E, weights, rootseen=[], tabs=0, costSoFar=0, r=float("inf")):

#    print ''.join(' ' for a in xrange(tabs)), node, "seen", rootseen

    # keep this node
    keepcost = float("inf")
    keep_nodes_removed = []
    keep_seen = [s for s in rootseen] + [node]
    for e in E: # check for edges to seen nodes
        if (e[0] == node and e[1] in rootseen) or (e[1] == node and e[0] in rootseen):
            break
    else: # we can keep this one
#        print ''.join(' ' for a in xrange(tabs)), "keeping", node, "seen", keep_seen, "checking", [e for e in E if e[0] == node or e[1] == node]
        keepcost = 0
        for e in list(E): # the order of the edge list may change during the loop, so we'll operate on a copy
            if e[1] not in keep_seen:
                morecost, keep_seen, moreremoved, done = MVC(e[1], E, weights, keep_seen, tabs+1, costSoFar+keepcost, r)
            elif e[0] not in keep_seen:
                morecost, keep_seen, moreremoved, done = MVC(e[0], E, weights, keep_seen, tabs+1, costSoFar+keepcost, r)
            else:
                continue
            keep_nodes_removed.extend(moreremoved)
            keepcost += morecost
            # BOUND
            if done == 2:
                return keepcost, keep_seen, keep_nodes_removed, 2 # done!
            if done == 0: # impossible (this can't be kept)
                keepcost = float("inf")
                break

    # remove this node
    remove_seen = [s for s in rootseen] + [node]
#    print ''.join(' ' for a in xrange(tabs)), "removing", node, "costs", weights[node], "seen", remove_seen
    removecost = weights[node]
    removed = []
    remove_nodes_removed = [node]
    i = 0
    while i < len(E): # remove edges
        if E[i][0] == node or E[i][1] == node:
            removed.append(E.pop(i))
        else:
            i += 1

    if len(E) == 0: # no more edges, great!
        # BOUND
        if costSoFar+removecost <= r: # under our threshold, great!
            return removecost, remove_seen, remove_nodes_removed, 2 # 2 means done!

    elif costSoFar+removecost < r: # only need to check this if it's going to be possible, we have to be able to remove at least one more SNP
#        print ''.join(' ' for a in xrange(tabs)), "removed to check", removed
        for e in removed: # recurse through all edges (although hereafter removed)
            if e[1] not in remove_seen: # node is also in seen
                morecost, remove_seen, moreremoved, done = MVC(e[1], E, weights, remove_seen, tabs+1, costSoFar+removecost, r)
            elif e[0] not in remove_seen:
                morecost, remove_seen, moreremoved, done = MVC(e[0], E, weights, remove_seen, tabs+1, costSoFar+removecost, r)
            else:
                continue
            remove_nodes_removed.extend(moreremoved)
            removecost += morecost
            if done == 2: # BOUND
                return removecost, remove_seen, remove_nodes_removed, 2 # done!
            if done == 0: # impossible
                removecost = float("inf")
                break

    else: # there are edges and the total cost is >= r, so we can't remove this one
        removecost = float("inf")

    E.extend(removed) # add edges back

#    print ''.join(' ' for a in xrange(tabs)), seen

    if min(costSoFar+keepcost, costSoFar+removecost) > r:
        status = 0 # impossible
    else:
        status = 1 # possible, but not yet definitive

    if keepcost <= removecost: # err on the side of keeping
#        print ''.join(' ' for a in xrange(tabs)), "keep", node, "(", weights[node], ") cost", keepcost
        return keepcost, keep_seen, keep_nodes_removed, status # cost of nodes removed, nodes seen under this thread (don't do them again)
    else:
#        print ''.join(' ' for a in xrange(tabs)), "remove", node, "(", weights[node], ") cost", removecost
        return removecost, remove_seen, remove_nodes_removed, status # cost of nodes removed, nodes seen under this thread (don't do them again)

# Utility:
# find all disconnected subgraphs
def subgraphs(E):
    groups = []
    group_vertices = []
    for e in E:
        first = None
        for g in xrange(len(groups)):
            if e[0] in group_vertices[g] or e[1] in group_vertices[g]:
                if first == None:
                    groups[g].append(e)
                    group_vertices[g].add(e[0])
                    group_vertices[g].add(e[1])
                    first = g
                else: # second and final
                    groups[first].extend(groups.pop(g))
                    group_vertices[first].update(group_vertices.pop(g))
                    break
        if first == None:
            groups.append([e])
            group_vertices.append(set([e[0], e[1]]))
    return groups, [list(g) for g in group_vertices]


# Test it with a couple graphs with nonobvious solutions and one with a cycle
if __name__ == "__main__":
    V = [0, 1, 2, 3, 4, 5, 6, 7] # not used

    E = [(0, 1), (1, 2),
         (3, 4), (3, 5), (4, 5), (5, 6)]
    print "All edges:", E

    weights = {0:1, 1:5, 2:1, 3:1, 4:2, 5:1, 6:10, 7:1}
    print "Weights:", weights

    # should split into 2 groups:
    # {0, 1, 2} {3, 4, 5, 6} (7 does not appear)
    # then remove {0, 2} at a cost of 2, and {3, 5} at a cost of 2

    print "Breaking into groups via subgraphs() ..."
    t0 = time.time()
    edge_groups, vertex_groups = subgraphs(E)
    t1 = time.time()
    print "Runtime: %.2f ms" % ((t1-t0) * 1000)

    score = 0
    for i in xrange(len(edge_groups)):
        print "\nGroup", i
        print "Vertices:", vertex_groups[i], "Edges:", edge_groups[i]
        t2 = time.time()
        result = MVC(vertex_groups[i][0], edge_groups[i], weights)
        t3 = time.time()
        print "MVC nodes to remove:", result[2], "(score: %i)" % result[0]
        print "Runtime: %.2f ms" % ((t3-t2) * 1000)
        score += result[0]

    print "\nFull graph MVC score: %i" % score
