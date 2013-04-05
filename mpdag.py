# Jeremy Wang
# Jan 15, 2012
#
# Maximum path through a directed acyclic graph
# - with the extra constraint that no more than 2 intervals may overlap at a point
#
#  compute the cover of single-overlap intervals with maximum total interval size
#  dynamic programming (compute distal to proximal intervals)

def maxSingleOverlap(invs):
    edges = [[] for i in invs] # edges of the form (score, sink[list], sink_start_pos)
    # IN THE FUTURE we don't need to keep the entire edge history, once you're two intervals removed, we can get rid of them
    edges[-1].append((invs[-1].end_index - invs[-1].start_index + 1, [], float("inf"))) # DP base case (last interval)
    for i in xrange(len(invs)-2, -1, -1):
        # find subsequent overlapping intervals
        j = i+1
        while j < len(invs) and invs[j].start_index <= invs[i].end_index:
            # look through edges leaving the overlapping interval
            # making sure we don't have more than 2x overlap
            viable_edges = [(e[0], e[1]) for e in edges[j] if invs[i].end_index < e[2]]
            if len(viable_edges) == 0:
                # no valid cover using this edge
                score = 0
                sink = [None]
            else:
                score, sink = max(viable_edges, key=lambda a: a[0])
            edges[i].append((score+(invs[i].end_index - invs[i].start_index + 1), [j]+sink, invs[j].start_index))
            j += 1
    # edges[0] is all edges out of 0
    # we get the absolute best if we test edges[0] and edges[n] for all n overlapping 0
    best = None
    i = 0
    while i < len(invs) and invs[i].start_index <= invs[0].end_index:
        for e in edges[i]:
            if best == None or e[0] > best[0]:
                best = (e[0], [i]+e[1]) # score, inv list
        i += 1
    return best # (score, inv list)

class Interval:
    def __init__(self, s, e):
        self.start_index = s
        self.end_index = e

if __name__ == "__main__":
    invs = [Interval(0, 5), Interval(2, 20), Interval(5, 21), Interval(18, 25), Interval(20, 30), Interval(22, 35)]
    print maxSingleOverlap(invs)
