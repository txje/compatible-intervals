
# ideally, find the maximum total 1:1 overlap between the two sets of intervals
# in reality, use only those overlapping regions which are the maximum overlap for both intervals
# I think it is currently O(n) in the number of intervals on one side since intervals may only overlap two others,
# ... so it is at most double overlap on both sides or 2*2 = 4n
def compareIntervals(invs0, invs1):
    j = 0
    score = 0
    for i in xrange(len(invs0)):
        maxoverlap = None
        jmax = None
        while j < len(invs1) and invs1[j][1] < invs0[i][0]:
            j += 1
        l = j
        while l < len(invs1) and invs1[l][0] <= invs0[i][1]:
            overlap = min(invs0[i][1], invs1[l][1]) - max(invs0[i][0], invs1[l][0])
            if maxoverlap == None or overlap > maxoverlap: # we keep track of only the FIRST maximum overlap occurrence - this is always best
                maxoverlap = overlap
                jmax = l
            l += 1
        if jmax != None: # at least one overlapping interval
            k = i + 1
            while k < len(invs0) and invs0[k][0] <= invs1[jmax][1]:
                overlap = min(invs0[k][1], invs1[jmax][1]) - max(invs0[k][0], invs1[jmax][0])
                if overlap > maxoverlap:
                    break
                k += 1
            else: # maximum overlap for both
                score += maxoverlap
                j = jmax + 1 # we have to move ahead so this doesn't get used again
    return score
