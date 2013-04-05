# Jeremy Wang
# Jan 10, 2012
#
# Compute compatible intervals allowing a certain number of 
# characters (SNPs) to be removed when computing incompatibility
#
# accepts a single file of pre-parsed haplotype of the form:
# pos,strn0,strn1,strn2...
# 100,0,1,0
# with haplotypes consisting of 0, 1, [N, H, D, V]

import sys
import math
import time
import matplotlib
matplotlib.use('Agg')
import pylab

# my own local library/script(s)
from mvc import * # minimum vertex cover (for removing fewest SNPs)
from mpdag import * # maximum path through a DAG (for finding maximum cover)

        
# ----------------------------------------------------------------
# Data structures
# ----------------------------------------------------------------

class SNP:
    def __init__(self, position, sdp_index, alleles):
        self.position = position
        self.sdp_index = sdp_index
        self.alleles = alleles
        
    def __repr__(self):
        return str(self.position)
    
    def __cmp__(self, other):
        return cmp(self.position, other.position)

class SDP:
    def __init__(self, p, length):
        self.prob = p # float 0.0 -> 1.0 -- in this case, 0 is minority allele, 1 in majority, and 0.5 is N, there is nothing in between
        self.length = length

    def __repr__(self):
        return self.toarray()
    
    def __str__(self):
        return ''.join(self.toarray())
    
    def toarray(self):
        sdparray = ['%.1f' % p for p in self.prob]
        return sdparray
    
    def __eq__(self, other):
        return self.prob == other.prob
    
    def __len__(self):
        return self.length
    
    def __getitem__(self, key):
        return self.prob[key]
    
    def __iter__(self):
        return (p for p in self.prob)

class Interval:
    def __init__(self, start_index, end_index):
        #the snps at both start and end indices are part of the interval
        self.start_index = start_index
        self.end_index = end_index
        
    def __repr__(self):
        return "(" + str(self.start_index) + "," + str(self.end_index) + ")"
    
    def __cmp__(self, other):
        return self.start_index - other.start_index;
    
    def size(self):
        return self.end_index - self.start_index + 1

class Chromosome:
    def __init__(self):
        self.chrName = ''
        self.strains = []
        self.SNPs = []
        self.SDPs = []
        
        self.left = []
        self.right = []
        self.cores = []
        self.uber = []
        self.maxk = []

    def sortSNPs(self):
        self.SNPs.sort()
    
    def lrScan(self, removable, inv_horizon=[], verbose=False):
        self.left = HorizonAllScan(self, 0, len(self.SNPs)-1, uber=False, horizon=inv_horizon, verbose=verbose)
        
    def rlScan(self, removable, inv_horizon=[], verbose=False):
        scan = HorizonAllScan(self, len(self.SNPs)-1, 0, uber=False, horizon=inv_horizon, verbose=verbose)
        reversed = []
        for i in xrange(len(scan) - 1, -1, -1):
            temp = scan[i].start_index
            scan[i].start_index = scan[i].end_index
            scan[i].end_index = temp
            reversed.append(scan[i])
        self.right = reversed
    
    def uberScan(self, removable, inv_horizon=[], verbose=False):
        print "Uber scan,", removable, "removable"
        if removable > 0:
            self.uber, self.uber_removed = CRAllScan(self, 0, len(self.SNPs)-1, uber=True, removable=removable, verbose=verbose)
        else:
            self.uber = HorizonAllScan(self, 0, len(self.SNPs)-1, uber=True, horizon=inv_horizon, verbose=verbose)
        
    def coreScan(self):
        if len(self.left) == 0 or not len(self.left) == len(self.right):
            print "Can't core scan because left != right or left = []"
            print "left: " + str(len(self.left)) + ", right: " + str(len(self.right))
            print "LR: [" + str(self.left[:10]) + " ... " + str(self.left[-10:]) + "]"
            print "RL: [" + str(self.right[:10]) + " ... " + str(self.right[-10:]) + "]"
            #return
        invs = []
        for i in xrange(min(len(self.left), len(self.right))):
            invs.append(Interval(self.left[i].start_index, self.right[i].end_index))
        self.cores = invs
        if len(self.left) == 0 or not len(self.left) == len(self.right):
            print "Cores: [" + str(self.cores[:10]) + " ... " + str(self.cores[-10:]) + "]"
    
    def maxkScan(self, removable, inv_horizon=[], verbose=False):
        if len(self.uber) == 0:
            t0 = time.time()
            self.uberScan(removable, inv_horizon, verbose=verbose)
            t1 = time.time()
            uber = t1-t0
            if verbose:
                print "  Uber scan: %.2f" % (t1-t0)
        if len(self.cores) == 0:
            if len(self.left) == 0:
                t0 = time.time()
                self.lrScan(removable, inv_horizon, verbose=verbose)
                t1 = time.time()
                lr = t1-t0
                if verbose:
                    print "  LR scan: %.2f" % (t1-t0)
            if len(self.right) == 0:
                t0 = time.time()
                self.rlScan(removable, inv_horizon, verbose=verbose)
                t1 = time.time()
                rl = t1-t0
                if verbose:
                    print "  RL scan: %.2f" % (t1-t0)
            t0 = time.time()
            self.coreScan()
            t1 = time.time()
            if verbose:
                print "  Core scan: %.2f" % (t1-t0)
        t0 = time.time()
        self.maxk = maxkScan(self)
        t1 = time.time()
        if verbose:
            print "  Max-K scan: %.2f" % (t1-t0)
        return lr, rl, uber

    def invFromType(self, type):
        if type == "left":
            return self.left
        if type == "right":
            return self.right
        if type == "cores":
            return self.cores
        if type == "uber":
            return self.uber
        if type == "maxk":
            return self.maxk

class Node:
    def __init__(self, strains, label="", description=""):
        self.strains = strains
        self.label = label
        self.description = description
        self.edges = []

class Edge:
    def __init__(self, direction, weight, end_node):
        self.direction = direction
        self.weight = weight
        self.sink = end_node

# Convert to binary tree with no labelled internal nodes
def printTree(t, strains):
    if len(t.edges) > 0:
        out = [printTree(e.sink, strains) + [e.weight] for e in t.edges]
        if len(t.strains) > 0:
            out.append(['|'.join(strains[s] for s in t.strains), 0]) # zero-weight edge
        # make bifurcating
        while(len(out) > 2): # sub1, sub2
            out = [[out[0], out[1], 0]] + out[2:]
        return out
    else:
        return ['|'.join(strains[s] for s in t.strains)]

# ----------------------------------------------------------------
# MaxK data structures and computation
# ----------------------------------------------------------------

class DPNode:
    def __init__(self, a, b):
        self.next = []
        self.a = a
        self.b = b
        self.back = None
        self.total = 0
        self.back_count = 0
        
    def AddNext(self, node):
        overlap = self.b + 1 - node.a
        if overlap >= 0:
            self.next.append(node)
            newTotal = self.total + overlap
            if newTotal == node.total:
                node.back = self
                self.back_count += 1
            if newTotal > node.total:
                node.total = newTotal
                node.back = self
                self.back_count = 1
        
class DPArray:
    def __init__(self):
        self.current = []
        self.last = None
        self.start = self.current
    
    def NextLevel(self):
        if len(self.current) > 0:
            self.last = self.current
            self.current = []
        
    def AddNode(self, a, b):
        newNode = DPNode(a, b)
        self.current.append(newNode)
        if self.last != None:
            for node in self.last:
                node.AddNext(newNode)
                
    def LastNode(self):
        max = -1
        lastNode = None
        for node in self.current:
            if node.total >= max:
                max = node.total
                lastNode = node
        
        return lastNode
    

def maxkScan(chrom):
    cores = chrom.cores
    uber = chrom.uber
    minDPArray = DPArray()
    last = 0
    for i in xrange(len(cores)):
        minDPArray.NextLevel()
        for j in xrange(last, len(uber)):
            if uber[j].start_index > cores[i].start_index:
                break
            elif (uber[j].end_index >= cores[i].end_index):
                minDPArray.AddNode(uber[j].start_index, uber[j].end_index)
        last = j
    invs = []
    node = minDPArray.LastNode()
    while node != None:
        invs.append(Interval(node.a, node.b))
        node = node.back
    invs.reverse()
    print "in maxkScan: [" + str(invs[:10]) + " ... " + str(invs[-10:]) + "]"
    return invs

# ----------------------------------------------------------------
# LR, RL, and Uber Scan computation
# ----------------------------------------------------------------

# Uber scan is going to be a little different for the CR version --
# Instead of going to the end of the current interval (the first
# incompatibility) and scanning backward, we'll slide the starting
# edge of the interval up just past the first incompatibility, but
# with the interval still containing all the rest.
# If we allow 5 SNPs removed, intervals should overlap up to 5 deep
# I think

# -------------- Character Removal with threshold ----------------
# ignore a certain number of SNPs entirely
#
# removable: number of SNPs we are allowed to ignore/remove
# overlap: number of removed SNPs adjacent intervals should overlap (should be <= removals)
# removable = 0, overlap = 0 should be a normal scan
# horizon: a set of intervals under which everything should be considered compatible

def CRAllScan(chrom, start, end, uber=False, removable=0, verbose=False):
    SNPs = chrom.SNPs
    SDPs = chrom.SDPs
    inc = 1
    if start > end:
        inc = -1
    invs = []
    removed_snps = []
    bounds = [start, start] # interval bounds
    
    sdp_hash = [0 for s in SDPs] # True/False list of seen SDPs
    sdp_list = [] # accumulating list of SDPs in this interval
    incompatibilities = set() # collection of found incompatibilities for each SDP
        
    removed = []
    h = 0 # index into the horizon list
    
    while not bounds[1] == end + inc:
        t0 = time.time()
        
        # slide starting point forward
        # we already have MORE THAN the limit of incompatibilities, so just slide until we get back to the limit
        if uber and bounds[1] > bounds[0]:
            # move starting point forward
            while True:
                # remove bounds[0] SDP
                sdp = SNPs[bounds[0]].sdp_index
                # move bounds[0] up (toward end)
                bounds[0] += inc
                # test to see if we can stop
                sdp_hash[sdp] -= 1
                if sdp_hash[sdp] == 0: # we don't need to change these until we get the last instance of an SDP
                    incompatibilities = set(i for i in incompatibilities if i[0] != sdp and i[1] != sdp) # easy
                    sdp_list.remove(sdp)
                    #print "  removing", sdp
                # we still need to check if this interval is achievable
                edge_groups, vertex_groups = subgraphs(incompatibilities)
                mincovers = [MVC(vertex_groups[i][0], edge_groups[i], sdp_hash, r=removable) for i in xrange(len(edge_groups))] # r=removable bounds each subgroup search to the threshold r
                score = sum([m[0] for m in mincovers]) # scores
                removed = set([a for m in mincovers for a in m[2]]) # SDP indices removed (have to collapse unattached subsets)
                if score <= removable: # gone far enough
                    #print '0', score, removed, [sdp_hash[a] for a in set([a for m in mincovers for a in m[2]])]
                    break
        else:
            # start from right after where the last one ended
            print "bounds are equal"
            bounds[0] = bounds[1]
            incompatibilities = set()
            sdp_hash = [0 for s in SDPs]
            sdp_list = []
        
        while bounds[1] != end + inc:
            sdp = SNPs[bounds[1]].sdp_index
            
            # there is a very tricky subtlety as to when this needs to happen
            # right after this action, bounds[1] = (last compatible SNP) + 2
            bounds[1] += inc
            
            if sdp_hash[sdp] == 0:
                for s in sdp_list:
                    if not compatible(SDPs[sdp], SDPs[s]):
                        incompatibilities.add((s,sdp))
                sdp_list.append(sdp)
                #print "  adding", sdp
            
            sdp_hash[sdp] += 1
            
            # First check to see if old removals will still work
            # if this SDP has already been seen and does not need to be removed, we can ignore it and move on
            if sdp_hash[sdp] == 1 or sdp in removed: # we need to find a new removal set, if possible
            
                # find a set of SDPs we can remove under our bound which eliminates ALL incompatibilities
                # construct a graph where each node represents an SDP with a value equal to the number of corresonding SNPs
                # a node has an edge to another node if they are incompatible
                # to solve, find the lowest total value nodes to remove such that there are no more edges
                # (minimum vertex cover problem)
                
                edge_groups, vertex_groups = subgraphs(incompatibilities)
#                print "----------------------------------------------------------------------------"
                mincovers = [MVC(vertex_groups[i][0], edge_groups[i], sdp_hash, r=removable) for i in xrange(len(edge_groups))]
                score = sum([m[0] for m in mincovers]) # scores
            
                if score > removable:
                    #print '1', score, set([a for m in mincovers for a in m[2]]), [sdp_hash[a] for a in set([a for m in mincovers for a in m[2]])]
                    #print incompatibilities
                    #for v in [v for g in vertex_groups for v in g]:
                    #    print '  ', v, ':', sdp_hash[v]
                    break
                
                removed = set([a for m in mincovers for a in m[2]]) # SDP indices removed (have to collapse unattached subsets)
    
        # convert removed SDP indices into a list of SNP indices
        removed = [i for i in xrange(bounds[0], bounds[1]-inc) if SNPs[i].sdp_index in removed]
        if len(removed) > removable:
            print "We removed too many SNPs -- this shouldn't happen"
    
        # Finish up this interval
        t1 = time.time()
        if verbose:
            print bounds[0], '-', (bounds[1]-inc*2), (t1-t0), '(', ((bounds[1]-inc*2-bounds[0])/(t1-t0)), 'snps/sec)', len(incompatibilities)
        invs.append(Interval(bounds[0], bounds[1] - inc*2))
        removed_snps.append(removed)
    
    return invs, removed_snps

# I swear there is a shorter, more efficient version of this somewhere else, perhaps grab
# that one and modify it for the horizon some time
def HorizonAllScan(chrom, start, end, uber=False, horizon=[], verbose=False):
    SNPs = chrom.SNPs
    SDPs = chrom.SDPs
    inc = 1
    if start > end:
        inc = -1
    invs = []
    i = start
    invStart = start
    scan_back = False
    sdp_hash = {}
    
    if inc < 0:
        # just reverse the horizon order
        horizon = [h for h in reversed(horizon)]
    h = 0 # horizon index
                
    while not i == end + inc:
        if scan_back:
        #scan backward
            while not invStart == start - inc:
                checkPos = i
                
                # seek through horizon until we find one which encompases invStart and i
                under_horizon = False
                while h > 0 and ((inc > 0 and horizon[h].start_index > invStart) or (inc < 0 and horizon[h].end_index < invStart)):
                    h -= 1
                if h < len(horizon) and horizon[h].end_index >= i and horizon[h].end_index >= invStart and horizon[h].start_index <= i and horizon[h].start_index <= invStart:
                    under_horizon = True # will skip the entire checking process later
                    
                string = str(SDPs[SNPs[invStart].sdp_index])
                if not sdp_hash.has_key(string) and not under_horizon: # do not re-check duplicate SDPs, this gives us linear time scan
                    while (not checkPos == invStart) and compatible(SDPs[SNPs[invStart].sdp_index], SDPs[SNPs[checkPos].sdp_index]):
                        checkPos -= inc
                    if not checkPos == invStart:
                        break
                    sdp_hash[string] = 1
                invStart -= inc
            invStart += inc
        if uber:
            scan_back = True
        #scan forward
        i += inc
        while not i == end + inc:
            checkPos = invStart
            
            # seek through horizon until we find one which encompases invStart and i
            under_horizon = False
            while h < len(horizon)-1 and ((inc > 0 and horizon[h].end_index < i) or (inc < 0 and horizon[h].start_index > i)):
                h += 1
            if h < len(horizon) and horizon[h].end_index >= i and horizon[h].end_index >= invStart and horizon[h].start_index <= i and horizon[h].start_index <= invStart:
                under_horizon = True # will skip the entire checking process later
            
            string = str(SDPs[SNPs[i].sdp_index])
            if not sdp_hash.has_key(string) and not under_horizon: # do not re-check duplicate SDPs, this gives us linear time scan:
                while (not checkPos == i) and compatible(SDPs[SNPs[i].sdp_index], SDPs[SNPs[checkPos].sdp_index]):
                    checkPos += inc
                if not checkPos == i:
                    if verbose:
                        print invStart, "-", (i-inc)
                    invs.append(Interval(invStart, i - inc))
                    if uber:
                        invStart = i - inc
                    else:
                        invStart = i
                    sdp_hash = {}
                    break
                sdp_hash[string] = 1
            i += inc
    invs.append(Interval(invStart, end))
    return invs


def compatible(sdp1, sdp2):
    # ---------------------------------------
    # return value:
    # True - compatible or False - incompatible
    # every allele value should be rounded, 0.5 should be ignored
    # ---------------------------------------
    g = [[False, False], [False, False]] # gametes
    for i in xrange(len(sdp1)):
        if sdp1[i] == 0.5 or sdp2[i] == 0.5:
            continue
        g[int(round(sdp1[i]))][int(round(sdp2[i]))] = True
    return not(g[0][0] and g[0][1] and g[1][0] and g[1][1])


def parse(filename):
    data = [line.strip().split(',') for line in open(filename, 'r').read().strip().split('\n')]
    strains = data[0][1:]
    
    sdps = {}
    sdp_count = 0
    snps = []
    length = len(strains)
    lessthantwotons = 0 # noninformative or singletons
    for d in data[1:]:#[:10000]: # skip header
        position = int(d[0])
        snp = d[1:]
        maf = min(snp.count('0'), snp.count('1')) # minor allele frequency
        if maf <= 1:
            lessthantwotons += 1
            continue # this is the way we do it in the original and it helps. that is all. (removes noninformative and singletons)
        alleles = []
        for i in xrange(length):
            letter = snp[i].upper()
            if letter == '1':
                a = 1
            elif letter == 'N' or letter == 'H' or letter == 'V' or letter == 'D':
                a = 0.5
            else:
                a = 0
            alleles.append(a)
        alleles = tuple(alleles)
        sdp_index = -1
        if not sdps.has_key(alleles):
            sdps[alleles] = [sdp_count, 1]
            sdp_index = sdp_count
            sdp_count += 1
        else:
            sdp_index = sdps[alleles][0]
            sdps[alleles][1] += 1
        snps.append(SNP(position, sdp_index, len(sdps.keys())))
    sdps = [SDP(k, len(k)) for k,v in sorted(sdps.iteritems(), key=lambda a:a[1][0])]
    print "Non-informative and singletons: %i" % lessthantwotons
    print "SDPs: ", len(sdps)
    chrom = Chromosome()
    chrom.SNPs = snps
    chrom.SDPs = sdps
    chrom.strains = strains
    print "SNPs: " + str(len(data)-1)
    return chrom

# read true 552 detectable simulated intervals
def readTruth(fname):
    data = [line.strip().split(',') for line in open(fname).read().strip().split('\n')]
    intervals = [(int(d[0]), int(d[1])) for d in data]
    return intervals

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

def main(outfile, fnames, removable):
    for fname in fnames:
        print "Reading %s..." % fname
        chrom = parse(fname)
        # SNPs should already be in ascending order by position
        
        #chrom.maxkScan(removable, verbose=True)
        #invs = chrom.invFromType('maxk')
        if removable > 0:
            chrom.uberScan(removable, verbose=True)
            invs = chrom.invFromType('uber')
            removed_snps = chrom.uber_removed # parallel to the interval list
            
            # check maximum interval overlap
            maxoverlap = 0
            trailing = 0
            count = 0
            for i in xrange(1, len(invs)):
                while invs[trailing].end_index < invs[i].start_index:
                    trailing += 1
                if i - trailing + 1 > maxoverlap:
                    maxoverlap = (i - trailing)
                    count = 0
                elif i - trailing + 1 == maxoverlap:
                    count += 1
                if invs[i].start_index <= invs[i-1].start_index or invs[i].end_index <= invs[i-1].end_index:
                    print "not strictly in order", invs[i-1].start_index, invs[i-1].end_index, '-', invs[i].start_index, invs[i].end_index
            print "Maximum overlap:", maxoverlap, "count", count
            
            score, keepers = maxSingleOverlap(invs)
            print score, keepers
            maxc = [invs[k] for k in keepers]
            maxc_removed = [removed_snps[k] for k in keepers]
            print "MaxC highest cover score:", score
            print len(maxc), "intervals"
            for i in xrange(2, len(maxc)):
                if maxc[i].start_index <= maxc[i-2].end_index:
                    print "extra overlap violation"
                if maxc[i].start_index > maxc[i-1].end_index:
                    print "no overlap violation"
            
            # now allow no removals, but make everything under the "horizon" of these intervals officially compatible
            chrom.uber = [] # reset the Uber scan
            chrom.maxkScan(0, invs, verbose=True)
        else:
            maxc = []
            chrom.maxkScan(0, verbose=True)
            
        print "left:", len(chrom.left)
        #print chrom.left
        print "right:", len(chrom.right)
        #print chrom.right
        print "uber:", len(chrom.uber)
        #print chrom.uber
        print "core:", len(chrom.cores)
        #print chrom.cores
        print "maxk:", len(chrom.maxk)
        #print chrom.maxk
        maxk = chrom.invFromType('maxk')
        
        # compare intervals to truth
        try:
            truth = readTruth('intervals_cliques/sim_maxk.csv')
            inv_total = sum([a[1]-a[0] for a in truth])
            tocompare0 = [(chrom.SNPs[i.start_index].position, chrom.SNPs[i.end_index].position) for i in maxc]
            tocompare1 = [(chrom.SNPs[i.start_index].position, chrom.SNPs[i.end_index].position) for i in maxk]
            score = compareIntervals(tocompare0, truth)
            print " Maxc overlap (SNPs): %i of %i (%.2f%%)" % (score, inv_total, float(score)/inv_total*100)
            score = compareIntervals(tocompare1, truth)
            print " Maxk overlap (SNPs): %i of %i (%.2f%%)" % (score, inv_total, float(score)/inv_total*100)
        except:
            print "No validation intervals."
            
        print "Writing intervals..."
        fout = open(outfile, 'w')
        # start,end,removed_snp,removed_snp, ...
        if "maxc" in outfile:
            fout.write('\n'.join(('%i,%i,' % (chrom.SNPs[maxc[i].start_index].position, chrom.SNPs[maxc[i].end_index].position) + ','.join(str(a) for a in maxc_removed[i])) for i in xrange(len(maxc))))
        elif "maxk" in outfile:
            fout.write('\n'.join(('%i,%i' % (chrom.SNPs[maxk[i].start_index].position, chrom.SNPs[maxk[i].end_index].position)) for i in xrange(len(maxk))))
        else:
            print "Output file does not appear to be MaxC or MaxK, so we aren't writing anything"
        fout.close()

        if removable > 0:
            print "Writing removals..."
            fout = open(outfile+'.rem', 'w')
            fout.write('\n'.join(','.join([str(r) for r in rem]) for rem in removed_snps))
            fout.close()

        
if __name__ == "__main__":
    print "usage: python cr_compatibility.py <outputFile> <removableSNPs> <filename> [<filename>, <filename>, ...]"
    if len(sys.argv) < 4:
        sys.exit(1)
    outfile = sys.argv[1]
    removable = int(sys.argv[2])
    fnames = sys.argv[3:]
    main(outfile, fnames, removable)
