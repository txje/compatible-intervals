# ----------------------------------------------------------------
# Jeremy Wang
# June 1, 2010
#
# Compute standard four-gamete compatible intervals including
# LR, LR, Uber, and Max-k using homozygous genotypes
#
# accepts a single file of pre-parsed haplotype of the form:
# pos,strn0,strn1,strn2...
# 100,0,1,0
# with haplotypes consisting of 0, 1, [N, H, D, V]
# -- or -- (binary where the first 4 bytes indicate the position)
# [byte][byte][byte][byte][0|1][0|1][0|1]...
# ----------------------------------------------------------------

import sys

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
    def __init__(self, bin, n, h, v, d, length):
        # each of bin, n, h, v, and d are bit vectors (integers) whose bits represent presence of the particular call
        # each vector stores the call for each strain on this SDP in the order they appear in the Chromosome.strains array from most to least significant bit
        # i.e. strains[0] call is vector >> length-1 (highest bit) and strains[length] call is vector & 1 or vector % 2 (lowest bit)
        self.bin = bin # 1 indicates minority allele (1)
        self.n = n # 1 indicates no call (N)
        self.h = h # 1 indicates heterozygous (H)
        self.v = v # 1 indicates VINO (V)
        self.d = d # 1 indicates deletion (D)
        # from process of elimination, a positions with 0 in all bit vectors is the majority allele (0)
        self.length = length

    def __repr__(self):
        return self.toarray()
    
    def __str__(self):
        return ''.join(self.toarray())
    
    def toarray(self):
        sdparray = []
        for i in xrange(self.length):
            if self.n >> self.length-1-i & 1 == 1:
                sdparray.append('N')
            elif self.h >> self.length-1-i & 1 == 1:
                sdparray.append('H')
            elif self.v >> self.length-1-i & 1 == 1:
                sdparray.append('V')
            elif self.d >> self.length-1-i & 1 == 1:
                sdparray.append('D')
            elif self.bin >> self.length-1-i & 1 == 1:
                sdparray.append('1')
            else:
                sdparray.append('0')
        return sdparray
    
    def __eq__(self, other):
        if self.toarray() == other.toarray():
            return True
        return False
        
    def singleton(self):
        j = 0
        while self.bin > 2**j:
            j += 1
        if self.bin == 2**j:
            return True
        return False
    
    def identical(self, other):
        if self.bin == other.bin and self.n == other.n and self.h == other.h and self.v == other.v and self.d == other.d and self.length == other.length:
            return True
        return False

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
    
    def lrScan(self):
        self.left = allScan(self, 0, len(self.SNPs)-1, uber=False)
        
    def rlScan(self):
        scan = allScan(self, len(self.SNPs)-1, 0, uber=False)
        reversed = []
        for i in xrange(len(scan) - 1, -1, -1):
            temp = scan[i].start_index
            scan[i].start_index = scan[i].end_index
            scan[i].end_index = temp
            reversed.append(scan[i])
        self.right = reversed
    
    def uberScan(self):
        self.uber = allScan(self, 0, len(self.SNPs)-1, uber=True)
        
    def coreScan(self):
        if len(self.left) == 0 or not len(self.left) == len(self.right):
            print "Can't core scan because left != right or left = []"
            print "left: " + str(len(self.left)) + ", right: " + str(len(self.right))
            return
        invs = []
        for i in xrange(len(self.left)):
            invs.append(Interval(self.left[i].start_index, self.right[i].end_index))
        self.cores = invs
    
    def maxkScan(self):
        if len(self.uber) == 0:
            self.uberScan()
        if len(self.cores) == 0:
            if len(self.left) == 0:
                self.lrScan()
            if len(self.right) == 0:
                self.rlScan()
            self.coreScan()
        self.maxk = maxkScan(self)

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
# LR, RL, and Uber Scan computation
# ----------------------------------------------------------------

# allScan performs any of left->right, right->left (if start > end), and uber scan (if uber flag is set)
# 8/2010 - added sdp_hash to prevent checking of duplicate SDPs a la the time proof
def allScan(chrom, start, end, uber=False):
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
    while not i == end + inc:
        if scan_back:
        #scan backward
            while not invStart == start - inc:
                checkPos = i
                string = str(SDPs[SNPs[invStart].sdp_index])
                if not sdp_hash.has_key(string): # do not re-check duplicate SDPs, this gives us linear time scan
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
            string = str(SDPs[SNPs[i].sdp_index])
            if not sdp_hash.has_key(string): # do not re-check duplicate SDPs, this gives us linear time scan
                while (not checkPos == i) and compatible(SDPs[SNPs[i].sdp_index], SDPs[SNPs[checkPos].sdp_index]):
                    checkPos += inc
                if not checkPos == i:
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

# compatible tests 4-gamete compatibility of sdp1 and sdp2 with 2 alleles, n, h, d, and v
def compatible(sdp1, sdp2):
    l = sdp1.length
    mask = (sdp1.n^(2**l-1)) & (sdp2.n^(2**l-1)) & (sdp1.h^(2**l-1)) & (sdp2.h^(2**l-1)) & (sdp1.d^(2**l-1)) & (sdp2.d^(2**l-1)) & (sdp1.v^(2**l-1)) & (sdp2.v^(2**l-1))
    gametes = 0
    if sdp1.bin & sdp2.bin & mask > 0:
        gametes += 1
    if (sdp1.bin^(2**l-1)) & sdp2.bin & mask > 0:
        gametes += 1
    if sdp1.bin & (sdp2.bin^(2**l-1)) & mask > 0:
        gametes += 1
    if (sdp1.bin^(2**l-1)) & (sdp2.bin^(2**l-1)) & mask > 0:
        gametes += 1
    if gametes == 4:
        return False
    return True

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
    return invs

# ----------------------------------------------------------
# Tools
# ----------------------------------------------------------

def printIntervals(invs, outfile, strains, tree=False, sdp=False):
    fout = open(outfile, 'w')
    out = "middle_start,middle_end,start,end,op_start,op_end,sdp_count,haplotypes"
    if tree:
        out += ",tree"
    if sdp:
        out += ",sdp"
    fout.write(out + "\n")
    length = chrom.SDPs[chrom.SNPs[0].sdp_index].length
    for i in xrange(len(invs)):
        inv = invs[i]
        start_index = inv.start_index
        end_index = inv.end_index
        start_pos = chrom.SNPs[start_index].position
        end_pos = chrom.SNPs[end_index].position
        mid_start_pos = chrom.SNPs[(invs[i-1].end_index+1 if i > 0 else start_index)].position
        mid_end_pos = chrom.SNPs[(invs[i+1].start_index-1 if i < len(invs)-1 else end_index)].position
        outside_start = (chrom.SNPs[start_index-1].position + 1) if start_index > 0 else start_pos
        outside_end = (chrom.SNPs[end_index+1].position - 1) if end_index < len(chrom.SNPs)-1 else end_pos
        
        sdp_dict = {}
        for snp in chrom.SNPs[start_index: end_index+1]:
            sdp_dict[snp.sdp_index] = 1
        sdp_count = len(sdp_dict)
        
        # ------- compute number of haplotypes ignoring Ns ----------------
        haps = ['' for j in xrange(length)]
        # construct haplotypes one SDP at a time
        for j in xrange(start_index, end_index+1):
            sdp_str = str(chrom.SDPs[chrom.SNPs[j].sdp_index])
            for k in xrange(len(sdp_str)):
                haps[k] += sdp_str[k]
        haps = list(set(haps)) # remove duplicates
        h = 0
        while h < len(haps):
            h2 = h+1
            while h2 < len(haps):
                for j in xrange(len(haps[h])):
                    if haps[h][j] != 'N' and haps[h2][j] != 'N' and haps[h][j] != haps[h2][j]:
                        h2 += 1
                        break
                else:
                    haps[h] = ''.join([haps[h][j] if haps[h2][j] == 'N' else haps[h2][j] for j in xrange(len(haps[h]))])
                    haps.pop(h2)
            h += 1
        hap_count = len(haps)
        
        # str(start_index) + "," + str(end_index) + "," + 
        fout.write(str(mid_start_pos) + "," + str(mid_end_pos) + "," + str(start_pos) + "," + str(end_pos) + "," + str(outside_start) + "," + str(outside_end) + "," + str(sdp_count) + "," + str(hap_count))
        if tree or sdp:
            if tree:
                fout.write("," + str(printTree(makeTree([chrom.SDPs[chrom.SNPs[j].sdp_index].toarray() for j in xrange(start_index, end_index+1)]), strains)))
            if sdp:
                for bin in sdp_bins.keys():
                    fout.write("," + str(bin))
        fout.write("\n")
    fout.close()

# takes a list of SNPs (not SDPs in that there may be duplicates)
# each SNP should be a list of alleles {0, 1, N, V, D, N}
def makeTree(snps):
    # expand tree using SNPs in descending order of how many 1s they have
    snps.sort(lambda a,b: b.count('1')-a.count('1'))
    root = Node(range(len(snps[0])))
    
    # every SNP with N-ish values adds support to every consistent 0/1 SNP
    nsnps = {}
    s = 0
    while s < len(snps):
        if snps[s].count('0') + snps[s].count('1') < len(snps[s]):
            n = ''.join(snps.pop(s))
            nsnps[n] = nsnps.get(n, 0) + 1
        else:
            snps[s] = [snps[s], 1] # SNP, count
            s += 1
    for n,v in nsnps.iteritems():
        for s in xrange(len(snps)):
            for i in xrange(len(n)):
                if n[i] in ['0', '1'] and n[i] != snps[s][0][i]:
                    break
            else:
                snps[s][1] += 1
    
    for snp in snps:
        indices = [i for i in xrange(len(snp[0])) if snp[0][i] == '1']
        current = root
        # navigate to the correct node
        while not indices[0] in current.strains:
            for edge in current.edges:
                # if this edge exactly matches the given SNP, increment the weight
                if indices == edge.direction:
                    edge.weight += snp[1]
                if indices[0] in edge.direction:
                    current = edge.sink
                    break
        if indices != current.strains:
            # push out a new subtree
            for index in indices:
                current.strains.remove(index)
            newnode = Node([a for a in indices])
            current.edges.append(Edge([a for a in indices], snp[1], newnode))
        elif len(current.edges) > 0:
            # push out a leaf node
            current.strains = []
            newnode = Node([a for a in indices])
            current.edges.append(Edge([a for a in indices], snp[1], newnode))
    return root

def parse(filename):
    data = [line.strip().split(',') for line in open(filename, 'r').read().strip().split('\n')]
    strains = data[0][1:]
    
    sdps = {}
    sdp_count = 0
    snps = []
    length = len(strains)
    for d in data[1:]: # skip header
        position = int(d[0])
        snp = d[1:]
        #SNPList = subset(SNPList, indices)
        alleles = []
        bin = 0
        n = 0
        h = 0
        v = 0
        d = 0
        ns = 0
        for i in xrange(length):
            letter = snp[i].upper()
            if letter == '1':
                bin += 2**(length-1-i)
            elif letter == 'N' or letter == 'H' or letter == 'V' or letter == 'D':
                if letter == 'N':# or letter == 'H' or letter == 'D':
                    n += 2**(length-1-i)
                if letter == 'H':
                    h += 2**(length-1-i)
                if letter == 'V':
                    v += 2**(length-1-i)
                if letter == 'D':
                    d += 2**(length-1-i)
        sdpkey = bin*(2**(length*4)) + n*(2**(length*3)) + h*(2**(length*2)) + v*(2**(length)) + d
        sdp_index = -1
        if not sdps.has_key(sdpkey):
            sdps[sdpkey] = [sdp_count, 1]
            sdp_index = sdp_count
            sdp_count += 1
        else:
            sdp_index = sdps[sdpkey][0]
            sdps[sdpkey][1] += 1
        snps.append(SNP(position, sdp_index, alleles))
    sdps = [SDP(k>>length*4, (k>>length*3)%(2**length), (k>>length*2)%(2**length), (k>>length)%(2**length), k%(2**length), length) for k,v in sorted(sdps.iteritems(), key=lambda a:a[1][0])]
    chrom = Chromosome()
    chrom.SNPs = snps
    chrom.SDPs = sdps
    chrom.strains = strains
    print "SNPs: " + str(len(data)-1)
    return chrom

if __name__ == "__main__":
    print "Usage: python compatinv.py <input> <output> <interval type>"
    
    infile = sys.argv[1]
    output = sys.argv[2]
    inv_type = sys.argv[3]
    
    print "Reading file..."
    chrom = parse(infile)
    
    print "Processing..."
    chrom.maxkScan()
    invs = chrom.invFromType(inv_type)
    
    print "Writing output..."
    printIntervals(invs, output, chrom.strains, tree=True)
