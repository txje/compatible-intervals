# Jeremy Wang
# Jan 24, 2013
#
# Draw compatibility matrices with intervals
# - forked from the much older version ~/projects/compatible_intervals/drawMatrix.py
#

import sys
import Image, ImageDraw, ImageFont

ALLELES = ['A', 'C', 'G', 'T']

colors = {'gray':(180,180,180),
          'darkgray':(100,100,100),
          'pink':(255,160,160),
          'green':(0,255,0),
          'yellow':(180,180,0)
         }

invcolors = [colors['green'], colors['yellow'], colors['pink']] #colors['darkgray'], 
         
def box(draw, x, y, x2, y2, color):
    draw.line([(x,y),(x2,y)], fill=color)
    draw.line([(x2,y),(x2,y2)], fill=color)
    draw.line([(x2,y2),(x,y2)], fill=color)
    draw.line([(x,y2),(x,y)], fill=color)

def compMatrix(imfile, snps, invsets, blocksize=4):
    n = len(snps)
    m = len(snps[0])
    img = Image.new("RGBA", (blocksize*2*n, blocksize*2*m + blocksize*n), (255,255,255, 255))
    draw = ImageDraw.Draw(img)
    height = blocksize*m + blocksize*n
    ytop = height - blocksize*m
    
    #haplotypes/SNPs
    for i in xrange(n):
        snp = snps[i]
        for j in xrange(len(snp)):
            bit = snp[j]

            x = blocksize*2*i
            y = blocksize*2*j + blocksize*n
            if (bit == '0'):
                fillColor=(0,0,255,255)
            elif (bit == '1'):
                fillColor=(255,255,0,255)
            elif bit == '2':
                fillColor=(0,255,0,255)
            else: # probably a '3' for N
                fillColor=(255,255,255,255)
            draw.rectangle((x,y,x+int(blocksize*2*0.8),y+int(blocksize*2*0.8)), fill=fillColor)
    #compatibility grid/triangle
    for i in xrange(0,n-1):
        for j in xrange(i+1, n):
            x = (j-i)*blocksize + int(blocksize*0.8) + blocksize*2*i
            y = blocksize*n - (j-i)*blocksize
            if compatible(snps[i],snps[j]):
                fcolor = (200,200,200,255) #gray
            else:
                fcolor = (255,0,0,255) #red
            #apply shading for intervals
            invs_covering = []
            for k in xrange(len(invsets)):
                invs = invsets[k]
                for inv in invs:
                    if i >= inv[0] and j <= inv[1]:
                        invs_covering.append(k)
            #shade = 0.25 + 0.75/(2**inv_cover)
            # use only the last tint seen
            invs_covering = invs_covering[-1:]
            # average colors with tints
            fcolor = ((fcolor[0] + sum(invcolors[k][0] for k in invs_covering)) / (len(invs_covering)+1),
                      (fcolor[1] + sum(invcolors[k][1] for k in invs_covering)) / (len(invs_covering)+1),
                      (fcolor[2] + sum(invcolors[k][2] for k in invs_covering)) / (len(invs_covering)+1),
                      fcolor[3])
            
            draw.polygon([(x-int(blocksize*0.8), y), (x, y-int(blocksize*0.8)), (x+int(blocksize*0.8), y), (x, y+int(blocksize*0.8))], fill=fcolor) # diamond
    
    #draw interval boundaries
    for i in xrange(len(invsets)):
        invs = invsets[i]
        for inv in invs:
            fcolor = (0,0,0,255) #black
            inv_height = ((inv[1]+1)-inv[0])*blocksize
            draw.line([(inv[0]*blocksize*2,ytop),(inv[0]*blocksize*2+inv_height,ytop-inv_height)], fill=fcolor, width=int(blocksize*0.4))
            draw.line([((inv[1]+1)*blocksize*2,ytop),(inv[0]*blocksize*2+inv_height,ytop-inv_height)], fill=fcolor, width=int(blocksize*0.4))
    
    img.save(imfile, "PNG")


# draws a compatibility matrix where each SNP is only one pixel
# looks like this, where each character is a pixel
'''
c
ci
icc
ccci
iccic
000000
100110
011101
'''
def condensedMatrix(imfile, snps, invsets):
    n = len(snps)
    m = len(snps[0])
    img = Image.new("RGBA", (n, n+m), (255,255,255, 255))
    height = n+m
    ytop = height - m
    
    #haplotypes/SNPs
    for i in xrange(n):
        snp = snps[i]
        for j in xrange(len(snp)):
            bit = snp[j]

            x = i
            y = height - 1 - j
            if (bit == '0'):
                fillColor=(0,0,255,255)
            elif (bit == '1'):
                fillColor=(255,255,0,255)
            elif bit == '2':
                fillColor=(0,255,0,255)
            else: # probably a '3' for N
                fillColor=(255,255,255,255)
            img.putpixel((x,y), fillColor)
    #compatibility grid/triangle
    for i in xrange(0,n-1):
        for j in xrange(i+1, n):
            x = i
            y = ytop - j + i
            if compatible(snps[i],snps[j]):
                fcolor = (200,200,200,255) #gray
            else:
                fcolor = (255,0,0,255) #red
            #apply shading for intervals
            invs_covering = []
            for k in xrange(len(invsets)):
                invs = invsets[k]
                for inv in invs:
                    if i >= inv[0] and j <= inv[1]:
                        invs_covering.append(k)
            #shade = 0.25 + 0.75/(2**inv_cover)
            # average colors with tints
            fcolor = ((fcolor[0] + sum(invcolors[k][0] for k in invs_covering)) / (len(invs_covering)+1),
                      (fcolor[1] + sum(invcolors[k][1] for k in invs_covering)) / (len(invs_covering)+1),
                      (fcolor[2] + sum(invcolors[k][2] for k in invs_covering)) / (len(invs_covering)+1),
                      fcolor[3])
            
            img.putpixel((x,y), fcolor)
    
    # we can't draw interval boundaries    
    img.save(imfile, "PNG")
    
def compatible(sdp1, sdp2):
    gametes = set()
    for i in range(len(sdp1)):
        if sdp1[i] not in ['0', '1'] or sdp2[i] not in ['0', '1']:
            continue
        gametes.add(sdp1[i] + sdp2[i])
    if len(gametes) == 4:
        return False
    return True

def readSNPs(fname):
    data = [line.strip().split(',') for line in open(fname, 'r').read().strip().split('\n')]
    strains = data[0][1:]
    snps = [[int(d[0]), ''.join(d[1:])] for d in data[1:]]
    return snps

def readIntervals(fname):
    data = [line.strip().split(',') for line in open(fname, 'r').read().strip().split('\n')]
    invs = [(int(d[0]), int(d[1])) for d in data]
    return invs
    
if __name__ == "__main__":
    if len(sys.argv) < 6:
        print "usage: python drawMatrix.py <SNPFile> <startPosition> <endPosition> <outputImageFile> <intervalFile> [<intervalFile> ...]"
    snpfile = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    imagefile = sys.argv[4]
    invfiles = sys.argv[5:]
    
    REDUCE = True # flag to remove consecutive identical SNPs, other reductions
    
    # for some reason, we get some intervals with 1 as an end position, ignore these
    invsets = [[a for a in readIntervals(ifile) if ((a[0] >= start and a[0] < end) or (a[1] <= end and a[1] > start)) and a[1] != 1] for ifile in invfiles] # (start,end), (start,end), ...
    snps = [s for s in readSNPs(snpfile) if s[0] >= start and s[0] <= end] # (pos, alleles), (pos, alleles), ...
    # normalize SNPs to {'0', '1', '2', and '3'} for majority, minority, het, and N, respectively
    s = 0
    while s < len(snps):
        alleles = {}
        for a in snps[s][1]:
            alleles[a] = alleles.get(a, 0) + 1
        alleles = sorted([(k,v) for k,v in alleles.iteritems()], key = lambda a: a[1]) # sort by allele frequency
        val = 0
        for a in alleles:
            if a in ALLELES:
                snps[s][1].replace(a, str(val))
                val += 1
                if val > 1:
                    val = 3 # 3rd and more alleles and Ns become 3
        snps[s][1] = snps[s][1].replace('H', '2')
        snps[s][1] = snps[s][1].replace('N', '3')
        # remove consecutive identical SNPs and singletons
        if REDUCE and snps[s][1].count('1') <= 1 or snps[s][1].count('0') <= 1:
            snps.pop(s)
        elif REDUCE and s > 0 and snps[s][1] == snps[s-1][1]:
            snps.pop(s)
        else:
            s += 1
    # convert intervals from bp to SNP indices
    for i in xrange(len(invsets)):
        for j in xrange(len(invsets[i])):
            start = None
            end = len(snps)-1
            for s in xrange(len(snps)):
                if start == None and snps[s][0] >= invsets[i][j][0]:
                    start = s
                if snps[s][0] <= invsets[i][j][1]:
                    end = s
                else:
                    break
            invsets[i][j] = [start, end]
    
    # keep only allele information
    snps = [s[1] for s in snps]
    
    # build compatibility matrix
    matrix = [[compatible(snps[s], snps[t]) for s in xrange(len(snps))] for t in xrange(len(snps))]
    
    # remove SNPs which are not incompatible with at least one SNP in its interval or a neighboring interval
    # -------------------------------------------------------------------------------------------------------
    '''
    tokeep = set()
    s = 0
    while s < len(snps):
        incompatible = []
        for invs in invsets:
            for i in xrange(len(invs)):
                if invs[i][0] <= s and invs[i][1] >= s: # in this interval
                    # SNP bounds in which we want to check if there are incompatibilities
                    start = invs[i][0] - (1 if invs[i][0] > 0 else 0)
                    end = invs[i+1 if i < len(invs)-1 and invs[i+1][0] < s else i][1]
                    if end < len(matrix) - 1:
                        end += 1
                    for a in xrange(start, end+1):
                        if not matrix[s][a]: # incompatible
                            incompatible.append(a)
                    break
        if len(incompatible) > 0:
            tokeep.add(s)
            for i in incompatible:
                tokeep.add(i)
        s += 1
    
    toremove = set(s for s in xrange(len(snps))) - tokeep
    
    for r in sorted(list(toremove), key = lambda a: a*-1):
        snps.pop(r)
        # ... AND fix interval bounds
        for i in xrange(len(invsets)):
            for j in xrange(len(invsets[i])):
                if invsets[i][j][0] >= r:
                    invsets[i][j][0] -= 1
                    invsets[i][j][1] -= 1
                elif invsets[i][j][1] >= r:
                    invsets[i][j][1] -= 1
    '''
    # -------------------------------------------------------------------------------------------------------
    
    print "%i SNPs" % len(snps)
    print "%i Intervals" % sum(len(i) for i in invsets)
    
    print invsets
    
    #condensedMatrix(imagefile, snps, invsets)
    compMatrix(imagefile, snps, invsets, blocksize=5)
