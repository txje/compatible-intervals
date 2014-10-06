# Jeremy Wang
# Jan 24, 2013
#
# Draw compatibility matrices with intervals
#

import sys
import argparse
import Image, ImageDraw, ImageFont

ALLELES = ['A', 'C', 'G', 'T']

colors = {'gray':(180,180,180),
          'darkgray':(100,100,100),
          'pink':(255,160,160),
          'green':(0,255,0),
          'yellow':(180,180,0)
         }

invcolors = [(50,50,50), colors['green'], colors['pink']]

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
      #invs_covering = invs_covering[-1:]
      # use only the first
      invs_covering = [invs_covering[0]] if len(invs_covering) > 0 else []
      # average colors with tints
      fcolor = ((fcolor[0] + sum(invcolors[k][0] for k in invs_covering)) / (len(invs_covering)+1),
                (fcolor[1] + sum(invcolors[k][1] for k in invs_covering)) / (len(invs_covering)+1),
                (fcolor[2] + sum(invcolors[k][2] for k in invs_covering)) / (len(invs_covering)+1),
                fcolor[3])

      draw.polygon([(x-int(blocksize*0.8), y), (x, y-int(blocksize*0.8)), (x+int(blocksize*0.8), y), (x, y+int(blocksize*0.8))], fill=fcolor) # diamond

  #draw interval boundaries
  print " ----------- ALL INTERVALS ARE DARK GRAY ------------"
  print " ----------- NOT SHOWING BOUNDARY LINES FOR FIRST SET OF INTERVALS ------------"
  for i in xrange(1, len(invsets)):
    invs = invsets[i]
    for inv in invs:
      fcolor = invcolors[i] #(0,0,0,255)
      inv_height = ((inv[1]+1)-inv[0])*blocksize
      draw.line([(inv[0]*blocksize*2,ytop),(inv[0]*blocksize*2+inv_height,ytop-inv_height)], fill=fcolor, width=int(blocksize*0.2))
      draw.line([((inv[1]+1)*blocksize*2,ytop),(inv[0]*blocksize*2+inv_height,ytop-inv_height)], fill=fcolor, width=int(blocksize*0.2))

  img.save(imfile, "JPEG")


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

