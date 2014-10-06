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
#
# Had to stop using an SDP hash in HorizonScan to avoid duplicate checking
# a la the proof because it didn't cooperate with horizon checking

import sys
import math

# my own local library/script(s)
from mvc import * # minimum vertex cover (for removing fewest SNPs)
from mpdag import * # maximum path through a DAG (for finding maximum cover)


def lrScan(snps, removable, inv_horizon=[], verbose=False):
  return HorizonAllScan(snps, 0, len(snps)-1, uber=False, horizon=inv_horizon, verbose=verbose)

def rlScan(snps, removable, inv_horizon=[], verbose=False):
  scan = HorizonAllScan(snps, len(snps)-1, 0, uber=False, horizon=inv_horizon, verbose=verbose)
  reversed = [[s[1], s[0]] for s in scan[::-1]]
  return reversed

def uberScan(snps, removable, inv_horizon=[], verbose=False):
  return HorizonAllScan(snps, 0, len(snps)-1, uber=True, horizon=inv_horizon, verbose=verbose)

def CRUberScan(snps, removable, inv_horizon=[], verbose=False):
  return CRAllScan(snps, 0, len(snps)-1, uber=True, removable=removable, verbose=verbose)


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

def CRAllScan(snps, start, end, uber=False, removable=0, verbose=False):
  inc = 1
  if start > end:
      inc = -1
  invs = []
  removed_snps = []
  bounds = [start, start] # interval bounds

  sdp_hash = {''.join(s):0 for s in snps}
  sdp_list = [] # accumulating list of SDPs in this interval
  incompatibilities = set() # collection of found incompatibilities for each SDP

  removed = []
  h = 0 # index into the horizon list

  while not bounds[1] == end + inc:

    # slide starting point forward
    # we already have MORE THAN the limit of incompatibilities, so just slide until we get back to the limit
    if uber and bounds[1] > bounds[0]:
      # move starting point forward
      while True:
        # remove bounds[0] SDP
        sdp = ''.join(snps[bounds[0]])
        # move bounds[0] up (toward end)
        bounds[0] += inc
        # test to see if we can stop
        sdp_hash[sdp] -= 1
        if sdp_hash[sdp] == 0: # we don't need to change these until we get the last instance of an SDP
          incompatibilities = set(i for i in incompatibilities if i[0] != sdp and i[1] != sdp) # easy
          sdp_list.remove(sdp)
          if verbose:
            print "  removing", sdp
        # we still need to check if this interval is achievable
        edge_groups, vertex_groups = subgraphs(incompatibilities)
        mincovers = [MVC(vertex_groups[i][0], edge_groups[i], sdp_hash, r=removable) for i in xrange(len(edge_groups))] # r=removable bounds each subgroup search to the threshold r
        score = sum([m[0] for m in mincovers]) # scores
        removed = set([a for m in mincovers for a in m[2]]) # SDP indices removed (have to collapse unattached subsets)
        if score <= removable: # gone far enough
          break
    else:
      # start from right after where the last one ended
      if verbose:
        print "bounds are equal"
      bounds[0] = bounds[1]
      incompatibilities = set()
      sdp_hash = {''.join(s):0 for s in snps}
      sdp_list = []

    while bounds[1] != end + inc:
      sdp = ''.join(snps[bounds[1]])

      # there is a very tricky subtlety as to when this needs to happen
      # right after this action, bounds[1] = (last compatible SNP) + 2
      bounds[1] += inc

      if sdp_hash[sdp] == 0:
        for s in sdp_list:
          if not compatible(sdp, s):
            incompatibilities.add((s,sdp))
        sdp_list.append(sdp)
        if verbose:
          print "  adding", sdp

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
        if verbose:
          print "----------------------------------------------------------------------------"
        mincovers = [MVC(vertex_groups[i][0], edge_groups[i], sdp_hash, r=removable) for i in xrange(len(edge_groups))]
        score = sum([m[0] for m in mincovers]) # scores

        if score > removable:
          if verbose:
            print '1', score, set([a for m in mincovers for a in m[2]]), [sdp_hash[a] for a in set([a for m in mincovers for a in m[2]])]
            print incompatibilities
            for v in [v for g in vertex_groups for v in g]:
              print '  ', v, ':', sdp_hash[v]
          break

        removed = set([a for m in mincovers for a in m[2]]) # SDP indices removed (have to collapse unattached subsets)

    # convert removed SDP indices into a list of SNP indices
    removed = [i for i in xrange(bounds[0], bounds[1]-inc) if ''.join(snps[i]) in removed]
    if len(removed) > removable:
      print "We removed too many SNPs -- this shouldn't happen"

    # Finish up this interval
    if verbose:
      print "Finished interval SNPs %i - %i, %i incompatibilities" % (bounds[0], (bounds[1]-inc*2), len(incompatibilities))
    invs.append((bounds[0], bounds[1] - inc*2))
    removed_snps.append(removed)

  return invs, removed_snps

# I swear there is a shorter, more efficient version of this somewhere else, perhaps grab
# that one and modify it for the horizon some time
def HorizonAllScan(snps, start, end, uber=False, horizon=[], verbose=False):
  inc = 1
  if start > end:
    inc = -1
  invs = []
  i = start
  invStart = start
  scan_back = False

  if inc < 0:
    # just reverse the horizon order
    horizon = [h for h in reversed(horizon)]
  h = 0 # horizon index

  while not i == end + inc:
    if scan_back:
    #scan backward
      while not invStart == start - inc:

        string = ''.join(snps[invStart])

        # seek through horizon until we find one which encompases invStart
        while h > 0 and ((inc > 0 and horizon[h][0] > invStart) or (inc < 0 and horizon[h][1] < invStart)):
          h -= 1

        checkPos = i

        while (not checkPos == invStart) and ((h < len(horizon) and horizon[h][1] >= checkPos and horizon[h][1] >= invStart and horizon[h][0] <= checkPos and horizon[h][0] <= invStart) or compatible(snps[invStart], snps[checkPos])):
          checkPos -= inc

        if not checkPos == invStart:
          break

        invStart -= inc
      invStart += inc
    if uber:
      scan_back = True
    #scan forward
    i += inc
    while not i == end + inc:

      string = ''.join(snps[i])

      # seek through horizon until we find one which encompases i
      while h < len(horizon)-1 and ((inc > 0 and horizon[h][1] < i) or (inc < 0 and horizon[h][0] > i)):
        h += 1

      checkPos = invStart

      while (not checkPos == i) and ((h < len(horizon) and horizon[h][1] >= i and horizon[h][1] >= checkPos and horizon[h][0] <= i and horizon[h][0] <= checkPos) or compatible(snps[i], snps[checkPos])):
        checkPos += inc

      if not checkPos == i:
        if verbose:
          print invStart, "-", (i-inc)
        invs.append((invStart, i - inc))
        if uber:
          invStart = i - inc
        else:
          invStart = i
        break

      i += inc
  invs.append((invStart, end))
  return invs


# returns compatibility of sdp1 vs sdp2
# Ns *should* not cause any problems, they just won't add any incompatibility
def compatible(sdp1, sdp2):
  # Returns 0 if they are incompatible
  # 1 if they are compatible
  gametes = {}
  for i in range(len(sdp1)):
    if sdp1[i] == '0' and sdp2[i] == '0':
      gametes['00'] = True
    elif sdp1[i] == '1' and sdp2[i] == '0':
      gametes['10'] = True
    elif sdp1[i] == '1' and sdp2[i] == '1':
      gametes['11'] = True
    elif sdp1[i] == '0' and sdp2[i] == '1':
      gametes['01'] = True
  return (len(gametes) < 4)
