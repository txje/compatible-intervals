# ----------------------------------------------------------------
# Jeremy Wang
# Sep 12, 2014
#
# Modified from Summer 2010 version to dramatically simplify input and
# data structures required
#
# Compute standard four-gamete compatible intervals including
# LR, LR, Uber, and Max-k using homozygous genotypes
# ----------------------------------------------------------------

def lrScan(snps):
  return allScan(snps, 0, len(snps)-1, uber=False)

def rlScan(snps):
  scan = allScan(snps, len(snps)-1, 0, uber=False)
  reversed = [(s[1], s[0]) for s in scan[::-1]]
  return reversed

def uberScan(snps):
  return allScan(snps, 0, len(snps)-1, uber=True)

# ----------------------------------------------------------------
# LR, RL, and Uber Scan computation
# ----------------------------------------------------------------

# allScan performs any of left->right, right->left (if start > end), and uber scan (if uber flag is set)
# 8/2010 - added sdp_hash to prevent checking of duplicate SDPs a la the time proof
def allScan(snps, start, end, uber=False):
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
        string = ''.join(snps[invStart])
        if not sdp_hash.has_key(string): # do not re-check duplicate SDPs, this gives us linear time scan
          while (not checkPos == invStart) and compatible(snps[invStart], snps[checkPos]):
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
      string = ''.join(snps[i])
      if not sdp_hash.has_key(string): # do not re-check duplicate SDPs, this gives us linear time scan
        while (not checkPos == i) and compatible(snps[i], snps[checkPos]):
          checkPos += inc
        if not checkPos == i:
          invs.append((invStart, i - inc))
          if uber:
            invStart = i - inc
          else:
            invStart = i
          sdp_hash = {}
          break
        sdp_hash[string] = 1
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
