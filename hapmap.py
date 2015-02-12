# ------------------------------------------------------
# Utilities to handle HapMap data
# ------------------------------------------------------

import numpy
import re

# genotypes:
# rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode NA06984 NA06985

# phased haplotypes:
# rs# phys_position NAXXXXXX_A NAXXXXXX_B NAXXXXXX_A NAXXXXXX_B


# loads the SNP text format
# into an array of positions
def loadPositions(f, verbose=False, delimiter=' ', pos_idx=0, sample_idx=1):
  data = [re.split(delimiter, line.strip()) for line in open(f).read().strip().split('\n')]

  # look at header, sample names
  header = data[0]
  if verbose:
    print "%i SNPs" % (len(data)-1)

  positions = numpy.array([int(data[d][pos_idx]) for d in xrange(1, len(data))], dtype="u4")

  return positions


# loads the SNP text format
# into a list of strings consisting only of {0,1,2}, one for each sample
def loadSNPmatrix(f, n_cutoff=None, verbose=False, delimiter=' ', pos_idx=0, sample_idx=1, filter_positions=None):
  data = [re.split(delimiter, line.strip()) for line in open(f).read().strip().split('\n')]

  # look at header, sample names
  header = data[0]
  if type(sample_idx) == type([]):
    sample_idx = [header.index(s) for s in sample_idx]
  else:
    sample_idx = range(sample_idx, len(header))
  if verbose:
    print "Loading %i samples..." % (len(sample_idx))
    print [header[s] for s in sample_idx]
    print "%i SNPs" % (len(data)-1)

  # pivot data table, convert alleles to {0,1,2,N}
  snps = []
  if filter_positions is not None:
    f = 0
  for i in xrange(1, len(data)):
    # filter by positions, if necessary
    if filter_positions is not None:
      while f < len(filter_positions) and filter_positions[f] < int(data[i][pos_idx]):
        f += 1
      if f < len(filter_positions) and filter_positions[f] != int(data[i][pos_idx]):
        continue
    alleles = ''.join([data[i][s] for s in sample_idx])
    allele_counts = sorted([(a, alleles.count(a)) for a in "ACGTN"], key = lambda a: -1*a[1])
    if sum([a[1] for a in allele_counts]) < len(alleles):
      if verbose:
        print "Unknown allele (probably indel {D,I}) in '%s'" % alleles
      continue
    majority = allele_counts[0]
    if allele_counts[3][1] > 0 or (allele_counts[2][1] > 0 and (allele_counts[3][0] == 'N' or allele_counts[4][0] == 'N')):
      if verbose:
        print "Three or more alleles (probably known on chrX):", allele_counts
      continue
    alleles = []
    for s in sample_idx:
      allele = data[i][s]
      if 'N' in allele:
        allele = 'N'
      elif len(allele) > 1 and allele[0] != allele[-1]:
        allele = '2'
      elif allele[0] == majority[0]:
        allele = '0'
      else:
        allele = '1'
      alleles.append(allele)
    snps.append(alleles)

  # compute and filter by N fraction
  if n_cutoff is not None:
    s = 0
    while s < len(samples):
      n_fraction = float(samples[s].count('N'))/len(samples[s])
      if verbose:
        print "%s: %.4f N" % (header[sample_idx[s]], n_fraction)
      if n_fraction > n_cutoff:
        samples.pop(s)
      else:
        s += 1

  if verbose:
    print "matrix: %i samples x %i SNPs" % (len(snps[0]), len(snps))

  return snps
