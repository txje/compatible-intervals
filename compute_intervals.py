import genotype_intervals, simple_haplotype_intervals as haplotype_intervals, maxk, haplotype_character_removal_intervals as hap_cr_intervals
import hapmap
import argparse

def writeIntervals(intervals, outfile):
  fout = open(outfile, 'w')
  fout.write("start,end")
  for i in intervals:
    fout.write("\n%i,%i" % (i[0], i[1]))
  fout.close()

def main(input, pos, sample, delimiter, prefix, outdir, cr=None, verbose=False):

  snpList = hapmap.loadSNPmatrix(input, verbose=verbose, delimiter=delimiter, pos_idx=pos, sample_idx=sample, filter_positions=None)

  # are these heterozygous?
  geno = False
  for s in snpList:
    if 2 in s:
      geno = True
      break

  if cr is None or cr == 0:
    if geno:
      uber, groups = genotype_intervals.uberScan(snpList)
      writeIntervals(uber, prefix + ".uber_intervals.csv")

      lr, groups = genotype_intervals.lrScan(snpList)
      writeIntervals(lr, prefix + ".lr_intervals.csv")

      rl, groups = genotype_intervals.rlScan(snpList)
      writeIntervals(rl, prefix + ".rl_intervals.csv")

    else:
      uber = haplotype_intervals.uberScan(snpList)
      writeIntervals(uber, prefix + ".uber_intervals.csv")

      lr = haplotype_intervals.lrScan(snpList)
      writeIntervals(lr, prefix + ".lr_intervals.csv")

      rl = haplotype_intervals.rlScan(snpList)
      writeIntervals(rl, prefix + ".rl_intervals.csv")

    intervals = maxk.maxkScan(lr, rl, uber)
    print "%i intervals" % len(intervals)
    writeIntervals(intervals, prefix + ".maxk_intervals.csv")

  # ---------- character removal ---------------
  else:
    removable = cr

    # compute CR uber interval horizon, then use that skyline to do normal LR, RL, Uber intervals
    # - this is necessary so that maxK works correctly
    horizon, removed = hap_cr_intervals.CRUberScan(snpList, removable, verbose=verbose)
    fname = prefix + ".cr_horizon.csv"
    print "Writing CR horizon intervals to %s" % fname
    writeIntervals(horizon, fname)

    uber = hap_cr_intervals.uberScan(snpList, removable, inv_horizon=horizon, verbose=verbose)
    fname = prefix + ".uber_intervals.csv"
    print "Writing uber intervals to %s" % fname
    writeIntervals(uber, fname)

    lr = hap_cr_intervals.lrScan(snpList, removable, inv_horizon=horizon, verbose=verbose)
    fname = prefix + ".lr_intervals.csv"
    print "Writing LR intervals to %s" % fname
    writeIntervals(lr, fname)

    rl = hap_cr_intervals.rlScan(snpList, removable, inv_horizon=horizon, verbose=verbose)
    fname = prefix + ".rl_intervals.csv"
    print "Writing RL intervals to %s" % fname
    writeIntervals(rl, fname)

    intervals = maxk.maxkScan(lr, rl, uber)
    print "%i intervals" % len(intervals)
    fname = prefix + ".maxk_intervals.csv"
    print "Writing maxK intervals to %s" % fname
    writeIntervals(intervals, fname)


if __name__ == "__main__":
  parser = argparse.ArgumentParser("Compute compatible intervals using haplotypes, genotypes, and/or character removal")
  parser.add_argument("input", help="Delimited haplotype/genotype input file")
  parser.add_argument("pos", type=int, help="Field index of position in input file")
  parser.add_argument("sample", type=int, help="Field index of first sample in input file")
  parser.add_argument("delimiter", help="RE string to delimit lines in the input")
  parser.add_argument("prefix", help="Output prefix")
  parser.add_argument("--cr", type=int, help="Use character removal, with the given number removed")
  parser.add_argument("--verbose", action="store_true", help="Print verbose output", default=False)
  args = parser.parse_args()
  main(args.input, args.pos, args.sample, args.delimiter, args.prefix, args.cr, args.verbose)
