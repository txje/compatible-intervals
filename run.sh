input=<delimited text file, one line per SNP>
pos=<index of SNP position in line>
sample=<index of first sample name in header line>
preifx=<output prefix, including directory>

python compute_intervals.py $input $pos $sample '$delimiter' $prefix --verbose
