input=example/example_haplotypes.txt
pos=1
sample=2
delimiter='\s'
prefix=example/example

cd ..
python compute_intervals.py $input $pos $sample $delimiter $prefix --verbose
