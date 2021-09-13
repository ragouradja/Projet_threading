pdb="../data/1fxr.pdb"
fasta="../data/1fxd.fasta"

for runs in `seq 1 10`;do
python double.py -p $pdb -f $fasta --cpu 8 --blosum  --zscore 10 --log with_options.log
python double.py -p $pdb -f $fasta --cpu 8 --blosum --weight 1  --zscore 10 --log with_options.log
python double.py -p $pdb -f $fasta --cpu 8 --zscore 100 --log with_options.log
python double.py -p $pdb -f $fasta --cpu 8  --dssp  --zscore 10 --log with_options.log
python double.py -p $pdb -f $fasta --cpu 8 --log without_options.log
done