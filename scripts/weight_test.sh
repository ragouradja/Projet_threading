pdb="../data/1fxr.pdb"
fasta="../data/1fxd.fasta"

for runs in `seq 1 5`;do
for weight in `seq 1 4`;do
python double.py -p $pdb -f $fasta --cpu 8 --blosum --weight $weight --log weight.log
done
done