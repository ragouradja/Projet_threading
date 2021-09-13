path=../data/globulin/

for fasta in $path*.fasta;do
for pdb in $path*.pdb;do
python double.py -p $pdb -f $fasta --blosum --weight -1 --log globulin_blosum.log
done
done