path_glob=../data/globulin/
path_imm=../data/immunoglobulin/

for fasta in $path_glob*.fasta;do
for pdb in $path_imm*.pdb;do
python double.py -p $pdb -f $fasta --blosum --weight -1 --chain H --log globulin_imm_blosum.log
done
done