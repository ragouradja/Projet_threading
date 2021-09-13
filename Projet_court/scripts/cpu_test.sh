pdb=$1
fasta=$2
for cpu in `seq 1 8`;do
python double.py -p $pdb -f $fasta --cpu $cpu --log plt_time64.log
done
