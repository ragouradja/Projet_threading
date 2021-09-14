pdb=$1
fasta=$2
for cpu in `seq 1 8`;do
python main.py -p $pdb -f $fasta --cpu $cpu --log plt_time64.log
done
