export BLASTDB:=/local_uniref/uniref/uniref90

cd dataset
time for i in *.fasta
do
  echo "Running $i on PSI-BLAST..."
  time psiblast -query $i -evalue 0.01 -db uniref90.fasta -num_iterations 3  -out Psiblast_Output/$i.psiblast -out_ascii_pssm PSSM/$i.pssm -num_threads 8

done

echo "Done !"

############################################

#u2352@stud36:~/Desktop/dataset$ 

#cd Batch_1
time for i in *.fasta 
do   
  echo "Running $i on PSI-BLAST..."
  time psiblast -query $i -evalue 0.01 -db /local_uniref/uniref/uniref90/uniref90.db -num_iterations 3  -out Psiblast_Output/$i.psiblast -out_ascii_pssm PSSM/$i.pssm -num_threads 8

done 
echo "Done !"



##########################################

scp -r dataset u2352@stud36.dbb.su.se:~/Desktop


cp Batch*/Psiblast_Output/*.psiblast Psiblast_Output/

cp Batch*/PSSM/*.pssm PSSM/

scp -r  u2352@stud24.dbb.su.se:~/Desktop/dataset ~/Desktop/Bioinf_Pro_Summer

scp -r dataset u2352@stud36.dbb.su.se:~/Desktop

htop