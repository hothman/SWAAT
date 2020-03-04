#!/bin/bash -ue
cat UGT1A1.fa
python /home/houcemeddine/modules/PRODRES/PRODRES/PRODRES.py --input UGT1A1.fa 	                  --output ./  	                  --pfam-dir /home/houcemeddine/modules/PRODRES/db/pfam 	                  --pfamscan-script /home/houcemeddine/modules/PRODRES/PfamScan/pfam_scan.pl 	                  --pfamscan_bitscore 2 	                  --fallback-db-fasta /home/houcemeddine/modules/PRODRES/db/uniprot/uniref90.fasta 	                  --second-search psiblast 	                  --psiblast_e-val 0.001 	                  --psiblast_iter 3 	                  --prodres-db /home/houcemeddine/modules/PRODRES/db/prodres_db.nr100.sqlite3

   rootname=`head -n 1   UGT1A1.fa |awk -F  "|" '{print $2}'`
   outputfilename=$(echo "$rootname".pssm )
   mv $rootname/outputs/psiPSSM.txt ./$outputfilename
