JID1=`qsub -h pbs_surfjet_3dnp.sh`
JID2=`qsub -W depend=afterok:$JID1 pbs_surfjet_3dnp.sh`

qrls $JID1 # Release first job
