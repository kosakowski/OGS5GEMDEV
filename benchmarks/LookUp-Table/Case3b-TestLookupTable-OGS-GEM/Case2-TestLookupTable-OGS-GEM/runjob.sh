#cd /home/kosakowski/work/SMAN/OpaCemConstPorFineGrid
#cp /home/kosakowski/SOFT/OGSGEMS/sources/BUILDOGSGEMMPI/bin/ogs ./
#nohup ogsgems4r cc1p |grep "Time step\|failed" 1>run.log 2>&1 
#
#PID=5870
#
#while [[ ( -d /proc/$PID ) && ( -z `grep zombie /proc/$PID/status` ) ]]; do
#    sleep 1
#done
#
nohup nice -19 /home/kosakowski/GitMigration/sources/build/bin/ogs bench |grep "Time step\|failed this" 1>run.log 2>&1 

