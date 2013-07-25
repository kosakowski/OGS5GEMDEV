#$ -N xxNAMExx                 
#$ -pe openib xxNCPUSxx
#$ -l h_rt=xxHOURSxx:00:00

# Send mail at the beginning and end of the job
#$ -m be

# Who to send mail to at the end of the job
#$ -M "${MAIL}"

# Initialise environment modules (not sure what this does but it looks good)
. /etc/profile.d/modules.sh

# Use OpenMPI and GNU compilers for Infiniband network
module load openmpi/gcc/64/1.2.5/openib

# --mca btl openib,self specifies to use Infiniband network running MPI program
# Note: for resource reservation use: -R y
mpiexec --mca btl openib,self -np $NSLOTS xxEXExx xxPCSxx
