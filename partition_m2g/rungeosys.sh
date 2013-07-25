#!/bin/bash

function usage { echo "

rungeosys.sh 0.0.1
Usage: rungeosys.sh -in PROJECT [ options ]

  -in PROJECT   the path to the directory containing the GeoSys project
  [-n NAME]     the name to run the job under (defaults to 'geosys')
  [-out OUT]    the directory to put the output in (defaults to PROJECT)
  [-cpus NCPUS] the number of CPUs to use (defaults to 2)
  [-exe EXE]    the path to the GeoSys executable (defaults to
                '/exports/work/geos_ess/geosys/bin/rf4')
  [-t HOURS]    the maximum run time on the cluster, jobs will be terminated
                by the Sun Grid Engine at t > HOURS

"
}


# GET THE COMMANDLINE PARAMETERS ______________________________________________

while [ "$1" != "" ]; do
    case $1 in
        -n     | --name                )    shift
                                            NAME=$1
                                            ;;
        -cpus  | --cpus                )    shift
                                            NCPUS=$1
                                            ;;
        -exe   | --exe-path-to-geosys  )    shift
                                            EXE=$1
                                            ;;
        -in    | --input-project       )    shift
                                            PROJ=$1
                                            ;;
        -t     | --time-in-hours       )    shift
                                            HOURS=$1
                                            ;;
        -out   | --output-directory    )    shift
                                            OUT=$1
                                            ;;
        -h     | --help                ) usage
                                         exit
                                         ;;
        *                              ) usage
                                         exit 1
    esac
    shift
done


# MAKE UP ANY DEFAULT PARAMETERS NECESSARY ____________________________________

NCPUS=${NCPUS:-2}   # because most machines are dual-core
EXE=${EXE:-/exports/work/geos_ess/geosys/bin/rf4}     # the GeoSys executable
OUT=${OUT:-${PROJ}} # output directory defaults to input directory
HOURS=${HOURS:-1}   # defaults to one hour run time

# PARTITION THE MESH __________________________________________________________
# This should find any problems with the command line parameters

partition.sh -n ${NCPUS} -if ${PROJ} || exit 1

pcs=`ls ${PROJ}/*.pcs`   # assumes only one pcs file
PCS=${pcs%*.pcs}
bn=`basename ${PROJ}`    # NAME of job defaults to part of the path to the PROJ
NAME=${NAME:-${bn##*.}}


# SETUP A JOBSCRIPT ___________________________________________________________

jobscript=${OUT}/geosys_jobscript_${NAME}.sh

# substitute the parameters into the custom jobscript
function setupJobscript
{
  sed -i s#xxHOURSxx#${HOURS}#g ${OUT}/geosys_jobscript_${NAME}.sh
  sed -i s#xxMAILxx#${USER}#g ${OUT}/geosys_jobscript_${NAME}.sh
  sed -i s#xxPCSxx#${PCS}#g ${OUT}/geosys_jobscript_${NAME}.sh
  sed -i s#xxPROJxx#${PROJ}#g ${OUT}/geosys_jobscript_${NAME}.sh
  sed -i s#xxNAMExx#${NAME}#g ${OUT}/geosys_jobscript_${NAME}.sh
  sed -i s#xxNCPUSxx#${NCPUS}#g ${OUT}/geosys_jobscript_${NAME}.sh
  sed -i s#xxEXExx#${EXE}#g ${OUT}/geosys_jobscript_${NAME}.sh
  sed -i s#xxOUTxx#${OUT}#g ${OUT}/geosys_jobscript_${NAME}.sh
}

# copy the template jobscript for this project
cp /exports/work/geos_ess/geosys/bin/geosys_jobscript.sh ${jobscript} && \
setupJobscript


# RUN THE JOBSCRIPT ___________________________________________________________

qsub ${OUT}/geosys_jobscript_${NAME}.sh

exit 0


