#!/bin/bash
#
# Copyright Myles English Mon, 27 Oct 2008 15:03:52 +0000
#
# A wrapper script to automate finite element mesh partitioning.
#
# Used when preparing a project to be run by GeoSys on parallel architecture.
# Depends on the GeoSys utility m2g and METIS.

USAGE="
Usage: partition.sh -n nDOMAINS [OPTIONS]                                  

Partition a GeoSys format mesh FILE (*.msh) into nDOMAINS number of domains.

Mandatory arguments to long options are mandatory for short options too.

 -n DOMAINS, --n-domains          the number of domains to decompose mesh into

 [-if FILE, --input-file=FILE]    mesh to be partitioned, defaults to *.msh file
                                  in the current directory.  If FILE is a 
                                  directory then it will be checked for a mesh
                                  file.

 [-of FILE, --output-file=FILE]   output the domain decomposition to FILE,
                                  defaults to the same directory as the input
                                  file {inputname}.ddc

 [-p PROGRAM, --program=PROGRAM] Executable to use to partition mesh. Defaults
                                 to partdmesh of the METIS package.

 [--help]                        Show usage.

Output: A GeoSys format domain decomposition file (*.ddc).
"

# DEFAULTS
program="partdmesh"
infile=.  # defaults to current directory, the *.msh file is found later
outfile=
allgood=0

# SUPPORTED PARTITIONING PROGRAMS
PROGRAMS="partdmesh"                          # add new programs to this list

function usage {
    echo -e "${USAGE}"
}

while [ "$1" != "" ]; do
    case $1 in
        -n  | --n-domains  )    shift
                                domains=$1
                                ;;
        -if | --input-file )    shift
                                infile=$1
                                ;;
        -of | --output-file )   shift
                                outfile=$1
                                ;;
        -p  | --program )       shift
                                program=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

function checkDomains {
  if [ $domains -eq $domains 2> /dev/null ]; then
    return 0
  else
    echo -e "\nError: ${domains} DOMAINS is not an integer\n"
    allgood=1
    return 1
  fi
}

function checkInFileExists {
  if [ -f "${infile}" ]; then       # its a file
    return 0                        # its a directory
  elif [ -d "${infile}" ]; then
    let counter=0
    lof=`ls ${infile}/*.msh`        # get all the *.msh files
    for i in $lof; do
        let counter++
        if [ ${counter} -gt 1 ]; then
            echo -e "\nError: ${lof} - more than one .msh file."
            allgood=1
            return 1
        else infile=$i
        fi
    done
  else 
    echo -e "\nError: ${infile} is neither a mesh file nor a directory containing one."
  fi
}

function checkOutFileDoesNotExist {
  if [ -f "${outfile}" ]; then
   echo -e "\nError: ${outfile} already exists! Will not overwrite."
   allgood=1
   return 1
  fi
}

function checkFileIsMesh {
    ext=${infile##*.}
    if [ ${ext} != "msh" ]; then
      echo -e "\nError: ${infile} must have a GeoSys file extension of '.msh'."
      allgood=1
      return 1
    fi
}

function checkProgram {
    for i in $PROGRAMS; do
        if [[ ${program}=$i ]]; then return 0; fi
    done
    echo -e "\nError: ${program} not supported, must be one of ${PROGRAMS}."
    allgood=1
    return 1
}

function partitionWithXpartdmesh
{
    m2g $GEOSYS_MESH_BASE #2                   # makes GEOSYS_PROJ_NAM.mesh
    $program ${GEOSYS_MESH_BASE}.mesh $domains # makes .mesh.epart.{nDOMAINS}
                                               # ...and .mesh.npart.{nDOMAINS}

    m2g $GEOSYS_MESH_BASE -${domains} # makes GEOSYS_MESH_BASE.${domains}ddc
                                # and GEOSYS_MESH_BASE[0-(${domains}-1)].msh

    mv ${GEOSYS_MESH_BASE}.${domains}ddc ${outfile:-${GEOSYS_MESH_BASE}.ddc} &&
    echo -e "\nMesh decomposition file ${outfile:-${GEOSYS_MESH_BASE}.ddc} written."

    # clean up
    rm ${GEOSYS_MESH_BASE}*[0-9].msh 2> /dev/null
    rm ${GEOSYS_MESH_BASE}.mesh 2> /dev/null
    rm ${GEOSYS_MESH_BASE}.mesh.{e,n}part.${domains} 2> /dev/null
}

# MAIN ___________________________________________

if [[ ! `which m2g` ]]; then
  echo -e "\n The program m2g could not be found! (You may need to add it to \
           your PATH environment variable.)"
  exit 1
fi

# TEST THE INPUT PARAMETERS
checkDomains
checkProgram
checkInFileExists && checkFileIsMesh
GEOSYS_MESH_BASE=${infile%*.msh}
checkOutFileDoesNotExist
if (( ${allgood}==1 )); then usage; exit 1; fi

# PARTITION THE MESH
eval partitionWithX${program}  # done like this so that it can be extended with
                               # other third party mesh partitioning tools