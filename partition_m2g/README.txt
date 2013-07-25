
PACKAGE:

	partition_m2g


AUTHORS:  Wenqing Wang, Myles English
LICENCE: None and GPL
COPYRIGHT:  Wenqing Wang, Myles English


DESCRIPTION:

	Tools for partitioning a GeoSys mesh into domains.


USAGE:

	See https://geosys.ufz.de/trac/wiki/DomainDecomposition


DETAILS:

        Parallel processing of a simulation requires that the finite element
        mesh be decomposed into domains so that each domain can be sent to a
        separate cpu.

	partition_m2g contains the following tools:

	m2g - converts a GeoSys mesh into a METIS format mesh.

	partition.sh - a wrapper script for m2g and a third party partitioning
	executable

	A third party finite element mesh partitioning tool will need to be
	installed seperately (i.e. METIS).



DOCUMENTATION:

        See https://geosys.ufz.de/trac/wiki/DomainDecomposition

	This directory contains the source code for the mesh partitioning
	tool referred to in the attachment on the Wiki page (linked to from
	the front page by the text "How to compile source code and run
	simulation on Linux/UNIX", titled "A short description of how to
	compile GeoSys and run benchmarks under UNIX/Linux").

	The URL to the attachment "Parallel Simulation by
	Rockflow/GeoSys:Howto" is:
	https://geosys.ufz.de/trac/attachment/wiki/LinuxPage/Parallel%20simulation%20by%20GeoSys.pdf


SEE ALSO:

    ./rungeosys.sh - a wrapper around partition.sh and a Sun Grid Engine
    jobscript for parallel architectures


TODO:

	o  Add Windows install details
	o  m2g - default destination of output files should be the current
	   working directory, not the one where the .msh file starts from
	o  Pipes for input and output


Copyright Myles English, 2008-10-23
