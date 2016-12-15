This is a Subversion repository; use the 'svnadmin' tool to examine
it.  Do not add, delete, or modify files here unless you know how
to avoid corrupting the repository.

Visit http://subversion.apache.org/ for more information.

How to get and compile this version under linux:

You need boost and boost development libraries to be installed (especially boost-threads).

git clone https://github.com/kosakowski/OGS5GEMDEV

cd OGSGEMDEV/sources
mkdir build 
cd build
cmake .. -DOGS_FEM_GEMS=ON
make
make install

optional CMAKE arguments:
-DBOOST_ROOT=/path/to/boost/

