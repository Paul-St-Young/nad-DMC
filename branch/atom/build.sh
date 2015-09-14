rm -rf build
mkdir build
cd build
export CC=mpicc
export CXX=mpicxx
cmake -Wno-dev ..
sed -i 's/\/usr\/bin\/svn//' CMakeCache.txt
cmake -Wno-dev ..
make -j16
