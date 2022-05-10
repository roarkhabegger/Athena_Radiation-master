python3 configure.py \
 --prob=ParkerInst \
 --nghost=2 \
 -b \
 -cr \
 -hdf5 \
 -mpi \
 --nscalars=2 \
 --mpiccmd=h5pcc \
 --cflag='-DH5_HAVE_PARALLEL -lstdc++' \
 --cxx='g++-simd'
