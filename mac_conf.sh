python3 configure.py \
 --prob=ParkerInst \
 --nghost=2 \
 --nscalars=2 \
 -b \
 -cr \
 -hdf5 \
 -mpi \
 -h5double   \
 --mpiccmd='/opt/usr/local/bin/mpicxx' \
 --include='/opt/usr/local' \
 --cflag='-DH5_HAVE_PARALLEL'
