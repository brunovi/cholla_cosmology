#!/bin/sh -e

myprefix=$MEMBERWORK/ast125/pfft
PFFT_VERSION=git
FFTW_VERSION=3.3.4.11
INSTDIR=$myprefix/pfft-$PFFT_VERSION
FFTWDIR=/opt/cray/fftw/3.3.4.11/interlagos
TMP="${PWD}/tmp-pfft-$PFFT_VERSION"
LOGFILE="${TMP}/build.log"

bash check if directory exists
if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
	read answer
	if [ ${answer} = "y" ]; then
		rm -rf $TMP
	else
		echo "Program aborted."
		exit 1
	fi
fi

mkdir $TMP && cd $TMP
# cd $TMP
wget https://github.com/mpip/pfft/archive/master.tar.gz
# wget http://www.tu-chemnitz.de/~potts/workgroup/pippig/software/pfft-$PFFT_VERSION.tar.gz
gzip -dc master | tar xvf -
# cd pfft-$PFFT_VERSION
# cp ./tmp-pfft-1.0.8-alpha/pfft-1.0.8-alpha/configure ./tmp-pfft-git/pfft-master/
cd pfft-master
./bootstrap.sh
./configure --prefix=$INSTDIR --disable-shared --enable-openmp \
  CPPFLAGS="-I$FFTWDIR/include" \
  LDFLAGS="-L$FFTWDIR/lib" \
  FC=ftn CC=cc MPICC=CC MPIFC=ftn 2>&1 | tee $LOGFILE

make -j 4 2>&1 | tee -a $LOGFILE
make install 2>&1 | tee -a $LOGFILE
