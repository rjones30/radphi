#!/bin/bash
#
# omega_sim.bash - condor job script to simulate omega photoproduction events 
#                  in the radphi detector and run the standard reconstruction
#                  on the results.
#
# author: richard.t.jones at uconn.edu
# version: june 10, 2014

if [[ $# = 0 ]]; then
    echo ""
    echo "Usage: omega_sim.bash <process> <offset>"
    echo " where <process> is the sequence number of the run to process"
    echo " and <offset> is added to it to compute the sequence to be"
    echo " generated."
    echo ""
    exit 1
fi

# give a third argument to simulate omega,Delta instead of omega,P
bgsim=$3

# point to input and output data staging areas
export SRM_URL=srm://grinch.phys.uconn.edu:8443/srm/managerv2?SFN=/Gluex/radphi/omega_sim-9-2014
export X509_USER_PROXY=`pwd`/x509_user_proxy

# define a few radphi-specific environment variables
export TARGET=`uname`
export TARGET_ARCH=`uname -i`
export ROOTSYS=/usr/local/root5
export RADPHI_HOME=/home/radphi/cvs
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${RADPHI_HOME}/lib.${TARGET}/${TARGET_ARCH}:$ROOTSYS/lib
export RADPHI_MAKETARGET=linux
export RADPHI_BIN=$RADPHI_HOME/bin.$RADPHI_MAKETARGET/$TARGET_ARCH
export RADPHI_CONFIG=$RADPHI_HOME/config.`uname -m`
export RADPHI_MAKE=$RADPHI_HOME/makefiles
export RADPHI_LOGDIR=$RADPHI_HOME/logfiles
export TAPEIO_MAXVBLOCK=102400
export LOGD_CONFIG_FILE=$RADPHI_CONFIG/logd.conf
export PATH=$PATH:$RADPHI_BIN:$RADPHI_HOME/scripts:$ROOTSYS/bin

# function for error exit

function error_exit () {
    exitcode=$1
    echo ERROR - exiting with error code $exitcode
    rm -f *.hbook,* *.hbook *.root *.itape *.hst
    rm -f decayt.dat events.in
    exit $exitcode
}

# print a header on stdout
echo Running on worker node `hostname -f`

# look up the sequence to be processed
seqno=`expr $1 + $2 + 1`
echo "Processing sequence $seqno"
omega_sim=omega${bgsim}_sim

# generate
rm -f decayt.dat || error_exit $?
if [[ "$bgsim" == "2p" ]]; then
    ln -s decayt-signal.dat decayt.dat || error_exit $?
elif [[ "$bgsim" == "3D+p" ]]; then
    ln -s decayt-background.dat decayt.dat || error_exit $?
elif [[ "$bgsim" == "2D+p" ]]; then
    ln -s decayt-background2.dat decayt.dat || error_exit $?
elif [[ "$bgsim" == "2D+n" ]]; then
    ln -s decayt-background3.dat decayt.dat || error_exit $?
elif [[ "$bgsim" == "2D0p" ]]; then
    ln -s decayt-background4.dat decayt.dat || error_exit $?
else
    echo "Unknown simulation type $bgsim"
    exit 2
fi
./mcwrap.$TARGET_ARCH -l4.390 -u5.390 -comega || error_exit $?

# simulate (change to true if you want to run gradphi)
if false; then
    echo omega.itape >events.in
    sed -i 's/^RNDM .*/RNDM '`expr \( $seqno % 215 \) + 1`'/' control.in
    ./gradphi.$TARGET_ARCH || error_exit $?
    simitape=${omega_sim}_$seqno.itape
    mv simData.itape $simitape || error_exit $?
else
    omega_sim=omega${bgsim}_gen
    simitape=${omega_sim}_$seqno.itape
    mv omega.itape $simitape
fi

# save the simulation results
echo saving $simitape to $SRM_URL
srmcp -2 -overwrite_mode=ALWAYS file:///`pwd`/$simitape $SRM_URL/$simitape || error_exit $?

# reconstruct
./mctuplet.$TARGET_ARCH --count=1000000000 $simitape || error_exit $?

# convert reconstructed hbook output files to root
n=0
while [[ -r mctuple.hbook,$n ]]; do
    ./h2root.$TARGET_ARCH mctuple.hbook,$n ${omega_sim}_${seqno}-$n.root 1 1 65536 || error_exit $?
    n=`expr $n + 1`
done
./h2root.$TARGET_ARCH mctuple.hbook ${omega_sim}-${seqno}-$n.root 1 1 65536 || error_exit $?

# save the output root files to srm
for rootfile in `ls *.root`; do
    echo saving $rootfile to $SRM_URL
    srmcp -2 -overwrite_mode=ALWAYS file:///`pwd`/$rootfile $SRM_URL/$rootfile || error_exit $?
done

# cleanup and exit
rm -f *.hbook,* *.hbook *.root *.itape *.hst
rm -f decayt.dat events.in
