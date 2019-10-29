#!/bin/bash
#
echo $1 $2
echo "DES StampProd"
echo `hostname`
echo `uname -a`

jobid=$1
objlist=$2

echo "JOBNAME="$jobid "OBJLIST="$objlist

HOSTNAME=`/bin/hostname`
BASEDIR=`pwd`
DATE=`/bin/date +%H%M%S`
#
echo "BASEDIR="$BASEDIR 
source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh 
setup python 2.7.9+0
setup easyaccess
setup swarp 2.36.2+3
setup psycopg2 2.4.6+7
setup despyastro 0.3.5+0
setup pil 1.1.7+13
setup pyfits 3.3+3
setup cfitsio 3.370+0
setup sextractor
setup psfex
setup pypsfex 
setup numpy
setup scipy
setup fitsio 0.9.8rc1+3
#setup extralibs

set -x

#

FITSCUTOUTS_DIR=${BASEDIR}/CoaddCutouts
mkdir ${jobid}
cd ${jobid}
echo "Workdir="`pwd`  >> ${jobid}.log

mkdir data
mkdir workdir

STORAGE=${BASEDIR}/${jobid}/
set +x
#setup FitsCutouts
#setup easyaccess
set -x
ls -l
echo `which python` >> ${jobid}.log 2>&1
#
echo "PRODUCT DIR="$FITSCUTOUTS_DIR  >> ${jobid}.log 2>&1
#
echo $HOME >> ${jobid}.log 2>&1
export DES_SERVICES=${HOME}/.desservices.ini
mkdir etc
cp -r ${FITSCUTOUTS_DIR}/etc/* etc/.
#

echo "HOSTNAME="`/bin/hostname` >> ${jobid}.log
echo $DES_SERVICES >> ${jobid}.log

echo "START at " `/bin/date +%H%M%S` >>${jobid}_start.log
#
# copy objlist from the base directory to workdir
#
cp ../$objlist .
ls -l >> ${jobid}.log 2>&1

${FITSCUTOUTS_DIR}/MakeFitsCutoutsY6.py -w . -d dessci -l ${objlist} -s 14.0 -r 1 >> ${jobid}.log 2>&1
#
echo "NOW CHECK RESULTS" 
# now check results
nc=`ls -l *_cutouts.fits |wc -l`
if [ $nc -gt 0 ]
   then
     tar -cvzf ${jobid}_cutouts.tgz *_cutouts.fits
     echo "FINISH at " `/bin/date +%H%M%S` >>${jobid}_stop.log
else
     echo "FAILED at " `/bin/date +%H%M%S` >>${jobid}_stop.log
fi
cd ../
exit 0
#
