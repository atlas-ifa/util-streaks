#!/bin/bash
# Run strker on a difference image and find sources.
# Syntax: strker.sh basename 

source $ATLAS_HOME/bin/common.bashrc
verify_atlas_env

# Function to report elapsed time since start of script
elapsed() {
  date +%s.%N | awk -v start=$startime '{printf "t= %.3f (strker) \n", $1-start}'
  return 0
}
startime=`date +%s.%N`

umask 0002
set -e

obs=$1
shift

# Bust up the observation into the pieces
split_imparts $obs

# Reduction directory
reddir=""

# Input image file
image=""

# Difference directory
difdir=""

# Input difference file
diff=""

# Output directory
outdir=/tmp

# Output extension
outext=strk

# clean up
CLEAN=2

# how much blabber do you want?  0: very little, 1: a bit, 2: a lot.
VERBOSE=0

# What SNR to search?
sig=20

# Half length after binning: len=bin*(2*HLEN+1)
hlen=8

# Override defaults by the command line
eval $@

echo $(elapsed) $obs strker started on $(hostname) at $(date)

echo $(elapsed) $obs args evaluated

if [[ $VERBOSE -gt 0 ]] ; then set -x ; fi

if [[ -z $reddir ]] ; then reddir=$ATLAS_HOME/red/$IMSITECAM/$IMNN ; fi
if [[ -z $difdir ]] ; then difdir=$ATLAS_HOME/diff/$IMSITECAM/$IMNN ; fi

if [[ -z $image ]] ; then image=$reddir/$obs.fits.fz ; fi
if [[ -z $diff ]] ;  then difdir=$difdir/$obs.diff.fz ; fi

if [[ ! -e $image ]] ; then
    echo Image file $image does not exist
    exit 1
fi

if [[ ! -e $diff ]] ; then
    echo Difference image $diff does not exist
    exit 1
fi

if [[ ! -e $difdir/$obs.ddc ]] ; then
    echo Difference ddc file $difdir/$obs.ddc does not exist
    exit 1
fi

# Get a list of known stars
obsrefcat.sh $obs mlim=14 border=0 > /tmp/$obs.rc2

echo $(elapsed) $obs star list collected

# Run strker
strker $diff $image -strk $outdir/$obs.$outext -srchsig $sig -hlen $hlen -star /tmp/$obs.rc2 -verbose $VERBOSE

echo $(elapsed) $obs strker run

# Create a list of RA,Dec
tail -n +2 $outdir/$obs.$outext | pix2sky -TSKpix $diff - > /tmp/$obs.tmprd

# Create a list of magnitudes
zp=$(fitshdr $diff -v MAGZPT)
etime=$(fitshdr $diff -v EXPTIME)

awk -v zp=$zp -v dt=$etime 'NR>1{peak=$3; sky=$5; len=$8; wid=$9; \
   area=sqrt(len*len+wid*wid)*wid; flux=peak*area/dt; \
   m=99.999; if(flux > 0) m=zp-2.5/log(10)*log(flux); \
   printf "%7.3f\n", m}' $outdir/$obs.$outext > /tmp/$obs.tmpm

# Create a ddc file header
grep "#" $difdir/$obs.ddc > $outdir/$obs.$outext.ddc

# reformat to ddc format with RA, Dec, m, etc.
tail -n +2 $outdir/$obs.$outext | paste -d ' ' /tmp/$obs.tmprd /tmp/$obs.tmpm - | \
   awk '{ dmag=9.99; if($6>0) dmag=sqrt($7/$6*$7/$6+0.001); det=4; \
    printf "%9.5f %9.5f %7.3f %6.3f %8.2f %8.2f %6.2f %6.2f %6.2f %2d %6.2f %3d %3d %3d %3d %3d %3d %3d %3d %3d %2d %7.1f %5.1f\n", \
    $1,$2,$3,dmag,$4,$5,$11,$12,$13,det,$17,$19,$18,$23,0,0,$21,$20,0,$22,0,0,0}' >> $outdir/$obs.$outext.ddc



echo $(elapsed) $obs ddc file created








# Final clean up
if [[ $CLEAN -gt 0 ]] ; then 
  rm -f /tmp/$obs.rc2 /tmp/$obs.tmprd /tmp/$obs.tmpm
  echo $(elapsed) $obs results cleaned
fi

echo $(elapsed) $obs done.

exit 0
