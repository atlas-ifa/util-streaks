#!/bin/bash

rootdir=/tmp

mkdir -p $rootdir/STRK
cd $rootdir/STRK

cat > strk.pro <<EOF
! Use monsta to make a streak-sigma image
! Syntax: monsta strk.pro phi bin len
! E.g. monsta strk.pro 156 2 10
! Reads   r,d,c,s.phi, writes sig.phi

  set eadu=2 bin={arg3} len={arg4} rtn=bin*bin*len

  rd 1 r.{arg2}
  rd 2 d.{arg2}
  rd 3 c.{arg2}
  rd 4 s.{arg2}
  cop 7 3
  mi 7 7
  cop 6 4
  mi 6 6
  ai 6 7
  sqrt 6
  cop 7  1
  mc 7 eadu/rtn
  sqrt 7
  cop 5  2
  si 5  6
  di 5  7

  wd 5 sig.{arg2}
q
EOF

len=10  minormax=6  dphi=$(( 120 / len))
for ((phi=0; phi<180; phi+=dphi)) ; do
  echo $phi
  p3=$(printf "%03d" $phi)
  strkonv /atlas/diff/02a/58071/02a58071o0503c.diff.fz d.$p3 -phi $phi -len $len -nlam 1 -dcs 0
  strkonv /atlas/diff/02a/58071/02a58071o0503c.diff.fz c.$p3 -phi $phi -len $len -nlam 1 -dcs 1
  strkonv /atlas/diff/02a/58071/02a58071o0503c.diff.fz s.$p3 -phi $phi -len $len -nlam 1 -dcs 2
  strkonv /atlas/red/02a/58071/02a58071o0503c.fits.fz r.$p3 -bkg -clip -1000,100000 -phi $phi -len 10 -nlam 1 -dcs 0

  echo monsta
  time monsta strk.pro $p3 2 $len
  echo tphot
  time tphot sig.$p3 -eadu 1 -rad $len -sig 5 -min 5 -move $len -fwmax 200 -chin 30 -out det.tmp
  awk -v MM=$minormax 'NR==1{print $0} $11<MM{print $0}' det.tmp > det.$p3
done

cat det.[0-9]?? | grep v sky > det.all

exit 0
