#!/bin/bash
tag=MAD
for i in me1A me1B me1C me1D me1E me1F me1G me1L me1M me1N me1O me1P 
#for i in me5A me6A me6B me6C me6D me6E me6F me6G me6H me6I me6J
do
#data 
  find -L /pnfs/minerva/persistent/users/hsu/MAD/data/${i}/ -name *.root -type f -exec realpath {} + | pnfs2xrootd.sh > playlist_ccqenue_data${i}_nx_${tag}.txt
#mc
  find -L /pnfs/minerva/persistent/users/hsu/MAD/mc/${i}/ -name *.root -type f -exec realpath {} + | pnfs2xrootd.sh > playlist_ccqenue_mc${i}_nx_${tag}.txt
#ncdif
  find /pnfs/minerva/persistent/users/hsu/MAD/mc/${i}-NCDIF/ -name *.root -exec realpath {} + | pnfs2xrootd.sh > playlist_ccqenue_mc${i}_NCDIF_${tag}.txt
#Big NuE
  find /pnfs/minerva/persistent/users/hsu/MAD/mc/BigNuE_off/${i}/ -name "*.root" -exec realpath {} + | pnfs2xrootd.sh > playlist_ccqenue_mc${i}_BigNuE_${tag}.txt
#ext 2p2h
  find /pnfs/minerva/persistent/users/hsu/MAD/mc/2p2h/${i}/ -name "*.root" -exec realpath {} + | pnfs2xrootd.sh > playlist_ccqenue_mc${i}_ext_2p2h_${tag}.txt
done

