#!/bin/bash
#input="/atlasgpfs01/usatlas/data/cher97/mc16_5TeV.txt"
input="../GetStuff/$2_root.txt"
#input="mc16_5TeV_short.txt"

cd ../FTFAX

#indexline=$1
linenumber=0
while IFS= read -r line
do
  if [ $1 -eq $linenumber ]
  	then
  	#python test.py "$line"
  	#echo $line
  	xrdcp root://dcgftp.usatlas.bnl.gov:1096/$line root://rftpexp.usatlas.bnl.gov:/atlasgpfs01/usatlas/data/cher97/
  fi
  linenumber=$((linenumber+1))
done < "$input"
