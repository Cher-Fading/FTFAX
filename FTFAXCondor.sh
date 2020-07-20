#!/bin/bash
#input="/atlasgpfs01/usatlas/data/cher97/mc16_5TeV.txt"
input="../GetStuff/$2_root_pnfs.txt"
#input="mc16_5TeV_short.txt"

cd ~/FTFAX

#indexline=$1
linenumber=0
while IFS= read -r line
do
  if [ $1 -eq $linenumber ]
  	then
	mkdir '/atlasgpfs01/usatlas/data/cher97/tempin'$linenumber
	xrdcp 'root://dcgftp.usatlas.bnl.gov:1096/'$line '/atlasgpfs01/usatlas/data/cher97/tempin'$linenumber
	filename=/atlasgpfs01/usatlas/data/cher97/tempin$linenumber/$(ls /atlasgpfs01/usatlas/data/cher97/tempin$linenumber)
  	root -b -q -l 'bTagJF_condor.C("'$filename'","'$2'",'$3')'
	sleep 2
	rm -rf '/atlasgpfs01/usatlas/data/cher97/tempin'$linenumber
  fi
  linenumber=$((linenumber+1))
done < "$input"


