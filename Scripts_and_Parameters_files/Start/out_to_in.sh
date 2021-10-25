#!/bin/bash

for  filein in confout_* 
do
	fileout=`echo $filein |  sed s/out/in/`
	echo "mv  $filein to  $fileout"
	mv $filein  $fileout
done
