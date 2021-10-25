for  filein in confout* 
do
	export  fileout=`echo $filein |  sed s/out/in/`
	echo "mv  $filein to  $fileout"
	mv $filein  $fileout
done
