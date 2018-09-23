export ANNAL=${DIR}"/Analysis/"

#export rep=`grep "n_skip" parameters`
#export n=`wc -l  Ener_scal | sed s/"Ener_scal"//`
#echo $n
#let n=2*n/3
#sed s/"$rep"/"n_skip  = $n"/ parameters > tmp
#mv tmp parameters


for filename in *_scal; do
    echo $filename
    export filename1=$filename"J"
    if [ "$filename1" -ot "$filename" ]; then
       cp $filename Var_scal
       $ANNAL/cov_scal.out
       mv "Var_scalJ"  $filename"J"
       for filename2 in Var_scal_Auto_*; do
	   if [ -e "${filename2}" ]; then
               NewName=`echo ${filename2} | sed s/Var_scal/${filename}/`
	       mv $filename2 $NewName
	   fi
       done
       rm Var_scal
    fi
done

for filename in *_eq; do
    echo $filename
    export filename1=$filename"JK"
    if [ "$filename1" -ot "$filename" ]; then
       ln $filename ineq
       $ANNAL/cov_eq.out
       mv "equalJ"  $filename"JK"
       mv "equalJR" $filename"JR"
       for filename2 in Var_eq_Auto_Tr*; do
	   if [ -e "${filename2}" ]; then
	      NewName=`echo ${filename2} | sed s/Var_eq/${filename}/`
	      mv $filename2 $NewName
	   fi
       done
       rm ineq
    fi
done


for filename in *_tau; do
    echo $filename
    ln $filename intau
    $ANNAL/cov_tau.out
    rm intau
    mv SuscepJ $filename"JK"
    for filename1 in g_*; do
#	echo $filename1
        export Name=`echo ${filename} | sed s/_tau//`	
        export Dir=`echo ${filename1} | sed s/g_/${Name}_/`
	echo $Dir
	if [ ! -e $Dir ]; then 
           mkdir $Dir
        fi
	cd  $Dir 
	mv ../$filename1 .
	cd  ..
    done
done

