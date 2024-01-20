#!/bin/bash
shopt -s nocasematch 
ALF_DIR=/Users/fassaad/Programs/ALF_Projects/Kondo_Impurities
Analysis=1 
if [ "$Analysis" = "1" ]; then
	if [ -f "Spin_suscept.dat" ]; then
		rm "Spin_suscept.dat"
	fi
fi
Variable=(Run L_Bath Ham_U Ham_Jk  Beta  Dtaux NW NS NB CPU)
Name=(Run L_Bath_R Ham_U_R Ham_Jk_R  Beta_R Dtau_R NW_R NS_R NB_R CPU_R)
#      0   1       2       3         4      5      6     7   8    9
length=${#Variable[@]}
{
read -a Variable
echo "${Variable[@]}"
echo "${Variable[0]}"
while [ ! "${Variable[0]}" = "stop" ];   do
    if [ "${Variable[0]}" = "Y" ]; then
		string=${Variable[4]}; B_R_dir="${string//\.0/}"
		echo "B_R_dir"
	
		export Dir="L${Variable[1]}_U${Variable[2]}_J${Variable[3]}_B${B_R_dir}_Dt${Variable[5]}"
		echo "$Dir" 
		if [ "$Analysis" = "1" ]; then
			cd $Dir ||  exit
	 		"$ALF_DIR"/Analysis/ana.out *
	 		{
				Line=()
				read -a Line 
				echo "${Line[0]} ${Line[1]} ${Line[2]} ${Line[3]}" 
				echo  "${Variable[4]} ${Line[2]} ${Line[3]}"   >> ../Spin_suscept.dat
			}<SpinZ_tauJK		
			cd ..
		else
			if [ ! -e "$Dir" ]; then
	 		  	mkdir "$Dir"
	  			cd  "$Dir" || exit
	   			cp  ../Start/* . 
	    		cd ..
			else
	 		   	cd "$Dir" || exit
	  			cp ../Start/parameters .
	   	 		bash out_to_in.sh
	    		cd ..
			fi
			cd "$Dir" || exit
			((i=1))
			while [  $i -lt "$length" ]; do
	 		  	sed s/"${Name[$i]}"/"${Variable[$i]}"/    parameters  > tmp
	  			mv tmp parameters
	   	 		((i=i+1)) 
			done
			echo "$Dir"
			echo  "nohup mpirun  -np 6 /home/debian/Kondo_Impurities/Prog/ALF.out & " > run.sh  
			bash run.sh 
			cd ..
   		fi
	fi
    read -a Variable
done
}<Sims

