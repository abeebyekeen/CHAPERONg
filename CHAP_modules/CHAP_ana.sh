#! /bin/bash

#CHAP_ana - The trajectory/end-point analysis module of CHAPERONg
#CHAPERONg - An automation program for GROMACS md simulation
#Author: Abeeb A. Yekeen
#Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date: 2022.02.11

set -e
set -o pipefail

#set version
# CHAPERONg_version="beta3"

#demA=$'\n\n'"#**********************************CHAPERONg**********************************#"$'\n'
demA=$'\n\n'"#================================= CHAPERONg =================================#"$'\n'
demB=$'\n'"#=============================================================================#"$'\n\n'

if [[ $initiator2 != 'avail' ]] ; then
	echo "$demA"$'Do not run modules independently!\nLaunch CHAPERONg with run_CHAPERONg-<version>!!'"$demB"	
	exit 1
fi	

filesuffx=''

#import module with collected parameters
#. "$CHAPERONg_PATH/CHAP_modules/CHAP_colPar-4.02.sh"

#call module with defined fxns
. "$CHAPERONg_PATH/CHAP_modules/CHAP_deffxn.sh"


if [[ "$sysType" == 1 ]]; then
	sysType="protein_only"
	wraplabel="noPBC"
elif [[ "$sysType" == 2 ]]; then
	sysType="protein_lig"
	wraplabel="center"
elif [[ "$sysType" == 3 ]]; then
	sysType="protein_dna"
	wraplabel="center"
fi

#if [[ "$PBCcorrectType" != '' ]] && [[ "$PBCcorrectType" == 'noPBC' || "$PBCcorrectType" == 'nojump' || "$PBCcorrectType" == 'fit' || "$PBCcorrectType" == 'center' ]]
#then wraplabel="$PBCcorrectType"
#fi

#if [[ -f "pbcmol" ]] ; then wraplabel="noPBC" ; rm pbcmol ; fi
#if [[ -f "pbcjump" ]] ; then wraplabel="nojump" ; rm pbcjump ; fi
#if [[ -f "pbcfit" ]] ; then wraplabel="fit" ; rm pbcfit ; fi
#if [[ -f "pbccenter" ]] ; then wraplabel="center" ; rm pbccenter ; fi

#Initialize default parameters
Analysis()
{
#check the type of corrected trajectory file to use for analysis
if [[ -f "pbcmol" ]] && [[ ! -f "pbcjump" && ! -f "pbcfit" &&  ! -f "pbccenter" &&  ! -f "pbccombo" ]]
	then wraplabel="noPBC" ; rm pbcmol
elif [[ -f "pbcjump" ]] && [[ ! -f "pbcmol" && ! -f "pbcfit" &&  ! -f "pbccenter" &&  ! -f "pbccombo" ]]
	then wraplabel="nojump" ; rm pbcjump
elif [[ -f "pbcfit" ]] && [[ ! -f "pbcmol" && ! -f "pbcjump" &&  ! -f "pbccenter" &&  ! -f "pbccombo" ]]
	then wraplabel="fit" ; rm pbcfit
elif [[ -f "pbccenter" ]] && [[ ! -f "pbcmol" && ! -f "pbcjump" &&  ! -f "pbcfit" &&  ! -f "pbccombo" ]]
	then wraplabel="center" ; rm pbccenter
elif [[ -f "pbccombo" ]] && [[ ! -f "pbcmol" && ! -f "pbcjump" &&  ! -f "pbcfit" &&  ! -f "pbccenter" ]]
	then wraplabel="combo" ; rm pbccombo
fi

cat << AnalysisList

Select your choice(s) from the options listed below:
Option  Analysis
  0     Recenter, rewrap & correct molecules for pbc (trjconv)
  1     Quality assurance analyses (thermodynamic properties)
  2     Root mean square deviation (RMSD)
  3     Root mean square fluctuation (RMSF)
  4     Radius of gyration (Rg)
  5     Hydrogen bonding analysis (hbond)
  6     Solvent accessible surface area (SASA)
  7     Principal component analysis (PCA)
  8     Secondary structure analysis
  9     Clustering analysis
  10    Kernel density estimation
  11    Make a movie of the simulation
  12    Free energy calculations using the MMPBSA method (g_mmpbsa)
  13    Free energy surface (FES) with gmx sham
  14    Free energy surface using the CHAPERONg FES scripts
  15    Interactive 3D plot of the FES (using md-davis)
  16    Interactive hydrogen bond matrix (using md-davis)
  17    Extract frames from the trajectory
  18    Make index groups (make_ndx)
  19    All analyses but 16 and 17
  20    All analyses but 0, 16 and 17
  21    All analyses but 0, 9, 16 and 17
  
AnalysisList

# read -p '*Enter one or more combinations of the options here (separated by a space): ' analyse

# analysis=" $analyse "

# while [[ "$analysis" != *" 0 "* && "$analysis" != *" 1 "* && "$analysis" != *" 2 "* && \
# 	"$analysis" != *" 3 "* && "$analysis" != *" 4 "* && "$analysis" != *" 5 "* && \
# 	"$analysis" != *" 6 "* && "$analysis" != *" 7 "* && "$analysis" != *" 8 "* && \
# 	"$analysis" != *" 9 "* && "$analysis" != *" 10 "* && "$analysis" != *" 11 "* && \
# 	"$analysis" != *" 12 "* && "$analysis" != *" 13 "* && "$analysis" != *" 14 "* && \
# 	"$analysis" != *" 15 "* && "$analysis" != *" 16 "* && "$analysis" != *" 17 "* && \
# 	"$analysis" != *" 18 "* ]] ; do
# 		echo $'\nYou entered: '"$analyse"$'\n'
# 		echo $'Please enter a valid number!!\n'
# 		read -p '*Enter one or more combinations of the options here (separated by a space): ' analyse
# 		analysis=" $analyse "
# done

read -p '*Enter one or more combinations of the options here (separated by a space): ' analyse

# create a bash array listing valid numbers
valid_numbers=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19)

while ! [[ "$analyse" =~ ^([[:space:]]*[0-9][[:space:]]*)+$ ]] && \
	! [[ "$analyse" =~ (^|[[:space:]])("${valid_numbers[@]}")([[:space:]]|$) ]]
do
	echo $'\n You entered: '"$analyse"$'\n'
	echo $' Please enter a valid number!!\n'
	read -p '*Enter one or more combinations of the options here (separated by a space): ' analyse
done

analysis=" $analyse "

if [[ "$analyse" == "10" ]]; then
	analysis="$analyse"
fi

if [[ "$coordinates" == '' ]]; then
	coordinates="$filenm"
fi

ScanTRAJ()
{
if [[ ! -f "trajectDetails.log" ]]; then
	echo "$demA"$' Checking the trajectory to extract info about the number of\n frames and simulation time'"$demB"
	sleep 2
	eval "$gmx_exe_path" check -f "${filenm}"_${wraplabel}.xtc |& tee trajectDetails.log
	No_of_frames=$(cat trajectDetails.log | grep "Last" | awk '{print $(NF-2)}')
	simDuratnps=$(cat trajectDetails.log | grep "Last" | awk '{print $NF}')
	simDuratnpsINT=$(echo ${simDuratnps%\.*})
	sim_timestep=$(cat trajectDetails.log | grep -A1 "Item" | awk '{print $NF}' | tail -n 1)
	#simDuratn_nsFloat=$(echo "${simDuratnps%\.*} / 1000" | bc -l)
	simDuratn_nsFloat=$(awk "BEGIN {print $simDuratnps / 1000}")
	simDuratnINTns=$(echo ${simDuratn_nsFloat%\.*})
	echo $simDuratnINTns > simulation_duration

	echo "$demA"$' Extract number of frames and simulation duration from trajectory...DONE'"$demB"
	sleep 2
else
	No_of_frames=$(cat trajectDetails.log | grep "Last" | awk '{print $(NF-2)}')
	simDuratnps=$(cat trajectDetails.log | grep "Last" | awk '{print $NF}')
	simDuratnpsINT=$(echo ${simDuratnps%\.*})
	sim_timestep=$(cat trajectDetails.log | grep -A1 "Item" | awk '{print $NF}' | tail -n 1)
	#simDuratn_nsFloat=$(echo "${simDuratnps%\.*} / 1000" | bc -l)
	simDuratn_nsFloat=$(awk "BEGIN {print $simDuratnps / 1000}")
	simDuratnINTns=$(echo ${simDuratn_nsFloat%\.*})
	echo $simDuratnINTns > simulation_duration
fi
}

notifyImgFail()
{
echo "$demA"$'CHAPERONg could not generate a finished image file from the .xvg output.'\
$'\nConfirm that your xmgrace/gracebat is functional!'"$demB"
sleep 2
}

createDIR()
{
	currentAnadir="$(pwd)""/$AnaName"
	nDir=1
	bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
	if [[ -d "$currentAnadir" ]]; then
		base_currentAnadir=$(basename "$currentAnadir")
		base_bkupAnadir=$(basename "$bkupAnadir")
		echo $'\n'"$base_currentAnadir"$' folder exists,\n'"backing it up as $base_bkupAnadir"
		while [[ -d "$bkupAnadir" ]]; do
		nDir=$(( nDir + 1 )); bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
		done
		mv "$currentAnadir" "$bkupAnadir" && mkdir ./$AnaName
		echo $'\n'"Backing up the last $AnaName folder and its contents as $base_bkupAnadir"
		sleep 1
		mv ${filenm}*"$filesuffx".png ${filenm}*"$filesuffx".xvg ${filenm}_"$filesuffx".png ./$AnaName || true
		mv ${filenm}_"$filesuffx".xvg ${filenm}_*"$filesuffx".png ./$AnaName || true
		mv ${filenm}_*"$filesuffx".xvg ${filenm}*"$filesuffx"*.xvg ${filenm}*"$filesuffx"*.png ./$AnaName || true
		mv *"$filesuffx"*.png *"$filesuffx".png *"$filesuffx".xvg *"$filesuffx"*.xvg ./$AnaName || true
		mv "${filenm}"*"$filesuffx"*".png" "${filenm}"*"$filesuffx"".png" "${filenm}"*"$filesuffx"".xvg" ./$AnaName || true
		mv ${filenm}_"$filesuffx"*.png ${filenm}_"$filesuffx"*.xvg ./$AnaName || true
		
	elif [[ ! -d "$currentAnadir" ]]; then mkdir ./$AnaName
		mv ${filenm}*"$filesuffx".png ${filenm}*"$filesuffx".xvg ${filenm}_"$filesuffx".png ./$AnaName || true
		mv ${filenm}_"$filesuffx".xvg ${filenm}_*"$filesuffx".png ./$AnaName || true
		mv ${filenm}_*"$filesuffx".xvg ${filenm}*"$filesuffx"*.xvg ${filenm}*"$filesuffx"*.png ./$AnaName || true
		mv *"$filesuffx"*.png *"$filesuffx".png *"$filesuffx".xvg *"$filesuffx"*.xvg ./$AnaName || true
		mv "${filenm}"*"$filesuffx"*".png" "${filenm}"*"$filesuffx"".png" "${filenm}"*"$filesuffx"".xvg" ./$AnaName || true
		mv ${filenm}_"$filesuffx"*.png ${filenm}_"$filesuffx"*.xvg ./$AnaName || true
		#mv ${filenm}_Rg_ns.png ${filenm}_Rg_ns.xvg ${filenm}_Rg.xvg ./$AnaName || true
	fi
}
		
DNAwrapAlt()
{
echo "$demA""CHAPERONg could not find any index file. Centering on protein instead of Protein_DNA!""$demB"
sleep 2
echo "Protein" 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_${wraplabel}.xtc -center -pbc mol -ur compact
}
analyser0()
{	
echo "$demA"$' Now recentering the protein and rewrapping molecules within the unit cell...\n'
if [[ $automode == "full" && $sysType == "protein_only" ]]; then
	if [[ "$PBCcorrectType" != '' && "$wraplabel" == 'noPBC' ]] ; then
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"noPBC".xtc -pbc mol -center
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'nojump' ]] ; then
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'combo' ]] ; then
		echo 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump
		echo 4 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_"nojump".xtc -o "${filenm}"_"nojump_fitTrans".xtc -fit translation
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_"nojump_fitTrans".xtc -o "${filenm}"_"combo".xtc -pbc mol -center
		rm "${filenm}"_"nojump".xtc "${filenm}"_"nojump_fitTrans".xtc
	else
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"noPBC".xtc -pbc mol -center
		echo "$demA"$' Now removing possible jumps in the trajectory...\n'
		sleep 1
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
	fi
				
elif [[ $automode != "full" && $sysType == "protein_only" ]]; then
	if [[ "$PBCcorrectType" != '' && "$wraplabel" == 'noPBC' ]] ; then
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		sleep 1
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"noPBC".xtc -pbc mol -center
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'nojump' ]] ; then
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		sleep 1
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'combo' ]] ; then
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		sleep 2
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
		echo $'**Choose (4) for centering and "System" (0) for output when prompted\n'
		sleep 2
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_"nojump".xtc -o "${filenm}"_"nojump_fitTrans".xtc -fit translation
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		sleep 2
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_"nojump_fitTrans".xtc -o "${filenm}"_"combo".xtc -pbc mol -center
	else
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"noPBC".xtc -pbc mol -center
		echo "$demA"$' Now removing possible jumps in the trajectory...\n'
		sleep 1
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		sleep 1
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
	fi
		
elif [[ $automode == "full" && $sysType == "protein_lig" ]]; then
	if [[ "$PBCcorrectType" != '' && "$wraplabel" == 'center' ]] ; then
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'fit' ]] ; then
		echo 4 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"fit".xtc -fit rot+trans
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'nojump' ]] ; then
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'combo' ]] ; then
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump
		echo 4 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_"nojump".xtc -o "${filenm}"_"nojump_fitTrans".xtc -fit translation
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_"nojump_fitTrans".xtc -o "${filenm}"_"combo".xtc -center -pbc mol -ur compact
		rm "${filenm}"_"nojump".xtc "${filenm}"_"nojump_fitTrans".xtc
	else
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
		echo "$demA"$' Now performing rotational and translational fitting...\n'
		echo 4 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o "${filenm}"_fit.xtc -fit rot+trans
		echo "$demA"$' Now removing possible jumps in the trajectory...\n'
		sleep 1
		echo 1 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact
	fi

elif [[ $automode != "full" && $sysType == "protein_lig" ]]; then
	if [[ "$PBCcorrectType" != '' && "$wraplabel" == 'center' ]] ; then
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		sleep 1
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'nojump' ]] ; then
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		sleep 1
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'fit' ]] ; then	
		echo $'**Choose (4) for centering and "System" (0) for output when prompted\n'
		sleep 2
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"fit".xtc -fit rot+trans
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'combo' ]] ; then
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		sleep 2
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
		echo $'**Choose (4) for centering and "System" (0) for output when prompted\n'
		sleep 2
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_"nojump".xtc -o "${filenm}"_"nojump_fitTrans".xtc -fit translation
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		sleep 2
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_"nojump_fitTrans".xtc -o "${filenm}"_"combo".xtc -center -pbc mol -ur compact
	else
		echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
		echo "$demA"$' Now performing rotational and translational fitting...\n'
		echo $'**Choose "Backbone" (4) to perform lsq fitting to protein backbone, and "System" (0) for output when prompted\n'
		eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o "${filenm}"_fit.xtc -fit rot+trans
	fi

elif [[ $automode == "full" && $sysType == "protein_dna" ]]; then
	if [[ "$PBCcorrectType" != '' && "$wraplabel" == 'center' ]] ; then
		echo "Protein_DNA" 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"center".xtc -center -pbc mol -ur compact || DNAwrapAlt
	elif [[ "$PBCcorrectType" != '' && "$wraplabel" == 'nojump' ]] ; then
		echo "Protein_DNA" 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact || DNAwrapAlt
	else
		echo "Protein_DNA" 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"center".xtc -center -pbc mol -ur compact || DNAwrapAlt
		echo "$demA"$' Now removing possible jumps in the trajectory...\n'
		sleep 1
		echo "Protein_DNA" 0 | eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact || DNAwrapAlt
	fi
elif [[ $automode != "full" ]] && [[ $sysType == "protein_dna" ]]; then
	echo $'**Choose "Protein_DNA" for centering and "System" (0) for output when prompted\n'
	eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	eval "$gmx_exe_path" trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact || DNAwrapAlt
fi
echo "$demA"$' Recenter the protein and rewrap molecules within the unit cell...DONE'"$demB"
sleep 2
}
if [[ "$analysis" == *" 0 "* ]]; then analyser0; fi

analyser1()
{
	echo "$demA"$' Now calculating post-MD thermodynamic parameters...\n\n'
	sleep 2

	echo "Temperature" | eval "$gmx_exe_path" energy -f "${filenm}".edr -o postMD_Temperature.xvg
	gracebat postMD_Temperature.xvg -hdevice PNG -autoscale xy -printfile postMD_Temperature.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "$demA"$' Calculate Temperature progression...DONE\n\n' ; sleep 2

	echo "Pressure" | eval "$gmx_exe_path" energy -f "${filenm}".edr -o postMD_Pressure.xvg
	gracebat postMD_Pressure.xvg -hdevice PNG -autoscale xy -printfile postMD_Pressure.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "$demA"$' Calculate Pressure progression...DONE\n\n' ; sleep 2

	echo "Density" | eval "$gmx_exe_path" energy -f "${filenm}".edr -o postMD_Density.xvg
	gracebat postMD_Density.xvg -hdevice PNG -autoscale xy -printfile postMD_Density.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "$demA"$' Calculate Density progression...DONE\n\n' ; sleep 2

	echo "Total-Energy" | eval "$gmx_exe_path" energy -f "${filenm}".edr -o postMD_TotalEnergy.xvg
	gracebat postMD_TotalEnergy.xvg -hdevice PNG -autoscale xy -printfile postMD_TotalEnergy.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "$demA"$' Calculate Total energy...DONE\n\n' ; sleep 2

	echo "Potential" | eval "$gmx_exe_path" energy -f "${filenm}".edr -o postMD_Potential.xvg
	gracebat postMD_Potential.xvg -hdevice PNG -autoscale xy -printfile postMD_Potential.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "$demA"$' Calculate Potential energy...DONE\n\n' ; sleep 2

	echo "Kinetic-En." | eval "$gmx_exe_path" energy -f "${filenm}".edr -o postMD_KineticEn.xvg
	gracebat postMD_KineticEn.xvg -hdevice PNG -autoscale xy -printfile postMD_KineticEn.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "$demA"$' Calculate Kinetic energy...DONE\n\n' ; sleep 2

	currentMDthermodyndir="$(pwd)""/postMD_thermodynamics"
	nMDtherm=1
	bkupMDtherm="$(pwd)""/#postMD_thermodynamics"".""backup.""$nMDtherm"
	base_bkupMDtherm=$(basename "$bkupMDtherm")
	if [[ -d "$currentMDthermodyndir" ]]; then
		echo $'\n'"$currentMDthermodyndir"$' folder exists,\n'"backing it up as $base_bkupMDtherm"
		sleep 1
		while [[ -d "$bkupMDtherm" ]]; do
			nMDtherm=$(( nMDtherm + 1 ))
			bkupMDtherm="$(pwd)""/#postMD_thermodynamics"".""backup.""$nMDtherm"
			base_bkupMDtherm=$(basename "$bkupMDtherm")
		done
		mv "$currentMDthermodyndir" "$bkupMDtherm" && mkdir ./postMD_thermodynamics || true
		echo $'\n'"Backing up the last postMD_thermodynamics folder and its contents as $base_bkupMDtherm"
		sleep 1
	elif [[ ! -d "$currentMDthermodyndir" ]]; then mkdir postMD_thermodynamics
	fi
	gracebat postMD_Potential.xvg postMD_KineticEn.xvg postMD_TotalEnergy.xvg -hdevice PNG \
	-autoscale xy -printfile postMD_Energies.png -fixed 7500 4000 -legend load || notifyImgFail
	mv postMD_Temperature.xvg postMD_KineticEn.xvg postMD_Potential.xvg ./postMD_thermodynamics || true
	mv postMD_Density.xvg postMD_Pressure.xvg postMD_TotalEnergy.xvg ./postMD_thermodynamics || true
	mv postMD_Temperature.png postMD_KineticEn.png postMD_Potential.png postMD_Density.png ./postMD_thermodynamics || true
	mv postMD_Pressure.png postMD_TotalEnergy.png postMD_Energies.png ./postMD_thermodynamics || true

	echo "$demA"$' Calculate post-MD thermodynamics parameters...DONE'"$demB"
	sleep 2

}
if [[ "$analysis" == *" 1 "* ]]; then analyser1; fi

altRMSD()
{
	echo "$demA"$'There are multiple groups identified as '"$ligname"\
		$'.\nCHAPERONg will try to guess the appropriate group to be used for '"$ligname"" RMSD calculations""$demB"
	sleep 2
	echo "$demA""CHAPERONg: Selecting group 13 for ""$ligname"\
		$'.\nIf this is wrong, terminate and re-run RMSD analysis without the automation flag!'"$demB"
	sleep 2

	echo 13 13 | eval "$gmx_exe_path" rms -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o ${filenm}_"$ligname"-rmsd.xvg -tu ns
}

analyser2()
{
	echo "$demA"$' Now calculating RMSD...\n'
	sleep 2
	if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]] && [[ $automode == "full" ]] ; then
		echo "Backbone" "Backbone" | eval "$gmx_exe_path" rms -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o ${filenm}_BB-rmsd.xvg -tu ns
			
		gracebat ${filenm}_BB-rmsd.xvg -hdevice PNG -autoscale xy -printfile ${filenm}_BB-rmsd.png \
		-fixed 7500 4000 -legend load || notifyImgFail
		
	elif [[ $automode == "full" ]] && [[ $sysType == "protein_lig" ]]; then
		echo 4 4 | eval "$gmx_exe_path" rms -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o ${filenm}_BB-rmsd.xvg -tu ns
			
		gracebat ${filenm}_BB-rmsd.xvg -hdevice PNG -autoscale xy -printfile ${filenm}_BB-rmsd.png \
		-fixed 7500 4000 -legend load || notifyImgFail
		echo "$demA"$'Protein RMSD calculation... DONE\n  Now calculating ligand RMSD...\n'
		sleep 2
			
		echo "$ligname" "$ligname" | eval "$gmx_exe_path" rms -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o \
		"$filenm"_"$ligname"-rmsd.xvg -tu ns || altRMSD
					
		gracebat "$filenm"_"$ligname"-rmsd.xvg -hdevice PNG -autoscale xy -printfile \
		"$filenm"_"$ligname"-rmsd.png -fixed 7500 4000 -legend load || notifyImgFail
			
	else
		eval "$gmx_exe_path" rms -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o ${filenm}_rmsd.xvg -tu ns
			
		gracebat ${filenm}_rmsd.xvg -hdevice PNG -autoscale xy -printfile ${filenm}_rmsd.png \
		-fixed 7500 4000 -legend load || notifyImgFail
	fi
	echo "$demA"$' Compute RMSD... DONE'"$demB"
	sleep 2
	AnaName="RMSD"
	filesuffx="rmsd"
	createDIR
	echo "$demA"$' Generate a finished figure of the RMSD plot... DONE'"$demB"
	sleep 2
}

if [[ "$analysis" == *" 2 "* ]]; then analyser2 ; fi

analyser3()
{
echo "$demA"$' Now calculating RMSF...\n'
if [[ $automode == "full" ]]; then
	echo "Backbone" "Backbone" | eval "$gmx_exe_path" rmsf -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o ${filenm}_BB-rmsf.xvg -res
	echo "$demA"$' RMSF with backbone lsq fitting and calculation...DONE'"$demB"
	sleep 2
	echo "C-alpha" "C-alpha" | eval "$gmx_exe_path" rmsf -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o ${filenm}_Calpha-rmsf.xvg -res
	echo "$demA"$' RMSF with Calpha lsq fitting and calculation...DONE'"$demB"
	sleep 2
	gracebat ${filenm}_BB-rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${filenm}_BB-rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat ${filenm}_Calpha-rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${filenm}_Calpha-rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat ${filenm}_BB-rmsf.xvg ${filenm}_Calpha-rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${filenm}_BB-Calpha-rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
else
	eval "$gmx_exe_path" rmsf -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o ${filenm}_rmsf.xvg -res
	echo "$demA"$' Compute RMSF...DONE'"$demB"
	sleep 2
	gracebat ${filenm}_rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${filenm}_rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
fi
	
AnaName="RMSF"
filesuffx="rmsf"
createDIR
	
echo "$demA"$' Generate finished figure(s) of the RMSF plot(s)... DONE'"$demB"
sleep 2
}
if [[ "$analysis" == *" 3 "* ]]; then analyser3 ; fi
	
analyser4()
{
echo "$demA"$' Now calculating Rg...\n'
if [[ $automode == "full" ]]; then
	echo "Protein" | eval "$gmx_exe_path" gyrate -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o ${filenm}_Rg.xvg
	echo "$demA"$' Compute radius of gyration...DONE'"$demB"
	sleep 2
	echo "$demA"$' Now converting Rg plot to ns format...\n'
	sleep 2
	grep "^[@#]" ${filenm}_Rg.xvg | sed "s/ps/ns/g" > ${filenm}_Rg_ns.xvg
	grep -v "^[@#]" ${filenm}_Rg.xvg | \
	awk '{print $1/1000"      "$2"      "$3"      "$4"     "$5}' >> ${filenm}_Rg_ns.xvg
else
	echo $'**In the following step, CHOOSE Protein (1) for Rg analysis\n\n'
	eval "$gmx_exe_path" gyrate -s "${filenm}".tpr -f "${filenm}"_${wraplabel}.xtc -o ${filenm}_Rg.xvg
	echo "$demA"$' Compute radius of gyration...DONE'"$demB"
	sleep 2
	echo "$demA"$' Now converting Rg plot to ns format...\n'
	sleep 2
	grep "^[@#]" ${filenm}_Rg.xvg | sed "s/ps/ns/g" > ${filenm}_Rg_ns.xvg
	grep -v "^[@#]" ${filenm}_Rg.xvg | \
	awk '{print $1/1000"      "$2"      "$3"      "$4"     "$5}' >> ${filenm}_Rg_ns.xvg
fi
echo "$demA"$' Generate ns and ps Rg plots...DONE'"$demB"
sleep 2
gracebat ${filenm}_Rg_ns.xvg -hdevice PNG -autoscale xy -printfile ${filenm}_Rg_ns.png -fixed 7500 4000 -legend load || notifyImgFail
	
AnaName="Rg"
filesuffx="Rg"
createDIR
echo "$demA"$' Generate a finished figure of the Rg plot... DONE'"$demB"
sleep 2
}

if [[ "$analysis" == *" 4 "* ]]; then analyser4 ; fi

altHBOND()
{
echo "$demA""There are multiple groups identified as ""$ligname""."\
$'\nCHAPERONg will try to guess the appropriate group to be used for protein-'"$ligname"" hbond calculations""$demB"

sleep 2

echo "$demA""CHAPERONg: Selecting group 13 for ""$ligname""."\
$'\nIf this is wrong, terminate and re-run hbond analysis without the automation flag!'"$demB"

sleep 2
echo 1 13 | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num hbnum_ProLig_${filenm}.xvg \
	-hbm hb_matrix_ProLig_${filenm}.xpm -hbn hb_index_ProLig_${filenm}.ndx -tu ns $hbthread
}

hbond_DNA1()
{
	echo "$demA"$' Now executing Intra-protein hydrogen bonding analysis...\n'
	sleep 2

	echo "Protein" "Protein" | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s ${filenm}.tpr -n index.ndx \
		-num hbnum_Pro_${filenm}.xvg -hbm hb_matrix_Pro_${filenm}.xpm -hbn hb_index_Pro_${filenm}.ndx -tu ns $hbthread

	echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
	sleep 2

	echo "$demA"$' Now executing Intra-DNA hydrogen bonding analysis...\n'
	sleep 2

	echo "DNA" "DNA" | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx \
		-num hbnum_DNA_${filenm}.xvg -hbm hb_matrix_DNA_${filenm}.xpm -hbn hb_index_DNA_${filenm}.ndx -tu ns $hbthread

	echo "$demA"$' Intra-DNA hydrogen bonding analysis...DONE'"$demB"
	sleep 2

	echo "$demA"$' Now executing Protein-DNA hydrogen bonding analysis...\n'
	echo "Protein" "DNA" | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -num \
		hbnum_Pro_DNA_${filenm}.xvg -hbm hb_matrix_Pro_DNA_${filenm}.xpm -hbn hb_index_Pro_DNA_${filenm}.ndx -tu ns $hbthread

	echo "$demA"$' Protein-DNA hydrogen bonding analysis... DONE'"$demB"
	sleep 2
}
hbond_DNA2()
{
echo "$demA"$' Now executing Intra-protein hydrogen bonding analysis...\n'
sleep 2

echo "Protein" "Protein" | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num \
	hbnum_Pro_${filenm}.xvg -hbm hb_matrix_Pro_${filenm}.xpm -hbn hb_index_Pro_${filenm}.ndx -tu ns $hbthread
gracebat hbnum_Pro_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
hbnum_Pro_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Intra-DNA hydrogen bonding analysis...\n'
sleep 2

echo "DNA" "DNA" | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num \
	hbnum_DNA_${filenm}.xvg -hbm hb_matrix_DNA_${filenm}.xpm -hbn hb_index_DNA_${filenm}.ndx -tu ns $hbthread
gracebat hbnum_DNA_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
hbnum_DNA_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
echo "$demA"$' Intra-DNA hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Protein-DNA hydrogen bonding analysis...\n'
echo "Protein" "DNA" | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num hbnum_Pro_DNA_${filenm}.xvg \
	-hbm hb_matrix_Pro_DNA_${filenm}.xpm -hbn hb_index_Pro_DNA_${filenm}.ndx -tu ns $hbthread
gracebat hbnum_Pro_DNA_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
hbnum_Pro_DNA_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
echo "$demA"$' Protein-DNA hydrogen bonding analysis... DONE'"$demB"
sleep 2
}

analyser5()
{
echo "$demA"$' Now executing H-bond analysis...\n'
if [[ $automode == "full" && $sysType == "protein_only" ]]; then
	echo 1 1 | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num hbnum_intraPro_${filenm}.xvg \
	-hbm hb_matrix_intraPro_${filenm}.xpm -hbn hb_index_intraPro_${filenm}.ndx -tu ns $hbthread
	echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	echo 1 "SOL" | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num hbnum_Pro-SOL_${filenm}.xvg \
	-hbm hb_matrix_Pro-SOL_${filenm}.xpm -hbn hb_index_Pro-SOL_${filenm}.ndx -tu ns $hbthread || \
	echo "$demA"$' There are multiple groups with the name SOL. Skipping...'
	echo "$demA"$' Protein-SOL hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_intraPro_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_intraPro_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat hbnum_Pro-SOL_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_Pro-SOL_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat hbnum_intraPro_${filenm}.xvg hbnum_Pro-SOL_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_intraPro_Pro-SOL_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail

elif [[ $automode == "full" && $sysType == "protein_lig" ]]; then
	echo 1 "$ligname" | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num hbnum_ProLig_${filenm}.xvg \
	-hbm hb_matrix_ProLig_${filenm}.xpm -hbn hb_index_ProLig_${filenm}.ndx -tu ns $hbthread || altHBOND
	echo "$demA"$' Protein-ligand hydrogen bonding analysis...DONE'"$demB"
	sleep 2

	echo 1 1 | eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num hbnum_intraPro_${filenm}.xvg \
	-hbm hb_matrix_intraPro_${filenm}.xpm -hbn hb_index_intraPro_${filenm}.ndx -tu ns $hbthread
	echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_ProLig_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_ProLig_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat hbnum_intraPro_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_intraPro_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == "protein_only" && $automode == "semi" ]] ; then
	eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num hbnum_${filenm}.xvg \
	-hbm hb_matrix_${filenm}.xpm -hbn hb_index_${filenm}.ndx -tu ns $hbthread
	echo "$demA"$' Hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail	
elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode == "semi" ]] ; then
	eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -num hbnum_${filenm}.xvg \
	-hbm hb_matrix_${filenm}.xpm -hbn hb_index_${filenm}.ndx -tu ns $hbthread || \
	eval "$gmx_exe_path" hbond -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -num hbnum_${filenm}.xvg \
	-hbm hb_matrix_${filenm}.xpm -hbn hb_index_${filenm}.ndx -tu ns $hbthread
	echo "$demA"$' Hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == "protein_dna" ]] && [[ $automode == "full" ]] ; then
	hbond_DNA1 || hbond_DNA2
	echo "$demA"$' Hydrogen bonding analysis...DONE'"$demB"
	sleep 2
fi
	
currenthbonddir="$(pwd)""/hbond"
nhbond=1
bkuphbonddir="$(pwd)""/#hbond"".""backup.""$nhbond"
base_bkuphbonddir=$(basename "$bkuphbonddir")
if [[ -d "$currenthbonddir" ]]; then
	echo $'\n'"$currenthbonddir"$' folder exists,\n'"backing it up as $base_bkuphbonddir"
	sleep 1
	while [[ -d "$bkuphbonddir" ]]; do
		nhbond=$(( nhbond + 1 ))
		bkuphbonddir="$(pwd)""/#hbond"".""backup.""$nhbond"
		base_bkuphbonddir=$(basename "$bkuphbonddir")
	done
	mv "$currenthbonddir" "$bkuphbonddir" && mkdir ./hbond || true
	echo $'\n'"Backing up the last hbond folder and its contents as $base_bkuphbonddir"
	sleep 1
	mv hbnum_*.png hbnum_*.xvg hbnum*.png hb_matrix_* hb_index_* ./hbond || true
elif [[ ! -d "$currenthbonddir" ]]; then
	mkdir ./hbond; mv hbnum_*.png hbnum_*.xvg hbnum*.png hb_matrix_* hb_index_* ./hbond || true
fi
echo "$demA"$' Generate finished figure(s) of the hbond plot(s)... DONE'"$demB"
}

if [[ "$analysis" == *" 5 "* ]]; then analyser5 ; fi

analyser6()
{
	echo "$demA"$' Now calculating solvent accessible surface area (SASA)...\n'
	if [[ $automode == "full" && $sysType == "protein_only" ]]; then
		echo 1 | eval "$gmx_exe_path" sasa -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o sasa_${filenm}.xvg -tu ns
		gracebat sasa_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
		sasa_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail		
	elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode == "semi" ]] ; then
		eval "$gmx_exe_path" sasa -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -o sasa_${filenm}.xvg -tu ns
		gracebat sasa_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
		sasa_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
	elif [[ $sysType == "protein_lig" && $automode == "full" ]]; then
		echo 1 | eval "$gmx_exe_path" sasa -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -o \
		sasa_Pro_${filenm}.xvg -tu ns
		gracebat sasa_Pro_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
		sasa_Pro_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
	elif [[ $sysType == "protein_dna" && $automode == "full" ]]; then
		echo "Protein" | eval "$gmx_exe_path" sasa -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -o \
		sasa_Pro_${filenm}.xvg -tu ns
		gracebat sasa_Pro_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
		sasa_Pro_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail	
		echo "$demA"$'Compute solvent accessible surface area (SASA) for DNA only...DONE'"$demB"	
		echo "DNA" | eval "$gmx_exe_path" sasa -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -o \
		sasa_DNA_${filenm}.xvg -tu ns
		gracebat sasa_Pro_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
		sasa_DNA_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail	
		echo "$demA"$'Compute solvent accessible surface area (SASA) for DNA only...DONE'"$demB"	
		echo "$demA"$'Now calculating solvent accessible surface area (SASA) for Protein-DNA complex...\n'	
		echo "Protein_DNA" | eval "$gmx_exe_path" sasa -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -o \
		sasa_Pro_DNA_${filenm}.xvg -tu ns
		gracebat sasa_Pro_DNA_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
		sasa_Pro_DNA_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail	
		echo "$demA"$'Now calculating solvent accessible surface area (SASA) for Protein-DNA complex...DONE'"$demB"	
	elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
		eval "$gmx_exe_path" sasa -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o sasa_${filenm}.xvg -tu ns
		gracebat sasa_${filenm}.xvg -hdevice PNG -autoscale xy -printfile \
		sasa_${filenm}.png -fixed 7500 4000 -legend load || notifyImgFail
	fi
	echo "$demA"$' Compute solvent accessible surface area (SASA)...DONE'"$demB"
	sleep 2
	currentSASAdir="$(pwd)""/SASA"
	nSASA=1
	bkupSASAdir="$(pwd)""/#SASA"".""backup.""$nSASA"
	base_bkupSASAdir=$(basename "$bkupSASAdir")
	if [[ -d "$currentSASAdir" ]]; then
		echo $'\n'"$currentSASAdir"$' folder exists,\n'"backing it up as $base_bkupSASAdir"
		sleep 1
		while [[ -d "$bkupSASAdir" ]]; do
			nSASA=$(( nSASA + 1 ))
			bkupSASAdir="$(pwd)""/#SASA"".""backup.""$nSASA"
			base_bkupSASAdir=$(basename "$bkupSASAdir")
		done
		mv "$currentSASAdir" "$bkupSASAdir" && mkdir ./SASA || true
		echo $'\n'"Backing up the last SASA folder and its contents as $base_bkupSASAdir"
		sleep 1
	elif [[ ! -d "$currentSASAdir" ]]; then mkdir ./SASA
	fi
	mv sasa*${filenm}.png sasa*${filenm}.xvg ./SASA || true
	echo "$demA"$' Generate a finished figure of the SASA plot... DONE'"$demB"
	sleep 2
}

if [[ "$analysis" == *" 6 "* ]]; then analyser6 ; fi

analyser7()
{
	echo "$demA"$' Now running principal component analysis (PCA)...\n'
	if [[ $automode == "full" && $sysType == "protein_only" ]]; then
		echo 4 4 | eval "$gmx_exe_path" covar -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
		-o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
		echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
		echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'
		echo 4 4 | eval "$gmx_exe_path" anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_${wraplabel}.xtc -eig \
		"${filenm}"_eigenval.xvg -s "${filenm}".tpr -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg
			
		gracebat PCA_2dproj_"${filenm}".xvg -hdevice PNG -autoscale xy -printfile \
		PCA_2dproj_"${filenm}".png -fixed 7500 4000 -legend load || notifyImgFail
	elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
		echo $'**Choose "Backbone" (4) twice when prompted\n'
		eval "$gmx_exe_path" covar -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
		echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
		echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'\
		$'**Choose "Backbone" (4) twice when prompted\n'
		eval "$gmx_exe_path" anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_${wraplabel}.xtc -eig "${filenm}"_eigenval.xvg \
		-s "${filenm}".tpr -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg
		gracebat PCA_2dproj_"${filenm}".xvg -hdevice PNG -autoscale xy -printfile \
		PCA_2dproj_"${filenm}".png -fixed 7500 4000 -legend load || notifyImgFail
	elif [[ $sysType == "protein_lig" ]] || [[ $sysType == "protein_dna" ]] && [[ $automode == "full" ]] ; then
		echo "Backbone" "Backbone" | eval "$gmx_exe_path" covar -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx \
		-o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
		echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
		echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'
		echo "Backbone" "Backbone" | eval "$gmx_exe_path" anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_${wraplabel}.xtc -eig \
		"${filenm}"_eigenval.xvg -s "${filenm}".tpr -n index.ndx -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg
		gracebat PCA_2dproj_"${filenm}".xvg -hdevice PNG -autoscale xy -printfile \
		PCA_2dproj_"${filenm}".png -fixed 7500 4000 -legend load || notifyImgFail
	elif [[ $sysType == "protein_lig" ]] || [[ $sysType == "protein_dna" ]] && [[ $automode == "semi" ]] ; then
		echo $'**Choose "Backbone" (4) twice when prompted\n'
		eval "$gmx_exe_path" covar -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
		-n index.ndx -o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
		echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
		echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'\
		$'**Choose "Backbone" (4) twice when prompted\n'
		eval "$gmx_exe_path" anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_${wraplabel}.xtc -eig \
		"${filenm}"_eigenval.xvg -s "${filenm}".tpr -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg	
	fi
	echo "$demA"$' Principal component analysis (PCA)...DONE'"$demB"
	sleep 2
	currentPCAdir="$(pwd)""/PCA"
	nPCA=1
	bkupPCAdir="$(pwd)""/#PCA"".""backup.""$nPCA"
	base_bkupPCAdir=$(basename "$bkupPCAdir")
	if [[ -d "$currentPCAdir" ]]; then
		echo $'\n'"$currentPCAdir"$' folder exists,\n'"backing it up as $base_bkupPCAdir"
		sleep 1
		while [[ -d "$bkupPCAdir" ]]; do
			nPCA=$(( nPCA + 1 )); bkupPCAdir="$(pwd)""/#PCA"".""backup.""$nPCA"
			base_bkupPCAdir=$(basename "$bkupPCAdir")
		done
		mv "$currentPCAdir" "$bkupPCAdir" && mkdir ./PCA
		echo $'\n'"Backing up the last PCA folder and its contents as $base_bkupPCAdir"
		sleep 1
	elif [[ ! -d "$currentPCAdir" ]]; then mkdir ./PCA
	fi
	mv PCA_2dproj_*.png *eigenval.xvg PCA_2dproj_*.xvg *_eigenvec.trr covar.log average.pdb dd?????? ./PCA || true
	echo "$demA"$' Generate finished figures of the PCA plots... DONE'"$demB"
	sleep 2
}

if [[ "$analysis" == *" 7 "* ]]; then analyser7 ; fi
	
dsspCheck="Avail"
DSSPfail()
{
	echo "$demA"$' DSSP could not be configured for gromacs!\n Skipping secondary structure analysis...\n'
	sleep 2
	dsspCheck="notAvail"
}
useCHAPdssp()
{
	echo "$demA"$' DSSP not detected on your machine.'\
	$'\nDo you want use the DSSP executable packaged with CHAPERONg?\n'
	read -p ' Enter a response here (yes or no): ' configDSSP
		while [[ "$configDSSP" != "yes" && "$movieLeng" != "y" \
		&& "$movieLeng" != "no" && "$movieLeng" != "n" ]]; do
			echo $'\nYou entered: '"$configDSSP"
			echo $'Please enter a valid response (yes/no or y/n)!!\n'
			read -p ' Enter a response here (y/n): ' configDSSP
		done
	if [[ $configDSSP == "yes" || $configDSSP == "y" ]] ; then
		export DSSP="$(echo $CHAPERONg_PATH)/CHAP_utilities/dssp-x64"
		alias DSSP="$(echo $CHAPERONg_PATH)/CHAP_utilities/dssp-x64"
	elif [[ $configDSSP == "no" || $configDSSP == "n" ]] ; then echo ""
	fi
}

useCHAPdsspGMX()
{
	echo "$demA"$' DSSP not detected on your machine.'\
	$'Either it is not installed or the\n environment has not been set (properly).'
	sleep 2
	echo $'\n Do you want configure the DSSP executable packaged with'\
	$'CHAPERONg for use\n by GMX?\n'
	sleep 2
	read -p ' Enter a response here (yes or no): ' configDSSP
		while [[ "$configDSSP" != "yes" && "$configDSSP" != "y" \
		&& "$configDSSP" != "no" && "$configDSSP" != "n" ]]; do
			echo $'\nYou entered: '"$configDSSP"
			echo $'Please enter a valid response (yes/no or y/n)!!\n'
			read -p ' Enter a response here (y/n): ' configDSSP
		done
	if [[ $configDSSP == "yes" || $configDSSP == "y" ]] ; then
		export DSSP="$(echo $CHAPERONg_PATH)/CHAP_utilities/dssp-x64"
		alias DSSP="$(echo $CHAPERONg_PATH)/CHAP_utilities/dssp-x64"

		echo "$demA"$' DSSP configured...\n'
		sleep 1
		echo $' Now attempting to run secondary analysis again...\n\n'
		sleep 2
		echo "MainChain" | eval "$gmx_exe_path" do_dssp -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
		-o ss_"${filenm}".xpm -tu ns -dt ${dt_dssp} || DSSPfail

	elif [[ $configDSSP == "no" || $configDSSP == "n" ]] ; then DSSPfail
	fi
}

analyser8()
{
	# echo "$demA"$' Checking DSSP availability and configuration...\n'
	# sleep 2
	# catchDSSP1error=''
	# DSSP -h &> tempdssp.temp1 || true
	# catchDSSP1error=$(cat tempdssp.temp1 | grep "not found")
	# cat tempdssp.temp1
	# echo $catchDSSP1error
	# sleep 3
	# cat tempdssp.temp1 | grep "not found"
	# sleep 3
	# echo "Pass 1"
	# #rm tempdssp.temp1

	# if [[ $catchDSSP1error == *"not found" ]]; then
	# 	echo "Pass 2"
	# 	catchdssp2error=''
	# 	dssp -h &> tempdssp.temp2 || true
	# 	catchdssp2error=$(cat tempdssp.temp2 | grep "not found")
	# 	cat tempdssp.temp2
	# 	echo $catchdssp2error
	# 	sleep 3
	# 	cat tempdssp.temp2 | grep "not found"
	# 	sleep 3
	# 	#rm tempdssp.temp2
	# 	echo "Pass 3"
	# 	if [[ $catchdssp2error == *"not found" ]]; then
	# 		echo "Pass 4"
	# 		useCHAPdssp
	# 		echo "Pass 5"
	# 	fi	
	# fi

	echo "$demA"$' Now computing secondary structure with DSSP...\n'
	sleep 2
	#CollectDSSPdt()

	#echo $'Enter a number below to set the time interval for frames to be used (dt).'\
	#$'\nCHAPERONg recommends 0.1 for 100ns mds, 0.2 for 200ns etc.\nYou can also enter 0 instead, for gmx default estimation.\n'

	#read -p 'Value of dt: ' dt_dssp
	#echo $'\nYou entered: '"$dt_dssp"$'\n'

	#sleep 2

	#if [[ $simDuratnINTns == 1000 ]]; then dt_dssp=
	#elif [[ $simDuratnINTns == 500 ]]; then dt_dssp="0.5"
	#elif [[ $simDuratnINTns == 200 ]]; then dt_dssp="0.2"
	#fi

	#dt_dssp=$(echo "scale=3; ${simDuratnINTns} / 1000" | bc -l)

	dt_dssp_alt()
	{
		simDuratnINTns=$(cat simulation_duration)
		dt_dssp=$(awk "BEGIN {print $simDuratnINTns / 1000}")
	}

	dt_dssp=$(awk "BEGIN {print $simDuratnINTns / 1000}") || dt_dssp_alt

	echo "MainChain" | eval "$gmx_exe_path" do_dssp -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
	-o ss_"${filenm}".xpm -tu ns -dt ${dt_dssp} || useCHAPdsspGMX
	if [[ "$dsspCheck" == "Avail" ]] ; then
		echo "$demA"$' Compute secondary structure...DONE'"$demB"
		sleep 1
		echo "$demA"$' Detecting colour codes assigned in the eps file...\n'
		sleep 2
		while IFS= read -r line; do
			scanSSname=$(echo "$line" | awk '{print $6}')
			if [[ "$scanSSname" == '"Coil"' ]] ; then coilCODEcoil=$(echo "$line" | awk '{print $3}')
			elif [[ "$scanSSname" == '"B-Sheet"' ]] ; then coilCODEbeta=$(echo "$line" | awk '{print $3}')
			elif [[ "$scanSSname" == '"A-Helix"' ]] ; then coilCODEhelix=$(echo "$line" | awk '{print $3}')
			fi
		done < ss_"${filenm}".xpm
		
		echo $' Reducing colour coding to helix-sheet-turn-coil...\n'
		sleep 2
		
		while IFS= read -r line; do
			scanSSname=$(echo "$line" | awk '{print $6}')
			if [[ "$scanSSname" == '"B-Bridge"' ]] ; then 
				echo "$line" | awk '{print $1 $2 $coilCODEbeta $4 $5 $7}' >> ss_"${filenm}"_HETC.xpm
			elif [[ "$scanSSname" == '"Bend"' ]] ; then
				echo "$line" | awk '{print $1 $2 $coilCODEcoil $4 $5 $7}' >> ss_"${filenm}"_HETC.xpm
			elif [[ "$scanSSname" == '"5-Helix"' ]] ; then
				echo "$line" | awk '{print $1 $2 $coilCODEhelix $4 $5 $7}' >> ss_"${filenm}"_HETC.xpm
			elif [[ "$scanSSname" == '"3-Helix"' ]] ; then
				echo "$line" | awk '{print $1 $2 $coilCODEhelix $4 $5 $7}' >> ss_"${filenm}"_HETC.xpm
			else echo "$line" >> ss_"${filenm}"_HETC.xpm
			fi
		done < ss_"${filenm}".xpm
		
		echo $' Converting output ss_xpm to an eps file...\n\n'
		sleep 2
		#eval "$gmx_exe_path" xpm2ps -f ss_"${filenm}"_HETC.xpm -o ss_"${filenm}"_colortype2.eps -rainbow blue || true
		eval "$gmx_exe_path" xpm2ps -f ss_"${filenm}"_HETC.xpm -o ss_"${filenm}"_colortype1.eps || true
		echo "$demA"$' Converting eps to pdf...\n'
		sleep 1
		ps2pdf ss_"${filenm}"_colortype1.eps ss_"${filenm}"_colortype1size2.pdf || true
		# ps2pdf -sPAPERSIZE=ledger ss_"${filenm}"_colortype1.eps ss_"${filenm}"_colortype1size2.pdf || true
		#ps2pdf -sPAPERSIZE=ledger ss_"${filenm}"_colortype2.eps ss_"${filenm}"_colortype2size2.pdf || true
		ps2pdf ss_"${filenm}"_colortype1.eps ss_"${filenm}"_colortype1size1.pdf || true
		#ps2pdf ss_"${filenm}"_colortype2.eps ss_"${filenm}"_colortype2size1.pdf || true

		echo $' Converting ss_eps to png...\n'
		sleep 1
		convert ss_"${filenm}"_colortype1.eps -trim -bordercolor ss_"${filenm}"_colortype1.png || echo "$demA"$' The program "convert" not found!\n'

		convert ss_"${filenm}"_colortype1.eps -trim -bordercolor white -units pixelsperinch -density 600 ss_"${filenm}"_colortype1.png || echo "$demA"$' The program "convert" not found!\n'
		
		currentSecStrdir="$(pwd)""/Secondary_structure"
		nSecStr=1
		bkupSecStrdir="$(pwd)""/#Secondary_structure"".""backup.""$nSecStr"
		base_bkupSecStrdir=$(basename "$bkupSecStrdir")
		if [[ -d "$currentSecStrdir" ]]; then
			echo $'\n'"$currentSecStrdir"$' folder exists,\n'"backing it up as $base_bkupSecStrdir"
			sleep 1
			while [[ -d "$bkupSecStrdir" ]]; do
				nSecStr=$(( nSecStr + 1 ))
				bkupSecStrdir="$(pwd)""/#Secondary_structure"".""backup.""$nSecStr"
				base_bkupSecStrdir=$(basename "$bkupSecStrdir")
			done
			mv "$currentSecStrdir" "$bkupSecStrdir" && mkdir ./Secondary_structure || true
			echo $'\n'"Backing up the last Secondary_structure folder and its contents as $base_bkupSecStrdir"
			sleep 1
			mv scount.xvg ss_*.xpm ss_*.eps ss_*.pdf ss_*.png ./Secondary_structure || true
		elif [[ ! -d "$currentSecStrdir" ]]; then
			mkdir Secondary_structure; mv scount.xvg ss_*.xpm ss_*.eps ss_*.pdf ss_*.png ./Secondary_structure || true
		fi
		echo "$demA"$' Secondary structure analysis...DONE'"$demB"
		sleep 2
	fi
}

if [[ "$analysis" == *" 8 "* ]]; then ScanTRAJ; analyser8 ; fi

analyser9()
{
	if [[ $frame_b == 0 && $frame_e == 0 ]] ; then clustr_range="-b 0"
	elif [[ $frame_b == 0 && $frame_e != 0 ]] ; then clustr_range="-b 0 -e $frame_e"
	elif [[ $frame_b != 0 && $frame_e == 0 ]] ; then clustr_range="-b $frame_b"
	elif [[ $frame_b != 0 && $frame_e != 0 ]] ; then clustr_range="-b $frame_b -e $frame_e"
	fi

	echo "$demA"$' Preparing to cluster frames from the trajectory...\n\n\n'
	sleep 2
	if [[ $automode == "full" && $sysType == "protein_only" ]]; then
		echo "Backbone" "Protein" | eval "$gmx_exe_path" cluster -f "${filenm}"_${wraplabel}.xtc \
		-s ${filenm}.tpr -method $method_clust -cutoff $cut_cl -g clustering_details.log \
		-cl clusters_representatives.pdb -dist clusters_rmsd_distribution.xvg $clustr_range \
		-clid cluster_id.xvg -clndx clusters_index.ndx -sz cluster_size.xvg -dt $dt
	elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
		eval "$gmx_exe_path" cluster -f "${filenm}"_${wraplabel}.xtc -s ${filenm}.tpr \
		-method $method_clust -cutoff $cut_cl $clustr_range -g clustering_details.log \
		-cl clusters_representatives.pdb -dt $dt -dist clusters_rmsd_distribution.xvg \
		-clid cluster_id.xvg -clndx clusters_index.ndx -sz cluster_size.xvg
	elif [[ $automode == "full" && $sysType == "protein_lig" ]] ; then
		echo "Backbone" "Protein_$ligname" | eval "$gmx_exe_path" cluster -f "${filenm}"_${wraplabel}.xtc \
		-s ${filenm}.tpr -method $method_clust -cutoff $cut_cl $clustr_range -g clustering_details.log \
		-cl clusters_representatives.pdb -dt $dt -dist clusters_rmsd_distribution.xvg -clid cluster_id.xvg \
		-clndx clusters_index.ndx -sz cluster_size.xvg -n index.ndx
	elif [[ $automode == "semi" && $sysType == "protein_lig" ]] ; then
		eval "$gmx_exe_path" cluster -f "${filenm}"_${wraplabel}.xtc -s ${filenm}.tpr \
		-method $method_clust -dt $dt -cutoff $cut_cl $clustr_range -g clustering_details.log -cl clusters_representatives.pdb -n index.ndx \
		-dist clusters_rmsd_distribution.xvg -clid cluster_id.xvg -clndx clusters_index.ndx -sz cluster_size.xvg
	elif [[ $automode == "full" && $sysType == "protein_dna" ]] ; then
		echo "Backbone" "Protein_DNA" | eval "$gmx_exe_path" cluster -f "${filenm}"_${wraplabel}.xtc \
		-s ${filenm}.tpr -dt $dt -method $method_clust -cutoff $cut_cl $clustr_range -g clustering_details.log \
		-cl clusters_representatives.pdb -dist clusters_rmsd_distribution.xvg -clid cluster_id.xvg -clndx clusters_index.ndx -sz cluster_size.xvg -n index.ndx	
	elif [[ $automode == "semi" && $sysType == "protein_dna" ]] ; then
		eval "$gmx_exe_path" cluster -f "${filenm}"_${wraplabel}.xtc -s ${filenm}.tpr -method $method_clust -dt $dt \
		-cutoff $cut_cl $clustr_range -g clustering_details.log -cl clusters_representatives.pdb -n index.ndx \
		-dist clusters_rmsd_distribution.xvg -clid cluster_id.xvg -clndx clusters_index.ndx -sz cluster_size.xvg
	fi

	currentClusteringdir="$(pwd)""/Clustering"
	nClustering=1
	bkupClusteringdir="$(pwd)""/#Clustering"".""backup.""$nClustering"
	base_bkupClusteringdir=$(basename "$bkupClusteringdir")
	if [[ -d "$currentClusteringdir" ]]; then
		echo $'\n'"$currentClusteringdir"$' folder exists,\n'"backing it up as $base_bkupClusteringdir"
		sleep 1
		while [[ -d "$bkupClusteringdir" ]]; do
			nClustering=$(( nClustering + 1 ))
			bkupClusteringdir="$(pwd)""/#Clustering"".""backup.""$nClustering"
			base_bkupClusteringdir=$(basename "$bkupClusteringdir")
		done
		mv "$currentClusteringdir" "$bkupClusteringdir" && mkdir ./Clustering || true
		echo $'\n'"Backing up the last Clustering folder and its contents as $base_bkupClusteringdir"
		sleep 1
	elif [[ ! -d "$currentClusteringdir" ]]; then mkdir ./Clustering
	fi
	eval "$gmx_exe_path" xpm2ps -f rmsd-clust.xpm -o rmsd-clust.eps
	ps2pdf rmsd-clust.eps rmsd-clust.pdf || true
	mv cluster*.log cluster*.pdb cluster*.xvg cluster*.ndx rmsd-clust.* ./Clustering || true
	echo "$demA"$' Cluster frames from the trajectory...DONE'"$demB"
	sleep 2
}

if [[ "$analysis" == *" 9 "* ]]; then analyser9 ; fi


analyser10()
{
	printf "$demA Preparing to estimate probability density function using KDE...\n\n"
	sleep 2
	echo " Select the data for PDF estimation"
	sleep 1
cat << AnalysisList

Select your choice(s) from the options listed below:
Option  Data
  1     Root mean square deviation (RMSD)
  2     Radius of gyration (Rg)
  3     Hydrogen bonds (hbond)
  4     Solvent accessible surface area (SASA)
  
AnalysisList

read -p ' Enter one or more options here (separated by a space): ' data_kde

# create a bash array listing valid numbers
valid_numbers=(1 2 3 4)

while ! [[ "$data_kde" =~ ^([[:space:]]*[0-9][[:space:]]*)+$ ]] && \
	! [[ "$data_kde" =~ (^|[[:space:]])("${valid_numbers[@]}")([[:space:]]|$) ]]
do
	printf "\n You entered: ${data_kde}\n\n"
	printf " Please enter a valid number!!\n\n"
	read -p ' Enter one or more options here (separated by a space): ' data_kde
done

data_kde_ext=("$data_kde")
count_data_in=0

printf "$demA Generating the input files for KDE\n\n"
for i in ${data_kde_ext[*]} ; do
	if (( $count_data_in == 0 )) ; then
		if [[ "$i" == 1 ]] ; then
			echo "RMSD" > CHAP_kde_dataset_list.dat
		elif [[ "$i" == 2 ]] ; then
			echo "Rg" > CHAP_kde_dataset_list.dat
		elif [[ "$i" == 3 ]] ; then
			echo "Hbond" > CHAP_kde_dataset_list.dat
		elif [[ "$i" == 4 ]] ; then
			echo "SASA" > CHAP_kde_dataset_list.dat				
		fi
	elif (( $count_data_in > 0 )) ; then
		if [[ "$i" == 1 ]] ; then
			echo "RMSD" >> CHAP_kde_dataset_list.dat
		elif [[ "$i" == 2 ]] ; then
			echo "Rg" >> CHAP_kde_dataset_list.dat
		elif [[ "$i" == 3 ]] ; then
			echo "Hbond" >> CHAP_kde_dataset_list.dat
		elif [[ "$i" == 4 ]] ; then
			echo "SASA" >> CHAP_kde_dataset_list.dat
		fi
	fi
	count_data_in=$(( count_data_in + 1 ))
done


}

if [[ "$analysis" == *" 10 "* ]]; then analyser10 ; fi

makeMoviePy1()
{
echo $'load PyMOLsession_allSet.pse\nmovie.produce dynamics_moviePy.mpg, quality 100'\
$'\nquit' > make1_movie_Pyscript.pml
pymol make1_movie_Pyscript.pml
}
	
makeMoviePy2()
{
echo $'load PyMOLsession_allSet.pse\nmovie.produce dynamics_moviePy.mpg, ray, quality=100'\
$'\nquit' > make2_movie_Pyscript.pml
pymol make2_movie_Pyscript.pml
}

##function specific for movie update
makeMoviePyx()
{
echo "cd ${movieDIRECORY}"$'\nload PyMOLsession_allSet.pse'\
$'\nmovie.produce dynamics_moviePy.mpg, quality 100\nquit' > make1_movie_Pyscript.pml
pymol make1_movie_Pyscript.pml
}
##function specific for movie update	
makeMoviePyy()
{
echo "cd ${movieDIRECORY}"$'\nload PyMOLsession_allSet.pse'\
$'\nmovie.produce dynamics_moviePy.mpg, ray, quality=100\nquit' > make2_movie_Pyscript.pml
pymol make2_movie_Pyscript.pml
}

# define function for variables_for_regMD_Movie
variables_for_regMD_Movie()
{
	message_Movie="Preparing to make a summary movie of the trajectory"
	trajectlog="trajectDetails.log"
	simulationcontext="simulation"
	xtcFileMovie="$filenm""_${wraplabel}"
	tprFileMovie="${filenm}"
	outXTCmovie="${filenm}""_trjEvery""$skimov""skipForMovie"
	movieDIRECORY="MOVIE"
}

analyser11()
{
	echo "$demA $message_Movie"

	if [[ $customframeNo == '' ]]; then
cat << askMovielength

Do you want to proceed to making a movie summarized into 200 frames?

  1) Yes.
  2) No, I want a different number of frames for the movie.

askMovielength

	read -p ' Enter 1 or 2 here: ' movieLeng
		while [[ "$movieLeng" != 1 && "$movieLeng" != 2 ]]; do
			printf "\n You entered: $movieLeng \n"
			printf ' Please enter a valid number (1 or 2)!!\n\n'
			read -p ' Enter 1 or 2 here: ' movieLeng
		done
	if [[ $movieLeng == 1 ]] ; then
		customframeNo_int=200
	elif [[ $movieLeng == 2 ]] ; then
		cat "$trajectlog" | grep -v "GROMACS reminds"
		sleep 2

		echo "$demA Above is a summary of your $simulationcontext trajectory."
		sleep 1
		printf ' You may find the info useful to provide a response to the prompt below.\n'
		sleep 2
		printf "$demA How many frames do you want the movie to be composed of?\n\n"
		sleep 1
		read -p ' *Please enter a value here: ' customframeNo
		printf " You entered: $customframeNo \n\n"
		customframeNo_int=$(echo ${customframeNo%\.*})
	fi 
elif [[ "$customframeNo" != '' ]]; then
	customframeNo_int=$(echo ${customframeNo%\.*})
fi

if (( $No_of_frames >= "$customframeNo_int" )) ; then
	skimov_raw=$(awk "BEGIN {print $No_of_frames / $customframeNo_int}")
	skimov=$(echo ${skimov_raw%\.*})
elif (( $No_of_frames < "$customframeNo_int" )) ; then 
	skimov=1
	echo "$demA"" Number of frames in the trajectory: ${No_of_frames}"\
	$'\n'" Total number of frames in the trajectory is less than $customframeNo!"\
	$'\n'" Using ${No_of_frames} frames directly.""$demB"
	customframeNo=$(echo ${No_of_frames})
	customframeNo_int=$(echo ${No_of_frames})
fi

echo "$demA"$' Will now extract frames to be used for the movie...\n\n'
sleep 2
#if [[ $sysType == 1 ]] || [[ $sysType == 2 ]] || [[ $sysType == 3 ]] && [[ $flw == 1 ]] ; then
echo 0 | eval "$gmx_exe_path" trjconv -f "$xtcFileMovie".xtc -s ${tprFileMovie}.tpr -o ${outXTCmovie}.xtc -skip $skimov
#elif [[ $sysType == 1 ]] || [[ $sysType == 2 ]] && [[ $flw == 0 ]] ; then
#gmx trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o "${outXTCmovie}.xtc" -skip $skimov
#fi
sleep 2
echo "$demA"$' Preparing to extract $customframeNo snapshots...\n'
sleep 2
if [[ $automode == "full" && $sysType == "protein_only" ]]; then
	echo 1 | eval "$gmx_exe_path" trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -o summaryForMovie.pdb
elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
	eval "$gmx_exe_path" trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -o summaryForMovie.pdb
elif [[ $automode == "semi" || $automode == "full" ]] && [[ $sysType == "protein_dna" ]]; then
	echo "Protein_DNA" | eval "$gmx_exe_path" trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -n index.ndx -o summaryForMovie.pdb
elif [[ $automode == "semi" && $sysType == "protein_lig" ]]; then
	eval "$gmx_exe_path" trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -n index.ndx -o summaryForMovie.pdb
elif [[ $automode == "full" && $sysType == "protein_lig" ]]; then
	echo "Protein_$ligname" | eval "$gmx_exe_path" trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -n index.ndx -o summaryForMovie.pdb
fi
currentMOVIEdir="$(pwd)""/${movieDIRECORY}"
nMOVIE=1
bkupMOVIEdir="$(pwd)""/#${movieDIRECORY}"".""backup.""$nMOVIE"
base_bkupMOVIEdir=$(basename "$bkupMOVIEdir")
if [[ -d "$currentMOVIEdir" ]]; then
	echo $'\n'"$currentMOVIEdir"$' folder exists,\n'"backing it up as $base_bkupMOVIEdir"
	sleep 1
	while [[ -d "$bkupMOVIEdir" ]]; do
		nMOVIE=$(( nMOVIE + 1 ))
		bkupMOVIEdir="$(pwd)""/#${movieDIRECORY}"".""backup.""$nMOVIE"
		base_bkupMOVIEdir=$(basename "$bkupMOVIEdir")
	done
	mv "$currentMOVIEdir" "$bkupMOVIEdir" && mkdir "${movieDIRECORY}" || true
	echo $'\n'"Backing up the last MOVIE folder and its contents as $base_bkupMOVIEdir"
	sleep 1
elif [[ ! -d "$currentMOVIEdir" ]]; then mkdir "${movieDIRECORY}"
fi	
	
echo $'load summaryForMovie.pdb\nsave PyMOLsession.pse\nintra_fit name ca+c+n+o\n'\
$'preset.pretty(selection='"'all'"$')\n'\
$'spectrum chain, green cyan orange magenta\ncolor atomic, (not elem C)\nbg white\n'\
$'set movie_loop, 0\nsmooth\norient\nviewport 760, 540\nzoom all, -10\n'\
$'set ray_trace_frames=1\nset ray_opaque_background, 0\nset cache_frames=0\n'\
$'mclear\n'"cd ${movieDIRECORY}"$'\nsave PyMOLsession_allSet.pse\nmpng frame_.png\nquit' > prep_movie_Pyscript.pml
	
echo "$demA"$'Now, PyMOL will do the job. You sit back and have a cup of tea...Cheers!'"$demB"
sleep 2
pyM=0
pymol prep_movie_Pyscript.pml || pyM=1
echo "$demA"$' Extract frames as images...DONE'"$demB""$demA"$' Now converting images to movie...\n'
cd ./${movieDIRECORY}
mov_make=0
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.gif || mov_make=1
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.mp4 || mov_make=2

if [[ "$mov_make" == 1 ]] && [[ "$pyM" == 0 ]]; then
	echo "$demA"$'The program "Convert/ImageMagick" could not be found.\nCHAPERONg detected'\
	$'"PyMOL" and will use it to make a movie which may,\n'\
	$'however, be of lesser quality'"$demB"
		
	makeMoviePy1
	makeMoviePy2
		
	echo "$demA"$'Movie (lesser quality) made with PyMOL...\n'
fi

if [[ "$mov_make" == 2 ]] && [[ "$pyM" == 0 ]]; then
	echo "$demA"$'A mp4 movie could not be made. This may be due to the program'\
	$'\n"Convert/ImageMagick" not being found, or some library is missing.'\
	$'\nCHAPERONg detected "PyMOL" and will use it to make a mp4 movie which may,\n'\
	$'however, be of lesser quality'"$demB"
		
	makeMoviePy1
	makeMoviePy2
		
	echo "$demA"$'Movie (lesser quality) made with PyMOL...\n'
fi

#mkdir ./frames
#mv frame_*.png ./frames || true ; mv ../summaryForMovie.pdb ./ || true; mv ../"${outXTCmovie}.xtc" ./ || true ; mv ../*.pse ./ || true ; rm ../*movie_Pyscript.pml || true
rm frame_*.png || true ; rm ../summaryForMovie.pdb || true
rm ../"${outXTCmovie}.xtc" || true 
mv ../*.pse ./ || true ; rm ../*_movie_Pyscript.pml || true
rm ./*_movie_Pyscript.pml || true
rm ./PyMOLsession.pse || true
#rm prep_movie_Pyscript.pml ../prep_movie_Pyscript.pml
#rm ../prep_movie_Pyscript.pml || true; cd ..
#rm *movie_Pyscript.pml prep_movie_Pyscript.pml || true
echo "$demA"$' Convert images to movie...DONE'"$demB"
cd ..
}

analyser11update()
{

echo "$demA"$'Preparing to make a summary movie from a preset PyMOL session\n'
sleep 2

currentMOVIEdir="$(pwd)""/${movieDIRECORY}"
if [[ ! -d "$currentMOVIEdir" ]]; then
	echo "No MOVIE directory from a previous run exists... Exiting"
	exit 1
fi
nMOVIE=1
currentMOVIEgif="$(pwd)""/${movieDIRECORY}/dynamics_movie.gif"
bkupMOVIEgif="$(pwd)""/${movieDIRECORY}/dynamics_movie_""backup""$nMOVIE"".gif"
if [[ -f "$currentMOVIEgif" ]]; then
	base_currentMOVIEgif=$(basename "$currentMOVIEgif")
	base_bkupMOVIEgif=$(basename "$bkupMOVIEgif")
	echo $'\n'"$base_currentMOVIEgif" "exists, backing it up as $base_bkupMOVIEgif"$'\n'
	sleep 1
	while [[ -f "$bkupMOVIEgif" ]]; do
	nMOVIE=$(( nMOVIE + 1 ))
	bkupMOVIEgif="$(pwd)""/${movieDIRECORY}/dynamics_movie_""backup""$nMOVIE"".gif"
	done
	mv "$currentMOVIEgif" "$bkupMOVIEgif" || true
	echo $'\n'"Backing up the last .gif MOVIE as $base_bkupMOVIEgif"
	sleep 1
fi	

nMOVIE=1
currentMOVIEmp4="$(pwd)""/${movieDIRECORY}/dynamics_movie.mp4"
bkupMOVIEmp4="$(pwd)""/${movieDIRECORY}/dynamics_movie_""backup""$nMOVIE"".mp4"
if [[ -f "$currentMOVIEmp4" ]]; then
	base_currentMOVIEmp4=$(basename "$currentMOVIEmp4")
	base_bkupMOVIEmp4=$(basename "$bkupMOVIEmp4")
	echo $'\n'"$base_currentMOVIEmp4"" exists, backing it up as $base_bkupMOVIEmp4"$'\n'
	sleep 1
	while [[ -f "$bkupMOVIEmp4" ]]; do
	nMOVIE=$(( nMOVIE + 1 ))
	bkupMOVIEmp4="$(pwd)""/${movieDIRECORY}/dynamics_movie_""backup""$nMOVIE"".mp4"
	done
	mv "$currentMOVIEmp4" "$bkupMOVIEmp4" || true
	echo $'\n'"Backing up the last .mp4 MOVIE as $base_bkupMOVIEmp4"
	sleep 1
fi	
	
echo "cd ${movieDIRECORY}"$'\nload PyMOLsession_allSet.pse\nmpng frame_.png\nquit' > prep_movie_Pyscript.pml
	
echo "$demA"$'Now, PyMOL will do the job. You sit back and have a cup of tea...Cheers!'"$demB"
sleep 2
pyM=0
pymol prep_movie_Pyscript.pml || pyM=1
echo "$demA"$'Extract frames as images...DONE'"$demB""$demA"$'Now converting images to movie...\n'
cd ./${movieDIRECORY}
mov_make=0
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.gif || mov_make=1
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.mp4 || mov_make=2

if [[ "$mov_make" == 1 ]] && [[ "$pyM" == 0 ]]; then
	echo "$demA"$'The program ''"'"Convert/ImageMagick "'"'"could not be found. CHAPERONg detected "'"'\
	$'PyMOL''"'" and will use it to make a movie which may, however, be of lesser quality""$demB"
		
	makeMoviePyx
	makeMoviePyy
		
	echo "$demA"$'Movie (lesser quality) made with PyMOL...\n'
fi

if [[ "$mov_make" == 2 ]] && [[ "$pyM" == 0 ]]; then
	echo "$demA"$'A mp4 movie could not be made. This may be due to the program'\
	$'"Convert/ImageMagick" not being found, or some library is missing.'\
	$'\nCHAPERONg detected "PyMOL" and will use it to make a mp4 movie which may,'\
	$'however, be of lesser quality'"$demB"
		
	makeMoviePyx
	makeMoviePyy
		
	echo "$demA"$'Movie (lesser quality) made with PyMOL...\n'
fi

rm frame_*.png || true 
# rm ../*movie_Pyscript.pml || true
rm ./*movie_Pyscript.pml || true
echo "$demA"$' Convert images to movie...DONE'"$demB"
cd ..
}

if [[ "$analyse" == "9" ]] && [[ -d "$movieDIRECORY" ]]; then
cat << MovChoic
$demA
Make a new movie or adjust (e.g. the orientation of) a previously prepared one?

  a     Make a new movie
  b     Adjust a previous one

MovChoic

	read -p '*Enter your choice here (a or b): ' moviechoic

	while [[ "$moviechoic" != "a" ]] && [[ "$moviechoic" != "b" ]] ; do
		echo $'\nYou entered: '"$moviechoic"$'\n'
		echo $'Please enter a valid letter!!\n'
		read -p '*Enter your choice here (a or b): ' moviechoic
	done
 
	if [[ "$moviechoic" == "a" ]]; then
		ScanTRAJ; variables_for_regMD_Movie; analyser11
	elif [[ "$moviechoic" == "b" ]]; then
		variables_for_regMD_Movie; analyser11update
	fi

elif [[ "$analyse" == "10" ]] && [[ ! -d "$movieDIRECORY" ]]; then
	ScanTRAJ; variables_for_regMD_Movie; analyser11
fi

if [[ "$analysis" == *" 11 "* ]]; then
	ScanTRAJ; variables_for_regMD_Movie; analyser11
fi
	
analyser12()
{
if [[ $mmGMXpath != '' ]] ; then
	indexer=''
	#if [[ $sysType == 1 ]]; then indexer=''
 	if [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] ; then indexer='-n index.ndx' ; fi
	echo "$demA"$' Preparing to generate input files for g_MMPBSA free energy calculations...\n'

	if [[ $trajFraction == '' ]]; then trajFraction=3 ; fi
	trajFraction=$(echo ${trajFraction%\.*})

	if (( $trajFraction == 1 )) ; then
		simDuratnps_lastFractn_beginINT=0
	elif (( $trajFraction >= 2 )) ; then
		trajFractionFactor=$(( trajFraction - 1 ))
		simDuratnps_lastFractn_begin=$(awk "BEGIN {print $simDuratnps * $trajFractionFactor / $trajFraction}")
		simDuratnps_lastFractn_beginINT=$(echo ${simDuratnps_lastFractn_begin%\.*})
	fi
	No_of_last_third_frames=$(awk "BEGIN {print $No_of_frames / $trajFraction}")
	No_of_last_third_framesINT=$(echo ${No_of_last_third_frames%\.*})

	if [[ $mmpbframesNo == '' ]]; then mmpbframesNo=100 ;fi	
	
	if (( $No_of_last_third_framesINT >= $mmpbframesNo )) ; then
		skipframegpsa=$(awk "BEGIN {print $No_of_last_third_frames / $mmpbframesNo }")
		skipframegpsaINT=$(echo ${skipframegpsa%\.*})
		echo "$demA"" Number of frames in the last third of the trajectory: ${No_of_last_third_frames}"\
		$'\n'" ${skipframegpsaINT} frames will be skipped at intervals to produce a total of ~100 frames for"\
		$'\n'" g_mmpbsa calculations.""$demB"
		
	elif (( $No_of_last_third_framesINT < $mmpbframesNo )) ; then
		skipframegpsaINT=1
		echo "$demA"" Number of frames in the last third of the trajectory: ${No_of_last_third_frames}"\
		$'\n'" Total number of frames in the trajectory is less than 200."\
		$'\n'" No frames will be skipped for g_mmpbsa calculations.""$demB"
	fi

	if [[ "$mmGMX" == "1" ]] ; then
		echo "$demA"$' Preparing to generate a compatible .tpr for g_MMPBSA...\n\n\n'
		sleep 2
		eval $mmGMXpath grompp -f md.mdp -c "${filenm}".gro -p topol.top -o \
		"${filenm}"_TPR_for_g_mmpbsa.tpr $indexer
		echo "$demA"$' Generate a compatible .tpr for g_MMPBSA...DONE'"$demB"
		sleep 2
		echo "$demA"$' Generating a compatible fraction of the trajectory for g_MMPBSA...\n\n\n'
		sleep 2
		if [[ $automode == "full" ]]; then
			echo 0 | eval $mmGMXpath trjconv -s "${filenm}"_TPR_for_g_mmpbsa.tpr -f "${filenm}"_${wraplabel}.xtc \
			-o "${filenm}"_lastFractntraj4_mmpbsa.xtc -b $simDuratnps_lastFractn_beginINT
		
			echo "$demA"$' Generate a compatible trajectory file for g_MMPBSA...DONE'"$demB"
			sleep 2
		
			echo "$demA Extracting $mmpbframesNo frames from the trajectory..."$'\n\n\n'
			echo 0 | eval $mmGMXpath trjconv -s "${filenm}"_TPR_for_g_mmpbsa.tpr -f "${filenm}"_lastFractntraj4_mmpbsa.xtc \
			-o "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -skip $skipframegpsaINT
		
			echo "$demA Extract $mmpbframesNo frames from the trajectory...DONE""$demB"
			sleep 2
		
		elif [[ $automode != "full" ]]; then
			eval $mmGMXpath trjconv -s "${filenm}"_TPR_for_g_mmpbsa.tpr -f "${filenm}"_${wraplabel}.xtc \
			-o "${filenm}"_lastFractntraj4_mmpbsa.xtc -b $simDuratnps_lastFractn_beginINT
		
			echo "$demA"$' Generate a compatible fraction of the trajectory for g_MMPBSA...DONE'"$demB"
			sleep 2
		
			echo "$demA Extracting $mmpbframesNo frames from the trajectory..."$'\n\n\n'
		
			eval $mmGMXpath trjconv -s "${filenm}"_TPR_for_g_mmpbsa.tpr -f "${filenm}"_lastFractntraj4_mmpbsa.xtc \
			-o "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -skip $skipframegpsaINT
		
			echo "$demA Extract $mmpbframesNo frames from the trajectory...DONE""$demB"
			sleep 2
		fi	
	
		echo $' Generate input files for g_MMPBSA free energy calculations...DONE'"$demB"
		sleep 2

		echo "$demA"$' Now preparing to run g_MMPBSA calculations...\n\n\n'
		if [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode == "full" ]]; then
			echo 1 "$ligname" | ${CHAPERONg_PATH}/CHAP_utilities/g_mmpbsa_pkg/g_mmpbsa -f \
			"${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -s \
			"${filenm}"_TPR_for_g_mmpbsa.tpr -n index.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp
		elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode != "full" ]] ; then
			${CHAPERONg_PATH}/CHAP_utilities/g_mmpbsa_pkg/g_mmpbsa -f \
			"${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -s \
			"${filenm}"_TPR_for_g_mmpbsa.tpr -n index.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp
		fi
		echo "$demA"$' Run g_MMPBSA calculations...DONE'"$demB"
		sleep 2
	elif [[ "$mmGMX" == '' ]] ; then
		echo "$demA"$' Now preparing to run g_MMPBSA calculations...\n\n\n'
		if [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode == "full" ]]; then
			echo 1 "$ligname" | g_mmpbsa -f "${filenm}"_${wraplabel}.xtc -s \
			"${filenm}".tpr -n index.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp || \
			echo "$demA"$' g_mmpbsa failed to run. Ensure your environments are properly set...\n'
		elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode != "full" ]] ; then
			g_mmpbsa -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -i pbsa.mdp -pdie 2 \
			-pbsa -decomp || echo "$demA"$' g_mmpbsa failed to run. Ensure your environments are properly set...\n'
		fi
	fi
	echo "$demA"$' Calculating average binding energy & contribution of residues...\n'
	python ${CHAPERONg_PATH}/CHAP_utilities/g_mmpbsa_pkg/MmPbSaStat.py \
	-m energy_MM.xvg -p polar.xvg -a apolar.xvg || \
	python3 ${CHAPERONg_PATH}/CHAP_utilities/g_mmpbsa_pkg/MmPbSaStatPy3.py \
	-m energy_MM.xvg -p polar.xvg -a apolar.xvg || true

	python ${CHAPERONg_PATH}/CHAP_utilities/g_mmpbsa_pkg/MmPbSaDecomp.py -bs \
	-nbs 2000 \-m contrib_MM.dat -p contrib_pol.dat -a contrib_apol.dat || \
	python3 ${CHAPERONg_PATH}/CHAP_utilities/g_mmpbsa_pkg/MmPbSaDecompPy3.py -bs \
	-nbs 2000 -m contrib_MM.dat -p contrib_pol.dat -a contrib_apol.dat || true

	if [[ "$mmGMX" == "1" ]] ; then
		if [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode == "full" ]]; then
			echo 1 "$ligname" | eval ${CHAPERONg_PATH}/CHAP_utilities/g_mmpbsa_pkg/energy2bfac -s \
			"${filenm}"_TPR_for_g_mmpbsa.tpr -i energyMapIn.dat
		elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode != "full" ]] ; then
			eval ${CHAPERONg_PATH}/CHAP_utilities/g_mmpbsa_pkg/energy2bfac -s \
			"${filenm}"_TPR_for_g_mmpbsa.tpr -i energyMapIn.dat
		fi
	elif [[ "$mmGMX" == '' ]] ; then
		if [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode == "full" ]]; then
			echo 1 "$ligname" | energy2bfac -s "${filenm}"_TPR_for_g_mmpbsa.tpr -i energyMapIn.dat || \
			echo "$demA"$'energy2bfac failed to run. Ensure your environment are properly set...\n'
		elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $automode != "full" ]] ; then
			energy2bfac -s "${filenm}"_TPR_for_g_mmpbsa.tpr -i energyMapIn.dat || \
			echo "$demA"$'energy2bfac failed to run. Ensure your environment are properly set...\n'
		fi
	fi

	echo "$demA"$' Calculate average binding energy & contribution of residues...DONE'"$demB"
	sleep 2

	AnaName="MMPBSA"
	currentAnadir="$(pwd)""/$AnaName"
	nDir=1
	bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
	if [[ -d "$currentAnadir" ]]; then
		base_currentAnadir=$(basename "$currentAnadir")
		base_bkupAnadir=$(basename "$bkupAnadir")
		echo $'\n'"$base_currentAnadir"$' folder exists,\n'"backing it up as $base_bkupAnadir"
		sleep 1
		while [[ -d "$bkupAnadir" ]]; do
			nDir=$(( nDir + 1 )); bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
		done
		mv "$currentAnadir" "$bkupAnadir" && mkdir ./$AnaName
		echo $'\n'"Backing up the last $AnaName folder and its contents as $base_bkupAnadir"
		sleep 1
	elif [[ ! -d "$currentAnadir" ]]; then mkdir ./$AnaName
	fi
	rm "${filenm}"_TPR_for_g_mmpbsa.tpr "${filenm}"_lastFractntraj4_mmpbsa.xtc \
	"${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc || true
	mv energy_MM.xvg polar.xvg apolar.xvg contrib_MM.dat contrib_pol.dat contrib_apol.dat ./$AnaName || true
	mv full_energy.dat summary_energy.dat final_contrib_energy.dat energyMapIn.dat ./$AnaName || true
	mv complex.pdb subunit_1.pdb subunit_2.pdb ./$AnaName || true

elif [[ $mmGMXpath == '' ]] ; then
	echo "$demA"$'GMX path for g_mmpbsa not set. Use the parFile option!\n'
fi

}

if [[ "$analysis" == *" 12 "* ]]; then ScanTRAJ; analyser12 ; fi

useFoundPCA_sham()
{
	echo "$demA"$' Preparing PCA-derived FES with gmx sham...\n'
	sleep 1
				
	eval "$gmx_exe_path" sham -f ./PCA/PCA_2dproj_$filenm.xvg -ls ./PCA/FEL_PCA_sham_$filenm.xpm -notime || true
			
	eval "$gmx_exe_path" xpm2ps -f ./PCA/FEL_PCA_sham_$filenm.xpm -o ./PCA/FEL_PCA_sham_$filenm.eps -rainbow red || true
		
	ps2pdf ./PCA/FEL_PCA_sham_$filenm.eps ./PCA/FEL_PCA_sham_${filenm}_landscape.pdf || true
	# ps2pdf -sPAPERSIZE=ledger ./PCA/FEL_PCA_sham_$filenm.eps ./PCA/FEL_PCA_sham_${filenm}_landscape.pdf || true
			
	ps2pdf ./PCA/FEL_PCA_sham_$filenm.eps ./PCA/FEL_PCA_sham_${filenm}_portrait.pdf || true

	pdf2ppm -png -r 600 ./PCA/FEL_PCA_sham_${filenm}_portrait.pdf ./PCA/FEL_PCA_sham_${filenm}_portrait.png || true

	convert ./PCA/FEL_PCA_sham_$filenm.eps -trim -bordercolor white ./PCA/FEL_PCA_sham_${filenm}_convertps2png.png || true

	convert ./PCA/FEL_PCA_sham_$filenm.eps -trim -bordercolor white -units pixelsperinch \
	-density 600 -resize 3000x5000 ./PCA/FEL_PCA_sham_${filenm}_convertps2pngfull.png || true
				
	currentFELPCAdir="$(pwd)""/PCA_FES_sham"
	nFELpca=1
	bkupFELpcadir="$(pwd)""/#PCA_FES_sham"".""backup.""$nFELpca"
	base_bkupFELpcadir=$(basename "$bkupFELpcadir")
	if [[ -d "$currentFELPCAdir" ]]; then
		echo $'\n'"$currentFELPCAdir"$' folder exists,\n'"backing it up as $base_bkupFELpcadir"
		sleep 1
		while [[ -d "$bkupFELpcadir" ]]; do
			nFELpca=$(( nFELpca + 1 ))
			bkupFELpcadir="$(pwd)""/#PCA_FES_sham"".""backup.""$nFELpca"
			base_bkupFELpcadir=$(basename "$bkupFELpcadir")
		done
		mv "$currentFELPCAdir" "$bkupFELpcadir" && mkdir ./PCA_FES_sham || true
		echo $'\n'"Backing up the last PCA_FES_sham folder and its contents as $base_bkupFELpcadir"
		sleep 1
	elif [[ ! -d "$currentFELPCAdir" ]]; then mkdir PCA_FES_sham
	fi
	mv ./PCA/FEL_PCA_sham_* enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./PCA_FES_sham || true
	echo "$demA"$' Prepare Gibbs FES with gmx sham...DONE'"$demB"
	sleep 2
	ana_folder="PCA_FES_sham"
}

useFoundRgRMSData_sham()
{
	echo "$demA"$' Preparing Rg Vs RMSD FES using gmx sham...\n\n'
	sleep 1
				
	eval "$gmx_exe_path" sham -f RgVsRMSD.xvg -ls FEL_sham_RgVsRMSD_$filenm.xpm -notime || true
			
	eval "$gmx_exe_path" xpm2ps -f FEL_sham_RgVsRMSD_$filenm.xpm -o FEL_sham_RgVsRMSD_$filenm.eps -rainbow red || true
		
	ps2pdf FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_landscape.pdf || true
	# ps2pdf -sPAPERSIZE=ledger FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_landscape.pdf || true
			
	ps2pdf FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_portrait.pdf || true

	pdf2ppm -png -r 600 FEL_sham_RgVsRMSD_${filenm}_portrait.pdf FEL_sham_RgVsRMSD_${filenm}_portrait.png || true

	convert FEL_sham_RgVsRMSD_$filenm.eps -trim -bordercolor white FEL_sham_RgVsRMSD_${filenm}_convertps2png.png || true

	convert FEL_sham_RgVsRMSD_$filenm.eps -trim -bordercolor white -units pixelsperinch \
	-density 600 -resize 3000x5000 FEL_sham_RgVsRMSD_${filenm}_convertps2pngfull.png || true
				
	currentFELshamRgVsRMSDdir="$(pwd)""/RgVsRMSD_FEL_sham"
	nFELRgVsRMSD=1
	bkupFELshamRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FEL_sham"".""backup.""$nFELRgVsRMSD"
	base_bkupFELshamRgVsRMSDdir=$(basename "$bkupFELshamRgVsRMSDdir")
	if [[ -d "$currentFELshamRgVsRMSDdir" ]]; then
		echo $'\n'"$currentFELshamRgVsRMSDdir"$' folder exists,\n'"backing it up as $base_bkupFELshamRgVsRMSDdir"
		sleep 1
		while [[ -d "$bkupFELshamRgVsRMSDdir" ]]; do
			nFELRgVsRMSD=$(( nFELRgVsRMSD + 1 ))
			bkupFELshamRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FEL_sham"".""backup.""$nFELRgVsRMSD"
			base_bkupFELshamRgVsRMSDdir=$(basename "$bkupFELshamRgVsRMSDdir")
		done
		mv "$currentFELshamRgVsRMSDdir" "$bkupFELshamRgVsRMSDdir" && mkdir ./RgVsRMSD_FEL_sham || true
		echo $'\n'"Backing up the last RgVsRMSD_FEL_sham folder and its contents as $base_bkupFELshamRgVsRMSDdir"
		sleep 1
	elif [[ ! -d "$currentFELshamRgVsRMSDdir" ]]; then mkdir RgVsRMSD_FEL_sham
	fi
	mv FEL_sham_RgVsRMSD_* RMSData.dat RgData.dat RgVsRMSD.xvg enthalpy.xpm \
	entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./RgVsRMSD_FEL_sham || true

	echo "$demA"$' Prepare Rg Vs RMSD FES with gmx sham...DONE'"$demB"
	sleep 2
	ana_folder="RgVsRMSD_FEL_sham"
}

order_parameters()
{
cat << orderPair

  Which order parameter pair do you want to use for the FES calculations?

    1) Principal components
    2) Rg versus RMSD
    3) Other user-provided pair of order parameters

orderPair

read -p '  Enter 1, 2 or 3 here: ' orderPair_choice
	while [[ $orderPair_choice != 1 && $orderPair_choice != 2 && $orderPair_choice != 3 ]]
	do
		echo $'\nYou entered: '"$orderPair_choice"
		echo $'Please enter a valid number (1, 2 or 3)!!\n'
		read -p '  Enter 1, 2 or 3 here: ' orderPair_choice
	done

if [[ $orderPair_choice == 1 ]] ; then
	echo "$demA"$' Checking the working directory for pre-calculated PCA_2d projection data...'"$demB"
	sleep 2
	exist2dPCA="$(pwd)""/PCA/PCA_2dproj_""$filenm"".xvg"
	if [[ -f "$exist2dPCA" ]] ; then
		if [[ $automode == "semi" ]] ; then
			echo $'CHAPERONg: Pre-calculated PCA_2d projection file found!\nFile found:'" $exist2dPCA"
			sleep 2
cat << askFELuseexist

Do you want to use this file for the 2D energetic landscape calculations?

  1) Yes, use the file above
  2) No, repeat PCA calculations and use the new output for FES plot
  3) No, I want to provide another file to be used

askFELuseexist

			read -p ' Enter 1, 2 or 3 here: ' PCFile
			while [[ "$PCFile" != 1 && "$PCFile" != 2 && "$PCFile" != 3 ]]; do
				echo $'\nYou entered: '"$PCFile"
				echo $'Please enter a valid number (1, 2 or 3)!!\n'
				read -p ' Enter 1, 2 or 3 here: ' PCFile
			done
		elif [[ $automode == "full" ]] ; then
			echo "$demA"$' Pre-calculated PCA_2d projection data found!\n File found:'" $exist2dPCA"\
			$'\n *CHAPERONg in auto mode\n'" $exist2dPCA will be used for FES plotting"
			sleep 2
		fi
	elif [[ ! -f "$exist2dPCA" ]] ; then analyser7		
	fi
	inputPair="$exist2dPCA"
#if order parameter pair is RgVsRMSD	
elif [[ $orderPair_choice == 2 ]] ; then
	existRg="$(pwd)""/Rg/""${filenm}_Rg_ns.xvg"
	if [[ -f "$existRg" ]] ; then
		if [[ $automode == "semi" ]] ; then
			echo "$demA"$' Pre-calculated Rg data found!\n File found:'" $existRg"$'\n'
			sleep 2
cat << askFELuseexist

Do you want to use this file for free energy surface calculations?

  1) Yes, use the file above
  2) No, repeat Rg calculations and use the new output for FES plot
  3) No, I want to provide another file to be used

askFELuseexist

			read -p ' Enter 1, 2 or 3 here: ' PCFile
			while [[ "$PCFile" != 1 && "$PCFile" != 2 && "$PCFile" != 3 ]]; do
				echo $'\nYou entered: '"$PCFile"
				echo $'Please enter a valid number (1, 2 or 3)!!\n'
				read -p ' Enter 1, 2 or 3 here: ' PCFile
			done
		elif [[ $automode == "full" ]] ; then
			echo "$demA"$' Pre-calculated Rg data found!\nFile found:'" $existRg"\
			$'\n *CHAPERONg in auto mode\n'" $existRg will be used for FES plotting"
			sleep 2
		fi
	elif [[ ! -f "$existRg" ]] ; then analyser4
	fi
	cat "$existRg" | grep -v "^[@#]" | awk '{print $2}' > RgData.dat

	inputRg_xvgData="$existRg"

	#check for pre-calculated RMSD data
	existRMSD="$(pwd)""/RMSD/""${filenm}_BB-rmsd.xvg"
	if [[ -f "$existRMSD" ]] ; then
		if [[ $automode == "semi" ]] ; then
			echo "$demA"$' Pre-calculated RMSD data found!\n File found:'" $existRMSD"$'\n'
			sleep 2
cat << askFELuseexist

Do you want to use this file for free energy landscape calculations?

  1) Yes, use the file above
  2) No, repeat RMSD calculations and use the new output for FES plot
  3) No, I want to provide another file to be used

askFELuseexist

			read -p ' Enter 1, 2 or 3 here: ' PCFile
			while [[ "$PCFile" != 1 && "$PCFile" != 2 && "$PCFile" != 3 ]]; do
				echo $'\nYou entered: '"$PCFile"
				echo $'Please enter a valid number (1, 2 or 3)!!\n'
				read -p ' Enter 1, 2 or 3 here: ' PCFile
			done
		elif [[ $automode == "full" ]] ; then
			echo "$demA"$' Pre-calculated RMSD data found!\n File found:'" $existRMSD"\
			$'\n *CHAPERONg in auto mode\n'" $existRMSD will be used for FES plotting"
			sleep 2
		fi
	elif [[ ! -f "$existRMSD" ]] ; then analyser2		
	fi
	cat "$existRMSD" | grep -v "^[@#]" | awk '{print $2}' > RMSData.dat
	echo "# This file contains the RMSD and Rg values extracted by CHAPERONg" > RgVsRMSD.xvg
	echo "# from the data generated by GROMACS..." >> RgVsRMSD.xvg
	echo "#" >> RgVsRMSD.xvg
	echo "@    title "$'"Plot of Rg against RMSD"' >> RgVsRMSD.xvg
	echo "@    xaxis  label "$'"'"RMSD (nm)"$'"' >> RgVsRMSD.xvg
	echo "@    yaxis  label "$'"'"Rg (nm)"$'"' >> RgVsRMSD.xvg
	echo "@TYPE xy" >> RgVsRMSD.xvg
	paste -d "        " RMSData.dat RgData.dat >> RgVsRMSD.xvg

	inputRMSD_xvgData="$existRMSD"
	inputPair="RgVsRMSD.xvg"
fi
}

analyser13()
{
	echo "$demA"$' Constructing free energy surface with gmx sham...\n'
	order_parameters
	if [[ $orderPair_choice == 1 && $automode == "semi" ]] ; then
		if [[ "$PCFile" == 1 ]]; then useFoundPCA_sham
		elif [[ "$PCFile" == 2 ]]; then analyser7; useFoundPCA_sham	
		elif [[ "$PCFile" == 3 ]]; then echo ""
			read -p ' Provide the path to the pre-calculated 2d_PCA projection file: ' precalcPCfile
	
			echo "$demA"$'Preparing PCA-derived FES with user-provided 2d_PCA file...\n\n'
			sleep 1

			felcal=0			

			eval "$gmx_exe_path" sham -f $precalcPCfile -ls FEL_PCA_sham_$filenm.xpm -notime || felcal=1
			if [[ "$felcal" == 0 ]] ; then
				
				eval "$gmx_exe_path" xpm2ps -f FEL_PCA_sham_$filenm.xpm -o FEL_PCA_sham_$filenm.eps -rainbow red || true

				ps2pdf FEL_PCA_sham_$filenm.eps FEL_PCA_sham_${filenm}_landscape.pdf || true
				# ps2pdf -sPAPERSIZE=ledger FEL_PCA_sham_$filenm.eps FEL_PCA_sham_${filenm}_landscape.pdf || true

				ps2pdf FEL_PCA_sham_$filenm.eps FEL_PCA_sham_${filenm}_portrait.pdf || true
			
				currentFELPCAdir="$(pwd)""/PCA_FES_sham"
				nFELpca=1
				bkupFELpcadir="$(pwd)""/#PCA_FES_sham"".""backup.""$nFELpca"
				base_bkupFELpcadir=$(basename "$bkupFELpcadir")
				if [[ -d "$currentFELPCAdir" ]]; then
					echo $'\n'"$currentFELPCAdir"$' folder exists,\n'"backing it up as $base_bkupFELpcadir"
					sleep 1
					while [[ -d "$bkupFELpcadir" ]]; do
						nFELpca=$(( nFELpca + 1 ))
						bkupFELpcadir="$(pwd)""/#PCA_FES_sham"".""backup.""$nFELpca"
						base_bkupFELpcadir=$(basename "$bkupFELpcadir")
					done
					mv "$currentFELPCAdir" "$bkupFELpcadir" && mkdir ./PCA_FES_sham || true
					echo $'\n'"Backing up the last PCA_FES_sham folder and its contents as $base_bkupFELpcadir"
					sleep 1
					mv FEL_PCA_sham_*.xpm enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg FEL_PCA_sham_*.eps FEL_PCA_sham_*.pdf ./PCA_FES_sham || true
				elif [[ ! -d "$currentFELPCAdir" ]]; then mkdir PCA_FES_sham
					mv FEL_PCA_sham* enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./PCA_FES_sham || true
				fi
				echo "$demA"$' Prepare PCA-based 2D energetic landscape using gmx sham...DONE'"$demB"
				sleep 2
				ana_folder="PCA_FES_sham"
			
			elif [[ "$felcal" == 1 ]] ; then
				echo $'Calculation failed.\n'" Please confirm that you have entered the right path/file as input!"
				sleep 1
			fi
		fi			
	elif [[ $orderPair_choice == 1 && $automode == "full" ]] ; then useFoundPCA_sham
	elif [[ $orderPair_choice == 2 && $automode == "semi" ]] ; then
		if [[ "$PCFile" == 1 ]]; then useFoundRgRMSData_sham
		elif [[ "$PCFile" == 2 ]]; then analyser2; analyser4; useFoundRgRMSData_sham	
		elif [[ "$PCFile" == 3 ]]; then echo ""
cat << inputFormat

Do you want to provide individual order parameter or a pre-combined pair?
  1) Individual Rg and RMSD data
  2) Pre-combined Rg Vs RMSD pair

inputFormat

			read -p ' Enter 1 or 2 here: ' inFormat
			while [[ "$inFormat" != 1 && "$inFormat" != 2 ]]; do
				echo $'\nYou entered: '"$inFormat"
				echo $'Please enter a valid number (1 or 2)!!\n'
				read -p ' Enter 1 or 2 here: ' inFormat
			done

			if [[ "$inFormat" == 1 ]]; then
				read -p ' Provide the path to the pre-calculated Rg.xvg data: ' precalcRg
				inputRg_xvgData="$precalcRg"
				echo ""

				read -p ' Provide the path to the pre-calculated RMSD.xvg data: ' precalcRMSD
				inputRMSD_xvgData="$precalcRMSD"

				echo "$demA"$'Preprocessing user-provided Rg and RMSD data files...\n\n'
				sleep 1
				cat "$precalcRg" | grep -v "^[@#]" | awk '{print $2}' > RgData.dat
				
				cat "$precalcRMSD" | grep -v "^[@#]" | awk '{print $2}' > RMSData.dat
				echo "# This file contains the RMSD and Rg values extracted by CHAPERONg" > RgVsRMSD.xvg
				echo "# from the data generated by GROMACS..." >> RgVsRMSD.xvg
				echo "#"  >> RgVsRMSD.xvg
				echo "@    title "$'"Plot of Rg against RMSD"' >> RgVsRMSD.xvg
				echo "@    xaxis  label "$'"'"RMSD (nm)"$'"' >> RgVsRMSD.xvg
				echo "@    yaxis  label "$'"'"Rg (nm)"$'"' >> RgVsRMSD.xvg
				echo "@TYPE xy" >> RgVsRMSD.xvg
				paste -d "        " RMSData.dat RgData.dat >> RgVsRMSD.xvg
				
				precalcRgRMS="./RgVsRMSD.xvg"
			elif [[ "$inFormat" == 2 ]]; then
				read -p ' Provide the path to the pre-calculated RgVsRMSD.xvg: ' precalcRgRMS
			fi

			felcal=0
			eval "$gmx_exe_path" sham -f $precalcRgRMS -ls FEL_sham_RgVsRMSD_$filenm.xpm -notime || felcal=1
		
			if [[ "$felcal" == 0 ]] ; then
				eval "$gmx_exe_path" xpm2ps -f FEL_sham_RgVsRMSD_$filenm.xpm -o FEL_sham_RgVsRMSD_$filenm.eps -rainbow red || true
	
				ps2pdf FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_landscape.pdf || true
				# ps2pdf -sPAPERSIZE=ledger FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_landscape.pdf || true
		
				ps2pdf FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_portrait.pdf || true

				pdf2ppm -png -r 600 FEL_sham_RgVsRMSD_${filenm}_portrait.pdf FEL_sham_RgVsRMSD_${filenm}_portrait.png || true

				convert FEL_sham_RgVsRMSD_$filenm.eps -trim -bordercolor white FEL_sham_RgVsRMSD_${filenm}_convertps2png.png || true

				convert FEL_sham_RgVsRMSD_$filenm.eps -trim -bordercolor white -units pixelsperinch \
				-density 600 -resize 3000x5000 FEL_sham_RgVsRMSD_${filenm}_convertps2pngfull.png || true
			
				currentFELshamRgVsRMSDdir="$(pwd)""/RgVsRMSD_FEL_sham"
				nFELRgVsRMSD=1
				bkupFELshamRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FEL_sham"".""backup.""$nFELRgVsRMSD"
				base_bkupFELshamRgVsRMSDdir=$(basename "$bkupFELshamRgVsRMSDdir")
				if [[ -d "$currentFELshamRgVsRMSDdir" ]]; then
					echo $'\n'"$currentFELshamRgVsRMSDdir"$' folder exists,\n'"backing it up as $base_bkupFELshamRgVsRMSDdir"
					sleep 1
					while [[ -d "$bkupFELshamRgVsRMSDdir" ]]; do
						nFELRgVsRMSD=$(( nFELRgVsRMSD + 1 ))
						bkupFELshamRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FEL_sham"".""backup.""$nFELRgVsRMSD"
						base_bkupFELshamRgVsRMSDdir=$(basename "$bkupFELshamRgVsRMSDdir")
					done
					mv "$currentFELshamRgVsRMSDdir" "$bkupFELshamRgVsRMSDdir" && \
					mkdir ./RgVsRMSD_FEL_sham || true
					echo $'\n'"Backing up the last RgVsRMSD_FEL_sham folder and its contents as $base_bkupFELshamRgVsRMSDdir"
					sleep 1
					mv FEL_sham_RgVsRMSD_* RgData.dat RMSData.dat RgVsRMSD.xvg enthalpy.xpm \
					entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./RgVsRMSD_FEL_sham || true
				elif [[ ! -d "$currentFELshamRgVsRMSDdir" ]]; then mkdir RgVsRMSD_FEL_sham
					mv FEL_sham_RgVsRMSD_* RgData.dat RMSData.dat RgVsRMSD.xvg enthalpy.xpm \
					entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./RgVsRMSD_FEL_sham || true
				fi
				echo "$demA"$' Prepare Rg Vs RMSD 2D energetic landscape with gmx sham...DONE'"$demB"
				sleep 2
				ana_folder="RgVsRMSD_FEL_sham"

			elif [[ "$felcal" == 1 ]] ; then
				echo $'Calculation failed.\n'" Please confirm that you have entered the right path/file as input!"
				sleep 1
			fi
		fi
 	elif [[ $orderPair_choice == 2 && $automode == "full" ]] ; then useFoundRgRMSData_sham

	elif [[ $orderPair_choice == 3 ]] ; then
cat << inputFormat

Do you want to provide individual order parameter or a pre-combined pair?
  1) Individual order parameters
  2) Pre-combined order parameter pair

inputFormat

		read -p ' Enter 1 or 2 here: ' inFormat
		while [[ "$inFormat" != 1 && "$inFormat" != 2 ]]; do
			echo $'\nYou entered: '"$inFormat"
			echo $'Please enter a valid number (1 or 2)!!\n'
			read -p ' Enter 1 or 2 here: ' inFormat
		done

		if [[ "$inFormat" == 1 ]]; then
			read -p ' Provide the path to the 1st order_parameter.xvg data: ' precalcOrderPar1
			echo ""
			read -p ' Provide the path to the 2nd order_parameter.xvg data: ' precalcOrderPar2
				
			echo "$demA"$'Preprocessing user-provided order parameter data files...\n\n'
			sleep 1
			cat "$precalcOrderPar1" | grep -v "^[@#]" | awk '{print $2}' > precalcOrderPar1.dat
			cat "$precalcOrderPar2" | grep -v "^[@#]" | awk '{print $2}' > precalcOrderPar2.dat
			echo "# This file contains the order parameters 1 and 2 data extracted by CHAPERONg" > OrderParameterPair.xvg
			echo "# from the data files provided by user..." >> OrderParameterPair.xvg
			echo "#" >> OrderParameterPair.xvg
			echo "@    title "$'"Plot of Order Parameter 1 against Order Parameter 2"' >> OrderParameterPair.xvg
			echo "@    xaxis  label "$'"Order Parameter 1"' >> OrderParameterPair.xvg
			echo "@    yaxis  label "$'"Order Parameter 2"' >> OrderParameterPair.xvg
			echo "@TYPE xy" >> OrderParameterPair.xvg
			paste -d "        " precalcOrderPar1.dat precalcOrderPar2.dat >> OrderParameterPair.xvg
				
			precalcOrderParPair="./OrderParameterPair.xvg"
		elif [[ "$inFormat" == 2 ]]; then
			read -p ' Provide the path to the pre-calculated order_parameter_pair.xvg: ' precalcOrderParPair
		fi

		felcal=0
		eval "$gmx_exe_path" sham -f $precalcOrderParPair -ls FEL_sham_OrderParameterPair_${filenm}.xpm -notime || felcal=1
		
		if [[ "$felcal" == 0 ]] ; then
			eval "$gmx_exe_path" xpm2ps -f FEL_sham_OrderParameterPair_$filenm.xpm -o \
			FEL_sham_OrderParameterPair_${filenm}.eps -rainbow red || true

			ps2pdf FEL_sham_OrderParameterPair_$filenm.eps \
			FEL_sham_OrderParameterPair_${filenm}_landscape.pdf || true
			# ps2pdf -sPAPERSIZE=ledger FEL_sham_OrderParameterPair_$filenm.eps \
			# FEL_sham_OrderParameterPair_${filenm}_landscape.pdf || true
	
			ps2pdf FEL_sham_OrderParameterPair_$filenm.eps \
			FEL_sham_OrderParameterPair_${filenm}_portrait.pdf || true

			pdf2ppm -png -r 600 FEL_sham_OrderParameterPair_${filenm}_portrait.pdf \
			FEL_sham_OrderParameterPair_${filenm}_portrait.png || true

			convert FEL_sham_OrderParameterPair_$filenm.eps -trim -bordercolor white \
			FEL_sham_OrderParameterPair_${filenm}_convertps2png.png || true

			convert FEL_sham_OrderParameterPair_$filenm.eps -trim -bordercolor white -units pixelsperinch \
			-density 600 -resize 3000x5000 FEL_sham_OrderParameterPair_${filenm}_convertps2pngfull.png || true
			
			currentFELshamOrderParameterPairdir="$(pwd)""/OrderParameterPair_FEL_sham"
			nFELOrderParameterPair=1
			bkupFELshamOrderParameterPairdir="$(pwd)""/#OrderParameterPair_FEL_sham"".""backup.""$nFELOrderParameterPair"
			base_bkupFELshamOrderParameterPairdir=$(basename "$bkupFELshamOrderParameterPairdir")
			if [[ -d "$currentFELshamOrderParameterPairdir" ]]; then
				echo $'\n'"$currentFELshamOrderParameterPairdir"$' folder exists,\n'"backing it up as $base_bkupFELshamOrderParameterPairdir"
				sleep 1
				while [[ -d "$bkupFELshamOrderParameterPairdir" ]]; do
					nFELOrderParameterPair=$(( nFELOrderParameterPair + 1 ))
					bkupFELshamOrderParameterPairdir="$(pwd)""/#OrderParameterPair_FEL_sham"".""backup.""$nFELOrderParameterPair"
					base_bkupFELshamOrderParameterPairdir=$(basename "$bkupFELshamOrderParameterPairdir")
				done
				mv "$currentFELshamOrderParameterPairdir" "$bkupFELshamOrderParameterPairdir" && mkdir ./OrderParameterPair_FEL_sham || true
				echo $'\n'"Backing up the last OrderParameterPair_FEL_sham folder and its contents as $base_bkupFELshamOrderParameterPairdir"
				sleep 1
				mv FEL_sham_OrderParameterPair_* precalcOrderPar1.dat precalcOrderPar2.dat OrderParameterPair.xvg \
				enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./OrderParameterPair_FEL_sham || true
			elif [[ ! -d "$currentFELshamOrderParameterPairdir" ]]; then mkdir OrderParameterPair_FEL_sham
				mv FEL_sham_OrderParameterPair_* precalcOrderPar1.dat precalcOrderPar2.dat OrderParameterPair.xvg \
				enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./OrderParameterPair_FEL_sham || true
			fi
			echo "$demA"$' Prepare 2D energetic landscape with gmx sham...DONE'"$demB"
			sleep 2
			ana_folder="OrderParameterPair_FEL_sham"
		elif [[ "$felcal" == 1 ]] ; then
			echo $'Calculation failed.\n'" Please confirm that you have entered the right path/file as input!"
			sleep 1
		fi
	fi

	# extract lowest free energy structures
	echo "$demA Identifying the lowest energy bins and frames"$'\n'
	sleep 2
	min0_bin_index=$(grep -F '0.000' ./${ana_folder}/shamlog.log | tail -n 1 | awk '{print $5}')
	echo " The bin with index $min0_bin_index contains the structures with the lowest energy"
	sleep 2
	echo $'\n Three representative structures will be extracted from this bin\n'
	sleep 1
	# min0_index_spaced=" $min0_index "
	min0_struct1_frame=$(grep -A1 "\[ $min0_bin_index \]" ./${ana_folder}/bindex.ndx | tail -n 1)
	min0_struct2_frame=$(grep -A2 "\[ $min0_bin_index \]" ./${ana_folder}/bindex.ndx | tail -n 1)
	min0_struct3_frame=$(grep -A3 "\[ $min0_bin_index \]" ./${ana_folder}/bindex.ndx | tail -n 1)
	# min0_struct1=(grep -FA1 \["$min0_index_spaced"\] ./${ana_folder}/bindex.ndx | tail -n 1)
	ScanTRAJ
	# sim_timestep
	echo $' Identifying the corresponding times for the lowest energy structures...'
	sleep 2
	min0_struct1_time=$(awk "BEGIN {print $sim_timestep * $min0_struct1_frame}")
	min0_struct2_time=$(awk "BEGIN {print $sim_timestep * $min0_struct2_frame}")
	min0_struct3_time=$(awk "BEGIN {print $sim_timestep * $min0_struct3_frame}")

	echo "$demA"$' Extracting lowest energy structures from the trajectory...\n\n\n'
	sleep 2
	structure1="${filenm}"_LowestEnergyBin_structure1_frame"$min0_struct1_frame".pdb
	structure2="${filenm}"_LowestEnergyBin_structure2_frame"$min0_struct2_frame".pdb
	structure3="${filenm}"_LowestEnergyBin_structure3_frame"$min0_struct3_frame".pdb

	if [[ $automode == "full" && $sysType == "protein_only" ]]; then
		echo 1 | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -o "$structure1" -dump "$min0_struct1_time"
		echo 1 | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -o "$structure2" -dump "$min0_struct2_time"
		echo 1 | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -o "$structure3" -dump "$min0_struct3_time"
	elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
		eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -o "$structure1" -dump "$min0_struct1_time"
		eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -o "$structure2" -dump "$min0_struct2_time"
		eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -o "$structure3" -dump "$min0_struct3_time"
	elif [[ $automode == "semi" || $automode == "full" ]] && [[ $sysType == "protein_dna" ]]; then
		echo "Protein_DNA" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -n index.ndx -o "$structure1" -dump "$min0_struct1_time"
		echo "Protein_DNA" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -n index.ndx -o "$structure2" -dump "$min0_struct2_time"
		echo "Protein_DNA" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -n index.ndx -o "$structure3" -dump "$min0_struct3_time"
	elif [[ $automode == "semi" && $sysType == "protein_lig" ]]; then
		eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s \
		"${filenm}".tpr -n index.ndx -o "$structure1" -dump "$min0_struct1_time"
		eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s \
		"${filenm}".tpr -n index.ndx -o "$structure2" -dump "$min0_struct2_time"
		eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s \
		"${filenm}".tpr -n index.ndx -o "$structure3" -dump "$min0_struct3_time"
	elif [[ $automode == "full" && $sysType == "protein_lig" ]]; then
		echo "Protein_$ligname" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -n index.ndx -o "$structure1" -dump "$min0_struct1_time"
		echo "Protein_$ligname" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -n index.ndx -o "$structure2" -dump "$min0_struct2_time"
		echo "Protein_$ligname" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
		-s "${filenm}".tpr -n index.ndx -o "$structure3" -dump "$min0_struct3_time"
	fi
	echo "$demA"$' Extract lowest energy structures from the trajectory...DONE\n\n'
	sleep 2
	mv "$structure1" "$structure2" "$structure3" ./${ana_folder}/ || true

cat << inform
 All outputs from the free energy surface calculations have been saved to
 the folder ${ana_folder}.

 In the prompt below, you may check the shamlog.log file to find the index of
 the bin of interest you may wish to extract from, and the bindex.ndx file to
 identify the frames in the bin.

inform

	get_more_structs=1
	while [[ "$get_more_structs" == 1 ]]
	do
		
cat << extractMoreStructs
 Do you want extract additional structures from the trajectory?

  1) Yes
  2) No

extractMoreStructs

		read -p ' Enter a response here (1 or 2): ' get_more_structs
		
		while [[ "$get_more_structs" != 1 && "$get_more_structs" != 2 ]]
		do
			echo $' \nPlease enter the appropriate response (1 or 2)!!\n'
			echo $' Extract additional structures from the trajectory??\n  1) Yes\n  2) No\n'
			read -p ' Enter 1 or 2 here: ' get_more_structs
		done
		
		if [[ "$get_more_structs" == 2 ]] ; then
			echo ""
		elif [[ "$get_more_structs" == 1 ]] ; then
			echo ""
			read -p ' Specify the frame number of the structure to extract: ' frame_no

			echo $'\n Identifying the corresponding time for the specified frame...'
			sleep 2
			spec_struct_time=$(awk "BEGIN {print $sim_timestep * $frame_no}")
			spec_struct="${filenm}"_structure_at_frame"$frame_no".pdb

			echo "$demA"$' Extracting the specified structure from the trajectory...\n\n\n'
			sleep 2

			if [[ $automode == "full" && $sysType == "protein_only" ]]; then
				echo 1 | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
				-s "${filenm}".tpr -o "$spec_struct" -dump "$spec_struct_time"
			elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
				eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
				-s "${filenm}".tpr -o "$spec_struct" -dump "$spec_struct_time"
			elif [[ $automode == "semi" || $automode == "full" ]] && [[ $sysType == "protein_dna" ]]; then
				echo "Protein_DNA" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
				-s "${filenm}".tpr -n index.ndx -o "$spec_struct" -dump "$spec_struct_time"
			elif [[ $automode == "semi" && $sysType == "protein_lig" ]]; then
				eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s \
				"${filenm}".tpr -n index.ndx -o "$spec_struct" -dump "$spec_struct_time"
			elif [[ $automode == "full" && $sysType == "protein_lig" ]]; then
				echo "Protein_$ligname" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc \
				-s "${filenm}".tpr -n index.ndx -o "$spec_struct" -dump "$spec_struct_time"
			fi
			echo "$demA"$' Extract the specified structure from the trajectory...DONE\n\n'
			sleep 2
			mv "$spec_struct" ./${ana_folder}/ || true
		fi
	done

	echo "$demA"$' Construct free energy landscape with gmx sham...DONE'"$demB"
	sleep 2
}

if [[ "$analysis" == *" 13 "* ]]; then analyser13 ; fi

useFoundPCA_FESPy()
{
	echo "$demA"$' Extracting principal components from 2d_projection data...\n\n'
	sleep 2
	cat $exist2dPCA | grep -v "^[@#]" | awk '{print $1}' > PC1.dat
	cat $exist2dPCA | grep -v "^[@#]" | awk '{print $2}' > PC2.dat
	paste -d "," PC1.dat PC2.dat > OrderParameterPair.dat
	cat PC1.dat | sort -n > sorted_PC1.dat
	cat PC2.dat | sort -n > sorted_PC2.dat

	echo $' Extract principal components from 2d_projection data...DONE'"$demB"
	sleep 2

	echo "$demA"$' Preparing parameters for FES calculations...\n'
	sleep 2

	echo $' Determining minimal and maximal data points...\n'
	sleep 1
	minPC1=$(head -1 sorted_PC1.dat)
	maxPC1=$(tail -1 sorted_PC1.dat)
	minPC2=$(head -1 sorted_PC2.dat)
	maxPC2=$(tail -1 sorted_PC2.dat)

	rm sorted_PC1.dat sorted_PC2.dat

	ScanTRAJ
	echo "minPar1,$minPC1"$'\n'"maxPar1,$maxPC1"$'\n'"minPar2,$minPC2" > CHAP_fes_Par.in
	echo "maxPar2,$maxPC2"$'\n'"no_of_frames,$No_of_frames" >> CHAP_fes_Par.in
	echo $'XaxisL,PC1\nYaxisL,PC2\n'"Temp,$Temp" >> CHAP_fes_Par.in
	echo $'outFilename,PCA_FES\nplotTitle,PCA-derived' >> CHAP_fes_Par.in
	echo $'x_bin_count,100\ny_bin_count,100' >> CHAP_fes_Par.in

	echo $' The input parameters for FES calculations have been prepared\n'\
	$'These parameters have been written to file (CHAP_fes_Par.in).\n'
	echo $'  Do you want to proceed?\n\n   (1) Yes\n   (2) No\n'
	read -p '  Enter a response here (1 or 2): ' para_set

	while [[ "$para_set" != 1 && "$para_set" != 2 ]]; do
			echo $'\n You entered: '"$para_set"
			echo $' Please enter a valid number (1 or 2)!!\n'
			read -p '  Enter 1 or 2 here: ' para_set
		done

	if [[ "$para_set" == 2 ]]; then
		echo ' Skipping FES calculation'
		sleep 2
	elif [[ "$para_set" == 1 ]]; then
		echo "$demA"$' Now running construct_free_en_surface.py to construct FES...\n'
		sleep 2
		python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_construct_free_en_surface.py || \
		python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_construct_free_en_surface.py

		echo $'\n Run construct_free_en_surface.py...DONE'
		sleep 2
		echo $'\n Cleaning up...'"$demB"
		#echo "$demA"$' Cleaning up...\n'
		sleep 1

		currentFESchapPCAdir="$(pwd)""/PCA_FES_chap"
		nFESPCA=1
		bkupFESchapPCAdir="$(pwd)""/#PCA_FES_chap"".""backup.""$nFESPCA"
		base_bkupFESchapPCAdir=$(basename "$bkupFESchapPCAdir")
		if [[ -d "$currentFESchapPCAdir" ]]; then
			echo $'\n'"$currentFESchapPCAdir"$' folder exists,\n'"backing it up as $base_bkupFESchapPCAdir"
			sleep 1
			while [[ -d "$bkupFESchapPCAdir" ]]; do
				nFESPCA=$(( nFESPCA + 1 )); bkupFESchapPCAdir="$(pwd)""/#PCA_FES_chap"".""backup.""$nFESPCA"
				base_bkupFESchapPCAdir=$(basename "$bkupFESchapPCAdir")
			done
			echo $'\n'"Backing up the last PCA_FES_chap folder and its contents as $base_bkupFESchapPCAdir"
			sleep 1
			mv "$currentFESchapPCAdir" "$bkupFESchapPCAdir" && mkdir ./PCA_FES_chap || true
		elif [[ ! -d "$currentFESchapPCAdir" ]]; then mkdir PCA_FES_chap
		fi

		OrderParameter1="PC1.dat"
		OrderParameter2="PC2.dat"
		fesFigure="PCA_FES.png"
		results_folder="PCA_FES_chap"
		# cat "$exist2dPCA" | grep -v "^[@#]" | awk '{print $1}' > SimTime.dat
		# fetch simulation time from the trajectory using the ScanTRAJ fxn
		checksimtime="SimTime.dat"
		if [[ ! -f "$checksimtime" ]] ; then
			echo "$demA"$' Extracting the simulation time-points from the trajectory...\n'
			ScanTRAJ
			increment_factor=$(awk "BEGIN {print $simDuratnINTns / $No_of_frames}")
			simtimeRecorded=0
			echo "$simtimeRecorded" > SimTime.dat
			while [[ "$simtimeRecorded" != "$simDuratnINTns" ]]; do
				simtimeRecorded=$(awk "BEGIN {print $simtimeRecorded + $increment_factor}")
				echo "$simtimeRecorded" >> SimTime.dat
				if [[ "$simtimeRecorded" == "$simDuratnINTns" ]]; then
					break
				fi
			done
			echo $' Extract simulation time-points from the trajectory...DONE'"$demB"
			sleep 2
		fi
	fi
}

useFoundRgRMSData_FESPy()
{
	echo "$demA"$' Preparing parameters for FES calculations...\n'
	sleep 2
	
	# cat RMSData.dat | sort -n > sorted_RMSData.dat
	# cat RgData.dat | sort -n > sorted_RgData.dat

	paste -d "," RMSData.dat RgData.dat > OrderParameterPair.dat

	echo $' Determining minimal and maximal data points...'
	sleep 1
	# minRMSD=$(head -1 sorted_RMSData.dat)
	# maxRMSD=$(tail -1 sorted_RMSData.dat)
	# minRg=$(head -1 sorted_RgData.dat)
	# maxRg=$(tail -1 sorted_RgData.dat)

	# rm sorted_RMSData.dat sorted_RgData.dat

	ScanTRAJ
	echo "no_of_frames,$No_of_frames" > CHAP_fes_Par.in
	echo $'XaxisL,RMSD (nm)\nYaxisL,Rg (nm)\n'"Temp,$Temp" >> CHAP_fes_Par.in
	echo $'outFilename,RgVsRMSD_FES\nplotTitle,Rg Vs RMSD' >> CHAP_fes_Par.in

	echo "$demA"$' Now running construct_free_en_surface.py to construct FES...\n'
	sleep 2
	python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_construct_free_en_surface.py || \
	python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_construct_free_en_surface.py

	echo $'\n Run construct_free_en_surface.py...DONE'"$demB"
	sleep 2
	echo $'\n Cleaning up...'"$demB"
	#echo "$demA"$' Cleaning up...\n'
	sleep 1
	
	currentFESchapRgVsRMSDdir="$(pwd)""/RgVsRMSD_FES_chap"
	nFESRgVsRMSD=1
	bkupFESchapRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FES_chap"".""backup.""$nFESRgVsRMSD"
	base_bkupFESchapRgVsRMSDdir=$(basename "$bkupFESchapRgVsRMSDdir")
	if [[ -d "$currentFESchapRgVsRMSDdir" ]]; then
		echo $'\n'"$currentFESchapRgVsRMSDdir"$' folder exists,\n'"backing it up as $base_bkupFESchapRgVsRMSDdir"
		sleep 1
		while [[ -d "$bkupFESchapRgVsRMSDdir" ]]; do
			nFESRgVsRMSD=$(( nFESRgVsRMSD + 1 ))
			bkupFESchapRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FES_chap"".""backup.""$nFESRgVsRMSD"
			base_bkupFESchapRgVsRMSDdir=$(basename "$bkupFESchapRgVsRMSDdir")
		done
		echo $'\n'"Backing up the last RgVsRMSD_FES_chap folder and its contents as $base_bkupFESchapRgVsRMSDdir"
		sleep 1
		mv "$currentFESchapRgVsRMSDdir" "$bkupFESchapRgVsRMSDdir" && mkdir ./RgVsRMSD_FES_chap || true
	elif [[ ! -d "$currentFESchapRgVsRMSDdir" ]]; then mkdir RgVsRMSD_FES_chap
	fi

	OrderParameter1="RMSData.dat"
	OrderParameter2="RgData.dat"
	fesFigure="RgVsRMSD_FES.png"
	results_folder="RgVsRMSD_FES_chap"
	#fetch simulation time from one of the .xvg files
	cat "$inputRg_xvgData" | grep -v "^[@#]" | awk '{print $1}' > SimTime.dat
	mv RgVsRMSD.xvg ./"$results_folder" || true

}

analyser14()
{	
	echo "$demA"$' Constructing FES using CHAPERONg energetic landscape scripts...\n'
	order_parameters
 	if [[ $orderPair_choice == 1 && $automode == "semi" ]] ; then
		if [[ "$PCFile" == 1 ]]; then useFoundPCA_FESPy
		elif [[ "$PCFile" == 2 ]]; then analyser7; useFoundPCA_FESPy	
		elif [[ "$PCFile" == 3 ]]; then echo ""
			read -p ' Provide the path to the pre-calculated 2d_PCA projection file: ' precalcPCfile
			exist2dPCA="$precalcPCfile"
			useFoundPCA_FESPy
		fi			
 	elif [[ $orderPair_choice == 1 && $automode == "full" ]] ; then useFoundPCA_FESPy
	elif [[ $orderPair_choice == 2 && $automode == "semi" ]] ; then
		if [[ "$PCFile" == 1 ]]; then useFoundRgRMSData_FESPy
		elif [[ "$PCFile" == 2 ]]; then analyser2; analyser4; useFoundRgRMSData_FESPy	
		elif [[ "$PCFile" == 3 ]]; then echo ""
			read -p ' Provide the path to the pre-calculated Rg.xvg data: ' precalcRg
			read -p ' Provide the path to the pre-calculated RMSD.xvg data: ' precalcRMSD
			
			echo "$demA"$' Pre-processing user-provided Rg Vs RMSD data files...\n\n'
			sleep 1
			cat "$precalcRg" | grep -v "^[@#]" | awk '{print $2}' > RgData.dat
			RgData="$precalcRg"
			cat "$precalcRMSD" | grep -v "^[@#]" | awk '{print $2}' > RMSData.dat
			echo "# This file contains the RMSD and Rg values extracted by CHAPERONg"
			echo "# from the data generated by GROMACS..."
			echo "#"
			echo "@    title "$'"Plot of Rg against RMSD"' > RgVsRMSD.xvg
			echo "@    xaxis  label "$'"'"RMSD (nm)"$'"' >> RgVsRMSD.xvg
			echo "@    yaxis  label "$'"'"Rg (nm)"$'"' >> RgVsRMSD.xvg
			echo "@TYPE xy" >> RgVsRMSD.xvg
			paste -d "        " RMSData.dat RgData.dat >> RgVsRMSD.xvg
			
			RMSData="$precalcRMSD"
			precalcRgRMS="./RgVsRMSD.xvg"
			useFoundRgRMSData_FESPy
		fi

	elif [[ $orderPair_choice == 2 && $automode == "full" ]] ; then
		useFoundRgRMSData_FESPy
		# mv RgVsRMSD.xvg ./"$results_folder" || true
	elif [[ $orderPair_choice == 3 ]] ; then
cat << inputFormat

Do you want to provide individual order parameter or a pre-combined pair?
  1) Individual order parameters
  2) Pre-combined order parameter pair

inputFormat

		read -p ' Enter 1 or 2 here: ' inFormat
		while [[ "$inFormat" != 1 && "$inFormat" != 2 ]]; do
			echo $'\nYou entered: '"$inFormat"
			echo $'Please enter a valid number (1 or 2)!!\n'
			read -p ' Enter 1 or 2 here: ' inFormat
		done

		if [[ "$inFormat" == 1 ]]; then
			read -p ' Provide the path to the 1st order_parameter.xvg data: ' precalcOrderPar1
			inputprecalcOrderPar1="$precalcOrderPar1"

			echo ""

			read -p ' Provide the path to the 2nd order_parameter.xvg data: ' precalcOrderPar2
			inputprecalcOrderPar2="$precalcOrderPar2"
				
			echo "$demA"$'Pre-processing user-provided order parameter data files...\n\n'
			sleep 1
			cat "$precalcOrderPar1" | grep -v "^[@#]" | awk '{print $2}' > precalcOrderPar1.dat
			cat "$precalcOrderPar2" | grep -v "^[@#]" | awk '{print $2}' > precalcOrderPar2.dat
			echo "# This file contains the order parameters 1 and 2 data extracted by CHAPERONg" > OrderParameterPair.xvg
			echo "# from the data files provided by user..." >> OrderParameterPair.xvg
			echo "#" >> OrderParameterPair.xvg
			echo "@    title "$'"Plot of Order Parameter 1 against Order Parameter 2"' >> OrderParameterPair.xvg
			echo "@    xaxis  label "$'"Order Parameter 1"' >> OrderParameterPair.xvg
			echo "@    yaxis  label "$'"Order Parameter 2"' >> OrderParameterPair.xvg
			echo "@TYPE xy" >> OrderParameterPair.xvg
			paste -d "        " precalcOrderPar1.dat precalcOrderPar2.dat >> OrderParameterPair.xvg
				
		elif [[ "$inFormat" == 2 ]]; then
			read -p ' Provide the path to the pre-calculated order_parameter_pair.xvg: ' precalcOrderParPair
			cat "$precalcOrderParPair" | grep -v "^[@#]" | awk '{print $1}' > precalcOrderPar1.dat
			cat "$precalcOrderParPair" | grep -v "^[@#]" | awk '{print $2}' > precalcOrderPar2.dat
		fi

		echo "$demA"$' Preparing parameters for FES calculations...\n'
		sleep 2
		cat precalcOrderPar1.dat | sort -n > sorted_precalcOrderPar1.dat
		cat precalcOrderPar2.dat | sort -n > sorted_precalcOrderPar2.dat

		echo $' Determining minimal and maximal data points...'
		sleep 1
		minParam1=$(head -1 sorted_precalcOrderPar1.dat)
		maxParam1=$(tail -1 sorted_precalcOrderPar1.dat)
		minParam2=$(head -1 sorted_precalcOrderPar2.dat)
		maxParam2=$(tail -1 sorted_precalcOrderPar2.dat)

		rm sorted_precalcOrderPar1.dat sorted_precalcOrderPar2.dat

		ScanTRAJ
		echo "minPar1,$minParam1"$'\n'"maxPar1,$maxParam1"$'\n'"minPar2,$minParam2" > CHAP_fes_Par.in
		echo "maxPar2,$maxParam2"$'\n'"no_of_frames,$No_of_frames"$'\n'"Temp,$Temp" >> CHAP_fes_Par.in
		echo $'XaxisL,PC1\nYaxisL,PC2\noutFilename,FES\nplotTitle,' >> CHAP_fes_Par.in

		echo "$demA"$' Now running construct_free_en_surface.py to construct FES...\n'
		sleep 2
		python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_construct_free_en_surface.py || \
		python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_construct_free_en_surface.py

		echo $'\n Run construct_free_en_surface.py...DONE'
		sleep 2
		echo $'\n Cleaning up...'"$demB"
		#echo "$demA"$' Cleaning up...\n'
		sleep 1

		currentFESchapdir="$(pwd)""/FES_chap"
		nFES=1
		bkupFESchapdir="$(pwd)""/#FES_chap"".""backup.""$nFES"
		base_bkupFESchapdir=$(basename "$bkupFESchapdir")
		if [[ -d "$currentFESchapdir" ]]; then
			echo $'\n'"$currentFESchapdir"$' folder exists,\n'"backing it up as $base_bkupFESchapdir"
			sleep 1
			while [[ -d "$bkupFESchapdir" ]]; do
				nFES=$(( nFES + 1 )); bkupFESchapdir="$(pwd)""/#FES_chap"".""backup.""$nFES"
				base_bkupFESchapdir=$(basename "$bkupFESchapdir")
			done
			echo $'\n'"Backing up the last FES_chap folder and its contents as $base_bkupFESchapdir"
			sleep 1
			mv "$currentFESchapdir" "$bkupFESchapdir" && mkdir ./FES_chap || true
			
		elif [[ ! -d "$currentFESchapdir" ]]; then mkdir FES_chap
		fi
		OrderParameter1="precalcOrderPar1.dat"
		OrderParameter2="precalcOrderPar2.dat"
		fesFigure="FES.png"
		results_folder="FES_chap"
		#fetch simulation time from one of the .xvg files
		cat "$inputprecalcOrderPar1" | grep -v "^[@#]" | awk '{print $1}' > SimTime.dat
		mv OrderParameterPair.xvg ./"$results_folder" || true
	fi

	cat OrderParameters1_2_dG.dat | grep --color=none "\S" > OrderParameters1_2_dG_nogap.dat
	# paste SimTime.dat OrderParameters1_2_dG_nogap.dat > SimTime_OrderParameters1_2_dG.dat
	#sort based on dG values
	# cat SimTime_OrderParameters1_2_dG.dat | sort -k 4,4 -n > SimTime_OrderParameters1_2_dG-sorted.dat
	cat OrderParameters1_2_dG_nogap.dat | sort -k 3,3 -n > OrderParameters1_2_dG_nogap-sorted.dat
	mv "$fesFigure" OrderParameterPair.dat CHAP_fes_Par.in ./"$results_folder" || true
	mv OrderParameters1_2_dG.dat OrderParameters1_2_dG_nogap.dat ./"$results_folder" || true
	
	bin_prob=0
	raise_bin_problem()
	{
		bin_prob=1
		echo "$demA"$' **NOTE:\n  The choice of the bin size for the 2D histogram is not suitable.'\
		$'\n  Lowest energy structures will not be extracted.'\
		$'\n  Landscape data points will not be mapped to simulation time.'\
		$'\n  Try experimenting with different bin counts. To do this,'\
		$'\n  you can alter the FES input parameter file "CHAP_fes_Par.in".'\
		$'\n  Please know that the problem could also be that'\
		$'\n  the number of frames in the trajectory is small.'
	}
	
	
	mv binning_summary.dat ./"$results_folder" || raise_bin_problem

	echo "$demA"$' Construct free energy surface...DONE'"$demB"
	sleep 2
	
	if (( "$bin_prob" == 0 )) ; then
		echo "$demA"$'\n Proceed to identify, extract the lowest energy structure from the landscape?\n'
		read -p ' Enter a response here (yes or no): ' getEnergyMin
		
		while [[ "$getEnergyMin" != "yes" && "$getEnergyMin" != "no" && \
			"$getEnergyMin" != '"yes"' && "$getEnergyMin" != '"no"' && \
			"$getEnergyMin" != "y" && "$getEnergyMin" != "n" && \
			"$getEnergyMin" != '"y"' && "$getEnergyMin" != '"n"' ]]
		do
			echo $' \nPlease enter the appropriate response (a "yes" or a "no")!!\n'
			echo $' Identify and extract the lowest energy structure from the landscape??\n'
			read -p ' Enter yes or no here: ' getEnergyMin
		done
		if [[ "$getEnergyMin" == "yes" || "$getEnergyMin" == '"yes"' || \
			"$getEnergyMin" == '"y"' || "$getEnergyMin" == "y" ]] ; then
			if [[ ! -d "collect_mappings" ]] ; then mkdir collect_mappings
			elif [[ -d "collect_mappings" ]] ; then
				rm -r collect_mappings
				mkdir collect_mappings
			fi
			paste SimTime.dat $OrderParameter1 $OrderParameter2 > SimTime_OrderParameters1_2.dat
			echo "$demA"$' Mapping landscape data point to simulation time...\n'
			sleep 2
			python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_map_fes_parameter_to_simTime.py || \
			python ${CHAPERONg_PATH}/CHAP_utilities/CHAP_map_fes_parameter_to_simTime.py

			echo $' Map landscape data point to simulation time...DONE\n'
			sleep 2

			tail -n +2 ./collect_mappings/1.txt | sort -k2,2n -k3,3n > ./collect_mappings/sorted_1.txt

			head -n 2 ./collect_mappings/1.txt | tail -n 1 > ./collect_mappings/EnergyMinim.txt

			echo $' Identifying the corresponding time for the lowest energy structure...\n'
			sleep 2
			python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_get_lowest_en_datapoint.py || \
			python ${CHAPERONg_PATH}/CHAP_utilities/CHAP_get_lowest_en_datapoint.py
			lowEn_time=$(tail -n 1 ./collect_mappings/lowest_energy_datapoints_timed.dat | awk '{print $1}')
			lowEn_time_ps=$(awk "BEGIN {print $lowEn_time * 1000}")
			rm ./collect_mappings/1.txt ./collect_mappings/sorted_1.txt

			echo $' Identify the corresponding time for the lowest energy structure...DONE'
			sleep 1
			
			echo "$demA"$' Extracting the lowest energy structure from the trajectory...\n\n\n'
			sleep 2

			if [[ $automode == "full" && $sysType == "protein_only" ]]; then
				echo 1 | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
				-o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
			elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
				eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
				-o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
			elif [[ $automode == "semi" || $automode == "full" ]] && [[ $sysType == "protein_dna" ]]; then
				echo "Protein_DNA" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
				-n index.ndx -o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
			elif [[ $automode == "semi" && $sysType == "protein_lig" ]]; then
				eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n \
				index.ndx -o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
			elif [[ $automode == "full" && $sysType == "protein_lig" ]]; then
				echo "Protein_$ligname" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s \
				"${filenm}".tpr -n index.ndx -o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
			fi
			echo "$demA"$' Extract lowest energy structure from the trajectory...DONE'"$demB"
			sleep 2

			mv "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb ./collect_mappings/ || true
			mv SimTime.dat OrderParameters1_2_dG_nogap-sorted.dat collect_mappings SimTime_OrderParameters1_2.dat ./"$results_folder" || true
			rm $OrderParameter1 $OrderParameter2
		
cat << extractMoreStructs
 Do you want to map all data points from the FES to simulation time and,
 optionally, extract additional structures from the trajectory based on
 free energy and order parameters?

  1) Yes
  2) No

extractMoreStructs

			read -p ' Enter a response here (1 or 2): ' getMoreStructs
			
			while [[ "$getMoreStructs" != 1 && "$getMoreStructs" != 2 ]]
			do
				echo $' \nPlease enter the appropriate response (1 or 2)!!\n'
				echo $' Proceed to mapping the FES data points to simulation time??\n  1) Yes\n  2) No\n'
				read -p ' Enter 1 or 2 here: ' getMoreStructs
			done
			
			if [[ "$getMoreStructs" == 1 ]] ; then
				echo "$demA"$' Mapping all data points from the Free Energy Surface to simulation time...\n'
				sleep 2
				cp ./"$results_folder"/SimTime_OrderParameters1_2.dat . || true
				cp ./"$results_folder"/OrderParameters1_2_dG_nogap-sorted.dat . || true
				if [[ ! -d "collect_mappings_extra" ]] ; then mkdir collect_mappings_extra
				elif [[ -d "collect_mappings_extra" ]] ; then
					rm -r collect_mappings_extra
					mkdir collect_mappings_extra
				fi
				python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_map_all_dataPoint_to_simTime.py || \
				python ${CHAPERONg_PATH}/CHAP_utilities/CHAP_map_all_dataPoint_to_simTime.py

				echo "$demA"$' Collecting approximate simulation entries for mapped data points...\n'
				sleep 2
				echo $' Sorting and pre-processing FES data points...\n'
				cd ./collect_mappings_extra
				dataRec_counter=0
				for dataRec in DataPt_*.txt
				do
					dataRec_counter=$(( dataRec_counter + 1 ))
					if (( "$dataRec_counter" == 1 )) ; then
						tail -n +2 "${dataRec}" | sort -k2,2n -k3,3n > sorted_"${dataRec}"
						mapped_fes_Data=$(head -n 2 "${dataRec}" | tail -n 1)
						echo "$dataRec"$'\t'"$mapped_fes_Data" > mapped_dataPt.txt
						rm "${dataRec}"
					elif (( "$dataRec_counter" > 1 )) ; then
						tail -n +2 "${dataRec}" | sort -k2,2n -k3,3n >> sorted_"${dataRec}"
						mapped_fes_Data=$(head -n 2 "${dataRec}" | tail -n 1)
						echo "$dataRec"$'\t'"$mapped_fes_Data" >> mapped_dataPt.txt
						rm "${dataRec}"
					fi
				done
				cd ../
				echo $' Generating approximated simulation times mapped with free energy...\n'
				python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_approx_dataPoint_simTime_for_freeEn.py || \
				python ${CHAPERONg_PATH}/CHAP_utilities/CHAP_approx_dataPoint_simTime_for_freeEn.py

				(head -n 1 ./collect_mappings_extra/mappedFESdataPoints_timed.dat && tail -n +2 \
				./collect_mappings_extra/mappedFESdataPoints_timed.dat | sort -k1,1n -k4,4n) > \
				./collect_mappings_extra/time-sorted_mappedFESdataPoints_timed.dat
				
				(head -n 1 ./collect_mappings_extra/mappedFESdataPoints_timed.dat && tail -n +2 \
				./collect_mappings_extra/mappedFESdataPoints_timed.dat | sort -k4,4n -k2,2n) > \
				./collect_mappings_extra/energy-sorted_mappedFESdataPoints_timed.dat
				
				mv ./collect_mappings_extra/mappedFESdataPoints_timed.dat ./"$results_folder"/collect_mappings/
				mv ./collect_mappings_extra/time-sorted_mappedFESdataPoints_timed.dat ./"$results_folder"/collect_mappings/
				mv ./collect_mappings_extra/energy-sorted_mappedFESdataPoints_timed.dat ./"$results_folder"/collect_mappings/
				rm -r ./collect_mappings_extra SimTime_OrderParameters1_2.dat OrderParameters1_2_dG_nogap-sorted.dat

				echo $'\n Collect approximate simulation entries for mapped data points...DONE\n'"$demB"
				sleep 2
				echo "$demA"$'\n **NOTE: An output file named mappedFESdataPoints_timed.dat and copies of it'\
				$'\n (sorted by time or energy) all containing the mapped simulation time, order'\
				$'\n parameters and the corresponding free energy have been generated and saved'\
				$'\n into the folder '"$results_folder""/collect_mappings."\
				$'\n\n **You may use this file to identify the sumulation time(s) of the'\
				$'\n structure(s) you may want to extract from the 2D representation of the free'\
				$'\n energy landscape.'"$demB"
				sleep 2

cat << extractMoreStructs
 Do you want to extract a structure from the trajectory based on free energy?

 1) Yes
 2) No

extractMoreStructs

				read -p ' Enter a response here (1 or 2): ' getMoreStructure
				
				while [[ "$getMoreStructure" != 1 && "$getMoreStructure" != 2 ]]
				do
					echo $' \nPlease enter the appropriate response (1 or 2)!!\n'
					echo $' Extract a structure from the trajectory based on free energy??\n  1) Yes\n  2) No\n'
					read -p ' Enter 1 or 2 here: ' getMoreStructure
				done
			
				while [[ "$getMoreStructure" == 1 ]]
				do
					echo ""
					read -p ' Specify the time (in ns) at which to extract a structure: ' frameTime

					frameTime_ps=$(awk "BEGIN {print $frameTime * 1000}")

					echo "$demA"$' Extracting structure from the trajectory at the specified time...\n\n'
					sleep 2

					if [[ $automode == "full" && $sysType == "protein_only" ]]; then
						echo 1 | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
						-o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
					elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
						eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
						-o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
					elif [[ $automode == "semi" || $automode == "full" ]] && [[ $sysType == "protein_dna" ]]; then
						echo "Protein_DNA" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr \
						-n index.ndx -o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
					elif [[ $automode == "semi" && $sysType == "protein_lig" ]]; then
						eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n \
						index.ndx -o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
					elif [[ $automode == "full" && $sysType == "protein_lig" ]]; then
						echo "Protein_$ligname" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s \
						"${filenm}".tpr -n index.ndx -o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
					fi
					echo "$demA"$' Extract structure from the trajectory at the specified time...DONE'"$demB"
					sleep 2

					mv "${filenm}"_Structure_at_Time"$frameTime_ps".pdb ./"$results_folder"/collect_mappings/ || true

					echo "$demA"$'\n The structure has been saved to the folder '"./$results_folder""/collect_mappings"

cat << extractMoreStructs

 Do you want to extract an additional structure from the trajectory based
 on free energy?

 1) Yes
 2) No

extractMoreStructs

					read -p ' Enter a response here (1 or 2): ' getMoreStructure
					
					while [[ "$getMoreStructure" != 1 && "$getMoreStructure" != 2 ]]
					do
						echo $' \nPlease enter the appropriate response (1 or 2)!!\n'
						echo $' Extract a structure from the trajectory based on free energy??\n  1) Yes\n  2) No\n'
						read -p ' Enter 1 or 2 here: ' getMoreStructure
					done
				done
			elif [[ "$getMoreStructs" != 1 ]] ; then echo ""
			fi
		elif [[ "$getEnergyMin" == "no" || "$getEnergyMin" == '"no"' ]] ; then
			mv SimTime.dat OrderParameters1_2_dG_nogap-sorted.dat ./"$results_folder" || true
		fi
	fi
}

if [[ "$analysis" == *" 14 "* ]]; then analyser14 ; fi

useFoundPCA_mdDavis()
{
	echo "$demA"$' Extracting principal components from 2d_projection data...\n\n'
	sleep 2
	
	cat $exist2dPCA | grep -v "^[@#]" | awk '{print $1}' > PC1_datapoints.dat
	cat $exist2dPCA | grep -v "^[@#]" | awk '{print $2}' > PC2_datapoints.dat

	echo $' Extract principal components from 2d_projection data...DONE'"$demB"
	sleep 2

	echo "$demA"$' Preparing parameters for free energy landscape calculations...\n'
	sleep 2

	checksimtime="SimTime.dat"
	if [[ ! -f "$checksimtime" ]] ; then
		# fetch simulation time from the trajectory using the ScanTRAJ fxn
		echo $' Extracting the simulation time-points from the trajectory...\n'
		ScanTRAJ
		increment_factor=$(awk "BEGIN {print $simDuratnINTns / $No_of_frames}")
		simtimeRecorded=0
		echo "$simtimeRecorded" > SimTime.dat
		while [[ "$simtimeRecorded" != "$simDuratnINTns" ]]; do
			simtimeRecorded=$(awk "BEGIN {print $simtimeRecorded + $increment_factor}")
			echo "$simtimeRecorded" >> SimTime.dat
			if [[ "$simtimeRecorded" == "$simDuratnINTns" ]]; then
				break
			fi
		done
		echo $' Extract simulation time-points from the trajectory...DONE'"$demB" ; sleep 2
	fi

	currentFELmdDavisPCAdir="$(pwd)""/PCA_3D-FEL_mdDavis"
	nFELPCA=1
	bkupFELmdDavisPCAdir="$(pwd)""/#PCA_3D-FEL_mdDavis"".""backup.""$nFELPCA"
	base_bkupFELmdDavisPCAdir=$(basename "$bkupFELmdDavisPCAdir")
	if [[ -d "$currentFELmdDavisPCAdir" ]]; then
		echo $'\n'"$currentFELmdDavisPCAdir"$' folder exists,\n'"backing it up as $base_bkupFELmdDavisPCAdir" ; sleep 1
		while [[ -d "$bkupFELmdDavisPCAdir" ]]; do
			nFELPCA=$(( nFELPCA + 1 )); bkupFELmdDavisPCAdir="$(pwd)""/#PCA_3D-FEL_mdDavis"".""backup.""$nFELPCA"
			base_bkupFELmdDavisPCAdir=$(basename "$bkupFELmdDavisPCAdir")
		done
		echo $'\n'"Backing up the last PCA_3D-FEL_mdDavis folder and its contents as $base_bkupFELmdDavisPCAdir" ; sleep 1
		mv "$currentFELmdDavisPCAdir" "$bkupFELmdDavisPCAdir" && mkdir ./PCA_3D-FEL_mdDavis || true
	elif [[ ! -d "$currentFELmdDavisPCAdir" ]]; then mkdir PCA_3D-FEL_mdDavis
	fi

	# fesFigure="PCA_FES.png"
	results_folder="PCA_3D-FEL_mdDavis"
	outName="PCA_"

	# cat $exist2dPCA | grep "^[@#]" > ./"$results_folder"/PC1.xvg
	# cat $exist2dPCA | grep "^[@#]" > ./"$results_folder"/PC2.xvg

	paste SimTime.dat PC1_datapoints.dat > PC1.dat
	paste SimTime.dat PC2_datapoints.dat > PC2.dat

	cat PC1.dat >> ./"$results_folder"/PC1.xvg
	cat PC2.dat >> ./"$results_folder"/PC2.xvg

	rm PC1_datapoints.dat PC2_datapoints.dat PC1.dat PC2.dat
	OrderParameter1="${results_folder}/PC1.xvg"
	OrderParameter2="${results_folder}/PC2.xvg"
}

useFoundRgRMSData_mdDavis()
{
	echo "$demA"$' Preparing parameters for 3D free energy landscape calculations...\n'
	sleep 2

	currentFELmdDavisRgVsRMSDdir="$(pwd)""/RgVsRMSD_3D-FEL_mdDavis"
	nFELRgVsRMSD=1
	bkupFELmdDavisRgVsRMSDdir="$(pwd)""/#RgVsRMSD_3D-FEL_mdDavis"".""backup.""$nFELRgVsRMSD"
	base_bkupFELmdDavisRgVsRMSDdir=$(basename "$bkupFELmdDavisRgVsRMSDdir")
	if [[ -d "$currentFELmdDavisRgVsRMSDdir" ]]; then
		echo $'\n'"$currentFELmdDavisRgVsRMSDdir"$' folder exists,\n'"backing it up as $base_bkupFELmdDavisRgVsRMSDdir" ; sleep 1
		while [[ -d "$bkupFELmdDavisRgVsRMSDdir" ]]; do
			nFELRgVsRMSD=$(( nFELRgVsRMSD + 1 )); bkupFELmdDavisRgVsRMSDdir="$(pwd)""/#RgVsRMSD_3D-FEL_mdDavis"".""backup.""$nFELRgVsRMSD"
			base_bkupFELmdDavisRgVsRMSDdir=$(basename "$bkupFELmdDavisRgVsRMSDdir")
		done
		echo $'\n'"Backing up the last RgVsRMSD_3D-FEL_mdDavis folder and its contents as $base_bkupFELmdDavisRgVsRMSDdir" ; sleep 1
		mv "$currentFELmdDavisRgVsRMSDdir" "$bkupFELmdDavisRgVsRMSDdir" && mkdir ./RgVsRMSD_3D-FEL_mdDavis || true
	elif [[ ! -d "$currentFELmdDavisRgVsRMSDdir" ]]; then mkdir RgVsRMSD_3D-FEL_mdDavis
	fi

	OrderParameter1="$existRg"
	OrderParameter2="$existRMSD"
	results_folder="RgVsRMSD_3D-FEL_mdDavis"
	outName="RgVsRMSD_"

}

analyser15()
{	
	echo "$demA"$' Constructing a 3D plot of the FES using md-davis...\n'
	order_parameters
 	if [[ $orderPair_choice == 1 && $automode == "semi" ]] ; then
		if [[ "$PCFile" == 1 ]]; then useFoundPCA_mdDavis
		elif [[ "$PCFile" == 2 ]]; then analyser7; useFoundPCA_mdDavis	
		elif [[ "$PCFile" == 3 ]]; then echo ""
			read -p ' Provide the path to the pre-calculated 2d_PCA projection file: ' precalcPCfile
			exist2dPCA="$precalcPCfile"
			useFoundPCA_mdDavis
		fi			
 	elif [[ $orderPair_choice == 1 && $automode == "full" ]] ; then useFoundPCA_mdDavis
	elif [[ $orderPair_choice == 2 && $automode == "semi" ]] ; then
		if [[ "$PCFile" == 1 ]]; then useFoundRgRMSData_mdDavis
		elif [[ "$PCFile" == 2 ]]; then analyser2; analyser4; useFoundRgRMSData_mdDavis
		elif [[ "$PCFile" == 3 ]]; then echo ""

			read -p ' Provide the path to the pre-calculated Rg.xvg data: ' precalcRg
			echo ""
			read -p ' Provide the path to the pre-calculated RMSD.xvg data: ' precalcRMSD

			existRg="$precalcRg"
			existRMSD="$precalcRMSD"	
			useFoundRgRMSData_mdDavis
 		fi
	elif [[ $orderPair_choice == 2 && $automode == "full" ]] ; then
		useFoundRgRMSData_mdDavis
 	elif [[ $orderPair_choice == 3 ]] ; then
		read -p ' Provide the path to the 1st order_parameter.xvg data: ' OrderParameter1
		echo ""
		read -p ' Provide the path to the 2nd order_parameter.xvg data: ' OrderParameter2
		
		echo "$demA"$' Preparing parameters for 3D free energy landscape calculations...\n'
		sleep 2
		
		results_folder="3D-FEL_mdDavis"
		outName=""

		currentFELmdDavisdir="$(pwd)""/3D-FEL_mdDavis"
		nFEL=1
		bkupFELmdDavisdir="$(pwd)""/#3D-FEL_mdDavis"".""backup.""$nFEL"
		base_bkupFELmdDavisdir=$(basename "$bkupFELmdDavisdir")
		if [[ -d "$currentFELmdDavisdir" ]]; then
			echo $'\n'"$currentFELmdDavisdir"$' folder exists,\n'"backing it up as $base_bkupFELmdDavisdir"
			sleep 1
			while [[ -d "$bkupFELmdDavisdir" ]]; do
				nFEL=$(( nFEL + 1 )); bkupFELmdDavisdir="$(pwd)""/#3D-FEL_mdDavis"".""backup.""$nFEL"
				base_bkupFELmdDavisdir=$(basename "$bkupFELmdDavisdir")
			done
			echo $'\n'"Backing up the last 3D-FEL_mdDavis folder and its contents as $base_bkupFELmdDavisdir"
			sleep 1
			mv "$currentFELmdDavisdir" "$bkupFELmdDavisdir" && mkdir ./3D-FEL_mdDavis || true
		elif [[ ! -d "$currentFELmdDavisdir" ]]; then mkdir 3D-FEL_mdDavis
		fi
	fi

	echo "$demA"$' Constructing 3D plot of the md-davis free energy surface...\n'
	sleep 2

	md-davis landscape_xvg -T $Temp -x ${OrderParameter1} -y ${OrderParameter2} -n \
	"3D_FEL_for_${filenm}" -l "3D_FEL_for_${filenm}" -o ${results_folder}/${filenm}_3D_${outName}FEL.html

	if [[ $orderPair_choice == 2 ]] ; then
		rm RgVsRMSD.xvg RgData.dat RMSData.dat
	fi
	echo $'\n Construct 3D plot of the md-davis free energy surface...DONE'"$demB"
	sleep 2
}

if [[ "$analysis" == *" 15 "* ]]; then analyser15 ; fi

getHBdataforMatrix()
{
	if [[ $sysType == "protein_only" ]]; then
        existXVGin="$(pwd)""/hbond/hbnum_intraPro_${filenm}.xvg"
        existMATRIXin="$(pwd)""/hbond/hb_matrix_intraPro_${filenm}.xpm"
        existINDEXin="$(pwd)""/hbond/hb_index_intraPro_${filenm}.ndx"
    elif [[ $sysType == "protein_lig" ]]; then
        existXVGin="$(pwd)""/hbond/hbnum_ProLig_${filenm}.xvg"
        existMATRIXin="$(pwd)""/hbond/hb_matrix_ProLig_${filenm}.xpm"
        existINDEXin="$(pwd)""/hbond/hb_index_ProLig_${filenm}.ndx"
    elif [[ $sysType == "protein_dna" ]]; then
        existXVGin="$(pwd)""/hbond/hbnum_Pro_DNA_${filenm}.xvg"
        existMATRIXin="$(pwd)""/hbond/hb_matrix_Pro_DNA_${filenm}.xpm"
        existINDEXin="$(pwd)""/hbond/hb_index_Pro_DNA_${filenm}.ndx"
    fi
}

hbondMatrix_useFoundData()
{
	currentHBmatrixdir="$(pwd)""/hbond_matrix"
	nHBmatrix=1
	bkupHBmatrixdir="$(pwd)""/#hbond_matrix"".""backup.""$nHBmatrix"
	base_bkupHBmatrixdir=$(basename "$bkupHBmatrixdir")
	if [[ -d "$currentHBmatrixdir" ]]; then
		echo $'\n'"$currentHBmatrixdir"$' folder exists,\n'"backing it up as $base_bkupHBmatrixdir"
		sleep 1
		while [[ -d "$bkupHBmatrixdir" ]]; do
			nHBmatrix=$(( nHBmatrix + 1 ))
			bkupHBmatrixdir="$(pwd)""/#hbond_matrix"".""backup.""$nHBmatrix"
			base_bkupHBmatrixdir=$(basename "$bkupHBmatrixdir")
		done
		echo $'\n'"Backing up the last hbond_matrix folder and its contents as $base_bkupHBmatrixdir"
		sleep 1
		mv "$currentHBmatrixdir" "$bkupHBmatrixdir" && mkdir ./hbond_matrix || true
	elif [[ ! -d "$currentHBmatrixdir" ]]; then mkdir hbond_matrix
	fi

	echo "$demA"$' Extracting a reference structure from the trajectory...\n\n\n'
	sleep 2
	if [[ $automode == "full" && $sysType == "protein_only" ]]; then
		echo 1 | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o hbond_matrix/${filenm}_referenceStructure.pdb -dump 0
	elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
		eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o hbond_matrix/${filenm}_referenceStructure.pdb -dump 0
	elif [[ $automode == "semi" || $automode == "full" ]] && [[ $sysType == "protein_dna" ]]; then
		echo "Protein_DNA" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -o hbond_matrix/${filenm}_referenceStructure.pdb -dump 0
	elif [[ $automode == "semi" && $sysType == "protein_lig" ]]; then
		eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -o hbond_matrix/${filenm}_referenceStructure.pdb -dump 0
	elif [[ $automode == "full" && $sysType == "protein_lig" ]]; then
		echo "Protein_$ligname" | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -o hbond_matrix/${filenm}_referenceStructure.pdb -dump 0
	fi
	echo "$demA"$' Extract a reference structure from the trajectory...DONE\n'
	echo $' Preparing the H-bond matrix...\n'
	sleep 2
	echo $' Detecting the list of H-bonds from the gmx index file...\n'
	sleep 2
	hbondList=$(grep "hbond" ${existINDEXin} | awk '{print $2}')
	echo " Calculating and saving the H-bond data..."
	sleep 2
	md-davis hbond -x ${existMATRIXin} -i ${existINDEXin} -s hbond_matrix/${filenm}_referenceStructure.pdb --save_pickle hbond_matrix/hb_data.p -g ${hbondList} > hbond_matrix/hb_counts.dat

	if [[ $automode == "semi" ]]; then
cat << hbcutAsk

  Proceed to plotting the H-bond matrix with a 33 % occurrence cut-off?

    1) Yes, set 33 % as the cut-off for the H-bonds to be included in the plot
    2) No, I want to set a different percentage cut-off

hbcutAsk
	
	read -p '  Enter 1 or 2 here: ' HBcutOff_ask
		while [[ "$HBcutOff_ask" != 1 && "$HBcutOff_ask" != 2 ]]; do
			echo $'\nYou entered: '"$HBcutOff_ask"
			echo $'Please enter a valid number (1 or 2)!!\n'
			read -p '  Enter 1 or 2 here: ' HBcutOff_ask
		done
	elif [[ $automode == "full" ]]; then HBcutOff_ask=1
	fi	
	if [[ $HBcutOff_ask == 1 ]] ; then
		echo $'\n Plotting the hydrogen bonds matrix...\n'
		sleep 2

		md-davis plot_hbond --percent --total_frames 101 --cutoff 33 \
		-o hbond_matrix/hbond_matrix_${filenm}.html hbond_matrix/hb_data.p
	elif [[ $HBcutOff_ask == 2 ]] ; then
		read -p ' Enter the percentage occurrence cut-off to use: ' HBcutOff
		echo $'\n Plotting the hydrogen bonds matrix...\n'
		sleep 2

		md-davis plot_hbond --percent --total_frames 101 --cutoff ${HBcutOff} \
		-o hbond_matrix/hbond_matrix_${filenm}.html hbond_matrix/hb_data.p
	fi
	echo $' Preparing the H-bond matrix...DONE'"$demB"
}

analyser16()
{
	getHBdataforMatrix
	if [[ -f "$existXVGin" ]] && [[ -f "$existMATRIXin" ]] && [[ -f "$existINDEXin" ]]
		then
		if [[ $automode == "semi" ]] ; then
			echo "$demA"$' Pre-calculated H-bond data found!\n\nFiles found:'\
			$'\n'"$existXVGin"$'\n'"$existMATRIXin"$'\n'"$existINDEXin"
			sleep 2
cat << askFELuseexist

 Do you want to use these data for plotting H-bond matrix?

   1) Yes, use the data above
   2) No, repeat H-dond analysis with gmx and use the new output
   3) No, I want to provide some other files to be used

askFELuseexist

			read -p ' Enter 1, 2 or 3 here: ' PCFile
			while [[ "$PCFile" != 1 && "$PCFile" != 2 && "$PCFile" != 3 ]]; do
				echo $'\nYou entered: '"$PCFile"
				echo $'Please enter a valid number (1, 2 or 3)!!\n'
				read -p ' Enter 1, 2 or 3 here: ' PCFile
			done

			if [[ $PCFile == 1 ]] ; then echo ""
			elif [[ $PCFile == 2 ]] ; then analyser5 ; getHBdataforMatrix
			elif [[ $PCFile == 3 ]] ; then
				read -p ' Provide the path to the pre-calculated H-bond .xvg data file: ' existXVGin
				echo ""
				read -p ' Provide the path to the pre-calculated H-bond matrix file: ' existMATRIXin
				echo ""
				read -p ' Provide the path to the pre-calculated H-bond index file: ' existINDEXin
			fi
		elif [[ $automode == "full" ]] ; then
			echo $'CHAPERONg: Pre-calculated H-bond data found!\n\nFiles found:'\
			$'\n'"$existXVGin"$'\n'"$existMATRIXin"$'\n'"$existINDEXin"\
			$'\n\n *CHAPERONg in auto mode; the data will be used for plotting H-bond matrix'
			sleep 2
			getHBdataforMatrix
		fi
	elif [[ -f "$existXVGin" ]] || [[ -f "$existMATRIXin" ]] || [[ -f "$existINDEXin" ]]
		then analyser5
	fi
	hbondMatrix_useFoundData

}

if [[ "$analysis" == *" 16 "* ]]; then analyser16 ; fi

analyser16()
{	
read -p '*Please enter the number of frames to skip at intervals: ' ski
if [[ "$ski" != "0" ]]; then skp="-skip ""$ski"
fi
echo "You entered: $ski"$'\n'
echo "$demA"$' Now extracting chosen frames...\n'
if [[ $automode == "full" && $sysType == "protein_only" ]]; then
	echo 0 | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o "${filenm}""_trjSystem_Every""$ski""skip.xtc" $skp
	echo 1 | eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o "${filenm}""_trjProtein_Every""$ski""skip.xtc" $skp
elif [[ $automode == "semi" || $automode == "full" ]] && [[ $sysType == "protein_lig" ]]; then
	eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}""_trj_Every""$ski""skip.xtc" $skp
elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
	eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o "${filenm}""_trj_Every""$ski""skip.xtc" $skp
else
	eval "$gmx_exe_path" trjconv -f "${filenm}"_${wraplabel}.xtc -s "${filenm}".tpr -o "${filenm}""_trj_Every""$ski""skip.xtc" $skp
fi
echo "$demA"$' Extract frames...DONE'"$demB"
sleep 2
}

#defining a function for make_ndx

makeNDXGroup2()
{
echo "$demA"$' Will now make index group(s)'"$demB"
sleep 2

read -p '**Provide a filename for the index group to be made: ' nameForIndex
echo "*You entered: ${nameForIndex}"$'\n\n**Select the groups to be indexed in the coming prompts\n\n'
sleep 2

ndxNAME="$nameForIndex"".ndx"

eval "$gmx_exe_path" make_ndx -f em.gro -o $ndxNAME

echo "$demA"" Make index group ${nameForIndex}... DONE""$demB"
sleep 2
}

if [[ "$analysis" == *" 17 "* ]]; then analyser16 ; fi

if [[ "$analysis" == *" 18 "* ]]; then makeNDXGroup2 ; fi

if [[ "$analysis" == *" 19 "* ]]
then analyser0; ScanTRAJ; analyser1; analyser2;	analyser3; analyser4; analyser5
	analyser6; analyser7; analyser8; analyser9; variables_for_regMD_Movie
	analyser11; analyser12; analyser13; analyser14; analyser15; analyser16
elif [[ "$analysis" == *" 20 "* ]]
then ScanTRAJ; analyser1; analyser2; analyser3; analyser4; analyser5
	analyser6; analyser7; analyser8; analyser9; variables_for_regMD_Movie
	analyser11; analyser12; analyser13; analyser14; analyser15; analyser16
elif [[ "$analysis" == *" 21 "* ]]
then ScanTRAJ; analyser1; analyser2; analyser3; analyser4; analyser5
	analyser6; analyser7; analyser8; analyser9; analyser12; analyser13
	analyser14; analyser15; analyser16
fi
}
