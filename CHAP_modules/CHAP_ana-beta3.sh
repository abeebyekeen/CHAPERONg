#! /bin/bash

#CHAP_ana - The trajectory/end-point analysis module of CHAPERONg
#CHAPERONg - An automation program for GROMACS md simulation
#Author: Abeeb A. Yekeen
#Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date: 2022.02.11

set -e
set -o pipefail

#set version
CHAPERONg_version="beta3"

#demA=$'\n\n'"#**********************************CHAPERONg**********************************#"$'\n'
demA=$'\n\n'"#================================= CHAPERONg =================================#"$'\n'
demB=$'\n'"#=============================================================================#"$'\n\n'

if [[ $initiator2 != 'avail' ]] ; then
	echo "$demA"$'Do not run modules independently!\nLaunch CHAPERONg with run_CHAPERONg-<version>!!'"$demB"	
	exit 1
fi	

filesuffx=''

#import module with collected parameters
#. CHAP_modules/CHAP_colPar-4.02.sh

#call module with defined fxns
. CHAP_modules/CHAP_deffxn-${CHAPERONg_version}.sh


if [[ "$sysType" == 1 ]]; then wraplabel="noPBC"
elif [[ "$sysType" == 2 ]]; then wraplabel="center"
elif [[ "$sysType" == 3 ]]; then wraplabel="center"
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
if [[ -f "pbcmol" ]] && [[ ! -f "pbcjump" && ! -f "pbcfit" &&  ! -f "pbccenter" ]] ; then wraplabel="noPBC" ; rm pbcmol
elif [[ -f "pbcjump" ]] && [[ ! -f "pbcmol" && ! -f "pbcfit" &&  ! -f "pbccenter" ]] ; then wraplabel="nojump" ; rm pbcjump
elif [[ -f "pbcfit" ]] && [[ ! -f "pbcmol" && ! -f "pbcjump" &&  ! -f "pbccenter" ]] ; then wraplabel="fit" ; rm pbcfit
elif [[ -f "pbccenter" ]] && [[ ! -f "pbcmol" && ! -f "pbcjump" &&  ! -f "pbcfit" ]] ; then wraplabel="center" ; rm pbccenter
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
  9     Make a movie of the simulation
  10    Free energy calculations using the MMPBSA method (g_mmpbsa)
  11    Construct free energy landscape with gmx sham
  12    Construct free energy surface using an external free energy function
  13    Construct free energy landscape using md-davis
  14    Extract frames from the trajectory
  15    Make index groups (make_ndx)
  16    All but 14 and 15
  17    All but 0, 14 and 15
  
AnalysisList

read -p '*Enter one or more combinations of the options here (separated by a space): ' analyse

analysis=" $analyse "

while [[ "$analysis" != *" 0 "* && "$analysis" != *" 1 "* && "$analysis" != *" 2 "* && \
	"$analysis" != *" 3 "* && "$analysis" != *" 4 "* && "$analysis" != *" 5 "* && \
	"$analysis" != *" 6 "* && "$analysis" != *" 7 "* && "$analysis" != *" 8 "* && \
	"$analysis" != *" 9 "* && "$analysis" != *" 10 "* && "$analysis" != *" 11 "* && \
	"$analysis" != *" 12 "* && "$analysis" != *" 13 "* && "$analysis" != *" 14 "* && \
	"$analysis" != *" 15 "* && "$analysis" != *" 16 "* ]] ; do
		echo $'\nYou entered: '"$analyse"$'\n'
		echo $'Please enter a valid number!!\n'
		read -p '*Enter one or more combinations of the options here (separated by a space): ' analyse
		analysis=" $analyse "
done

if [[ "$analyse" == "9" ]]; then
	analysis="$analyse"
fi

if [[ "$coordinates" == '' ]]; then
	coordinates="$filenm"
fi

ScanTRAJ()
{
if [[ ! -f "trajectDetails.log" ]]; then
	echo "$demA"$' Checking the trajectory to extract info about\n number of frames and simulation time'"$demB"
	eval $gmx_exe_path check -f "${filenm}"_"${wraplabel}".xtc |& tee trajectDetails.log
	No_of_frames=$(cat trajectDetails.log | grep "Last" | awk '{print $(NF-2)}')
	simDuratnps=$(cat trajectDetails.log | grep "Last" | awk '{print $NF}')
	#simDuratn_nsFloat=$(echo "${simDuratnps%\.*} / 1000" | bc -l)

	simDuratn_nsFloat=$(awk "BEGIN {print $simDuratnps / 1000}")
	simDuratnINTns=$(echo ${simDuratn_nsFloat%\.*})
	echo $simDuratnINTns > simulation_duration

	echo "$demA"$'\nExtract number of frames and simulation duration from trajectory...DONE'"$demB"
else
	No_of_frames=$(cat trajectDetails.log | grep "Last" | awk '{print $(NF-2)}')
	simDuratnps=$(cat trajectDetails.log | grep "Last" | awk '{print $NF}')
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
			mv ${coordinates}*"$filesuffx".png ${coordinates}*"$filesuffx".xvg ${coordinates}_"$filesuffx".png ./$AnaName || true
		mv ${coordinates}_"$filesuffx".xvg ${coordinates}_*"$filesuffx".png ./$AnaName || true
		mv ${coordinates}_*"$filesuffx".xvg ${coordinates}*"$filesuffx"*.xvg ${coordinates}*"$filesuffx"*.png ./$AnaName || true
		mv *"$filesuffx"*.png *"$filesuffx".png *"$filesuffx".xvg *"$filesuffx"*.xvg ./$AnaName || true
		mv "${coordinates}"*"$filesuffx"*".png" "${coordinates}"*"$filesuffx"".png" "${coordinates}"*"$filesuffx"".xvg" ./$AnaName || true
		mv ${coordinates}_"$filesuffx"*.png ${coordinates}_"$filesuffx"*.xvg ./$AnaName || true
		
	elif [[ ! -d "$currentAnadir" ]]; then mkdir ./$AnaName
		mv ${coordinates}*"$filesuffx".png ${coordinates}*"$filesuffx".xvg ${coordinates}_"$filesuffx".png ./$AnaName || true
		mv ${coordinates}_"$filesuffx".xvg ${coordinates}_*"$filesuffx".png ./$AnaName || true
		mv ${coordinates}_*"$filesuffx".xvg ${coordinates}*"$filesuffx"*.xvg ${coordinates}*"$filesuffx"*.png ./$AnaName || true
		mv *"$filesuffx"*.png *"$filesuffx".png *"$filesuffx".xvg *"$filesuffx"*.xvg ./$AnaName || true
		mv "${coordinates}"*"$filesuffx"*".png" "${coordinates}"*"$filesuffx"".png" "${coordinates}"*"$filesuffx"".xvg" ./$AnaName || true
		mv ${coordinates}_"$filesuffx"*.png ${coordinates}_"$filesuffx"*.xvg ./$AnaName || true
		#mv ${filenm}_Rg_ns.png ${filenm}_Rg_ns.xvg ${filenm}_Rg.xvg ./$AnaName || true
	fi
}
		
DNAwrapAlt()
{
echo "$demA""CHAPERONg could not find any index file. Centering on protein instead of Protein_DNA!""$demB"
echo "Protein" 0 | eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"${wraplabel}".xtc -center -pbc mol -ur compact
}
analyser0()
{	
echo "$demA"$' Now recentering the protein and rewrapping molecules within the unit cell...\n'
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 1 0 | eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"noPBC".xtc -pbc mol -center
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	echo 1 0 | eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
				
elif [[ $flw != 1 ]] && [[ $sysType == 1 ]]; then
	echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
	eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"noPBC".xtc -pbc mol -center
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
	sleep 1
	eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
		
elif [[ $flw == 1 ]] && [[ $sysType == 2 ]]; then
	echo 1 0 | eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
	echo "$demA"$' Now performing rotational and translational fitting...\n'
	echo 4 0 | eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o "${filenm}"_fit.xtc -fit rot+trans
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	echo 1 0 | eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact

elif [[ $flw != 1 ]] && [[ $sysType == 2 ]]; then
	echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
	eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
	echo "$demA"$' Now performing rotational and translational fitting...\n'
	echo $'**Choose "Backbone" (4) to perform lsq fitting to protein backbone, and "System" (0) for output when prompted\n'
		eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o "${filenm}"_fit.xtc -fit rot+trans

elif [[ $flw == 1 ]] && [[ $sysType == 3 ]]; then
	echo "Protein_DNA" 0 | eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"center".xtc -center -pbc mol -ur compact || DNAwrapAlt
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	echo "Protein_DNA" 0 | eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact || DNAwrapAlt
elif [[ $flw != 1 ]] && [[ $sysType == 3 ]]; then
	echo $'**Choose "Protein_DNA" for centering and "System" (0) for output when prompted\n'
	eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	eval $gmx_exe_path trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact || DNAwrapAlt
fi
echo "$demA"$' Recenter the protein and rewrap molecules within the unit cell...DONE'"$demB"
sleep 2
}
if [[ "$analysis" == *" 0 "* ]]; then analyser0; fi

analyser1()
{
	echo $'CHAPERONg: Program still under development. Check for an update later!\n   Thanks!!'
}
if [[ "$analysis" == *" 1 "* ]]; then analyser1; fi

altRMSD()
{
echo "$demA"$'There are multiple groups identified as '"$ligname"$'.\nCHAPERONg will try to guess the appropriate group to be used for '"$ligname"" RMSD calculations""$demB"

sleep 2
echo "$demA""CHAPERONg: Selecting group 13 for ""$ligname"$'.\nIf this is wrong, terminate and re-run RMSD analysis without the automation flag!'"$demB"

sleep 2

echo 13 13 | eval $gmx_exe_path rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_"$ligname"-rmsd.xvg -tu ns
}

analyser2()
{
echo "$demA"$' Now calculating RMSD...\n'
if [[ $sysType == 1 ]] || [[ $sysType == 3 ]] && [[ $flw == 1 ]] ; then
	echo "Backbone" "Backbone" | eval $gmx_exe_path rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${filenm}_BB-rmsd.xvg -tu ns
		
	gracebat ${filenm}_BB-rmsd.xvg -hdevice PNG -autoscale xy -printfile ${filenm}_BB-rmsd.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	
elif [[ $flw == 1 ]] && [[ $sysType == 2 ]]; then
	echo 4 4 | eval $gmx_exe_path rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${filenm}_BB-rmsd.xvg -tu ns
		
	gracebat ${filenm}_BB-rmsd.xvg -hdevice PNG -autoscale xy -printfile ${filenm}_BB-rmsd.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "$demA"$'Protein RMSD calculation... DONE\n  Now calculating ligand RMSD...\n'
	sleep 2
		
	echo "$ligname" "$ligname" | eval $gmx_exe_path rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o \
	${coordinates}_"$ligname"-rmsd.xvg -tu ns || altRMSD
				
	gracebat ${coordinates}_"$ligname"-rmsd.xvg -hdevice PNG -autoscale xy -printfile \
	${coordinates}_"$ligname"-rmsd.png -fixed 7500 4000 -legend load || notifyImgFail
		
else
	eval $gmx_exe_path rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_rmsd.xvg -tu ns
		
	gracebat ${coordinates}_rmsd.xvg -hdevice PNG -autoscale xy -printfile ${coordinates}_rmsd.png \
	-fixed 7500 4000 -legend load || notifyImgFail
fi
echo "$demA"$' Compute RMSD... DONE'"$demB"

AnaName="RMSD"
filesuffx="rmsd"
createDIR
echo "$demA"$' Generate a finished figure of the RMSD plot... DONE'"$demB"
}

if [[ "$analysis" == *" 2 "* ]]; then analyser2 ; fi

analyser3()
{
echo "$demA"$' Now calculating RMSF...\n'
if [[ $flw == 1 ]]; then
	echo "Backbone" "Backbone" | eval $gmx_exe_path rmsf -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_BB-rmsf.xvg -res
	echo "$demA"$' RMSF with backbone lsq fitting and calculation...DONE'"$demB"
	echo "C-alpha" "C-alpha" | eval $gmx_exe_path rmsf -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_Calpha-rmsf.xvg -res
	echo "$demA"$' RMSF with Calpha lsq fitting and calculation...DONE'"$demB"
	gracebat ${coordinates}_BB-rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${coordinates}_BB-rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat ${coordinates}_Calpha-rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${coordinates}_Calpha-rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat ${coordinates}_BB-rmsf.xvg ${coordinates}_Calpha-rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${coordinates}_BB-Calpha-rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
else
	eval $gmx_exe_path rmsf -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_rmsf.xvg -res
	echo "$demA"$' Compute RMSF...DONE'"$demB"
	gracebat ${coordinates}_rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${coordinates}_rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
fi
	
AnaName="RMSF"
filesuffx="rmsf"
createDIR
	
echo "$demA"$'Generate finished figure(s) of the RMSF plot(s)... DONE'"$demB"
}
if [[ "$analysis" == *" 3 "* ]]; then analyser3 ; fi
	
analyser4()
{
echo "$demA"$' Now calculating Rg...\n'
if [[ $flw == 1 ]]; then
	echo "Protein" | eval $gmx_exe_path gyrate -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${filenm}_Rg.xvg
	echo "$demA"$' Compute radius of gyration...DONE'"$demB"
	echo "$demA"$' Now converting Rg plot to ns format...\n'
	sleep 2
	grep "^[@#]" ${filenm}_Rg.xvg | sed "s/ps/ns/g" > ${filenm}_Rg_ns.xvg
	grep -v "^[@#]" ${filenm}_Rg.xvg | \
	awk '{print $1/1000"      "$2"      "$3"      "$4"     "$5}' >> ${filenm}_Rg_ns.xvg
else
	echo $'**In the following step, CHOOSE Protein (1) for Rg analysis\n\n'
	eval $gmx_exe_path gyrate -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${filenm}_Rg.xvg
	echo "$demA"$' Compute radius of gyration...DONE'"$demB"
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
#currentAnadir="$(pwd)""/$AnaName"
#	nDir=1
#	bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
#	if [[ -d "$currentAnadir" ]]; then
#		echo "$currentAnadir" "exists, backing it up as $bkupAnadir"
#		while [[ -d "$bkupAnadir" ]]; do
#		nDir=$(( nDir + 1 )); bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
#		done
#		mv "$currentAnadir" "$bkupAnadir" && mkdir ./$AnaName
#		echo "Backing up the last $AnaName folder and its contents as $bkupAnadir"
#		mv ${coordinates}_"$filesuffx"*.png ${coordinates}_"$filesuffx".xvg ${coordinates}_"$filesuffx"*.xvg ./$AnaName || true
#	elif [[ ! -d "$currentAnadir" ]]; then mkdir ./$AnaName
#mv ${filenm}_Rg_ns.png ${filenm}_Rg_ns.xvg ${filenm}_Rg.xvg ./$AnaName || true
#fi
#mv *"_Rg_ns.png" *"_Rg_ns.xvg" ./$AnaName || true
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
echo 1 13 | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_ProLig_${coordinates}.xvg -tu ns $hbthread
}

hbond_DNA1()
{
echo "$demA"$' Now executing Intra-protein hydrogen bonding analysis...\n'
sleep 2

echo "Protein" "Protein" | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -num hbnum_Pro_${coordinates}.xvg -tu ns $hbthread

echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Intra-DNA hydrogen bonding analysis...\n'
sleep 2

echo "DNA" "DNA" | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -num hbnum_DNA_${coordinates}.xvg -tu ns $hbthread

echo "$demA"$' Intra-DNA hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Protein-DNA hydrogen bonding analysis...\n'
echo "Protein" "DNA" | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -num hbnum_Pro_DNA_${coordinates}.xvg -tu ns $hbthread

echo "$demA"$' Protein-DNA hydrogen bonding analysis... DONE'"$demB"
sleep 2
}
hbond_DNA2()
{
echo "$demA"$' Now executing Intra-protein hydrogen bonding analysis...\n'
sleep 2

echo "Protein" "Protein" | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_Pro_${coordinates}.xvg -tu ns $hbthread
gracebat hbnum_Pro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
hbnum_Pro_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Intra-DNA hydrogen bonding analysis...\n'
sleep 2

echo "DNA" "DNA" | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_DNA_${coordinates}.xvg -tu ns $hbthread
gracebat hbnum_DNA_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
hbnum_DNA_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
echo "$demA"$' Intra-DNA hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Protein-DNA hydrogen bonding analysis...\n'
echo "Protein" "DNA" | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_Pro_DNA_${coordinates}.xvg -tu ns $hbthread
gracebat hbnum_Pro_DNA_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
hbnum_Pro_DNA_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
echo "$demA"$' Protein-DNA hydrogen bonding analysis... DONE'"$demB"
sleep 2
}

analyser5()
{
echo "$demA"$' Now executing H-bond analysis...\n'
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 1 1 | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_intraPro_${coordinates}.xvg -tu ns $hbthread
	echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	echo 1 "SOL" | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num \
	hbnum_Pro-SOL_${coordinates}.xvg -tu ns $hbthread || \
	echo " There are multiple groups with the name SOL. Skipping..."
	echo "$demA"$' Protein-SOL hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_intraPro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_intraPro_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat hbnum_Pro-SOL_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_Pro-SOL_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat hbnum_intraPro_${coordinates}.xvg hbnum_Pro-SOL_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_intraPro_Pro-SOL_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail

elif [[ $flw == 1 ]] && [[ $sysType == 2 ]]; then
	echo 1 "$ligname" | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_ProLig_${coordinates}.xvg -tu ns $hbthread || altHBOND
	echo "$demA"$' Protein-ligand hydrogen bonding analysis...DONE'"$demB"
	sleep 2

	echo 1 1 | eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_intraPro_${coordinates}.xvg -tu ns $hbthread
	echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_ProLig_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_ProLig_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat hbnum_intraPro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_intraPro_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == 1 ]] && [[ $flw == 0 ]] ; then
	eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_${coordinates}.xvg -tu ns $hbthread
	echo "$demA"$' Hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail	
elif [[ $sysType == 2 || "$sysType" == 3 ]] && [[ $flw == 0 ]] ; then
	eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -num hbnum_${coordinates}.xvg -tu ns $hbthread || \
	eval $gmx_exe_path hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_${coordinates}.xvg -tu ns $hbthread
	echo "$demA"$' Hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ "$sysType" == 3 ]] && [[ $flw == 1 ]] ; then
	hbond_DNA1 || hbond_DNA2
	echo "$demA"$' Hydrogen bonding analysis...DONE'"$demB"
	sleep 2
fi
	
currenthbonddir="$(pwd)""/hbond"
nhbond=1
bkuphbonddir="$(pwd)""/#hbond"".""backup.""$nhbond"
if [[ -d "$currenthbonddir" ]]; then
	echo $'\n'"$currenthbonddir"$' folder exists,\n'"backing it up as $bkuphbonddir"
	while [[ -d "$bkuphbonddir" ]]; do
	nhbond=$(( nhbond + 1 )); bkuphbonddir="$(pwd)""/#hbond"".""backup.""$nhbond"
	done
	mv "$currenthbonddir" "$bkuphbonddir" && mkdir ./hbond || true
	echo $'\n'"Backing up the last hbond folder and its contents as $bkuphbonddir"
	mv hbnum_*.png hbnum_*.xvg hbnum*.png ./hbond || true
elif [[ ! -d "$currenthbonddir" ]]; then
	mkdir ./hbond; mv hbnum_*.png hbnum_*.xvg hbnum*.png ./hbond || true
fi
echo "$demA"$'Generate finished figure(s) of the hbond plot(s)... DONE'"$demB"
}

if [[ "$analysis" == *" 5 "* ]]; then analyser5 ; fi

analyser6()
{
echo "$demA"$' Now calculating solvent accessible surface area (SASA)...\n'
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 1 | eval $gmx_exe_path sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o sasa_${coordinates}.xvg -tu ns
	gracebat sasa_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail		
elif [[ $sysType == 2 ]] || [[ "$sysType" == 3 ]] && [[ $flw == 0 ]] ; then
	eval $gmx_exe_path sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o sasa_${coordinates}.xvg -tu ns
	gracebat sasa_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == 2 ]] && [[ $flw == 1 ]]; then
	echo 1 | eval $gmx_exe_path sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o \
	sasa_Pro_${coordinates}.xvg -tu ns
	gracebat sasa_Pro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_Pro_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ "$sysType" == 3 ]] && [[ $flw == 1 ]]; then
	echo "Protein" | eval $gmx_exe_path sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o \
	sasa_Pro_${coordinates}.xvg -tu ns
	gracebat sasa_Pro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_Pro_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail	
	echo "$demA"$'Compute solvent accessible surface area (SASA) for DNA only...DONE'"$demB"	
	echo "DNA" | eval $gmx_exe_path sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o \
	sasa_DNA_${coordinates}.xvg -tu ns
	gracebat sasa_Pro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_DNA_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail	
	echo "$demA"$'Compute solvent accessible surface area (SASA) for DNA only...DONE'"$demB"	
	echo "$demA"$'Now calculating solvent accessible surface area (SASA) for Protein-DNA complex...\n'	
	echo "Protein_DNA" | eval $gmx_exe_path sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o \
	sasa_Pro_DNA_${coordinates}.xvg -tu ns
	gracebat sasa_Pro_DNA_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_Pro_DNA_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail	
	echo "$demA"$'Now calculating solvent accessible surface area (SASA) for Protein-DNA complex...DONE'"$demB"	
elif [[ $flw == 0 ]] && [[ $sysType == 1 ]]; then
	eval $gmx_exe_path sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o sasa_${coordinates}.xvg -tu ns
	gracebat sasa_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
fi
echo "$demA"$' Compute solvent accessible surface area (SASA)...DONE'"$demB"
currentSASAdir="$(pwd)""/SASA"
nSASA=1
bkupSASAdir="$(pwd)""/#SASA"".""backup.""$nSASA"
if [[ -d "$currentSASAdir" ]]; then
	echo $'\n'"$currentSASAdir"$' folder exists,\n'"backing it up as $bkupSASAdir"
	while [[ -d "$bkupSASAdir" ]]; do
	nSASA=$(( nSASA + 1 )); bkupSASAdir="$(pwd)""/#SASA"".""backup.""$nSASA"
	done
	mv "$currentSASAdir" "$bkupSASAdir" && mkdir ./SASA || true
	echo $'\n'"Backing up the last SASA folder and its contents as $bkupSASAdir"
	mv sasa*${coordinates}.png sasa*${coordinates}.xvg ./SASA || true
elif [[ ! -d "$currentSASAdir" ]]; then
	mkdir ./SASA; mv sasa*${coordinates}.png sasa*${coordinates}.xvg ./SASA || true
fi
	echo "$demA"$' Generate a finished figure of the SASA plot... DONE'"$demB"
}

if [[ "$analysis" == *" 6 "* ]]; then analyser6 ; fi

analyser7()
{
echo "$demA"$' Now running principal component analysis (PCA)...\n'
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 4 4 | eval $gmx_exe_path covar -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr \
	-o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
	echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
	echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'
	echo 4 4 | eval $gmx_exe_path anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_"${wraplabel}".xtc -eig \
	"${filenm}"_eigenval.xvg -s "${filenm}".tpr -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg
		
	gracebat PCA_2dproj_"${filenm}".xvg -hdevice PNG -autoscale xy -printfile \
	PCA_2dproj_"${filenm}".png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $flw == 0 ]] && [[ $sysType == 1 ]]; then
	echo $'**Choose "Backbone" (4) twice when prompted\n'
	eval $gmx_exe_path covar -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
	echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
	echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'\
	$'**Choose "Backbone" (4) twice when prompted\n'
	eval $gmx_exe_path anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_"${wraplabel}".xtc -eig "${filenm}"_eigenval.xvg \
	-s "${filenm}".tpr -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg
	gracebat PCA_2dproj_"${filenm}".xvg -hdevice PNG -autoscale xy -printfile \
	PCA_2dproj_"${filenm}".png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == 2 ]] || [[ "$sysType" == 3 ]] && [[ $flw == 1 ]] ; then
	echo "Backbone" "Backbone" | eval $gmx_exe_path covar -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx \
	-o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
	echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
	echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'
	echo "Backbone" "Backbone" | eval $gmx_exe_path anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_"${wraplabel}".xtc -eig \
	"${filenm}"_eigenval.xvg -s "${filenm}".tpr -n index.ndx -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg
	gracebat PCA_2dproj_"${filenm}".xvg -hdevice PNG -autoscale xy -printfile \
	PCA_2dproj_"${filenm}".png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == 2 ]] || [[ "$sysType" == 3 ]] && [[ $flw == 0 ]] ; then
	echo $'**Choose "Backbone" (4) twice when prompted\n'
	eval $gmx_exe_path covar -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}"_eigenval.xvg \
	-v "${filenm}"_eigenvec.trr
	echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
	echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'\
	$'**Choose "Backbone" (4) twice when prompted\n'
	eval $gmx_exe_path anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_"${wraplabel}".xtc -eig "${filenm}"_eigenval.xvg \
	-s "${filenm}".tpr -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg	
fi
echo "$demA"$' Principal component analysis (PCA)...DONE'"$demB"
currentPCAdir="$(pwd)""/PCA"
nPCA=1
bkupPCAdir="$(pwd)""/#PCA"".""backup.""$nPCA"
if [[ -d "$currentPCAdir" ]]; then
	echo $'\n'"$currentPCAdir"$' folder exists,\n'"backing it up as $bkupPCAdir"
	while [[ -d "$bkupPCAdir" ]]; do
	nPCA=$(( nPCA + 1 )); bkupPCAdir="$(pwd)""/#PCA"".""backup.""$nPCA"
	done
	mv "$currentPCAdir" "$bkupPCAdir" && mkdir ./PCA
	echo $'\n'"Backing up the last PCA folder and its contents as $bkupPCAdir"
	mv PCA_2dproj_*.png *eigenval.xvg PCA_2dproj_*.xvg *_eigenvec.trr covar.log average.pdb ./PCA || true
elif [[ ! -d "$currentPCAdir" ]]; then
	mkdir ./PCA
	mv PCA_2dproj_*.png *eigenval.xvg PCA_2dproj_*.xvg *_eigenvec.trr covar.log average.pdb dd?????? ./PCA || true
fi
echo "$demA"$'Generate finished figures of the PCA plots... DONE'"$demB"
}

if [[ "$analysis" == *" 7 "* ]]; then analyser7 ; fi
	
dsspCheck="Avail"
checkDSSP()
{
echo "$demA"$'DSSP not configured for gromacs!\nSkipping secondary structure analysis...\n'
sleep 2
dsspCheck="notAvail"
}

analyser8()
{

echo "$demA"$' Now computing secondary structure with DSSP...\n'

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

echo "MainChain" | eval $gmx_exe_path do_dssp -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o ss_"${filenm}".xpm -tu ns -dt ${dt_dssp} || checkDSSP
if [[ "$dsspCheck" == "Avail" ]] ; then
	echo "$demA"$' Compute secondary structure...DONE'"$demB"
	sleep 1
	
	echo "$demA"$' Detecting colour codes assigned in the eps file...\n'
	while IFS= read -r line; do
		scanSSname=$(echo "$line" | awk '{print $6}')
		if [[ "$scanSSname" == '"Coil"' ]] ; then coilCODEcoil=$(echo "$line" | awk '{print $3}')
		elif [[ "$scanSSname" == '"B-Sheet"' ]] ; then coilCODEbeta=$(echo "$line" | awk '{print $3}')
		elif [[ "$scanSSname" == '"A-Helix"' ]] ; then coilCODEhelix=$(echo "$line" | awk '{print $3}')
		fi
	done < ss_"${filenm}".xpm
	
	echo "$demA"$' Reducing colour coding to helix-sheet-turn-coil...\n'
	
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
	
	echo "$demA"$' Converting output ss_xpm to an eps file...\n'
	#eval $gmx_exe_path xpm2ps -f ss_"${filenm}"_HETC.xpm -o ss_"${filenm}"_colortype2.eps -rainbow blue || true
	eval $gmx_exe_path xpm2ps -f ss_"${filenm}"_HETC.xpm -o ss_"${filenm}"_colortype1.eps || true
	echo "$demA"$' Converting eps to pdf...\n' || true
	ps2pdf -sPAPERSIZE=ledger ss_"${filenm}"_colortype1.eps ss_"${filenm}"_colortype1size2.pdf || true
	#ps2pdf -sPAPERSIZE=ledger ss_"${filenm}"_colortype2.eps ss_"${filenm}"_colortype2size2.pdf || true
	ps2pdf ss_"${filenm}"_colortype1.eps ss_"${filenm}"_colortype1size1.pdf || true
	#ps2pdf ss_"${filenm}"_colortype2.eps ss_"${filenm}"_colortype2size1.pdf || true

	echo "$demA"$' Converting ss_eps to png...\n'
	
	convert ss_"${filenm}"_colortype1.eps -trim -bordercolor ss_"${filenm}"_colortype1.png || echo "$demA"$' The program "convert" not found!\n'

	convert ss_"${filenm}"_colortype1.eps -trim -bordercolor white -units pixelsperinch -density 600 ss_"${filenm}"_colortype1.png || echo "$demA"$' The program "convert" not found!\n'
	
	currentSecStrdir="$(pwd)""/Secondary_structure"
	nSecStr=1
	bkupSecStrdir="$(pwd)""/#Secondary_structure"".""backup.""$nSecStr"
	if [[ -d "$currentSecStrdir" ]]; then
		echo $'\n'"$currentSecStrdir"$' folder exists,\n'"backing it up as $bkupSecStrdir"
		while [[ -d "$bkupSecStrdir" ]]; do
		nSecStr=$(( nSecStr + 1 )); bkupSecStrdir="$(pwd)""/#Secondary_structure"".""backup.""$nSecStr"
		done
		mv "$currentSecStrdir" "$bkupSecStrdir" && mkdir ./Secondary_structure || true
		echo $'\n'"Backing up the last Secondary_structure folder and its contents as $bkupSecStrdir"
		mv scount.xvg ss_*.xpm ss_*.eps ss_*.pdf ss_*.png ./Secondary_structure || true
	elif [[ ! -d "$currentSecStrdir" ]]; then
		mkdir Secondary_structure; mv scount.xvg ss_*.xpm ss_*.eps ss_*.pdf ss_*.png ./Secondary_structure || true
	fi
	echo "$demA"$' Secondary structure analysis...DONE'"$demB"
fi
}

if [[ "$analysis" == *" 8 "* ]]; then ScanTRAJ; analyser8 ; fi

makeMoviePy1()
{
echo $'load PyMOLsession_allSet.pse\nmovie.produce dynamics_moviePy.mpg, quality 100\nquit' > make1_movie_Pyscript.pml
pymol make1_movie_Pyscript.pml
}
	
makeMoviePy2()
{
echo $'load PyMOLsession_allSet.pse\nmovie.produce dynamics_moviePy.mpg, ray, quality=100\nquit' > make2_movie_Pyscript.pml
pymol make2_movie_Pyscript.pml
}

##function specific for movie update
makeMoviePyx()
{
echo $'cd ./MOVIE\nload PyMOLsession_allSet.pse\nmovie.produce dynamics_moviePy.mpg, quality 100\nquit' > make1_movie_Pyscript.pml
pymol make1_movie_Pyscript.pml
}
##function specific for movie update	
makeMoviePyy()
{
echo $'cd ./MOVIE\nload PyMOLsession_allSet.pse\nmovie.produce dynamics_moviePy.mpg, ray, quality=100\nquit' > make2_movie_Pyscript.pml
pymol make2_movie_Pyscript.pml
}


analyser9()
{
echo "$demA"$'Preparing to make a summary movie of the trajectory'

cat << askMovielength

Do you want to proceed to making a movie summarized into 200 frames?

  1) Yes.
  2) No, I want a different number of frames for the movie.

askMovielength

read -p ' Enter 1 or 2 here: ' movieLeng
	while [[ "$movieLeng" != 1 && "$movieLeng" != 2 ]]; do
		echo $'\nYou entered: '"$movieLeng"
		echo $'Please enter a valid number (1 or 2)!!\n'
		read -p ' Enter 1 or 2 here: ' movieLeng
	done
	
if [[ $movieLeng == 2 ]] ; then
	cat trajectDetails.log | grep -v "GROMACS reminds"
	echo "$demA"$'Above is a summary of your simulation trajectory.\nUse the info to provide a response below.\n'
	read -p '*Please enter the number of frames to skip at intervals: ' skimov

	echo "You entered: $skimov"$'\n'

elif [[ $movieLeng == 1 ]] ; then
	if (( $No_of_frames >= 200 )) ; then 
		skimov_raw=$(awk "BEGIN {print $No_of_frames / 200}")
		skimov=$(echo ${skimov_raw%\.*})
		echo "$demA""Number of frames in the trajectory: ${No_of_frames}"$'\n'"${skimov} frames will be \
		skipped at intervals to produce a total of ~200 frames for movie""$demB"
		
	elif (( $No_of_frames < 200 )) ; then
		skimov=1
		echo "$demA""Number of frames in the trajectory: ${No_of_frames}"$'\n'"Total number of frames \
		in the trajectory is less than 200."$'\n'"No frames will be skipped.""$demB"
	fi
fi


echo "$demA"$' Will now extract frames to be used for the movie...'
sleep 2
#if [[ $sysType == 1 ]] || [[ $sysType == 2 ]] || [[ $sysType == 3 ]] && [[ $flw == 1 ]] ; then
echo 0 | eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -skip $skimov
#elif [[ $sysType == 1 ]] || [[ $sysType == 2 ]] && [[ $flw == 0 ]] ; then
#gmx trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -skip $skimov
#fi
sleep 2
echo "$demA"$' Preparing to extract 200 snapshots...\n'
sleep 2
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 1 | eval $gmx_exe_path trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -o summaryForMovie.pdb
elif [[ $flw == 0 ]] && [[ $sysType == 1 ]]; then
	eval $gmx_exe_path trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -o summaryForMovie.pdb
elif [[ $flw == 0 || $flw == 1 ]] && [[ $sysType == 3 ]]; then
	echo "Protein_DNA" | eval $gmx_exe_path trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -n index.ndx -o summaryForMovie.pdb
elif [[ $flw == 0 ]] && [[ $sysType == 2 ]]; then
	eval $gmx_exe_path trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -n index.ndx -o summaryForMovie.pdb
elif [[ $flw == 1 ]] && [[ $sysType == 2 ]]; then
	echo "Protein_$ligname" | eval $gmx_exe_path trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -n index.ndx -o summaryForMovie.pdb
fi
currentMOVIEdir="$(pwd)""/MOVIE"
nMOVIE=1
bkupMOVIEdir="$(pwd)""/#MOVIE"".""backup.""$nMOVIE"
if [[ -d "$currentMOVIEdir" ]]; then
	echo $'\n'"$currentMOVIEdir"$' folder exists,\n'"backing it up as $bkupMOVIEdir"
	while [[ -d "$bkupMOVIEdir" ]]; do
	nMOVIE=$(( nMOVIE + 1 )); bkupMOVIEdir="$(pwd)""/#MOVIE"".""backup.""$nMOVIE"
	done
	mv "$currentMOVIEdir" "$bkupMOVIEdir" && mkdir ./MOVIE || true
	echo $'\n'"Backing up the last MOVIE folder and its contents as $bkupMOVIEdir"
elif [[ ! -d "$currentMOVIEdir" ]]; then mkdir MOVIE
fi	
	
echo $'load summaryForMovie.pdb\nsave PyMOLsession.pse\nintra_fit name ca+c+n+o\n'\
$'preset.pretty(selection='"'""all""')"\
$'\nspectrum chain, green cyan orange magenta\ncolor atomic, (not elem C)\nbg white\n'\
$'set movie_loop, 0\nsmooth\norient\nviewport 760, 540\nzoom all, -10\n'\
$'set ray_trace_frames=1\nset ray_opaque_background, 0\nset cache_frames=0\n'\
$'mclear\ncd ./MOVIE\nsave PyMOLsession_allSet.pse\nmpng frame_.png\nquit' > prep_movie_Pyscript.pml
	
echo "$demA"$'Now, PyMOL will do the job. You sit back and have a cup of tea...Cheers!'"$demB"
sleep 2
pyM=0
pymol prep_movie_Pyscript.pml || pyM=1
echo "$demA"$'Extract frames as images...DONE'"$demB""$demA"$'Now converting images to movie...\n'
cd ./MOVIE
mov_make=0
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.gif || mov_make=1
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.mp4 || mov_make=2

if [[ "$mov_make" == 1 ]] && [[ "$pyM" == 0 ]]; then
	echo "$demA"$'The program "Convert/ImageMagick" could not be found.\nCHAPERONg detected'\
	$'"PyMOL" and will use it to make a movie which may, however, be of lesser quality'"$demB"
		
	makeMoviePy1
	makeMoviePy2
		
	echo "$demA"$'Movie (lesser quality) made with PyMOL...\n'
fi

if [[ "$mov_make" == 2 ]] && [[ "$pyM" == 0 ]]; then
	echo "$demA"$'A mp4 movie could not be made. This may be due to the program'\
	$'"Convert/ImageMagick" not being found, or some library is missing.'\
	$'\nCHAPERONg detected "PyMOL" and will use it to make a mp4 movie which may,'\
	$'however, be of lesser quality'"$demB"
		
	makeMoviePy1
	makeMoviePy2
		
	echo "$demA"$'Movie (lesser quality) made with PyMOL...\n'
fi

#mkdir ./frames
#mv frame_*.png ./frames || true ; mv ../summaryForMovie.pdb ./ || true; mv ../"${filenm}""_trjEvery""$skimov""skipForMovie.xtc" ./ || true ; mv ../*.pse ./ || true ; rm ../*movie_Pyscript.pml || true
rm frame_*.png || true ; rm ../summaryForMovie.pdb || true
rm ../"${filenm}""_trjEvery""$skimov""skipForMovie.xtc" || true 
mv ../*.pse ./ || true ; rm ../*movie_Pyscript.pml || true
rm ./PyMOLsession.pse || true
#rm prep_movie_Pyscript.pml ../prep_movie_Pyscript.pml
#rm ../prep_movie_Pyscript.pml || true; cd ..
#rm *movie_Pyscript.pml prep_movie_Pyscript.pml || true
echo "$demA"$' Convert images to movie...DONE'"$demB"
cd ..
}


analyser9update()
{

echo "$demA"$'Preparing to make a summary movie from a preset PyMOL session\n'
sleep 2

currentMOVIEdir="$(pwd)""/MOVIE"
if [[ ! -d "$currentMOVIEdir" ]]; then
	echo "No MOVIE directory from a previous run exists... Exiting"
	exit 1
fi
nMOVIE=1
currentMOVIEgif="$(pwd)""/MOVIE/dynamics_movie.gif"
bkupMOVIEgif="$(pwd)""/MOVIE/dynamics_movie_""backup""$nMOVIE"".gif"
if [[ -f "$currentMOVIEgif" ]]; then
	base_currentMOVIEgif=$(basename "$currentMOVIEgif")
	base_bkupMOVIEgif=$(basename "$bkupMOVIEgif")
	echo $'\n'"$base_currentMOVIEgif" "exists, backing it up as $base_bkupMOVIEgif"$'\n'
	while [[ -f "$bkupMOVIEgif" ]]; do
	nMOVIE=$(( nMOVIE + 1 )); bkupMOVIEgif="$(pwd)""/MOVIE/dynamics_movie_""backup""$nMOVIE"".gif"
	done
	mv "$currentMOVIEgif" "$bkupMOVIEgif" || true
	echo $'\n'"Backing up the last .gif MOVIE as $base_bkupMOVIEgif"
fi	

nMOVIE=1
currentMOVIEmp4="$(pwd)""/MOVIE/dynamics_movie.mp4"
bkupMOVIEmp4="$(pwd)""/MOVIE/dynamics_movie_""backup""$nMOVIE"".mp4"
if [[ -f "$currentMOVIEmp4" ]]; then
	base_currentMOVIEmp4=$(basename "$currentMOVIEmp4")
	base_bkupMOVIEmp4=$(basename "$bkupMOVIEmp4")
	echo $'\n'"$base_currentMOVIEmp4"" exists, backing it up as $base_bkupMOVIEmp4"$'\n'
	while [[ -f "$bkupMOVIEmp4" ]]; do
	nMOVIE=$(( nMOVIE + 1 )); bkupMOVIEmp4="$(pwd)""/MOVIE/dynamics_movie_""backup""$nMOVIE"".mp4"
	done
	mv "$currentMOVIEmp4" "$bkupMOVIEmp4" || true
	echo $'\n'"Backing up the last .mp4 MOVIE as $base_bkupMOVIEmp4"
fi	
	
echo $'cd ./MOVIE\nload PyMOLsession_allSet.pse\nmpng frame_.png\nquit' > prep_movie_Pyscript.pml
	
echo "$demA"$'Now, PyMOL will do the job. You sit back and have a cup of tea...Cheers!'"$demB"
sleep 2
pyM=0
pymol prep_movie_Pyscript.pml || pyM=1
echo "$demA"$'Extract frames as images...DONE'"$demB""$demA"$'Now converting images to movie...\n'
cd ./MOVIE
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
rm ../*movie_Pyscript.pml || true
echo "$demA"$' Convert images to movie...DONE'"$demB"
cd ..
}

if [[ "$analyse" == "9" ]]; then
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
 
if [[ "$moviechoic" == "a" ]]; then ScanTRAJ; analyser9
elif [[ "$moviechoic" == "b" ]]; then analyser9update; fi

fi

if [[ "$analysis" == *" 9 "* ]]; then ScanTRAJ; analyser9 ; fi
	
analyser10()
{
if [[ $mmGMXpath != '' ]] ; then	
	#if [[ $sysType == 1 ]]; then indexer=''
 	if [[ $sysType == 2 || $sysType == 3 ]] ; then indexer='-n index.ndx' ; fi
	echo "$demA"$'Preparing to generate input files for g_MMPBSA free energy calculations...\n'

	if [[ $trajFraction == '' ]]; then trajFraction=3 ; fi

	No_of_last_third_frames=$(awk "BEGIN {print $No_of_frames / $trajFraction}")
	No_of_last_third_framesINT=$(echo ${No_of_last_third_frames%\.*})

	trajFractionFactor=$(( trajFraction - 1 ))

	simDuratnps_lastFractn_begin=$(awk "BEGIN {print $simDuratnps * $trajFractionFactor / $trajFraction}")
	simDuratnps_lastFractn_beginINT=$(echo ${simDuratnps_lastFractn_begin%\.*})


	if [[ $mmpbframesNo == '' ]]; then mmpbframesNo=100 ;fi	
	
	if (( $No_of_last_third_framesINT >= $mmpbframesNo )) ; then
		skipframegpsa=$(awk "BEGIN {print $No_of_last_third_frames / $mmpbframesNo }")
		skipframegpsaINT=$(echo ${skipframegpsa%\.*})
		echo "$demA""Number of frames in the last third of the trajectory: ${No_of_last_third_frames}"$'\n'"${skipframegpsaINT} "\
		"frames will be skipped at intervals to produce a total of ~100 frames for g_mmpbsa calculations.""$demB"
		
	elif (( $No_of_last_third_framesINT < $mmpbframesNo )) ; then
		skipframegpsaINT=1
		echo "$demA""Number of frames in the last third of the trajectory: ${No_of_last_third_frames}"$'\n'"Total number of frames "\
		"in the trajectory is less than 200."$'\n'"No frames will be skipped for g_mmpbsa calculations.""$demB"
	fi

	if [[ "$mmGMX" == "1" ]] ; then
		echo "$demA"$' Preparing to generate a compatible .tpr for g_MMPBSA...\n'
		eval $mmGMXpath grompp -f md.mdp -c $coordinates_raw -p topol.top -o \
		"${filenm}"_TPR_for_g_mmpbsa.tpr $indexer
		echo "$demA"$'Generate a compatible .tpr for g_MMPBSA...DONE'"$demB"
		#fi
		echo "$demA"$'Generating a compatible fraction of the trajectory for g_MMPBSA...\n'
		if [[ $flw == 1 ]]; then
			echo 0 | eval $mmGMXpath trjconv -s "${filenm}"_TPR_for_g_mmpbsa.tpr -f "${filenm}"_"${wraplabel}".xtc \
			-o "${filenm}"_lastFractntraj4_mmpbsa.xtc -b $simDuratnps_lastFractn_beginINT
		
			echo "$demA"$'Generate a compatible last 3rd trajectory file for g_MMPBSA...DONE'"$demB"
		
			echo "$demA"$'Extracting 100 frames from the last 3rd of the trajectory...\n'
			echo 0 | eval $mmGMXpath trjconv -s "${filenm}"_TPR_for_g_mmpbsa.tpr -f "${filenm}"_lastFractntraj4_mmpbsa.xtc \
			-o "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -skip $skipframegpsaINT
		
			echo "$demA"$'Extract 100 frames from the last 3rd of the trajectory...DONE'"$demB"
		
		elif [[ $flw != 1 ]]; then
			eval $mmGMXpath trjconv -s "${filenm}"_TPR_for_g_mmpbsa.tpr -f "${filenm}"_"${wraplabel}".xtc \
			-o "${filenm}"_lastFractntraj4_mmpbsa.xtc -b $simDuratnps_lastFractn_beginINT
		
			echo "$demA"$'Generate a compatible fraction of the trajectory for g_MMPBSA...DONE'"$demB"
		
			echo "$demA"$'Extracting 100 frames from the last 3rd of the trajectory...\n'
		
			eval $mmGMXpath trjconv -s "${filenm}"_TPR_for_g_mmpbsa.tpr -f "${filenm}"_lastFractntraj4_mmpbsa.xtc \
			-o "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -skip $skipframegpsaINT
		
			echo "$demA"$'Extract 100 frames from the last 3rd of the trajectory...DONE'"$demB"
		
		fi	
	
		echo "$demA"$'Generate input files for g_MMPBSA free energy calculations...DONE'"$demB"

		echo "$demA"$' Now preparing to run g_MMPBSA calculations...\n'
		if [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw == 1 ]]; then
			echo 1 "$ligname" | ./utilities/g_mmpbsa_pkg/g_mmpbsa -f "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -s \
			"${filenm}"_TPR_for_g_mmpbsa.tpr -n index.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp
		elif [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw != 1 ]] ; then
			./utilities/g_mmpbsa_pkg/g_mmpbsa -f "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -s \
			"${filenm}"_TPR_for_g_mmpbsa.tpr -n index.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp
		fi
		echo "$demA"$'Run g_MMPBSA calculations...DONE'"$demB"
	elif [[ "$mmGMX" == '' ]] ; then
		echo "$demA"$' Now preparing to run g_MMPBSA calculations...\n'
		if [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw == 1 ]]; then
			echo 1 "$ligname" | g_mmpbsa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -i pbsa.mdp \
			-pdie 2 -pbsa -decomp || echo "$demA"$'g_mmpbsa failed to run. Ensure your environments are properly set...\n'
		elif [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw != 1 ]] ; then
			g_mmpbsa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -i pbsa.mdp \
			-pdie 2 -pbsa -decomp || echo "$demA"$'g_mmpbsa failed to run. Ensure your environments are properly set...\n'
		fi
	fi
	echo "$demA"$'Now calculating average binding energy & contribution of residues...\n'
	python ./utilities/g_mmpbsa_pkg/MmPbSaStat.py -m energy_MM.xvg -p polar.xvg -a apolar.xvg || \
	python3 ./utilities/g_mmpbsa_pkg/MmPbSaStatPy3.py -m energy_MM.xvg -p polar.xvg -a apolar.xvg || true

	python ./utilities/g_mmpbsa_pkg/MmPbSaDecomp.py -bs -nbs 2000 -m contrib_MM.dat -p contrib_pol.dat -a contrib_apol.dat || \
	python3 ./utilities/g_mmpbsa_pkg/MmPbSaDecompPy3.py -bs -nbs 2000 -m contrib_MM.dat -p contrib_pol.dat -a contrib_apol.dat || true

	if [[ "$mmGMX" == "1" ]] ; then
		if [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw == 1 ]]; then
			echo 1 "$ligname" | ./utilities/g_mmpbsa_pkg/energy2bfac -s "${filenm}"_TPR_for_g_mmpbsa.tpr -i energyMapIn.dat
		elif [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw != 1 ]] ; then
			./utilities/g_mmpbsa_pkg/energy2bfac -s "${filenm}"_TPR_for_g_mmpbsa.tpr -i energyMapIn.dat
		fi
	elif [[ "$mmGMX" == '' ]] ; then
		if [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw == 1 ]]; then
			echo 1 "$ligname" | energy2bfac -s "${filenm}"_TPR_for_g_mmpbsa.tpr -i energyMapIn.dat || \
			echo "$demA"$'energy2bfac failed to run. Ensure your environment are properly set...\n'
		elif [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw != 1 ]] ; then
			energy2bfac -s "${filenm}"_TPR_for_g_mmpbsa.tpr -i energyMapIn.dat || \
			echo "$demA"$'energy2bfac failed to run. Ensure your environment are properly set...\n'
		fi
	fi

	echo "$demA"$'Calculate average binding energy & contribution of residues...DONE'"$demB"

	AnaName="MMPBSA"
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
		rm "${filenm}"_TPR_for_g_mmpbsa.tpr "${filenm}"_lastFractntraj4_mmpbsa.xtc "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc || true
		mv energy_MM.xvg polar.xvg apolar.xvg contrib_MM.dat contrib_pol.dat contrib_apol.dat ./$AnaName || true
		mv full_energy.dat summary_energy.dat final_contrib_energy.dat energyMapIn.dat ./$AnaName || true
		mv complex.pdb subunit_1.pdb subunit_2.pdb ./$AnaName || true

	elif [[ ! -d "$currentAnadir" ]]; then mkdir ./$AnaName
		rm "${filenm}"_TPR_for_g_mmpbsa.tpr "${filenm}"_lastFractntraj4_mmpbsa.xtc "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc || true
		mv energy_MM.xvg polar.xvg apolar.xvg contrib_MM.dat contrib_pol.dat contrib_apol.dat ./$AnaName || true
		mv full_energy.dat summary_energy.dat final_contrib_energy.dat energyMapIn.dat ./$AnaName || true
		mv complex.pdb subunit_1.pdb subunit_2.pdb ./$AnaName || true
	fi

elif [[ $mmGMXpath == '' ]] ; then
	echo "$demA"$'GMX path for g_mmpbsa not set. Use the parFile option!\n'
fi

}

if [[ "$analysis" == *" 10 "* ]]; then ScanTRAJ; analyser10 ; fi

useFoundPCA_sham()
{
echo "$demA"$' Preparing PCA-derived FEL with gmx sham...\n'
sleep 1
			
eval $gmx_exe_path sham -f ./PCA/PCA_2dproj_$filenm.xvg -ls ./PCA/FEL_PCA_sham_$filenm.xpm -notime || true
		
eval $gmx_exe_path xpm2ps -f ./PCA/FEL_PCA_sham_$filenm.xpm -o ./PCA/FEL_PCA_sham_$filenm.eps -rainbow red || true
	
ps2pdf -sPAPERSIZE=ledger ./PCA/FEL_PCA_sham_$filenm.eps ./PCA/FEL_PCA_sham_${filenm}_landscape.pdf || true
		
ps2pdf ./PCA/FEL_PCA_sham_$filenm.eps ./PCA/FEL_PCA_sham_${filenm}_portrait.pdf || true

pdf2ppm -png -r 600 ./PCA/FEL_PCA_sham_${filenm}_portrait.pdf ./PCA/FEL_PCA_sham_${filenm}_portrait.png || true

convert ./PCA/FEL_PCA_sham_$filenm.eps -trim -bordercolor white ./PCA/FEL_PCA_sham_${filenm}_convertps2png.png || true

convert ./PCA/FEL_PCA_sham_$filenm.eps -trim -bordercolor white -units pixelsperinch -density 600 -resize 3000x5000 ./PCA/FEL_PCA_sham_${filenm}_convertps2pngfull.png || true
			
currentFELPCAdir="$(pwd)""/PCA_FEL_sham"
nFELpca=1
bkupFELpcadir="$(pwd)""/#PCA_FEL_sham"".""backup.""$nFELpca"
if [[ -d "$currentFELPCAdir" ]]; then
	echo $'\n'"$currentFELPCAdir"$' folder exists,\n'"backing it up as $bkupFELpcadir"
	while [[ -d "$bkupFELpcadir" ]]; do
	nFELpca=$(( nFELpca + 1 )); bkupFELpcadir="$(pwd)""/#PCA_FEL_sham"".""backup.""$nFELpca"
	done
	mv "$currentFELPCAdir" "$bkupFELpcadir" && mkdir ./PCA_FEL_sham || true
	echo $'\n'"Backing up the last PCA_FEL_sham folder and its contents as $bkupFELpcadir"
	mv ./PCA/FEL_PCA_sham_* enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./PCA_FEL_sham || true
elif [[ ! -d "$currentFELPCAdir" ]]; then mkdir PCA_FEL_sham
	mv ./PCA/FEL_PCA_sham_* enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./PCA_FEL_sham || true
fi
echo "$demA"$' Prepare Gibbs FEL with gmx sham...DONE'"$demB"
}

useFoundRgRMSData_sham()
{
echo "$demA"$' Preparing Rg Vs RMSD FEL using gmx sham...\n\n'
sleep 1
			
eval $gmx_exe_path sham -f RgVsRMSD.xvg -ls FEL_sham_RgVsRMSD_$filenm.xpm -notime || true
		
eval $gmx_exe_path xpm2ps -f FEL_sham_RgVsRMSD_$filenm.xpm -o FEL_sham_RgVsRMSD_$filenm.eps -rainbow red || true
	
ps2pdf -sPAPERSIZE=ledger FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_landscape.pdf || true
		
ps2pdf FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_portrait.pdf || true

pdf2ppm -png -r 600 FEL_sham_RgVsRMSD_${filenm}_portrait.pdf FEL_sham_RgVsRMSD_${filenm}_portrait.png || true

convert FEL_sham_RgVsRMSD_$filenm.eps -trim -bordercolor white FEL_sham_RgVsRMSD_${filenm}_convertps2png.png || true

convert FEL_sham_RgVsRMSD_$filenm.eps -trim -bordercolor white -units pixelsperinch -density 600 -resize 3000x5000 FEL_sham_RgVsRMSD_${filenm}_convertps2pngfull.png || true
			
currentFELshamRgVsRMSDdir="$(pwd)""/RgVsRMSD_FEL_sham"
nFELRgVsRMSD=1
bkupFELshamRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FEL_sham"".""backup.""$nFELRgVsRMSD"
if [[ -d "$currentFELshamRgVsRMSDdir" ]]; then
	echo $'\n'"$currentFELshamRgVsRMSDdir"$' folder exists,\n'"backing it up as $bkupFELshamRgVsRMSDdir"
	while [[ -d "$bkupFELshamRgVsRMSDdir" ]]; do
	nFELRgVsRMSD=$(( nFELRgVsRMSD + 1 )); bkupFELshamRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FEL_sham"".""backup.""$nFELRgVsRMSD"
	done
	mv "$currentFELshamRgVsRMSDdir" "$bkupFELshamRgVsRMSDdir" && mkdir ./RgVsRMSD_FEL_sham || true
	echo $'\n'"Backing up the last RgVsRMSD_FEL_sham folder and its contents as $bkupFELshamRgVsRMSDdir"
	mv FEL_sham_RgVsRMSD_* RMSData.dat RgData.dat RgVsRMSD.xvg enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./RgVsRMSD_FEL_sham || true
elif [[ ! -d "$currentFELshamRgVsRMSDdir" ]]; then mkdir RgVsRMSD_FEL_sham
	mv FEL_sham_RgVsRMSD_* RMSData.dat RgData.dat RgVsRMSD.xvg enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./RgVsRMSD_FEL_sham || true
fi
echo "$demA"$' Prepare Rg Vs RMSD FEL with gmx sham...DONE'"$demB"
}

order_parameters()
{
cat << orderPair

Which order parameter pair do you want to use for FEL calculations?

  1) Principal components
  2) Rg versus RMSD
  3) Other user-provided pair of order parameters

orderPair

read -p ' Enter 1, 2 or 3 here: ' orderPair_choice
	while [[ $orderPair_choice != 1 && $orderPair_choice != 2 && $orderPair_choice != 3 ]]
	do
		echo $'\nYou entered: '"$orderPair_choice"
		echo $'Please enter a valid number (1, 2 or 3)!!\n'
		read -p ' Enter 1, 2 or 3 here: ' orderPair_choice
	done

if [[ $orderPair_choice == 1 ]] ; then
	echo "$demA"$' Checking the working directory for pre-calculated PCA_2d projection data...'"$demB"
	sleep 2
	exist2dPCA="$(pwd)""/PCA/PCA_2dproj_""$filenm"".xvg"
	if [[ -f "$exist2dPCA" ]] ; then
		if [[ "$flw" == 0 ]] ; then
			echo $'CHAPERONg: Pre-calculated PCA_2d projection file found!\nFile found:'" $exist2dPCA"
			sleep 2
cat << askFELuseexist

Do you want to use this file for free energy landscape calculations?

  1) Yes, use the file above
  2) No, repeat PCA calculations and use the new output for FEL plot
  3) No, I want to provide another file to be used

askFELuseexist

			read -p ' Enter 1, 2 or 3 here: ' PCFile
			while [[ "$PCFile" != 1 && "$PCFile" != 2 && "$PCFile" != 3 ]]; do
				echo $'\nYou entered: '"$PCFile"
				echo $'Please enter a valid number (1, 2 or 3)!!\n'
				read -p ' Enter 1, 2 or 3 here: ' PCFile
			done
		elif [[ "$flw" == 1 ]] ; then
			echo "$demA"$' Pre-calculated PCA_2d projection data found!\n File found:'" $exist2dPCA"\
			$'\n *CHAPERONg in auto mode\n'" $exist2dPCA will be used for FEL plotting"
		fi
	elif [[ ! -f "$exist2dPCA" ]] ; then analyser7		
	fi
	inputPair="$exist2dPCA"
#if order parameter pair is RgVsRMSD	
elif [[ $orderPair_choice == 2 ]] ; then
	existRg="$(pwd)""/Rg/""${filenm}_Rg_ns.xvg"
	if [[ -f "$existRg" ]] ; then
		if [[ "$flw" == 0 ]] ; then
			echo "$demA"$' Pre-calculated Rg data found!\n File found:'" $existRg"$'\n'
			sleep 2
cat << askFELuseexist

Do you want to use this file for free energy landscape calculations?

  1) Yes, use the file above
  2) No, repeat Rg calculations and use the new output for FEL plot
  3) No, I want to provide another file to be used

askFELuseexist

			read -p ' Enter 1, 2 or 3 here: ' PCFile
			while [[ "$PCFile" != 1 && "$PCFile" != 2 && "$PCFile" != 3 ]]; do
				echo $'\nYou entered: '"$PCFile"
				echo $'Please enter a valid number (1, 2 or 3)!!\n'
				read -p ' Enter 1, 2 or 3 here: ' PCFile
			done
		elif [[ "$flw" == 1 ]] ; then
			echo "$demA"$' Pre-calculated Rg data found!\nFile found:'" $existRg"\
			$'\n *CHAPERONg in auto mode\n'" $existRg will be used for FEL plotting"
		fi
	elif [[ ! -f "$existRg" ]] ; then analyser4
	fi
	cat "$existRg" | grep -v "^[@#]" | awk '{print $2}' > RgData.dat

	inputRg_xvgData="$existRg"

	#check for pre-calculated RMSD data
	existRMSD="$(pwd)""/RMSD/""${filenm}_BB-rmsd.xvg"
	if [[ -f "$existRMSD" ]] ; then
		if [[ "$flw" == 0 ]] ; then
			echo "$demA"$' Pre-calculated RMSD data found!\n File found:'" $existRMSD"$'\n'
			sleep 2
cat << askFELuseexist

Do you want to use this file for free energy landscape calculations?

  1) Yes, use the file above
  2) No, repeat RMSD calculations and use the new output for FEL plot
  3) No, I want to provide another file to be used

askFELuseexist

			read -p ' Enter 1, 2 or 3 here: ' PCFile
			while [[ "$PCFile" != 1 && "$PCFile" != 2 && "$PCFile" != 3 ]]; do
				echo $'\nYou entered: '"$PCFile"
				echo $'Please enter a valid number (1, 2 or 3)!!\n'
				read -p ' Enter 1, 2 or 3 here: ' PCFile
			done
		elif [[ "$flw" == 1 ]] ; then
			echo "$demA"$' Pre-calculated RMSD data found!\n File found:'" $existRMSD"\
			$'\n *CHAPERONg in auto mode\n'" $existRMSD will be used for FEL plotting"
		fi
	elif [[ ! -f "$existRMSD" ]] ; then analyser2		
	fi
	cat "$existRMSD" | grep -v "^[@#]" | awk '{print $2}' > RMSData.dat
	echo "# This file contains the RMSD and Rg values extracted by CHAPERONg" > RgVsRMSD.xvg
	echo "# from the data generated by GROMACS..." >> RgVsRMSD.xvg
	echo "#" >> RgVsRMSD.xvg
	echo "@    title "$'"Plot of Rg against RMSD"' >> RgVsRMSD.xvg
	echo "@    xaxis  label "$'"RMSD (nm)"' >> RgVsRMSD.xvg
	echo "@    yaxis  label "$'"Rg (nm)"' >> RgVsRMSD.xvg
	echo "@TYPE xy" >> RgVsRMSD.xvg
	paste -d "        " RMSData.dat RgData.dat >> RgVsRMSD.xvg

	inputRMSD_xvgData="$existRMSD"
	inputPair="RgVsRMSD.xvg"
fi
}

analyser11()
{
	order_parameters
	if [[ $orderPair_choice == 1 && "$flw" == 0 ]] ; then
		if [[ "$PCFile" == 1 ]]; then useFoundPCA_sham
		elif [[ "$PCFile" == 2 ]]; then analyser7; useFoundPCA_sham	
		elif [[ "$PCFile" == 3 ]]; then echo ""
			read -p ' Provide the path to the pre-calculated 2d_PCA projection file: ' precalcPCfile
	
			echo "$demA"$'Preparing PCA-derived FEL with user-provided 2d_PCA file...\n\n'
			sleep 1

			felcal=0			

			eval $gmx_exe_path sham -f $precalcPCfile -ls FEL_PCA_sham_$filenm.xpm -notime || felcal=1
			if [[ "$felcal" == 0 ]] ; then
				
				eval $gmx_exe_path xpm2ps -f FEL_PCA_sham_$filenm.xpm -o FEL_PCA_sham_$filenm.eps -rainbow red || true

				ps2pdf -sPAPERSIZE=ledger FEL_PCA_sham_$filenm.eps FEL_PCA_sham_${filenm}_landscape.pdf || true

				ps2pdf FEL_PCA_sham_$filenm.eps FEL_PCA_sham_${filenm}_portrait.pdf || true
			
				currentFELPCAdir="$(pwd)""/PCA_FEL_sham"
				nFELpca=1
				bkupFELpcadir="$(pwd)""/#PCA_FEL_sham"".""backup.""$nFELpca"
				if [[ -d "$currentFELPCAdir" ]]; then
					echo $'\n'"$currentFELPCAdir"$' folder exists,\n'"backing it up as $bkupFELpcadir"
					while [[ -d "$bkupFELpcadir" ]]; do
						nFELpca=$(( nFELpca + 1 )); bkupFELpcadir="$(pwd)""/#PCA_FEL_sham"".""backup.""$nFELpca"
					done
					mv "$currentFELPCAdir" "$bkupFELpcadir" && mkdir ./PCA_FEL_sham || true
					echo $'\n'"Backing up the last PCA_FEL_sham folder and its contents as $bkupFELpcadir"
					mv FEL_PCA_sham_*.xpm enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg FEL_PCA_sham_*.eps FEL_PCA_sham_*.pdf ./PCA_FEL_sham || true
				elif [[ ! -d "$currentFELPCAdir" ]]; then mkdir PCA_FEL_sham
					mv FEL_PCA_sham* enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./PCA_FEL_sham || true
				fi
				echo "$demA"$' Prepare PCA-based FEL using gmx sham...DONE'"$demB"
			
			elif [[ "$felcal" == 1 ]] ; then
				echo $'Calculation failed.\n'" Please confirm that you have entered the right path/file as input!"
			fi
		fi			
	elif [[ $orderPair_choice == 1 && $flw == 1 ]] ; then useFoundPCA_sham
	elif [[ $orderPair_choice == 2 && $flw == 0 ]] ; then
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
				echo "@    xaxis  label "$'"RMSD (nm)"' >> RgVsRMSD.xvg
				echo "@    yaxis  label "$'"Rg (nm)"' >> RgVsRMSD.xvg
				echo "@TYPE xy" >> RgVsRMSD.xvg
				paste -d "        " RMSData.dat RgData.dat >> RgVsRMSD.xvg
				
				precalcRgRMS="./RgVsRMSD.xvg"
			elif [[ "$inFormat" == 2 ]]; then
				read -p ' Provide the path to the pre-calculated RgVsRMSD.xvg: ' precalcRgRMS
			fi

			felcal=0
			eval $gmx_exe_path sham -f $precalcRgRMS -ls FEL_sham_RgVsRMSD_$filenm.xpm -notime || felcal=1
		
			if [[ "$felcal" == 0 ]] ; then
				eval $gmx_exe_path xpm2ps -f FEL_sham_RgVsRMSD_$filenm.xpm -o FEL_sham_RgVsRMSD_$filenm.eps -rainbow red || true
	
				ps2pdf -sPAPERSIZE=ledger FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_landscape.pdf || true
		
				ps2pdf FEL_sham_RgVsRMSD_$filenm.eps FEL_sham_RgVsRMSD_${filenm}_portrait.pdf || true

				pdf2ppm -png -r 600 FEL_sham_RgVsRMSD_${filenm}_portrait.pdf FEL_sham_RgVsRMSD_${filenm}_portrait.png || true

				convert FEL_sham_RgVsRMSD_$filenm.eps -trim -bordercolor white FEL_sham_RgVsRMSD_${filenm}_convertps2png.png || true

				convert FEL_sham_RgVsRMSD_$filenm.eps -trim -bordercolor white -units pixelsperinch -density 600 -resize 3000x5000 FEL_sham_RgVsRMSD_${filenm}_convertps2pngfull.png || true
			
				currentFELshamRgVsRMSDdir="$(pwd)""/RgVsRMSD_FEL_sham"
				nFELRgVsRMSD=1
				bkupFELshamRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FEL_sham"".""backup.""$nFELRgVsRMSD"
				if [[ -d "$currentFELshamRgVsRMSDdir" ]]; then
					echo $'\n'"$currentFELshamRgVsRMSDdir"$' folder exists,\n'"backing it up as $bkupFELshamRgVsRMSDdir"
					while [[ -d "$bkupFELshamRgVsRMSDdir" ]]; do
						nFELRgVsRMSD=$(( nFELRgVsRMSD + 1 )); bkupFELshamRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FEL_sham"".""backup.""$nFELRgVsRMSD"
					done
					mv "$currentFELshamRgVsRMSDdir" "$bkupFELshamRgVsRMSDdir" && mkdir ./RgVsRMSD_FEL_sham || true
					echo $'\n'"Backing up the last RgVsRMSD_FEL_sham folder and its contents as $bkupFELshamRgVsRMSDdir"
					mv FEL_sham_RgVsRMSD_* RgData.dat RMSData.dat RgVsRMSD.xvg enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./RgVsRMSD_FEL_sham || true
				elif [[ ! -d "$currentFELshamRgVsRMSDdir" ]]; then mkdir RgVsRMSD_FEL_sham
					mv FEL_sham_RgVsRMSD_* RgData.dat RMSData.dat RgVsRMSD.xvg enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./RgVsRMSD_FEL_sham || true
				fi
				echo "$demA"$' Prepare Rg Vs RMSD FEL with gmx sham...DONE'"$demB"

			elif [[ "$felcal" == 1 ]] ; then
				echo $'Calculation failed.\n'" Please confirm that you have entered the right path/file as input!"
			fi
		fi
 	elif [[ $orderPair_choice == 2 && $flw == 1 ]] ; then useFoundRgRMSData_sham

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
		eval $gmx_exe_path sham -f $precalcOrderParPair -ls FEL_sham_OrderParameterPair_${filenm}.xpm -notime || felcal=1
		
		if [[ "$felcal" == 0 ]] ; then
			eval $gmx_exe_path xpm2ps -f FEL_sham_OrderParameterPair_$filenm.xpm -o FEL_sham_OrderParameterPair_${filenm}.eps -rainbow red || true

			ps2pdf -sPAPERSIZE=ledger FEL_sham_OrderParameterPair_$filenm.eps FEL_sham_OrderParameterPair_${filenm}_landscape.pdf || true
	
			ps2pdf FEL_sham_OrderParameterPair_$filenm.eps FEL_sham_OrderParameterPair_${filenm}_portrait.pdf || true

			pdf2ppm -png -r 600 FEL_sham_OrderParameterPair_${filenm}_portrait.pdf FEL_sham_OrderParameterPair_${filenm}_portrait.png || true

			convert FEL_sham_OrderParameterPair_$filenm.eps -trim -bordercolor white FEL_sham_OrderParameterPair_${filenm}_convertps2png.png || true

			convert FEL_sham_OrderParameterPair_$filenm.eps -trim -bordercolor white -units pixelsperinch -density 600 -resize 3000x5000 FEL_sham_OrderParameterPair_${filenm}_convertps2pngfull.png || true
			
			currentFELshamOrderParameterPairdir="$(pwd)""/OrderParameterPair_FEL_sham"
			nFELOrderParameterPair=1
			bkupFELshamOrderParameterPairdir="$(pwd)""/#OrderParameterPair_FEL_sham"".""backup.""$nFELOrderParameterPair"
			if [[ -d "$currentFELshamOrderParameterPairdir" ]]; then
				echo $'\n'"$currentFELshamOrderParameterPairdir"$' folder exists,\n'"backing it up as $bkupFELshamOrderParameterPairdir"
				while [[ -d "$bkupFELshamOrderParameterPairdir" ]]; do
					nFELOrderParameterPair=$(( nFELOrderParameterPair + 1 )); bkupFELshamOrderParameterPairdir="$(pwd)""/#OrderParameterPair_FEL_sham"".""backup.""$nFELOrderParameterPair"
				done
				mv "$currentFELshamOrderParameterPairdir" "$bkupFELshamOrderParameterPairdir" && mkdir ./OrderParameterPair_FEL_sham || true
				echo $'\n'"Backing up the last OrderParameterPair_FEL_sham folder and its contents as $bkupFELshamOrderParameterPairdir"
				mv FEL_sham_OrderParameterPair_* precalcOrderPar1.dat precalcOrderPar2.dat OrderParameterPair.xvg enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./OrderParameterPair_FEL_sham || true
			elif [[ ! -d "$currentFELshamOrderParameterPairdir" ]]; then mkdir OrderParameterPair_FEL_sham
				mv FEL_sham_OrderParameterPair_* precalcOrderPar1.dat precalcOrderPar2.dat OrderParameterPair.xvg enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./OrderParameterPair_FEL_sham || true
			fi
			echo "$demA"$' Prepare FEL with gmx sham...DONE'"$demB"
		elif [[ "$felcal" == 1 ]] ; then
			echo $'Calculation failed.\n'" Please confirm that you have entered the right path/file as input!"
		fi
	fi
	echo "$demA"$' Construct free energy landscape with gmx sham...DONE'"$demB"

}

if [[ "$analysis" == *" 11 "* ]]; then analyser11 ; fi

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

	echo $' Determining minimal and maximal data points...'
	sleep 1
	minPC1=$(head -1 sorted_PC1.dat)
	maxPC1=$(tail -1 sorted_PC1.dat)
	minPC2=$(head -1 sorted_PC2.dat)
	maxPC2=$(tail -1 sorted_PC2.dat)

	rm sorted_PC1.dat sorted_PC2.dat

	echo "minPar1,$minPC1"$'\n'"maxPar1,$maxPC1"$'\n'"minPar2,$minPC2"$'\n'"maxPar2,$maxPC2" > CHAP_fes_Par.in
	echo $'XaxisL,PC1\nYaxisL,PC2\nxbin,100\nybin,100\nTemp,300' >> CHAP_fes_Par.in
	echo $'outFilename,PCA_FES\nplotTitle,PCA-derived' >> CHAP_fes_Par.in

	echo "$demA"$' Now running construct_free_en_surface.py to construct FES...\n'
	sleep 2
	python3 ./utilities/CHAP_construct_free_en_surface.py || python3 ./utilities/CHAP_construct_free_en_surface.py

	echo $'\n Run construct_free_en_surface.py...DONE'
	sleep 2
	echo "$demA"$' Cleaning up...\n'

	currentFESchapPCAdir="$(pwd)""/PCA_FES_chap"
	nFESPCA=1
	bkupFESchapPCAdir="$(pwd)""/#PCA_FES_chap"".""backup.""$nFESPCA"
	if [[ -d "$currentFESchapPCAdir" ]]; then
		echo $'\n'"$currentFESchapPCAdir"$' folder exists,\n'"backing it up as $bkupFESchapPCAdir"
		while [[ -d "$bkupFESchapPCAdir" ]]; do
			nFESPCA=$(( nFESPCA + 1 )); bkupFESchapPCAdir="$(pwd)""/#PCA_FES_chap"".""backup.""$nFESPCA"
		done
		echo $'\n'"Backing up the last PCA_FES_chap folder and its contents as $bkupFESchapPCAdir"
		mv "$currentFESchapPCAdir" "$bkupFESchapPCAdir" && mkdir ./PCA_FES_chap || true
	elif [[ ! -d "$currentFESchapPCAdir" ]]; then mkdir PCA_FES_chap
	fi

	OrderParameter1="PC1.dat"
	OrderParameter2="PC2.dat"
	fesFigure="PCA_FES.png"
	results_folder="PCA_FES_chap"
	#fetch simulation time from the .xvg file
	cat "$exist2dPCA" | grep -v "^[@#]" | awk '{print $1}' > SimTime.dat

}

useFoundRgRMSData_FESPy()
{
	echo "$demA"$' Preparing parameters for FES calculations...\n'
	sleep 2
	
	cat RMSData.dat | sort -n > sorted_RMSData.dat
	cat RgData.dat | sort -n > sorted_RgData.dat

	paste -d "," RMSData.dat RgData.dat > OrderParameterPair.dat

	echo $' Determining minimal and maximal data points...'
	sleep 1
	minRMSD=$(head -1 sorted_RMSData.dat)
	maxRMSD=$(tail -1 sorted_RMSData.dat)
	minRg=$(head -1 sorted_RgData.dat)
	maxRg=$(tail -1 sorted_RgData.dat)

	rm sorted_RMSData.dat sorted_RgData.dat

	echo "minPar1,$minRMSD"$'\n'"maxPar1,$maxRMSD"$'\n'"minPar2,$minRg"$'\n'"maxPar2,$maxRg" > CHAP_fes_Par.in
	echo $'XaxisL,RMSD (nm)\nYaxisL,Rg (nm)\nxbin,100\nybin,100\nTemp,300' >> CHAP_fes_Par.in
	echo $'outFilename,RgVsRMSD_FES\nplotTitle,Rg Vs RMSD' >> CHAP_fes_Par.in

	echo "$demA"$' Now running construct_free_en_surface.py to construct FES...\n'
	sleep 2
	python3 ./utilities/CHAP_construct_free_en_surface.py || python3 ./utilities/CHAP_construct_free_en_surface.py

	echo $'\n Run construct_free_en_surface.py...DONE'"$demB"
	sleep 2
	
	echo "$demA"$' Cleaning up...\n'
	
	currentFESchapRgVsRMSDdir="$(pwd)""/RgVsRMSD_FES_chap"
	nFESRgVsRMSD=1
	bkupFESchapRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FES_chap"".""backup.""$nFESRgVsRMSD"
	if [[ -d "$currentFESchapRgVsRMSDdir" ]]; then
		echo $'\n'"$currentFESchapRgVsRMSDdir"$' folder exists,\n'"backing it up as $bkupFESchapRgVsRMSDdir"
		while [[ -d "$bkupFESchapRgVsRMSDdir" ]]; do
			nFESRgVsRMSD=$(( nFESRgVsRMSD + 1 )); bkupFESchapRgVsRMSDdir="$(pwd)""/#RgVsRMSD_FES_chap"".""backup.""$nFESRgVsRMSD"
		done
		echo $'\n'"Backing up the last RgVsRMSD_FES_chap folder and its contents as $bkupFESchapRgVsRMSDdir"
		mv "$currentFESchapRgVsRMSDdir" "$bkupFESchapRgVsRMSDdir" && mkdir ./RgVsRMSD_FES_chap || true
	elif [[ ! -d "$currentFESchapRgVsRMSDdir" ]]; then mkdir RgVsRMSD_FES_chap
	fi

	OrderParameter1="RMSData.dat"
	OrderParameter2="RgData.dat"
	fesFigure="RgVsRMSD_FES.png"
	results_folder="RgVsRMSD_FES_chap"
	#fetch simulation time from one of the .xvg files
	cat "$inputRg_xvgData" | grep -v "^[@#]" | awk '{print $1}' > SimTime.dat

}

analyser12()
{	
	order_parameters
 	if [[ $orderPair_choice == 1 && "$flw" == 0 ]] ; then
		if [[ "$PCFile" == 1 ]]; then useFoundPCA_FESPy
		elif [[ "$PCFile" == 2 ]]; then analyser7; useFoundPCA_FESPy	
		elif [[ "$PCFile" == 3 ]]; then echo ""
			read -p ' Provide the path to the pre-calculated 2d_PCA projection file: ' precalcPCfile
			exist2dPCA="$precalcPCfile"
			useFoundPCA_FESPy
		fi			
 	elif [[ $orderPair_choice == 1 && $flw == 1 ]] ; then useFoundPCA_FESPy
	elif [[ $orderPair_choice == 2 && $flw == 0 ]] ; then
		if [[ "$PCFile" == 1 ]]; then useFoundRgRMSData_FESPy
		elif [[ "$PCFile" == 2 ]]; then analyser2; analyser4; useFoundRgRMSData_FESPy	
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
				echo "@    xaxis  label "$'"RMSD (nm)"' >> RgVsRMSD.xvg
				echo "@    yaxis  label "$'"Rg (nm)"' >> RgVsRMSD.xvg
				echo "@TYPE xy" >> RgVsRMSD.xvg
				paste -d "        " RMSData.dat RgData.dat >> RgVsRMSD.xvg
				
				RMSData="$precalcRMSD"
				precalcRgRMS="./RgVsRMSD.xvg"
			elif [[ "$inFormat" == 2 ]]; then
				read -p ' Provide the absolute path to a pre-calculated RgVsRMSD.xvg: ' precalcRgRMS
			fi

			exist2dPCA="$precalcPCfile"
			useFoundRgRMSData_FESPy
		fi

	elif [[ $orderPair_choice == 2 && $flw == 1 ]] ; then useFoundRgRMSData_FESPy
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

		echo "minPar1,$minParam1"$'\n'"maxPar1,$maxParam1"$'\n'"minPar2,$minParam2"$'\n'"maxPar2,$maxParam2" > CHAP_fes_Par.in
		echo $'XaxisL,PC1\nYaxisL,PC2\nxbin,100\nybin,100\nTemp,300' >> CHAP_fes_Par.in
		echo $'outFilename,FES\nplotTitle,' >> CHAP_fes_Par.in

		echo "$demA"$' Now running construct_free_en_surface.py to construct FES...\n'
		sleep 2
		python3 ./utilities/CHAP_construct_free_en_surface.py || python3 ./utilities/CHAP_construct_free_en_surface.py

		echo $'\n Run construct_free_en_surface.py...DONE'
		sleep 2
		echo "$demA"$' Cleaning up...\n'

		currentFESchapdir="$(pwd)""/FES_chap"
		nFES=1
		bkupFESchapdir="$(pwd)""/#FES_chap"".""backup.""$nFES"
		if [[ -d "$currentFESchapdir" ]]; then
			echo $'\n'"$currentFESchapdir"$' folder exists,\n'"backing it up as $bkupFESchapdir"
			while [[ -d "$bkupFESchapdir" ]]; do
				nFES=$(( nFES + 1 )); bkupFESchapdir="$(pwd)""/#FES_chap"".""backup.""$nFES"
			done
			echo $'\n'"Backing up the last FES_chap folder and its contents as $bkupFESchapdir"
			mv "$currentFESchapdir" "$bkupFESchapdir" && mkdir ./FES_chap || true
			
		elif [[ ! -d "$currentFESchapdir" ]]; then mkdir FES_chap
		fi
		OrderParameter1="precalcOrderPar1.dat"
		OrderParameter2="precalcOrderPar2.dat"
		fesFigure="FES.png"
		results_folder="FES_chap"
		#fetch simulation time from one of the .xvg files
		cat "$inputprecalcOrderPar1" | grep -v "^[@#]" | awk '{print $1}' > SimTime.dat

	fi

	cat OrderParameters1_2_dG.dat | grep --color=none "\S" > OrderParameters1_2_dG_nogap.dat
	# paste SimTime.dat OrderParameters1_2_dG_nogap.dat > SimTime_OrderParameters1_2_dG.dat
	#sort based on dG values
	# cat SimTime_OrderParameters1_2_dG.dat | sort -k 4,4 -n > SimTime_OrderParameters1_2_dG-sorted.dat
	cat OrderParameters1_2_dG_nogap.dat | sort -k 3,3 -n > OrderParameters1_2_dG_nogap-sorted.dat
	mv "$fesFigure" OrderParameterPair.dat CHAP_fes_Par.in OrderParameters1_2_dG.dat ./"$results_folder" || true
	mv RgVsRMSD.xvg OrderParameters1_2_dG_nogap.dat ./"$results_folder" || true

	echo "$demA"$' Construct free energy surface...DONE'"$demB"
	
	echo "$demA"$'\n Proceed to identify, extract the lowest energy structure from the landscape?\n'
	read -p ' Enter a response here (yes or no): ' getEnergyMin
	
	while [[ "$getEnergyMin" != "yes" && "$getEnergyMin" != "no" && \
		"$getEnergyMin" != '"yes"' && "$getEnergyMin" != '"no"' ]]
	do
		echo $' \nPlease enter the appropriate response (a "yes" or a "no")!!\n'
		echo $' Identify and extract the lowest energy structure from the landscape??\n'
		read -p ' Enter yes or no here: ' getEnergyMin
	done
	if [[ "$getEnergyMin" == "yes" || "$getEnergyMin" == '"yes"' ]] ; then
		mkdir collect_mappings
		paste SimTime.dat $OrderParameter1 $OrderParameter2 > SimTime_OrderParameters1_2.dat
		echo "$demA"$' Mapping landscape data point to simulation time...\n'
		sleep 2
		python3 ./utilities/CHAP_map_fes_parameter_to_simTime.py || python ./utilities/CHAP_map_fes_parameter_to_simTime.py

		echo $' Map landscape data point to simulation time...DONE\n'
		sleep 2

		tail -n +2 ./collect_mappings/1.txt | sort -k2,2n -k3,3n > ./collect_mappings/sorted_1.txt

		head -n 2 ./collect_mappings/1.txt | tail -n 1 > ./collect_mappings/EnergyMinim.txt

		echo $' Identifying the corresponding time for the lowest energy structure...\n'
		sleep 2
		python3 ./utilities/CHAP_get_lowest_en_datapoint.py || python ./utilities/CHAP_get_lowest_en_datapoint.py
		lowEn_time=$(tail -n 1 ./collect_mappings/lowest_energy_datapoints_timed.dat | awk '{print $1}')
		lowEn_time_ps=$(awk "BEGIN {print $lowEn_time * 1000}")

		echo $' Identify the corresponding time for the lowest energy structure...DONE'
		sleep 1
		
		echo "$demA"$' Extracting lowest energy structure from the trajectory...\n\n'
		sleep 2

		if [[ $flw == 1 && $sysType == 1 ]]; then
			echo 1 | eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
		elif [[ $flw == 0 && $sysType == 1 ]]; then
			eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
		elif [[ $flw == 0 || $flw == 1 ]] && [[ $sysType == 3 ]]; then
			echo "Protein_DNA" | eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
		elif [[ $flw == 0 && $sysType == 2 ]]; then
			eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
		elif [[ $flw == 1 && $sysType == 2 ]]; then
			echo "Protein_$ligname" | eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}"_LowestEnergy_time"$lowEn_time_ps".pdb -dump "$lowEn_time_ps"
		fi
		echo "$demA"$' Extract lowest energy structure from the trajectory...DONE'"$demB"

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
		
		echo "$demA"$' Mapping all data points from the Free Enery Surface to simulation time...\n'
		sleep 2

		if [[ "$getMoreStructs" == 1 ]] ; then
			cp ./"$results_folder"/SimTime_OrderParameters1_2.dat . || true
			cp ./"$results_folder"/OrderParameters1_2_dG_nogap-sorted.dat . || true
			mkdir collect_mappings_extra
			python3 ./utilities/CHAP_map_all_dataPoint_to_simTime.py || \
			python ./utilities/CHAP_map_all_dataPoint_to_simTime.py

			echo "$demA"$' Collecting approximate simulation entries for mapped data points...\n'
			sleep 2
			echo $' Sorting and pre-processing FES data points...'
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
			python3 ./utilities/CHAP_approx_dataPoint_simTime_for_freeEn.py || \
			python ./utilities/CHAP_approx_dataPoint_simTime_for_freeEn.py

			sort ./collect_mappings_extra/mappedFESdataPoints_timed.dat -k1,1n -k4,4n > \
			./collect_mappings_extra/time-sorted_mappedFESdataPoints_timed.dat
			sort ./collect_mappings_extra/mappedFESdataPoints_timed.dat -k4,4n -k2,2n > \
			./collect_mappings_extra/energy-sorted_mappedFESdataPoints_timed.dat
			
			mv ./collect_mappings_extra/mappedFESdataPoints_timed.dat ./"$results_folder"/collect_mappings/
			mv ./collect_mappings_extra/time-sorted_mappedFESdataPoints_timed.dat ./"$results_folder"/collect_mappings/
			mv ./collect_mappings_extra/energy-sorted_mappedFESdataPoints_timed.dat ./"$results_folder"/collect_mappings/
			rm -r ./collect_mappings_extra SimTime_OrderParameters1_2.dat OrderParameters1_2_dG_nogap-sorted.dat

			echo $'\n Collect approximate simulation entries for mapped data points...DONE\n'"$demB"
			echo "$demA"$'\n **NOTE: An output file named mappedFESdataPoints_timed.dat and copies of it'\
			$'\n sorted by time and energy all containing the mapped simulation time, order'\
			$'\n parameters and the corresponding free energy has been generated and saved'\
			$'\n into the folder '"$results_folder""/collect_mappings."\
			$'\n\n **You may use this file to identify the sumulation time(s) of the'\
			$'\n structure(s) you may want to extract from the free energy landscape'"$demB"

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

				read -p ' Specify the time (in ns) at which to extract a structure: ' frameTime

				frameTime_ps=$(awk "BEGIN {print $frameTime * 1000}")

				echo "$demA"$' Extracting structure at the specified from the trajectory...\n\n'
				sleep 2

				if [[ $flw == 1 && $sysType == 1 ]]; then
					echo 1 | eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
				elif [[ $flw == 0 && $sysType == 1 ]]; then
					eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
				elif [[ $flw == 0 || $flw == 1 ]] && [[ $sysType == 3 ]]; then
					echo "Protein_DNA" | eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
				elif [[ $flw == 0 && $sysType == 2 ]]; then
					eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
				elif [[ $flw == 1 && $sysType == 2 ]]; then
					echo "Protein_$ligname" | eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}"_Structure_at_Time"$frameTime_ps".pdb -dump "$frameTime_ps"
				fi
				echo "$demA"$' Extract structure at the specified from the trajectory...DONE'"$demB"

				mv "${filenm}"_Structure_at_Time"$frameTime_ps".pdb ./"$results_folder"/collect_mappings/ || true

				echo "$demA"$'\n The structure has been saved to the folder '"$results_folder""/collect_mappings"

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
	elif [[ "$getEnergyMin" == "no" || "$getEnergyMin" == '"no"' ]] ; then echo ""
	fi
}

if [[ "$analysis" == *" 12 "* ]]; then analyser12 ; fi

analyser13()
{	
echo $'CHAPERONg: Program still under development. Check for an update later!\n   Thanks!!'	
}

if [[ "$analysis" == *" 13 "* ]]; then analyser13 ; fi


analyser14()
{	
read -p '*Please enter the number of frames to skip at intervals: ' ski
if [[ "$ski" != "0" ]]; then skp="-skip ""$ski"
fi
echo "You entered: $ski"$'\n'
echo "$demA"$' Now extracting chosen frames...\n'
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 0 | eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trjSystem_Every""$ski""skip.xtc" $skp
	echo 1 | eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trjProtein_Every""$ski""skip.xtc" $skp
elif [[ $flw == 0 ]] || [[ $flw == 1 ]] && [[ $sysType == 2 ]]; then
	eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}""_trj_Every""$ski""skip.xtc" $skp
elif [[ $flw == 0 ]] && [[ $sysType == 1 ]]; then
	eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trj_Every""$ski""skip.xtc" $skp
else
	eval $gmx_exe_path trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trj_Every""$ski""skip.xtc" $skp
fi
echo "$demA"$' Extract frames...DONE'"$demB"
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

eval $gmx_exe_path make_ndx -f em.gro -o $ndxNAME

echo "$demA"" Make index group ${nameForIndex}... DONE""$demB"
}


if [[ "$analysis" == *" 14 "* ]]; then analyser14 ; fi

if [[ "$analysis" == *" 15 "* ]]; then makeNDXGroup2 ; fi

if [[ "$analysis" == *" 16 "* ]]; then analyser0; ScanTRAJ; analyser1; analyser2; analyser3
	analyser4; analyser5; analyser6; analyser7; analyser8; analyser9; analyser10; analyser11; analyser12 
elif [[ "$analysis" == *" 17 "* ]]; then ScanTRAJ; analyser1; analyser2; analyser3
	analyser4; analyser5; analyser6; analyser7; analyser8; analyser9; analyser10; analyser11; analyser12 	
fi
}
