#! /bin/bash

#CHAP_ana -- The trajectory/end-point analysis module of CHAPERONg
#CHAPERONg -- An automation program for GROMACS md simulation
#Author -- Abeeb A. Yekeen
#Contact -- yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date -- 2022.02.11


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
  9     Make a movie of the simulation summarized into 200 frames
  10    Free energy calculations with the MMPBSA method
  11    Construct PCA-derived Gibbs free energy landscape (gmx sham)
  12    Construct free energy surface based on a pair of ordered parameters
  13    Extract frames from the trajectory
  14    Make index groups (make_ndx)
  15    All but 13 and 14
  16    All but 0, 13 and 14
  
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
	echo "$demA"$'Checking the trajectory to extract info about number of frames and simulation time'"$demB"
	gmx check -f "${filenm}"_"${wraplabel}".xtc |& tee trajectDetails.log
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
		echo "$base_currentAnadir" "exists, backing it up as $base_bkupAnadir"
		while [[ -d "$bkupAnadir" ]]; do
		nDir=$(( nDir + 1 )); bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
		done
		mv "$currentAnadir" "$bkupAnadir" && mkdir ./$AnaName
		echo "Backing up the last $AnaName folder and its contents as $base_bkupAnadir"
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
		#mv ${coordinates}_Rg_ns.png ${coordinates}_Rg_ns.xvg ${coordinates}_Rg.xvg ./$AnaName || true
	fi
}
		
DNAwrapAlt()
{
echo "$demA""CHAPERONg could not find any index file. Centering on protein instead of Protein_DNA!""$demB"
echo "Protein" 0 | gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"${wraplabel}".xtc -center -pbc mol -ur compact
}
analyser0()
{	
echo "$demA"$' Now recentering the protein and rewrapping molecules within the unit cell...\n'
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 1 0 | gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"noPBC".xtc -pbc mol -center
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	echo 1 0 | gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
				
elif [[ $flw != 1 ]] && [[ $sysType == 1 ]]; then
	echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
	gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"noPBC".xtc -pbc mol -center
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
	sleep 1
	gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -pbc nojump -center
		
elif [[ $flw == 1 ]] && [[ $sysType == 2 ]]; then
	echo 1 0 | gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
	echo "$demA"$' Now performing rotational and translational fitting...\n'
	echo 4 0 | gmx trjconv -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o "${filenm}"_fit.xtc -fit rot+trans
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	echo 1 0 | gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact

elif [[ $flw != 1 ]] && [[ $sysType == 2 ]]; then
	echo $'**Choose "Protein" (1) for centering and "System" (0) for output when prompted\n'
	gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
	echo "$demA"$' Now performing rotational and translational fitting...\n'
	echo $'**Choose "Backbone" (4) to perform lsq fitting to protein backbone, and "System" (0) for output when prompted\n'
		gmx trjconv -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o "${filenm}"_fit.xtc -fit rot+trans

elif [[ $flw == 1 ]] && [[ $sysType == 3 ]]; then
	echo "Protein_DNA" 0 | gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"center".xtc -center -pbc mol -ur compact || DNAwrapAlt
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	echo "Protein_DNA" 0 | gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact || DNAwrapAlt
elif [[ $flw != 1 ]] && [[ $sysType == 3 ]]; then
	echo $'**Choose "Protein_DNA" for centering and "System" (0) for output when prompted\n'
	gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"center".xtc -center -pbc mol -ur compact
	echo "$demA"$' Now removing possible jumps in the trajectory...\n'
	sleep 1
	gmx trjconv -s "${filenm}".tpr -f "${filenm}".xtc -n index.ndx -o "${filenm}"_"nojump".xtc -center -pbc nojump -ur compact || DNAwrapAlt
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

echo 13 13 | gmx rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_"$ligname"-rmsd.xvg -tu ns
}

analyser2()
{
echo "$demA"$' Now calculating RMSD...\n'
if [[ $sysType == 1 ]] || [[ $sysType == 3 ]] && [[ $flw == 1 ]]  ; then
	echo "Backbone" "Backbone" | gmx rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_BB-rmsd.xvg -tu ns
		
	gracebat ${coordinates}_BB-rmsd.xvg -hdevice PNG -autoscale xy -printfile ${coordinates}_BB-rmsd.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	
elif [[ $flw == 1 ]] && [[ $sysType == 2 ]]; then
	echo 4 4 | gmx rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_BB-rmsd.xvg -tu ns
		
	gracebat ${coordinates}_BB-rmsd.xvg -hdevice PNG -autoscale xy -printfile ${coordinates}_BB-rmsd.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "$demA"$'Protein RMSD calculation... DONE\n  Now calculating ligand RMSD...\n'
	sleep 2
		
	echo "$ligname" "$ligname" | gmx rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o \
	${coordinates}_"$ligname"-rmsd.xvg -tu ns || altRMSD
				
	gracebat ${coordinates}_"$ligname"-rmsd.xvg -hdevice PNG -autoscale xy -printfile \
	${coordinates}_"$ligname"-rmsd.png -fixed 7500 4000 -legend load || notifyImgFail
		
else
	gmx rms -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_rmsd.xvg -tu ns
		
	gracebat ${coordinates}_rmsd.xvg -hdevice PNG -autoscale xy -printfile ${coordinates}_rmsd.png \
	-fixed 7500 4000 -legend load || notifyImgFail
fi
echo "$demA"$' Compute RMSD... DONE'"$demB"

AnaName="RMSD"
filesuffx="rmsd"
createDIR
echo "$demA"$'Generate a finished figure of the RMSD plot... DONE'"$demB"
}

if [[ "$analysis" == *" 2 "* ]]; then analyser2 ; fi

analyser3()
{
echo "$demA"$' Now calculating RMSF...\n'
if [[ $flw == 1 ]]; then
	echo "Backbone" "Backbone" | gmx rmsf -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_BB-rmsf.xvg -res
	echo "$demA"$' RMSF with backbone lsq fitting and calculation...DONE'"$demB"
	echo "C-alpha" "C-alpha" | gmx rmsf -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_Calpha-rmsf.xvg -res
	echo "$demA"$' RMSF with Calpha lsq fitting and calculation...DONE'"$demB"
	gracebat ${coordinates}_BB-rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${coordinates}_BB-rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat ${coordinates}_Calpha-rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${coordinates}_Calpha-rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat ${coordinates}_BB-rmsf.xvg ${coordinates}_Calpha-rmsf.xvg -hdevice PNG -autoscale xy -printfile \
	${coordinates}_BB-Calpha-rmsf.png -fixed 7500 4000 -legend load || notifyImgFail
else
	gmx rmsf -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_rmsf.xvg -res
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
	echo "Protein" | gmx gyrate -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_Rg.xvg
	echo "$demA"$' Compute radius of gyration...DONE'"$demB"
	echo "$demA"$' Now converting Rg plot to ns format...\n'
	sleep 2
	grep "^[@#]" ${coordinates}_Rg.xvg | sed "s/ps/ns/g" > ${coordinates}_Rg_ns.xvg
	grep -v "^[@#]" ${coordinates}_Rg.xvg | \
	awk '{print $1/1000"      "$2"      "$3"      "$4"     "$5}' >> ${coordinates}_Rg_ns.xvg
else
	echo $'**In the following step, CHOOSE Protein (1) for Rg analysis\n\n'
	gmx gyrate -s "${filenm}".tpr -f "${filenm}"_"${wraplabel}".xtc -o ${coordinates}_Rg.xvg
	echo "$demA"$' Compute radius of gyration...DONE'"$demB"
	echo "$demA"$' Now converting Rg plot to ns format...\n'
	sleep 2
	grep "^[@#]" ${coordinates}_Rg.xvg | sed "s/ps/ns/g" > ${coordinates}_Rg_ns.xvg
	grep -v "^[@#]" ${coordinates}_Rg.xvg | \
	awk '{print $1/1000"      "$2"      "$3"      "$4"     "$5}' >> ${coordinates}_Rg_ns.xvg
fi
echo "$demA"$' Generate ns and ps Rg plots...DONE'"$demB"
sleep 2
gracebat ${coordinates}_Rg_ns.xvg -hdevice PNG -autoscale xy -printfile ${coordinates}_Rg_ns.png -fixed 7500 4000 -legend load || notifyImgFail
	
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
#mv ${coordinates}_Rg_ns.png ${coordinates}_Rg_ns.xvg ${coordinates}_Rg.xvg ./$AnaName || true
#fi
#mv *"_Rg_ns.png" *"_Rg_ns.xvg" ./$AnaName || true
echo "$demA"$'Generate a finished figure of the Rg plot... DONE'"$demB"
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
echo 1 13 | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_ProLig_${coordinates}.xvg -tu ns $hbthread
}

hbond_DNA1()
{
echo "$demA"$' Now executing Intra-protein hydrogen bonding analysis...\n'
sleep 2

echo "Protein" "Protein" | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -num hbnum_Pro_${coordinates}.xvg -tu ns $hbthread

echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Intra-DNA hydrogen bonding analysis...\n'
sleep 2

echo "DNA" "DNA" | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -num hbnum_DNA_${coordinates}.xvg -tu ns $hbthread

echo "$demA"$' Intra-DNA hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Protein-DNA hydrogen bonding analysis...\n'
echo "Protein" "DNA" | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -num hbnum_Pro_DNA_${coordinates}.xvg -tu ns $hbthread

echo "$demA"$' Protein-DNA hydrogen bonding analysis... DONE'"$demB"
sleep 2
}
hbond_DNA2()
{
echo "$demA"$' Now executing Intra-protein hydrogen bonding analysis...\n'
sleep 2

echo "Protein" "Protein" | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_Pro_${coordinates}.xvg -tu ns $hbthread
gracebat hbnum_Pro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
hbnum_Pro_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Intra-DNA hydrogen bonding analysis...\n'
sleep 2

echo "DNA" "DNA" | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_DNA_${coordinates}.xvg -tu ns $hbthread
gracebat hbnum_DNA_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
hbnum_DNA_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
echo "$demA"$' Intra-DNA hydrogen bonding analysis...DONE'"$demB"
sleep 2

echo "$demA"$' Now executing Protein-DNA hydrogen bonding analysis...\n'
echo "Protein" "DNA" | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_Pro_DNA_${coordinates}.xvg -tu ns $hbthread
gracebat hbnum_Pro_DNA_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
hbnum_Pro_DNA_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
echo "$demA"$' Protein-DNA hydrogen bonding analysis... DONE'"$demB"
sleep 2
}

analyser5()
{
echo "$demA"$' Now executing H-bond analysis...\n'
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 1 1 | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_intraPro_${coordinates}.xvg -tu ns $hbthread
	echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	echo 1 "SOL" | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num \
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
	echo 1 "$ligname" | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_ProLig_${coordinates}.xvg -tu ns $hbthread || altHBOND
	echo "$demA"$' Protein-ligand hydrogen bonding analysis...DONE'"$demB"
	sleep 2

	echo 1 1 | gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_intraPro_${coordinates}.xvg -tu ns $hbthread
	echo "$demA"$' Intra-protein hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_ProLig_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_ProLig_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
	gracebat hbnum_intraPro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_intraPro_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == 1 ]] && [[ $flw == 0 ]] ; then
	gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_${coordinates}.xvg -tu ns $hbthread
	echo "$demA"$' Hydrogen bonding analysis...DONE'"$demB"
	sleep 2
	gracebat hbnum_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	hbnum_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail	
elif [[ $sysType == 2 || "$sysType" == 3 ]] && [[ $flw == 0 ]] ; then
	gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -num hbnum_${coordinates}.xvg -tu ns $hbthread || \
	gmx hbond -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -num hbnum_${coordinates}.xvg -tu ns $hbthread
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
	echo "$currenthbonddir" "exists, backing it up as $bkuphbonddir"
	while [[ -d "$bkuphbonddir" ]]; do
	nhbond=$(( nhbond + 1 )); bkuphbonddir="$(pwd)""/#hbond"".""backup.""$nhbond"
	done
	mv "$currenthbonddir" "$bkuphbonddir" && mkdir ./hbond || true
	echo "Backing up the last hbond folder and its contents as $bkuphbonddir"
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
	echo 1 | gmx sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o sasa_${coordinates}.xvg -tu ns
	gracebat sasa_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail		
elif [[ $sysType == 2 ]] || [[ "$sysType" == 3 ]] && [[ $flw == 0 ]] ; then
	gmx sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o sasa_${coordinates}.xvg -tu ns
	gracebat sasa_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == 2 ]] && [[ $flw == 1 ]]; then
	echo 1 | gmx sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o \
	sasa_Pro_${coordinates}.xvg -tu ns
	gracebat sasa_Pro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_Pro_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ "$sysType" == 3 ]] && [[ $flw == 1 ]]; then
	echo "Protein" | gmx sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o \
	sasa_Pro_${coordinates}.xvg -tu ns
	gracebat sasa_Pro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_Pro_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail	
	echo "$demA"$'Compute solvent accessible surface area (SASA) for DNA only...DONE'"$demB"	
	echo "DNA" | gmx sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o \
	sasa_DNA_${coordinates}.xvg -tu ns
	gracebat sasa_Pro_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_DNA_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail	
	echo "$demA"$'Compute solvent accessible surface area (SASA) for DNA only...DONE'"$demB"	
	echo "$demA"$'Now calculating solvent accessible surface area (SASA) for Protein-DNA complex...\n'	
	echo "Protein_DNA" | gmx sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o \
	sasa_Pro_DNA_${coordinates}.xvg -tu ns
	gracebat sasa_Pro_DNA_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_Pro_DNA_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail	
	echo "$demA"$'Now calculating solvent accessible surface area (SASA) for Protein-DNA complex...DONE'"$demB"	
elif [[ $flw == 0 ]] && [[ $sysType == 1 ]]; then
	gmx sasa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o sasa_${coordinates}.xvg -tu ns
	gracebat sasa_${coordinates}.xvg -hdevice PNG -autoscale xy -printfile \
	sasa_${coordinates}.png -fixed 7500 4000 -legend load || notifyImgFail
fi
echo "$demA"$' Compute solvent accessible surface area (SASA)...DONE'"$demB"
currentSASAdir="$(pwd)""/SASA"
nSASA=1
bkupSASAdir="$(pwd)""/#SASA"".""backup.""$nSASA"
if [[ -d "$currentSASAdir" ]]; then
	echo "$currentSASAdir" "exists, backing it up as $bkupSASAdir"
	while [[ -d "$bkupSASAdir" ]]; do
	nSASA=$(( nSASA + 1 )); bkupSASAdir="$(pwd)""/#SASA"".""backup.""$nSASA"
	done
	mv "$currentSASAdir" "$bkupSASAdir" && mkdir ./SASA || true
	echo "Backing up the last SASA folder and its contents as $bkupSASAdir"
	mv sasa*${coordinates}.png sasa*${coordinates}.xvg ./SASA || true
elif [[ ! -d "$currentSASAdir" ]]; then
	mkdir ./SASA; mv sasa*${coordinates}.png sasa*${coordinates}.xvg ./SASA || true
fi
	echo "$demA"$'Generate a finished figure of the SASA plot... DONE'"$demB"
}

if [[ "$analysis" == *" 6 "* ]]; then analyser6 ; fi

analyser7()
{
echo "$demA"$' Now running principal component analysis (PCA)...\n'
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 4 4 | gmx covar -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr \
	-o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
	echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
	echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'
	echo 4 4 | gmx anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_"${wraplabel}".xtc -eig \
	"${filenm}"_eigenval.xvg -s "${filenm}".tpr -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg
		
	gracebat PCA_2dproj_"${filenm}".xvg -hdevice PNG -autoscale xy -printfile \
	PCA_2dproj_"${filenm}".png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $flw == 0 ]] && [[ $sysType == 1 ]]; then
	echo $'**Choose "Backbone" (4) twice when prompted\n'
	gmx covar -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
	echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
	echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'\
	$'**Choose "Backbone" (4) twice when prompted\n'
	gmx anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_"${wraplabel}".xtc -eig "${filenm}"_eigenval.xvg \
	-s "${filenm}".tpr -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg
	gracebat PCA_2dproj_"${filenm}".xvg -hdevice PNG -autoscale xy -printfile \
	PCA_2dproj_"${filenm}".png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == 2 ]] || [[ "$sysType" == 3 ]] && [[ $flw == 1 ]] ; then
	echo "Backbone" "Backbone" | gmx covar -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx \
	-o "${filenm}"_eigenval.xvg -v "${filenm}"_eigenvec.trr
	echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
	echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'
	echo "Backbone" "Backbone" | gmx anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_"${wraplabel}".xtc -eig \
	"${filenm}"_eigenval.xvg -s "${filenm}".tpr -n index.ndx -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg
	gracebat PCA_2dproj_"${filenm}".xvg -hdevice PNG -autoscale xy -printfile \
	PCA_2dproj_"${filenm}".png -fixed 7500 4000 -legend load || notifyImgFail
elif [[ $sysType == 2 ]] || [[ "$sysType" == 3 ]] && [[ $flw == 0 ]] ; then
	echo $'**Choose "Backbone" (4) twice when prompted\n'
	gmx covar -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}"_eigenval.xvg \
	-v "${filenm}"_eigenvec.trr
	echo "$demA"$' Compute and diagonalize covariance matrix...DONE'"$demB"
	echo "$demA"$' Now analyzing eigenvectors and calculating overlap between components...\n'\
	$'**Choose "Backbone" (4) twice when prompted\n'
	gmx anaeig -v "${filenm}"_eigenvec.trr -f "${filenm}"_"${wraplabel}".xtc -eig "${filenm}"_eigenval.xvg \
	-s "${filenm}".tpr -first 1 -last 2 -2d PCA_2dproj_"${filenm}".xvg	
fi
echo "$demA"$' Principal component analysis (PCA)...DONE'"$demB"
currentPCAdir="$(pwd)""/PCA"
nPCA=1
bkupPCAdir="$(pwd)""/#PCA"".""backup.""$nPCA"
if [[ -d "$currentPCAdir" ]]; then
	echo "$currentPCAdir" "exists, backing it up as $bkupPCAdir"
	while [[ -d "$bkupPCAdir" ]]; do
	nPCA=$(( nPCA + 1 )); bkupPCAdir="$(pwd)""/#PCA"".""backup.""$nPCA"
	done
	mv "$currentPCAdir" "$bkupPCAdir" && mkdir ./PCA
	echo "Backing up the last PCA folder and its contents as $bkupPCAdir"
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

echo "MainChain" | gmx do_dssp -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o ss_"${filenm}".xpm -tu ns -dt ${dt_dssp} || checkDSSP
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
	#gmx xpm2ps -f ss_"${filenm}"_HETC.xpm -o ss_"${filenm}"_colortype2.eps -rainbow blue || true
	gmx xpm2ps -f ss_"${filenm}"_HETC.xpm -o ss_"${filenm}"_colortype1.eps || true
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
		echo "$currentSecStrdir" "exists, backing it up as $bkupSecStrdir"
		while [[ -d "$bkupSecStrdir" ]]; do
		nSecStr=$(( nSecStr + 1 )); bkupSecStrdir="$(pwd)""/#Secondary_structure"".""backup.""$nSecStr"
		done
		mv "$currentSecStrdir" "$bkupSecStrdir" && mkdir ./Secondary_structure || true
		echo "Backing up the last Secondary_structure folder and its contents as $bkupSecStrdir"
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
echo "$demA"$'Now preparing to make a summary movie of the trajectory'

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
echo 0 | gmx trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -skip $skimov
#elif [[ $sysType == 1 ]] || [[ $sysType == 2 ]] && [[ $flw == 0 ]] ; then
#gmx trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -skip $skimov
#fi
sleep 2
echo "$demA"$' Preparing to extract 200 snapshots...\n'
sleep 2
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 1 | gmx trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -o summaryForMovie.pdb
elif [[ $flw == 0 ]] && [[ $sysType == 1 ]]; then
	gmx trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -o summaryForMovie.pdb
elif [[ $flw == 0 ]] || [[ $flw == 1 ]] && [[ $sysType == 3 ]]; then
	echo "Protein_DNA" | gmx trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -n index.ndx -o summaryForMovie.pdb
elif [[ $flw == 0 ]] && [[ $sysType == 2 ]]; then
	gmx trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -n index.ndx -o summaryForMovie.pdb
elif [[ $flw == 1 ]] && [[ $sysType == 2 ]]; then
	echo "Protein_$ligname" | gmx trjconv -f "${filenm}""_trjEvery""$skimov""skipForMovie.xtc" -s "${filenm}".tpr -n index.ndx -o summaryForMovie.pdb
fi
currentMOVIEdir="$(pwd)""/MOVIE"
nMOVIE=1
bkupMOVIEdir="$(pwd)""/#MOVIE"".""backup.""$nMOVIE"
if [[ -d "$currentMOVIEdir" ]]; then
	echo "$currentMOVIEdir" "exists, backing it up as $bkupMOVIEdir"
	while [[ -d "$bkupMOVIEdir" ]]; do
	nMOVIE=$(( nMOVIE + 1 )); bkupMOVIEdir="$(pwd)""/#MOVIE"".""backup.""$nMOVIE"
	done
	mv "$currentMOVIEdir" "$bkupMOVIEdir" && mkdir ./MOVIE || true
	echo "Backing up the last MOVIE folder and its contents as $bkupMOVIEdir"
elif [[ ! -d "$currentMOVIEdir" ]]; then mkdir MOVIE
fi	
	
echo $'load summaryForMovie.pdb\nsave PyMOLsession.pse\nintra_fit name ca+c+n+o\n'\
$'preset.pretty(selection='"'""all""')"\
$'\nspectrum chain, green cyan orange magenta\ncolor atomic, (not elem C)\nbg white\n'\
$'set movie_loop, 0\nsmooth\norient\nviewport 760, 540\nzoom all, -10\n'\
$'set ray_trace_frames=1\nset ray_opaque_background, 0\nset cache_frames=0\n'\
$'mclear\ncd ./MOVIE\nsave PyMOLsession_allSet.pse\nmpng frame_.png\nquit'  > prep_movie_Pyscript.pml
	
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

echo "$demA"$'Now preparing to make a summary movie from a preset PyMOL session\n'
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
	echo "$base_currentMOVIEgif" "exists, backing it up as $base_bkupMOVIEgif"$'\n'
	while [[ -f "$bkupMOVIEgif" ]]; do
	nMOVIE=$(( nMOVIE + 1 )); bkupMOVIEgif="$(pwd)""/MOVIE/dynamics_movie_""backup""$nMOVIE"".gif"
	done
	mv "$currentMOVIEgif" "$bkupMOVIEgif" || true
	echo "Backing up the last .gif MOVIE as $base_bkupMOVIEgif"
fi	

nMOVIE=1
currentMOVIEmp4="$(pwd)""/MOVIE/dynamics_movie.mp4"
bkupMOVIEmp4="$(pwd)""/MOVIE/dynamics_movie_""backup""$nMOVIE"".mp4"
if [[ -f "$currentMOVIEmp4" ]]; then
	base_currentMOVIEmp4=$(basename "$currentMOVIEmp4")
	base_bkupMOVIEmp4=$(basename "$bkupMOVIEmp4")
	echo "$base_currentMOVIEmp4" "exists, backing it up as $base_bkupMOVIEmp4"$'\n'
	while [[ -f "$bkupMOVIEmp4" ]]; do
	nMOVIE=$(( nMOVIE + 1 )); bkupMOVIEmp4="$(pwd)""/MOVIE/dynamics_movie_""backup""$nMOVIE"".mp4"
	done
	mv "$currentMOVIEmp4" "$bkupMOVIEmp4" || true
	echo "Backing up the last .mp4 MOVIE as $base_bkupMOVIEmp4"
fi	
	
echo $'cd ./MOVIE\nload PyMOLsession_allSet.pse\nmpng frame_.png\nquit'  > prep_movie_Pyscript.pml
	
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

#if [[ $sysType == 1 ]]; then indexer=''
if [[ $sysType == 2 || $sysType == 3 ]]; then indexer='-n index.ndx' ; fi
echo "$demA"$'Preparing to generate input files for g_MMPBSA free energy calculations...\n'

if [[ $trajFraction == '' ]]; then
	trajFraction=3
fi

No_of_last_third_frames=$(awk "BEGIN {print $No_of_frames / $trajFraction}")
No_of_last_third_framesINT=$(echo ${No_of_last_third_frames%\.*})

trajFractionFactor=$(( trajFraction - 1 ))

simDuratnps_lastFractn_begin=$(awk "BEGIN {print $simDuratnps * $trajFractionFactor / $trajFraction}")
simDuratnps_lastFractn_beginINT=$(echo ${simDuratnps_lastFractn_begin%\.*})


if [[ $mmpbframesNo == '' ]]; then
	mmpbframesNo=100
fi	
	
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

if [[ "$mmGMX" == "1" && "$mmGMXpath" != '' ]] ; then
	echo "$demA"$'Preparing to generate a compatible .tpr for g_MMPBSA...\n'
	eval $mmGMXpath grompp -f md.mdp -c $coordinates_raw -p topol.top -o \
	"${filenm}"_TPR_for_g_mmpbsa.tpr $indexer
	echo "$demA"$'Generate a compatible .tpr for g_MMPBSA...DONE'"$demB"
	#fi
	echo "$demA"$'Generating a compatible last 3rd of the trajectory for g_MMPBSA...\n'
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
		
		echo "$demA"$'Generate a compatible last 3rd trajectory file for g_MMPBSA...DONE'"$demB"
		
		echo "$demA"$'Extracting 100 frames from the last 3rd of the trajectory...\n'
		
		eval $mmGMXpath trjconv -s "${filenm}"_TPR_for_g_mmpbsa.tpr -f "${filenm}"_lastFractntraj4_mmpbsa.xtc \
		-o "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -skip $skipframegpsaINT
		
		echo "$demA"$'Extract 100 frames from the last 3rd of the trajectory...DONE'"$demB"
		
	fi	
	
	echo "$demA"$'Generate input files for g_MMPBSA free energy calculations...DONE'"$demB"

	echo "$demA"$'Now preparing to run g_MMPBSA calculations...\n'
	if [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw == 1 ]]; then
		echo 1 "$ligname" | ./utilities/g_mmpbsa_pkg/g_mmpbsa -f "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -s \
		"${filenm}"_TPR_for_g_mmpbsa.tpr -n index.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp
	elif [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw != 1 ]] ; then
		./utilities/g_mmpbsa_pkg/g_mmpbsa -f "${filenm}"_"$mmpbframesNo"frames_4_mmpbsa.xtc -s \
		"${filenm}"_TPR_for_g_mmpbsa.tpr -n index.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp
	fi
	echo "$demA"$'Run g_MMPBSA calculations...DONE'"$demB"
elif [[ "$mmGMX" == '' ]] ; then
	echo "$demA"$'Now preparing to run g_MMPBSA calculations...\n'
	if [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw == 1 ]]; then
		echo 1 "$ligname" | g_mmpbsa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -i pbsa.mdp \
		-pdie 2 -pbsa -decomp || echo "$demA"$'g_mmpbsa failed to run. Ensure your environment are properly set...\n'
	elif [[ $sysType == 2 || $sysType == 3 ]] && [[ $flw != 1 ]] ; then
		g_mmpbsa -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -i pbsa.mdp \
		-pdie 2 -pbsa -decomp || echo "$demA"$'g_mmpbsa failed to run. Ensure your environment are properly set...\n'
	fi
fi
echo "$demA"$'Now calculating average binding energy & contribution of residues...\n'
python ./utilities/g_mmpbsa_pkg/MmPbSaStat.py -m energy_MM.xvg -p polar.xvg -a apolar.xvg || \
python3 ./utilities/g_mmpbsa_pkg/MmPbSaStatPy3.py -m energy_MM.xvg -p polar.xvg -a apolar.xvg || true

python ./utilities/g_mmpbsa_pkg/MmPbSaDecomp.py -bs -nbs 2000 -m contrib_MM.dat -p contrib_pol.dat -a contrib_apol.dat || \
python3 ./utilities/g_mmpbsa_pkg/MmPbSaDecompPy3.py -bs -nbs 2000 -m contrib_MM.dat -p contrib_pol.dat -a contrib_apol.dat || true

if [[ "$mmGMX" == "1" && "$mmGMXpath" != '' ]] ; then
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
	echo "$base_currentAnadir" "exists, backing it up as $base_bkupAnadir"
	while [[ -d "$bkupAnadir" ]]; do
		nDir=$(( nDir + 1 )); bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
	done
	mv "$currentAnadir" "$bkupAnadir" && mkdir ./$AnaName
	echo "Backing up the last $AnaName folder and its contents as $base_bkupAnadir"
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

}

if [[ "$analysis" == *" 10 "* ]]; then ScanTRAJ; analyser10 ; fi

useFoundPCA()
{
echo "$demA"$'Now preparing Gibbs free energy landscape...\n'
sleep 1
			
gmx sham -f ./PCA/PCA_2dproj_$filenm.xvg -ls ./PCA/FEL_PCA_gibbs_$filenm.xpm -notime || true
		
gmx xpm2ps -f ./PCA/FEL_PCA_gibbs_$filenm.xpm -o ./PCA/FEL_PCA_gibbs_$filenm.eps -rainbow red || true
	
ps2pdf -sPAPERSIZE=ledger ./PCA/FEL_PCA_gibbs_$filenm.eps ./PCA/FEL_PCA_gibbs_${filenm}_landscape.pdf || true
		
ps2pdf ./PCA/FEL_PCA_gibbs_$filenm.eps ./PCA/FEL_PCA_gibbs_${filenm}_portrait.pdf || true

pdf2ppm -png -r 600 ./PCA/FEL_PCA_gibbs_${filenm}_portrait.pdf ./PCA/FEL_PCA_gibbs_${filenm}_portrait.png || true

convert ./PCA/FEL_PCA_gibbs_$filenm.eps -trim -bordercolor white ./PCA/FEL_PCA_gibbs_${filenm}_convertps2png.png || true

convert ./PCA/FEL_PCA_gibbs_$filenm.eps -trim -bordercolor white -units pixelsperinch -density 600 -resize 3000x5000  ./PCA/FEL_PCA_gibbs_${filenm}_convertps2pngfull.png || true
			
currentFELPCAdir="$(pwd)""/PCA_FreeEnergyLandscape"
nFELpca=1
bkupFELpcadir="$(pwd)""/#PCA_FreeEnergyLandscape"".""backup.""$nFELpca"
if [[ -d "$currentFELPCAdir" ]]; then
	echo "$currentFELPCAdir" "exists, backing it up as $bkupFELpcadir"
	while [[ -d "$bkupFELpcadir" ]]; do
	nFELpca=$(( nFELpca + 1 )); bkupFELpcadir="$(pwd)""/#PCA_FreeEnergyLandscape"".""backup.""$nFELpca"
	done
	mv "$currentFELPCAdir" "$bkupFELpcadir" && mkdir ./PCA_FreeEnergyLandscape || true
	echo "Backing up the last PCA_FreeEnergyLandscape folder and its contents as $bkupFELpcadir"
	mv ./PCA/FEL_PCA_gibbs_* enthalpy.xpm entropy.xpm prob.xpm shamlog.log bindex.ndx ener.xvg ./PCA_FreeEnergyLandscape || true
elif [[ ! -d "$currentFELPCAdir" ]]; then mkdir PCA_FreeEnergyLandscape
	mv ./PCA/FEL_PCA_gibbs_* ./PCA/enthalpy.xpm ./PCA/entropy.xpm ./PCA/prob.xpm ./PCA/shamlog.log ./PCA/bindex.ndx ./PCA/ener.xvg  ./PCA_FreeEnergyLandscape || true
fi
echo "$demA"$'Prepare Gibbs free energy landscape...DONE'"$demB"
}

analyser11()
{

echo "$demA"$'Checking the working directory for pre-calculated principal components...'"$demB"
sleep 2
exist2dPCA="$(pwd)""/PCA/PCA_2dproj_""$filenm"".xvg"
if [[ -f "$exist2dPCA" ]] ; then
	if [[ "$flw" == 0 ]] ; then
		echo $'CHAPERONg: Pre-calculated PCA_2d projection file found!\nFile found:'" $exist2dPCA"
		sleep 2
cat << askFELuseexist

Do you want to use this file for free energy landscape calculations?

  1) Yes, use the file above.
  2) No, repeat PCA calculations and use the new output for free energy lanscape plot.
  3) No, I want to provide another file to be used

askFELuseexist

		read -p ' Enter 1, 2 or 3 here: ' PCFile
		while [[ "$PCFile" != 1 ]] && [[ "$PCFile" != 2 ]] && [[ "$PCFile" != 3 ]]; do
			echo $'\nYou entered: '"$PCFile"
			echo $'Please enter a valid number (1, 2 or 3)!!\n'
			read -p ' Enter 1, 2 or 3 here: ' PCFile
		done

		if [[ "$PCFile" == 1 ]]; then useFoundPCA
		elif [[ "$PCFile" == 2 ]]; then analyser7; useFoundPCA	
		elif [[ "$PCFile" == 3 ]]; then echo ""
			read -p ' Please provide the path to the pre-calculated 2d_PCA projection file: ' precalcPCfile
		
			echo "$demA"$'Now preparing Gibbs free energy landscape with the user-provided 2d_PCA file...\n'
			sleep 1
			felcal=0			
			gmx sham -f $precalcPCfile -ls FEL_PCA_gibbs_$filenm.xpm -notime || felcal=1
			if [[ "$felcal" == 0 ]] ; then
				
				gmx xpm2ps -f FEL_PCA_gibbs_$filenm.xpm -o FEL_PCA_gibbs_$filenm.eps -rainbow red || true
			
				ps2pdf -sPAPERSIZE=ledger FEL_PCA_gibbs_$filenm.eps FEL_PCA_gibbs_${filenm}_landscape.pdf || true
			
				ps2pdf FEL_PCA_gibbs_$filenm.eps FEL_PCA_gibbs_${filenm}_portrait.pdf || true
			
				currentFELPCAdir="$(pwd)""/PCA_FreeEnergyLandscape"
				nFELpca=1
				bkupFELpcadir="$(pwd)""/#PCA_FreeEnergyLandscape"".""backup.""$nFELpca"
				if [[ -d "$currentFELPCAdir" ]]; then
					echo "$currentFELPCAdir" "exists, backing it up as $bkupFELpcadir"
					while [[ -d "$bkupFELpcadir" ]]; do
					nFELpca=$(( nFELpca + 1 )); bkupFELpcadir="$(pwd)""/#PCA_FreeEnergyLandscape"".""backup.""$nFELpca"
					done
					mv "$currentFELPCAdir" "$bkupFELpcadir" && mkdir ./PCA_FreeEnergyLandscape || true
					echo "Backing up the last PCA_FreeEnergyLandscape folder and its contents as $bkupFELpcadir"
					mv FEL_PCA_gibbs_*.xpm FEL_PCA_gibbs_*.eps FEL_PCA_gibbs_*.pdf ./PCA_FreeEnergyLandscape || true
				elif [[ ! -d "$currentFELPCAdir" ]]; then mkdir PCA_FreeEnergyLandscape; mv FEL_PCA_* ./PCA_FreeEnergyLandscape || true
				fi
				echo "$demA"$'Prepare Gibbs free energy landscape...DONE'"$demB"
			
			elif [[ "$felcal" == 1 ]] ; then
				echo $'Calculation failed.\n'" Please confirm that you have entered the right path/file as input!"
			fi
		fi
	
	elif [[ "$flw" == 1 ]] ; then 
	
		echo "$demA"$'Pre-calculated PCA_2d projection file found!\nFile found:'" $exist2dPCA"$'\nCHAPERONg in auto mode. Hence, '"$exist2dPCA will be used for free energy landscape creation"
		
		useFoundPCA
	fi
elif [[ ! -f "$exist2dPCA" ]] ; then analyser7; useFoundPCA		
fi 
	#echo $'CHAPERONg: Program still under development. Check for an update later!\n   Thanks!!'
}

if [[ "$analysis" == *" 11 "* ]]; then analyser11 ; fi

analyser12()
{	
echo $'CHAPERONg: Program still under development. Check for an update later!\n   Thanks!!'	
}

if [[ "$analysis" == *" 12 "* ]]; then analyser12 ; fi
	
analyser13()
{	
read -p '*Please enter the number of frames to skip at intervals: ' ski
if [[ "$ski" != "0" ]]; then skp="-skip ""$ski"
fi
echo "You entered: $ski"$'\n'
echo "$demA"$' Now extracting chosen frames...\n'
if [[ $flw == 1 ]] && [[ $sysType == 1 ]]; then
	echo 0 | gmx trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trjSystem_Every""$ski""skip.xtc" $skp
	echo 1 | gmx trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trjProtein_Every""$ski""skip.xtc" $skp
elif [[ $flw == 0 ]] || [[ $flw == 1 ]] && [[ $sysType == 2 ]]; then
	gmx trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -n index.ndx -o "${filenm}""_trj_Every""$ski""skip.xtc" $skp
elif [[ $flw == 0 ]] && [[ $sysType == 1 ]]; then
	gmx trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trj_Every""$ski""skip.xtc" $skp
else
	gmx trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${filenm}""_trj_Every""$ski""skip.xtc" $skp
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

gmx make_ndx -f em.gro -o $ndxNAME

echo "$demA"" Make index group ${nameForIndex}... DONE""$demB"
}


if [[ "$analysis" == *" 13 "* ]]; then analyser13 ; fi

if [[ "$analysis" == *" 14 "* ]]; then makeNDXGroup2 ; fi

if [[ "$analysis" == *" 15 "* ]]; then analyser0; ScanTRAJ; analyser1; analyser2; analyser3
	analyser4; analyser5; analyser6; analyser7; analyser8; analyser9; analyser10; analyser11; analyser12 
elif [[ "$analysis" == *" 16 "* ]]; then ScanTRAJ; analyser1; analyser2; analyser3
	analyser4; analyser5; analyser6; analyser7; analyser8; analyser9; analyser10; analyser11; analyser12 	
fi
}
