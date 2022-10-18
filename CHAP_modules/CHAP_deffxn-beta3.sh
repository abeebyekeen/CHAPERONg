#! /bin/bash

#CHAP_deffxn - The function definition module of CHAPERONg
#CHAPERONg - An automation program for GROMACS md simulation
#Author: Abeeb A. Yekeen
#Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date: 2022.02.11


set -e
set -o pipefail

#demA=$'\n\n'"#**********************************CHAPERONg**********************************#"$'\n'
demA=$'\n\n'"#================================= CHAPERONg =================================#"$'\n'
demB=$'\n'"#=============================================================================#"$'\n\n'
#demB=$'\n'"#*****************************************************************************#"$'\n\n'


if [[ $initiator1 != 'avail' ]] ; then
	echo "$demA"$' Do not run modules independently!\nLaunch CHAPERONg with run_CHAPERONg-<version>!!'"$demB"	
	exit 1
fi	

#import module with collected parameters
#. ./CHAP_colPar-4.00.sh


##Defining primary functions

#defining function for topol.top modification during ligand restrain step

modTop1()
{
catchLigTOPO="; Include ligand topology"
ligtopSt=0
while IFS= read -r line; do
	if [[ "$ligtopSt" == 1 ]] ; then
		ligtopSt=$(( ligtopSt + 1 ))
		lignamFile=$(echo "$line" | awk '{print $2}' | tr -d '"' | tr -d '.itp')
	elif [[ "$line" == "$catchLigTOPO" ]]; then
		ligtopSt=$(( ligtopSt + 1 ))
	fi
done < topol.top
}

modTop2()
{
catchLigITP="#include "$'"'"${lignamFile}.itp"$'"'
catchWaterTOP="; Include water topology"
lnum=0
ff_pos=0
#echo "$catchy"
while IFS= read -r line; do
	lnum=$(( lnum + 1 ))
	echo "$line" >> topol.top.temp
		if [[ "$line" == "$catchLigITP" ]]; then
			echo $'\n; Ligand position restraints\n#ifdef POSRES\n'\
$'#include '$'"'$'posre_'"$lignamFile"$'.itp'$'"\n#endif' >> topol.top.temp
		elif [[ "$ff_pos" == 1 ]]; then
			incl_ffWatMod_statement=$(echo "$line")
			ff_pos=$(( ff_pos + 1 ))
		elif [[ "$line" == "$catchWaterTOP" ]]; then
			ff_pos=$(( ff_pos + 1 ))
		fi
done < topol.top
ff_watMod_dir=$(echo "$incl_ffWatMod_statement" | awk '{print $2}')

mv topol.top topolPreResmod.top && mv topol.top.temp topol.top

}

modTop2forAcpype()
{
	catchLigITP="#include "$'"'"${lignamFile}.itp"$'"'
	catchMolType="[ moleculetype ]"
	lnum=0
	while IFS= read -r line; do
		lnum=$(( lnum + 1 ))
		echo "$line" >> topol.top.temp
			if [[ "$line" == "$catchLigITP" ]]; then
				echo $'\n; Ligand position restraints\n#ifdef POSRES\n'\
				$'#include '$'"'$'posre_'"$lignamFile"$'.itp'$'"\n#endif' >> topol.top.temp
			fi
	done < topol.top
	mv topol.top topolPreResmod.top && mv topol.top.temp topol.top
}

#defining a function for make_ndx

makeNDXGroup()
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


#defining functions for different steps of the md pipeline
s0GenTop()
{
	echo "$demA"$' Generating protein topology...'"$demB"
	sleep 2
	if [[ $coordinates_raw == '' ]]; then
		echo "$demA"" ERROR! No coordinates filename is given! Use the --input flag!!""$demB"
		Help; Credit; exit 1
	fi
	if [[ "$ffUse" == "wd" || "$ffUse" == '"wd"' ]] && [[ "$wat" != "" ]]; then
		echo "$demA"" GROMACS will use the force-field present in the working directory!""$demB"
		sleep 2
		echo 1 | eval $gmx_exe_path pdb2gmx -f $coordinates.pdb -o ${coordinates}_processed.gro $extr ${wmodel}
	elif [[ "$ffUse" != "" && "$ffUse" != "wd" && "$ffUse" != '"wd"' ]]; then
		echo "$demA"" gmx will use the specified $ffUse force-field""$demB"
		sleep 2
		eval $gmx_exe_path pdb2gmx -f $coordinates.pdb -o ${coordinates}_processed.gro $extr -ff $ffUse ${wmodel}
	elif [[ "$ffUse" == "" ]]; then
		eval $gmx_exe_path pdb2gmx -f $coordinates.pdb -o ${coordinates}_processed.gro $extr ${wmodel}
	fi
	echo "$demA"" Generate protein topology...DONE""$demB"
	sleep 2

	if [[ "$mdType" == 2 ]] ; then
		echo "$demA"$'\n Before you proceed, ensure you modify your topology to indicate the'\
		$'immobile\n reference for the pulling simulation.\n\n'\
		$'When you are done, enter "yes" below to proceed to unit cell definition...\n'

		read -p ' Do you want to proceed? (yes/no): ' procd2boxdef

		while [[ "$procd2boxdef" != "yes" && "$procd2boxdef" != "no" && \
			"$procd2boxdef" != '"yes"' && "$procd2boxdef" != '"no"' ]]
		do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p ' Do you want to proceed? (yes/no): ' procd2boxdef
		done
		if [[ "$procd2boxdef" == "yes" || "$procd2boxdef" == '"yes"' ]] ; then echo ""
		elif [[ "$procd2boxdef" == "no" || "$procd2boxdef" == '"no"' ]] ; then exit 0
		fi
	fi
}

s0CharmmSta()
{ ffcount=0
	for ffd in charmm*.ff; do
		if [[ "$ffd" == '' ]]; then
			echo "$demA"" Provide below the path/name of the charmm forcefield to be used"$'\n'
			read -p " Path (or name if in current directory) of forcefield: " forceffd
		else
			ffcount=$(( ffcount + 1))
			if [[ $ffcount == 2 ]] ; then 
				echo "$demA""Multiple forcefield folders found in the working directory! \
				Please keep only the one to be used and try again!!"
				exit 1
			fi
			forceffd="$ffd"
			echo "$demA"$'\nForcefield found in the current directory: '"$forceffd"$'\n'
			sleep 3
		fi
		break	
	done

	charm2gmxCount=0
	for cgenffscr in cgenff_charmm2gmx*.py; do
		if [[ "$cgenffscr" == '' ]]; then
			echo "$demA"" Please copy a suitable cgenff_charmm2gmx*.py script in the current \
			directory and repeat this step."$'\n'
			exit 1
		else
			charm2gmxCount=$(( charm2gmxCount + 1))
			if [[ $charm2gmxCount == 2 ]] ; then 
				echo "$demA"$'Multiple cgenff_charmm2gmx*.py scripts found in the working directory!'\
				$'Please keep only the appropriate one in the working directory.'\
				$'\nIf you are unsure, you may try each one at a time!!'
				exit 1
			fi
			cgenffscript="$cgenffscr"
			echo "charmm2gmx script found in the current directory: $cgenffscript"$'\n'
			sleep 3
		fi
		break	
	done

	ligcount=0
	for ligmol2 in *.mol2; do
		if [[ "$ligmol2" == '' ]]; then
			echo "$demA"" No mol2 ligand found in the current directory."
			exit 1
		else
			ligcount=$(( ligcount + 1))
			if [[ $ligcount == 2 ]] ; then 
				echo "$demA"$'Multiple mol2 ligands found in the working directory!'\
				$'Please keep only the required ligand and try again!!'
				exit 1
			fi
			ligandmol2="$ligmol2"
			echo "Ligand (mol2) found in the current directory: $ligandmol2"$'\n'
			sleep 3
		fi	
	done

	strcount=0
	for strLig in *.str; do
		if [[ "$strLig" == '' ]]; then
			echo "$demA"" No CHARMM stream file (.str) found in the current directory."
			exit 1
		else
		strcount=$(( strcount + 1))
			if [[ $strcount == 2 ]] ; then 
				echo "$demA"$'Multiple stream (.str) files found in the working directory!'\
				$'Please keep only the required one and try again!!'
				exit 1
			fi
		strLigand="$strLig"
		echo "Topology stream (str) file found in the current directory: $strLigand"$'\n'
		sleep 3
		fi	
	done

	echo "Ligand name (as will be written to your topol.top file): $ligname"\
	$'\n\nThe above files will be passed to charmm2gmx...'"$demB"
	sleep 2

	cgenff_charmm2gmxTry1()
	{
		python $cgenffscript $ligname $ligandmol2 $strLigand $forceffd
	}

	cgenff_charmm2gmxTry2()
	{
		echo $'\nCHAPERONg: cgenff_charmm2gmx failed. Trying again...\n'
		sleep 2
		python3 $cgenffscript $ligname $ligandmol2 $strLigand $forceffd
	}

	cgenff_charmm2gmxTry3()
	{
		echo $'\nCHAPERONg: cgenff_charmm2gmx failed. Trying again for the last time...\n'
		sleep 2
		python2 $cgenffscript $ligname $ligandmol2 $strLigand $forceffd
	}

	echo "$demA"" Will now convert CHARMM topology to GROMACS format...""$demB"

	cgenff_charmm2gmxTry1 || cgenff_charmm2gmxTry2 || cgenff_charmm2gmxTry3

	echo "$demA"" Convert CHARMM topology to GROMACS format... DONE""$demB"

}

s0CharmmStb()
{
LigINIcount=0
for LigINI in *_ini.pdb; do
	if [[ "$LigINI" == '' ]]; then
		echo "$demA"$'The current directory contains no ligand_ini.pdb file.'\
		$'\nPerhaps cgenff_charmm2gmx.py has not been run on your stream file (.str) yet.'
		exit 1
	else
	LigINIcount=$(( LigINIcount + 1))
		if [[ $LigINIcount == 2 ]] ; then 
			echo "$demA"$'Multiple ligand_ini.pdb files found in the working directory!'\
			$'\nPlease keep only the required one and try again!!'
			exit 1
		fi
	LigINIpdb="$LigINI"
	echo "$demA"" Ligand_ini.pdb file found in the current directory: $LigINIpdb""$demB"
	sleep 3
	fi	
done

ligPDB=$(basename "$LigINIpdb" _ini.pdb)

echo "$demA"" Now converting ligand_ini.pdb file to ligand.gro format...""$demB"

eval $gmx_exe_path editconf -f $LigINIpdb -o "$ligPDB"".gro"

echo "$demA"" Convert ligand_ini.pdb file to ligand.gro format...DONE"$'\n'
}

s0CharmmStc()
{
	if [[ $coordinates == *"_processed" ]]; then
		coordinates=$(basename "$coordinates" _processed) 
	fi
	LineCounterRecept1=0
	while IFS= read -r coordlineRec; do
		LineCounterRecept1=$(( LineCounterRecept1 + 1 ))
		if [[ $LineCounterRecept1 == 2 ]] ; then 
			atomCountRecept=$(echo "$coordlineRec" | awk '{print $1}')
		fi
	done < "$coordinates"_processed.gro

	echo "$demA"" Your ""$coordinates""_processed.gro file contains $atomCountRecept atoms"$'\n'

	LineCounterLig=0
	while IFS= read -r coordlineLig; do
		LineCounterLig=$(( LineCounterLig + 1 ))
		if [[ $LineCounterLig == 2 ]] ; then 
			atomCountLig=$(echo "$coordlineLig" | awk '{print $1}')
		fi
		if [[ "$LineCounterLig" != 2 && "$LineCounterLig" != 1 && "$LineCounterLig" != 0 ]]; then
			checkBox=$(echo "$coordlineLig" | awk '{print $1$2}')
			if [[ "$checkBox" == *.*.* || "$coordlineLig" == '' ]] ; then continue
			else
				echo "$coordlineLig" >> tempLigGro
			fi
		fi
	done < "$ligPDB"".gro"
	echo " Your ""$ligPDB"".gro"" file contains $atomCountLig atoms""$demB"

	TotalAtomCount=$(( atomCountRecept + atomCountLig ))

	echo "$demA"" Now preparing protein-ligand complex..."$'\n'

	LineCounterRecept2=0
	while IFS= read -r coordlineRec2; do
		LineCounterRecept2=$(( LineCounterRecept2 + 1 ))
		if [[ $LineCounterRecept2 == 1 ]] ; then 
			echo "$coordinates""-$ligname""_complex" >> "$coordinates""-$ligname"".gro"
			continue
		fi
		if [[ $LineCounterRecept2 == 2 ]] ; then 
			echo "$TotalAtomCount" >> "$coordinates""-$ligname"".gro"
			continue
		fi
		if [[ $LineCounterRecept2 != 0 && $LineCounterRecept2 != 1 && $LineCounterRecept2 != 2 ]]; then
			checkBox=$(echo "$coordlineRec2" | awk '{print $1$2}')
			if [[ "$checkBox" == *.*.* || "$coordlineRec2" == '' ]] ; then
				cat tempLigGro >> "$coordinates""-$ligname"".gro"
				echo "$coordlineRec2" >> "$coordinates""-$ligname"".gro"
				continue
			fi
			echo "$coordlineRec2" >> "$coordinates""-$ligname"".gro"
		fi
	done < "$coordinates"_processed.gro

	rm tempLigGro
	echo "$demA"" Prepare protein-ligand complex...DONE""$demB"

	coordinates="$coordinates""-$ligname"
}

s0CharmmStd()
{
	echo "$demA"" Preparing the topology for the complex..."

	Ligprmcount=0
	for ligprm in *.prm; do
		if [[ "$ligprm" == '' ]]; then
			echo "$demA""The current directory contains no ligand parameter (.prm) file. "\
			"Perhaps cgenff_charmm2gmx.py has not been run on your stream file (.str) yet."
			exit 1
		else
			Ligprmcount=$(( Ligprmcount + 1))
			if [[ $Ligprmcount == 2 ]] ; then 
				echo "$demA""Multiple ligand parameter (.prm) files found in the working "\
				"directory! Please keep only the required one and try again!!"
				exit 1
			fi
			LigprmTop="$ligprm"
			echo $'\n'" Ligand parameter file found in the current directory: $LigprmTop"$'\n\n'
			sleep 1
		fi	
	done

	for ligtop in *.itp; do
		if [[ "$ligtop" == '' ]]; then
			echo "$demA""The current directory contains no ligand topology (.itp) file. \
			Perhaps cgenff_charmm2gmx.py has not been run on your stream file (.str) yet."
			exit 1
		else
			if [[ "$ligtop" == 'posre.itp' ]] ; then continue ; fi
			Ligtopcount=$(( Ligtopcount + 1))
			if [[ $Ligtopcount == 2 ]] ; then 
				echo "$demA""Multiple ligand topology (.prm) files found in the working \
				directory! Please keep only the required one and try again!!"
				exit 1 ; fi
			LigitpTop="$ligtop"
			echo " Ligand topology file found in the current directory: $LigitpTop"
			sleep 1
		fi	
	done

	echo "$demA"" Now modifying topol.top file to include the ligand parameters..."$'\n\n'

	catchffpar="; Include forcefield parameters"
	catchPosRes="; Include Position restraint file"
	catchWaterTOP="; Include water topology"
	#initiate counters to track target lines that determine points of insertion or replacement
	lnumffpar=7; lnumposre=7; lnumMOLlist=7; lnumSysName=7; posres_pos=0; ff_pos=0
	while IFS= read -r line; do
		lnumffpar=$(( lnumffpar + 1 ))
		lnumposre=$(( lnumposre + 1 ))
		lnumMOLlist=$(( lnumMOLlist + 1 ))
		lnumSysName=$(( lnumSysName + 1 ))
		if [[ $lnumSysName == 3 ]] ; then
			echo "$coordinates""_complex" >> topol.top.temp ; continue; fi
		echo "$line" >> topol.top.temp
		if [[ $line == $catchffpar ]]; then	lnumffpar=1 ; fi
		if [[ $lnumffpar == 2 ]] ; then
			echo $'\n; Include ligand parameters\n#include "'"$LigprmTop"$'"' >> topol.top.temp ; fi
		if [[ $line == $catchPosRes ]]; then lnumposre=1 ; fi
		if [[ $lnumposre == 4 ]] ; then
			echo $'\n; Include ligand topology\n#include "'"$LigitpTop"$'"' >> topol.top.temp ; fi
		if [[ $line == "[ molecules ]" ]] ; then lnumMOLlist=1 ; fi
		if [[ $lnumMOLlist == 3 ]] ; then
			prochain="$line"
			echo "$ligname                 1" >> topol.top.temp ; fi
		if [[ "$ff_pos" == 1 ]]; then
			incl_ffWatMod_statement=$(echo "$line")
			ff_pos=$(( ff_pos + 1 )) ; fi
		if [[ $line == "$catchWaterTOP" ]]; then
			ff_pos=$(( ff_pos + 1 )) ; fi
		if [[ $line == "[ system ]" ]] ; then lnumSysName=1 ; fi

	done < topol.top

	ff_watMod_dir=$(echo "$incl_ffWatMod_statement" | awk '{print $2}')

	mv topol.top topolPreLigParmod.top && mv topol.top.temp topol.top

	echo " Modify topol.top file to include the ligand parameters...DONE"

	echo $'\n\n'" Prepare the topology for the complex...DONE""$demB"
	sleep 1
	echo "$demA"$'To prepare the topology for the complex, topol.top has been modified as follows:\n\n'\
	$'    **THE LINES\n\n'\
	$'; Include ligand parameters\n#include "'"$LigprmTop"$'"'\
	$'\n\n  **HAVE BEEN WRITTEN BETWEEN THE LINES\n\n'\
	$'; Include forcefield parameters\n#include **'".ff/forcefield.itp"$'\n\n'\
	$'    **AND THE LINES\n\n'\
	$''"[ moleculetype ]"$'\n; Name            nrexcl\n'

	sleep 3

	echo $'    **TO OBTAIN:\n\n'\
	$'; Include forcefield parameters\n#include **'".ff/forcefield.itp"$'\n'\
	$'\n; Include ligand parameters\n#include "'"$LigprmTop"$'"\n\n'\
	$''"[ moleculetype ]"$'\n; Name            nrexcl'

	sleep 2

	echo "$demA"$'\n    **ALSO, THE LINES\n\n'\
	$'; Include ligand topology\n#include "'"$LigitpTop"$'"\n'\
	$'\n  **HAVE BEEN WRITTEN BETWEEN THE LINES\n\n'\
	$'; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'\
	$'\n  **AND THE LINES\n\n'\
	$'; Include water topology\n#include '"$ff_watMod_dir"
	sleep 3

	echo $'\n    **TO OBTAIN:\n\n'\
	$'; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'\
	$'\n; Include ligand topology\n#include "'"$LigitpTop"$'"\n'\
	$'\n; Include water topology\n#include '"$ff_watMod_dir"
	sleep 2

	echo "$demA"$'\n    **ALSO, THE LINE\n\n'\
	$''"$ligname                 1"\
	$'\n\n    **HAS BEEN WRITTEN AFTER THE LINES\n\n'\
	$''"[ molecules ]"$'\n; Compound        #mols\n'"$prochain"$'\n'
	sleep 2

	echo $'    **TO OBTAIN:\n\n'\
	$''"[ molecules ]"$'\n; Compound        #mols\n'"$prochain"$'\n'\
	$''"$ligname                 1"
	sleep 2

	echo "$demA"$'\nBefore the modifications, topol.top was backed up to'\
	$'topolPreLigParmod.top.\nYou may want to check the new topol.top now to verify the'\
	$'modifications.\n\n WHEN YOU ARE DONE, you may proceed to the next step with a "yes" below...\n'

	read -p 'Do you want to proceed? (yes/no): ' procds

	while [[ "$procds" != "yes" && "$procds" != "no" && "$procds" != '"yes"' && "$procds" != '"no"' ]]
	do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p 'Do you want to proceed? (yes/no): ' procds
	done
	if [[ "$procds" == "yes" || "$procds" == '"yes"' ]] ; then echo ""
	elif [[ "$procds" == "no" || "$procds" == '"no"' ]] ; then exit 0
	fi
}


s0PrdrgSta()
{
	MolTypcatch=0; namcatch=0
	while IFS= read -r line; do
		findnam1=$(echo "$line" | awk '{print $1}')
		findnam2=$(echo "$line" | awk '{print $2}')
		findnam3=$(echo "$line" | awk '{print $3}')
		findnam="$findnam1"",""$findnam2"",""$findnam3"","
		#LineCounterRecept1=$(( LineCounterRecept1 + 1 ))
		if [[ "$line" == "[ moleculetype ]" ]] ; then MolTypcatch=1 ; continue; fi
		if [[ "$MolTypcatch" == 1 && "$findnam" == *"Name"*"nrexcl"* ]] ; then namcatch=1 ; continue ; fi
		if [[ "$namcatch" == 1 ]] ; then ligID="$findnam1" ; break ; fi
	done < ./prodrg/DRGGMX.ITP
	
	#rename PDB--residue name in the .gro file with the user-provided lig name
	#lignameDRGFIN="1""$ligname"
	cp ./prodrg/DRGFIN.GRO .
	sed -i "s|$ligID|$ligname|g" ./DRGFIN.GRO
	mv ./DRGFIN.GRO ./"$ligname"".gro"
	
	if [[ $coordinates == *"_processed" ]]; then
		coordinates=$(basename "$coordinates" _processed) 
	fi
	LineCounterRecept1=0
	while IFS= read -r coordlineRec; do
		LineCounterRecept1=$(( LineCounterRecept1 + 1 ))
		if [[ $LineCounterRecept1 == 2 ]] ; then 
			atomCountRecept=$(echo "$coordlineRec" | awk '{print $1}')
		fi
	done < "$coordinates"_processed.gro

	echo "$demA"" Your ""$coordinates""_processed.gro file contains $atomCountRecept atoms"$'\n'

	LineCounterLig=0
	while IFS= read -r coordlineLig; do
		LineCounterLig=$(( LineCounterLig + 1 ))
		if [[ $LineCounterLig == 2 ]] ; then 
			atomCountLig=$(echo "$coordlineLig" | awk '{print $1}')
		fi
		if [[ "$LineCounterLig" != 2 && "$LineCounterLig" != 1 && "$LineCounterLig" != 0 ]]; then
			checkBox=$(echo "$coordlineLig" | awk '{print $1$2}')
			if [[ "$checkBox" == *.*.* || "$coordlineLig" == '' ]] ; then continue
			else
				echo "$coordlineLig" >> tempLigGro
			fi
		fi
	done < "$ligname"".gro"

	echo " Your ""$ligname"".gro"" file contains $atomCountLig atoms""$demB"

	TotalAtomCount=$(( atomCountRecept + atomCountLig ))

	echo "$demA"" Now preparing protein-ligand complex..."$'\n'

	LineCounterRecept2=0
	while IFS= read -r coordlineRec2; do
		LineCounterRecept2=$(( LineCounterRecept2 + 1 ))
		if [[ $LineCounterRecept2 == 1 ]] ; then 
			echo "$coordinates""-$ligname""_complex" >> "$coordinates""-$ligname"".gro"
			continue
		fi
		if [[ $LineCounterRecept2 == 2 ]] ; then 
			echo "$TotalAtomCount" >> "$coordinates""-$ligname"".gro"
			continue
		fi
		if [[ $LineCounterRecept2 != 0 && $LineCounterRecept2 != 1 && $LineCounterRecept2 != 2 ]]; then
			checkBox=$(echo "$coordlineRec2" | awk '{print $1$2}')
			if [[ "$checkBox" == *.*.* || "$coordlineRec2" == '' ]] ; then
				cat tempLigGro >> "$coordinates""-$ligname"".gro"
				echo "$coordlineRec2" >> "$coordinates""-$ligname"".gro"
				continue
			fi
			echo "$coordlineRec2" >> "$coordinates""-$ligname"".gro"
		fi
	done < "$coordinates"_processed.gro

	rm tempLigGro
	echo " Prepare protein-ligand complex...DONE""$demB"

	coordinates="$coordinates""-$ligname"
}

s0PrdrgStb()
{
	MolTypcatch=0; namcatch=0
	while IFS= read -r line; do
		findnam1=$(echo "$line" | awk '{print $1}')
		findnam2=$(echo "$line" | awk '{print $2}')
		findnam3=$(echo "$line" | awk '{print $3}')
		findnam="$findnam1"",""$findnam2"",""$findnam3"","
		#LineCounterRecept1=$(( LineCounterRecept1 + 1 ))
		if [[ "$line" == "[ moleculetype ]" ]] ; then MolTypcatch=1 ; continue; fi
		if [[ "$MolTypcatch" == 1 && "$findnam" == *"Name"*"nrexcl"* ]] ; then namcatch=1 ; continue ; fi
		if [[ "$namcatch" == 1 ]] ; then ligID="$findnam1" ; break ; fi
	done < ./prodrg/DRGGMX.ITP
	
	#rename PDB--residue name in the .itp file with the user-provided lig name
	cp ./prodrg/DRGGMX.ITP .
	sed -i "s|$ligID|$ligname|g" ./DRGGMX.ITP
	mv ./DRGGMX.ITP ./"$ligname"".itp"
	
	echo "$demA"" Preparing the topology for the complex..."$'\n'

	Ligtopcount=0
	for ligtop in *.itp; do
		if [[ "$ligtop" == '' ]]; then
			echo "$demA""The current directory contains no ligand topology (.itp) file. Perhaps you have not "\
			"copied the prodrg folder into the working directory, or the folder has not been named accordingly!"
			exit 1
		else
			if [[ "$ligtop" == 'posre.itp' || "$ligtop" == "$ligitp" || "$ligtop" == 'posre'*'.itp' ]]
			then continue
			fi
			Ligtopcount=$(( Ligtopcount + 1))
			if [[ $Ligtopcount == 2 ]] ; then 
				echo "$demA"$'Multiple ligand topology (.itp) files found in the working '\
				$'directory! Please keep only the required one and try again!!'
				exit 1 
			fi
			if [[ "$ligtop" == "$ligname"".itp" ]] ; then LigitpTop="$ligtop"
				echo " Ligand topology file found in the current directory: $LigitpTop"
			fi
			sleep 1
		fi	
	done

	echo "$demA"" Now modifying topol.top file to include the ligand parameters and topology..."$'\n\n'

	#catchffpar="; Include forcefield parameters"
	catchPosRes="; Include Position restraint file"
	catchWaterTOP="; Include water topology"
	#initiate counters to track target lines that determine points of insertion or replacement
	#lnumffpar=7; 
	lnumposre=7; lnumMOLlist=7; lnumSysName=7; posres_pos=0; ff_pos=0
	while IFS= read -r line; do
		lnumffpar=$(( lnumffpar + 1 ))
		lnumposre=$(( lnumposre + 1 ))
		lnumMOLlist=$(( lnumMOLlist + 1 ))
		lnumSysName=$(( lnumSysName + 1 ))
		if [[ $lnumSysName == 3 ]] ; then
			echo "$coordinates""_complex" >> topol.top.temp ; continue; fi
		echo "$line" >> topol.top.temp
		#if [[ $line == $catchffpar ]]; then	lnumffpar=1 ; fi
		#if [[ $lnumffpar == 2 ]] ; then
		#	echo $'\n; Include ligand parameters\n#include "'"$LigprmTop"$'"' >> topol.top.temp ; fi
		if [[ $line == $catchPosRes ]]; then lnumposre=1 ; fi
		if [[ $lnumposre == 4 ]] ; then
			echo $'\n; Include ligand topology\n#include "'"$LigitpTop"$'"' >> topol.top.temp ; fi
		if [[ $line == "[ molecules ]" ]] ; then lnumMOLlist=1 ; fi
		if [[ $lnumMOLlist == 3 ]] ; then
			prochain="$line"
			echo "$ligname                 1" >> topol.top.temp ; fi
		if [[ "$ff_pos" == 1 ]]; then
			incl_ffWatMod_statement=$(echo "$line")
			ff_pos=$(( ff_pos + 1 )) ; fi
		if [[ $line == "$catchWaterTOP" ]]; then
			ff_pos=$(( ff_pos + 1 )) ; fi
		if [[ $line == "[ system ]" ]] ; then lnumSysName=1 ; fi

	done < topol.top

	ff_watMod_dir=$(echo "$incl_ffWatMod_statement" | awk '{print $2}')

	mv topol.top topolPreLigParmod.top && mv topol.top.temp topol.top

	echo " Modify topol.top file to include the ligand parameters and topology...DONE"

	echo $'\n\n'" Prepare the topology for the complex...DONE"
	sleep 1
	echo "$demA"$'To prepare the topology for the complex, '\
	$'topol.top has been modified as follows:\n\n'
	echo "$demA"$'\n    **THE LINES\n\n'\
	$'; Include ligand topology\n#include "'"$LigitpTop"$'"\n'\
	$'\n  **HAVE BEEN WRITTEN BETWEEN THE LINES\n\n'\
	$'; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'\
	$'\n  **AND THE LINES\n\n'\
	$'; Include water topology\n#include '"$ff_watMod_dir"
	sleep 3

	echo $'\n    **TO OBTAIN:\n\n'\
	$'; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'\
	$'\n; Include ligand topology\n#include "'"$LigitpTop"$'"\n'\
	$'\n; Include water topology\n#include '"$ff_watMod_dir"
	sleep 2

	echo "$demA"$'\n    **ALSO, THE LINE\n\n'\
	$''"$ligname                 1"\
	$'\n\n    **HAS BEEN WRITTEN AFTER THE LINES\n\n'\
	$''"[ molecules ]"$'\n; Compound        #mols\n'"$prochain"$'\n'
	sleep 2

	echo $'    **TO OBTAIN:\n\n'\
	$''"[ molecules ]"$'\n; Compound        #mols\n'"$prochain"$'\n'\
	$''"$ligname                 1"
	sleep 2

	echo "$demA"$'\nThe topol.top file before the modifications has been backed up to '\
	$'topolPreLigParmod.top.\nYou may want to check the new topol.top now to verify the '\
	$'modifications.\nWHEN YOU ARE DONE, you may proceed to the next step with a "yes" below...\n\n'

	read -p 'Do you want to proceed? (yes/no): ' procds

	while [[ "$procds" != "yes" && "$procds" != "no" && "$procds" != '"yes"' && "$procds" != '"no"' ]]
	do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p 'Do you want to proceed? (yes/no): ' procds
	done
	if [[ "$procds" == "yes" || "$procds" == '"yes"' ]] ; then echo ""
	elif [[ "$procds" == "no" || "$procds" == '"no"' ]] ; then exit 0
	fi

}

#Process topology for Acpype
s0AcpypeSta()
{
	#rename PDB--residue name in the .gro file with the user-provided lig name
	#lignameDRGFIN=" ""$ligname"
	for liggro in ./acpype/* ; do
		if [[ $liggro == *"_GMX.gro" ]] ; then cp $liggro ./"$ligname"".gro" ; break ; fi
	done

	#mv $liggro ./"$ligname"".gro"
	sed -i "s|UNL|$ligname|g" ./"$ligname"".gro"
	
	if [[ $coordinates == *"_processed" ]]; then
		coordinates=$(basename "$coordinates" _processed) 
	fi
	LineCounterRecept1=0
	while IFS= read -r coordlineRec; do
		LineCounterRecept1=$(( LineCounterRecept1 + 1 ))
		if [[ $LineCounterRecept1 == 2 ]] ; then 
			atomCountRecept=$(echo "$coordlineRec" | awk '{print $1}')
		fi
	done < "$coordinates"_processed.gro

	echo "$demA"" Your ""$coordinates""_processed.gro file contains $atomCountRecept atoms"$'\n'

	LineCounterLig=0
	while IFS= read -r coordlineLig; do
		LineCounterLig=$(( LineCounterLig + 1 ))
		if [[ $LineCounterLig == 2 ]] ; then 
			atomCountLig=$(echo "$coordlineLig" | awk '{print $1}')
		fi
		if [[ "$LineCounterLig" != 2 && "$LineCounterLig" != 1 && "$LineCounterLig" != 0 ]]; then
			checkBox=$(echo "$coordlineLig" | awk '{print $1$2}')
			if [[ "$checkBox" == *.*.* || "$coordlineLig" == '' ]] ; then continue
			else
				echo "$coordlineLig" >> tempLigGro
			fi
		fi
	done < "$ligname"".gro"

	echo " Your ""$ligname"".gro"" file contains $atomCountLig atoms""$demB"

	TotalAtomCount=$(( atomCountRecept + atomCountLig ))

	echo "$demA"" Now preparing protein-ligand complex..."$'\n'

	LineCounterRecept2=0
	while IFS= read -r coordlineRec2; do
		LineCounterRecept2=$(( LineCounterRecept2 + 1 ))
		if [[ $LineCounterRecept2 == 1 ]] ; then 
			echo "$coordinates""-$ligname""_complex" >> "$coordinates""-$ligname"".gro"
			continue
		fi
		if [[ $LineCounterRecept2 == 2 ]] ; then 
			echo "$TotalAtomCount" >> "$coordinates""-$ligname"".gro"
			continue
		fi
		if [[ $LineCounterRecept2 != 0 && $LineCounterRecept2 != 1 && $LineCounterRecept2 != 2 ]]; then
			checkBox=$(echo "$coordlineRec2" | awk '{print $1$2}')
			if [[ "$checkBox" == *.*.* || "$coordlineRec2" == '' ]] ; then
				cat tempLigGro >> "$coordinates""-$ligname"".gro"
				echo "$coordlineRec2" >> "$coordinates""-$ligname"".gro"
				continue
			fi
			echo "$coordlineRec2" >> "$coordinates""-$ligname"".gro"
		fi
	done < "$coordinates"_processed.gro

	rm tempLigGro
	echo " Prepare protein-ligand complex...DONE""$demB"

	coordinates="$coordinates""-$ligname"
}

#Prepare topology for acpype protein-lig complex
s0AcpypeStb()
{
	#rename PDB--residue name in the .itp file with the user-provided lig name
	#lignameDRGFIN=" ""$ligname"
	for ligitp in ./acpype/* ; do
		if [[ $ligitp == *"_GMX.itp" ]] ; then
		#cp $ligitp .
		cp $ligitp ./"$ligname"".itp.temp" ; break ; fi
	done
	
	#cp ./prodrg/DRGGMX.ITP .
	sed -i "s|UNL|$ligname|g" ./"$ligname"".itp.temp"
	
	echo "$demA"" Extracting ligand parameters from the Acpype-derived topology..."$'\n'
	#initiate counters to track target lines that determine points of insertion or replacement
	atmTypeBeg=0; atmTypeEnd=0
	while IFS= read -r line; do
		checkatmField=$(echo "$line" | awk '{print $1$2}')
		if [[ $line == *"created by"* ]] ; then
			echo "$line" >> "$ligname"".prm"
			echo "; Parameter extracted by CHAPERONg" >> "$ligname"".prm"
		fi	
		if [[ $line == "[ atomtypes ]" ]] ; then atmTypeBeg=1; fi
		if [[ $checkatmField == "" || $checkatmField == " " ]] && [[ $atmTypeEnd == 0 && $atmTypeBeg == 1 ]] ;
			then
			echo "" >> "$ligname"".prm"
			echo "$line" >> "$ligname"".prm"
			atmTypeBeg=2 ; atmTypeEnd=1
		fi
		if [[ $atmTypeBeg == 1 ]] ; then echo "$line" >> "$ligname"".prm" ; fi

		if [[ $atmTypeBeg == 0 || $atmTypeBeg == 2 ]] ; then echo "$line" >> "$ligname"".itp" ;	fi

	done < "$ligname"".itp.temp"
	rm "$ligname"".itp.temp"
	echo " Extract ligand parameters from the Acpype-derived topology...DONE""$demA"

	echo "$demA"" Preparing the topology for the complex..."$'\n'

	Ligprmcount=0
	for ligprm in *.prm; do
		if [[ "$ligprm" == '' ]]; then
			echo "$demA""The current directory contains no ligand parameter (.prm) file. "\
			"Perhaps CHAPERONg could not process your Acpype-derived topology ($ligname"_GMX.itp")."
			exit 1
		else
			Ligprmcount=$(( Ligprmcount + 1))
			if [[ $Ligprmcount == 2 ]] ; then 
				echo "$demA""Multiple ligand parameter (.prm) files found in the working "\
				"directory! Please keep only the required one and try again!!"
				exit 1
			fi
			LigprmTop="$ligprm"
			echo $'\n'" Ligand parameter file found in the current directory: $LigprmTop"$'\n\n'
			sleep 1
		fi	
	done

	Ligtopcount=0
	for ligtop in *.itp; do
		if [[ "$ligtop" == '' ]]; then
			echo "$demA""The current directory contains no ligand topology (.itp) file. Perhaps you have not "\
			"copied the prodrg folder into the working directory, or the folder has not been named accordingly!"
			exit 1
		else
			if [[ "$ligtop" == 'posre.itp' || "$ligtop" == "$ligitp" || "$ligtop" == 'posre'*'.itp' ]]
			then continue
			fi
			Ligtopcount=$(( Ligtopcount + 1))
			if [[ $Ligtopcount == 2 ]] ; then 
				echo "$demA"$'Multiple ligand topology (.itp) files found in the working '\
				$'directory! Please keep only the required one and try again!!'
				exit 1 
			fi
			if [[ "$ligtop" == "$ligname"".itp" ]] ; then LigitpTop="$ligtop"
				echo " Ligand topology file found in the current directory: $LigitpTop"
			fi
			sleep 1
		fi	
	done

	echo "$demA"" Now modifying topol.top file to include the ligand parameters and topology..."$'\n\n'

	catchffpar="; Include forcefield parameters"
	catchPosRes="; Include Position restraint file"
	catchWaterTOP="; Include water topology"
	#initiate counters to track target lines that determine points of insertion or replacement
	lnumffpar=7; 
	lnumposre=7; lnumMOLlist=7; lnumSysName=7; posres_pos=0; ff_pos=0
	while IFS= read -r line; do
		lnumffpar=$(( lnumffpar + 1 ))
		lnumposre=$(( lnumposre + 1 ))
		lnumMOLlist=$(( lnumMOLlist + 1 ))
		lnumSysName=$(( lnumSysName + 1 ))
		if [[ $lnumSysName == 3 ]] ; then
			echo "$coordinates""_complex" >> topol.top.temp ; continue; fi
		echo "$line" >> topol.top.temp
		if [[ $line == $catchffpar ]]; then	lnumffpar=1 ; fi
		if [[ $lnumffpar == 2 ]] ; then
			echo $'\n; Include ligand parameters\n#include "'"$LigprmTop"$'"' >> topol.top.temp ; fi
		if [[ $line == $catchPosRes ]]; then lnumposre=1 ; fi
		if [[ $lnumposre == 4 ]] ; then
			echo $'\n; Include ligand topology\n#include "'"$LigitpTop"$'"' >> topol.top.temp ; fi
		if [[ $line == "[ molecules ]" ]] ; then lnumMOLlist=1 ; fi
		if [[ $lnumMOLlist == 3 ]] ; then
			prochain="$line"
			echo "$ligname                 1" >> topol.top.temp ; fi
		if [[ "$ff_pos" == 1 ]]; then
			incl_ffWatMod_statement=$(echo "$line")
			ff_pos=$(( ff_pos + 1 )) ; fi
		if [[ $line == "$catchWaterTOP" ]]; then
			ff_pos=$(( ff_pos + 1 )) ; fi
		if [[ $line == "[ system ]" ]] ; then lnumSysName=1 ; fi

	done < topol.top

	ff_watMod_dir=$(echo "$incl_ffWatMod_statement" | awk '{print $2}')

	mv topol.top topolPreLigParmod.top && mv topol.top.temp topol.top

	echo " Modify topol.top file to include the ligand parameters and topology...DONE"

	echo $'\n\n'" Prepare the topology for the complex...DONE"
	sleep 1
	echo "$demA"$'To prepare the topology for the complex, topol.top has been modified as follows:\n\n'
	echo "$demA"$'\n    **THE LINES\n\n'\
	$'; Include ligand parameters\n#include "'"$LigprmTop"$'"\n'\
	$'\n  **HAVE BEEN WRITTEN BETWEEN THE LINES\n\n'\
	$'; Include forcefield parameters\n#include **'".ff/forcefield.itp"$'\n\n'\
	$'    **AND THE LINES\n\n'\
	$''"[ moleculetype ]"$'\n; Name            nrexcl\n'
	sleep 3

	echo $'    **TO OBTAIN:\n\n'\
	$'; Include forcefield parameters\n#include **'".ff/forcefield.itp"$'\n'\
	$'\n; Include ligand parameters\n#include "'"$LigprmTop"$'"\n'\
	$''"[ moleculetype ]"$'\n; Name            nrexcl'
	sleep 2

	echo "$demA"$'\n    **ALSO, THE LINES\n\n'\
	$'; Include ligand topology\n#include "'"$LigitpTop"$'"\n'\
	$'\n  **HAVE BEEN WRITTEN BETWEEN THE LINES\n\n'\
	$'; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'\
	$'\n  **AND THE LINES\n\n'\
	$'; Include water topology\n#include '"$ff_watMod_dir"
	sleep 3

	echo $'\n    **TO OBTAIN:\n\n'\
	$'; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'\
	$'\n; Include ligand topology\n#include "'"$LigitpTop"$'"\n'\
	$'\n; Include water topology\n#include '"$ff_watMod_dir"
	sleep 2

	echo "$demA"$'\n    **ALSO, THE LINE\n\n'"$ligname                 1"\
	$'\n\n    **HAS BEEN WRITTEN AFTER THE LINES\n\n'"[ molecules ]"\
	$'\n; Compound        #mols\n'"$prochain"$'\n'
	sleep 2

	echo $'    **TO OBTAIN:\n\n'"[ molecules ]"\
	$'\n; Compound        #mols\n'"$prochain"\
	$'\n'"$ligname                 1"
	sleep 2

	echo "$demA"$'\nThe topol.top file before the modifications has been backed up to '\
	$'topolPreLigParmod.top.\nYou may want to check the new topol.top now to verify the '\
	$'modifications.\nWHEN YOU ARE DONE, you may proceed to the next step with a "yes" below...\n\n'

	read -p 'Do you want to proceed? (yes/no): ' procds

	while [[ "$procds" != "yes" && "$procds" != "no" && "$procds" != '"yes"' && "$procds" != '"no"' ]]
	do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p 'Do you want to proceed? (yes/no): ' procds
	done
	if [[ "$procds" == "yes" || "$procds" == '"yes"' ]] ; then echo ""
	elif [[ "$procds" == "no" || "$procds" == '"no"' ]] ; then exit 0
	fi

}

#Process topology for LigParGen
s0ligpargenSta()
{
	#find topology file
	for ligLPGitp in ./ligpargen/tmp/* ; do
		if [[ $ligLPGitp == *".itp" ]] ; then LPGitp="$ligLPGitp" ; break ; fi
	done
	#find ligand ID in topology file
	MolTypcatch=0; namcatch=0
	while IFS= read -r line; do
		findnam1=$(echo "$line" | awk '{print $1}')
		findnam2=$(echo "$line" | awk '{print $2}')
		findnam3=$(echo "$line" | awk '{print $3}')
		findnam="$findnam1"",""$findnam2"",""$findnam3"","
		#LineCounterRecept1=$(( LineCounterRecept1 + 1 ))
		if [[ "$line" == "[ moleculetype ]" ]] ; then MolTypcatch=1 ; continue; fi
		if [[ "$MolTypcatch" == 1 && "$findnam" == *"Name"*"nrexcl"* ]] ; then namcatch=1 ; continue ; fi
		if [[ "$namcatch" == 1 ]] ; then ligID="$findnam1" ; break ; fi
	done < "$LPGitp"
	
	#rename PDB--residue name in the .gro file with the user-provided lig name
	#lignameDRGFIN=" ""$ligname"
	for liggro in ./ligpargen/tmp/* ; do
		if [[ $liggro == *".gro" ]] ; then cp $liggro ./"$ligname"".gro" ; break ; fi
	done

	#mv $liggro ./"$ligname"".gro"
	sed -i "s|$ligID|$ligname|g" ./"$ligname"".gro"
	
	if [[ $coordinates == *"_processed" ]]; then
		coordinates=$(basename "$coordinates" _processed) 
	fi
	LineCounterRecept1=0
	while IFS= read -r coordlineRec; do
		LineCounterRecept1=$(( LineCounterRecept1 + 1 ))
		if [[ $LineCounterRecept1 == 2 ]] ; then 
			atomCountRecept=$(echo "$coordlineRec" | awk '{print $1}')
		fi
	done < "$coordinates"_processed.gro

	echo "$demA"" Your ""$coordinates""_processed.gro file contains $atomCountRecept atoms"$'\n'

	LineCounterLig=0
	while IFS= read -r coordlineLig; do
		LineCounterLig=$(( LineCounterLig + 1 ))
		if [[ $LineCounterLig == 2 ]] ; then 
			atomCountLig=$(echo "$coordlineLig" | awk '{print $1}')
		fi
		if [[ "$LineCounterLig" != 2 && "$LineCounterLig" != 1 && "$LineCounterLig" != 0 ]]; then
			checkBox=$(echo "$coordlineLig" | awk '{print $1$2}')
			if [[ "$checkBox" == *.*.* || "$coordlineLig" == '' ]] ; then continue
			else
				echo "$coordlineLig" >> tempLigGro
			fi
		fi
	done < "$ligname"".gro"

	echo " Your ""$ligname"".gro"" file contains $atomCountLig atoms""$demB"

	TotalAtomCount=$(( atomCountRecept + atomCountLig ))

	echo "$demA"" Now preparing protein-ligand complex..."$'\n'

	LineCounterRecept2=0
	while IFS= read -r coordlineRec2; do
		LineCounterRecept2=$(( LineCounterRecept2 + 1 ))
		if [[ $LineCounterRecept2 == 1 ]] ; then 
			echo "$coordinates""-$ligname""_complex" >> "$coordinates""-$ligname"".gro"
			continue
		fi
		if [[ $LineCounterRecept2 == 2 ]] ; then 
			echo "$TotalAtomCount" >> "$coordinates""-$ligname"".gro"
			continue
		fi
		if [[ $LineCounterRecept2 != 0 && $LineCounterRecept2 != 1 && $LineCounterRecept2 != 2 ]]; then
			checkBox=$(echo "$coordlineRec2" | awk '{print $1$2}')
			if [[ "$checkBox" == *.*.* || "$coordlineRec2" == '' ]] ; then
				cat tempLigGro >> "$coordinates""-$ligname"".gro"
				echo "$coordlineRec2" >> "$coordinates""-$ligname"".gro"
				continue
			fi
			echo "$coordlineRec2" >> "$coordinates""-$ligname"".gro"
		fi
	done < "$coordinates"_processed.gro

	rm tempLigGro
	echo " Prepare protein-ligand complex...DONE""$demB"

	coordinates="$coordinates""-$ligname"
}

#Prepare topology for LigParGen protein-lig complex
s0ligpargenStb()
{
	#find topology file
	for ligitp in ./ligpargen/tmp/* ; do
		if [[ $ligitp == *".itp" ]] ; then
			LPGitp="$ligitp" 
			#cp $ligitp .
			cp $ligitp ./"$ligname"".itp.temp"
			break
		fi
	done
	#detect ligand ID in topology file and store into the variable ligID
	MolTypcatch=0; namcatch=0
	while IFS= read -r line; do
		findnam1=$(echo "$line" | awk '{print $1}')
		findnam2=$(echo "$line" | awk '{print $2}')
		findnam3=$(echo "$line" | awk '{print $3}')
		findnam="$findnam1"",""$findnam2"",""$findnam3"","
		#LineCounterRecept1=$(( LineCounterRecept1 + 1 ))
		if [[ "$line" == "[ moleculetype ]" ]] ; then MolTypcatch=1 ; continue; fi
		if [[ "$MolTypcatch" == 1 && "$findnam" == *"Name"*"nrexcl"* ]] ; then namcatch=1 ; continue ; fi
		if [[ "$namcatch" == 1 ]] ; then ligID="$findnam1" ; break ; fi
	done < "$LPGitp"
	
	#rename PDB--residue name in the .itp file with the user-provided lig name
	sed -i "s|$ligID|$ligname|g" ./"$ligname"".itp.temp"
	
	echo "$demA"" Extracting ligand parameters from the LigParGen-derived topology..."$'\n'
	#initiate counters to track target lines that determine points of insertion or replacement
	atmTypeBeg=0; atmTypeEnd=0
	while IFS= read -r line; do
		if [[ $line == *"GENERATED BY"* ]] ; then
			echo "$line" >> "$ligname"".prm"
			echo "; Parameters extracted by CHAPERONg" >> "$ligname"".prm"
		fi	
		if [[ $line == "[ atomtypes ]" ]] ; then atmTypeBeg=1; fi
		if [[ $line == "[ moleculetype ]" && $atmTypeEnd == 0 && $atmTypeBeg == 1 ]] ; then
			echo "" >> "$ligname"".prm" ; atmTypeBeg=2 ; atmTypeEnd=1
		fi
		if [[ $atmTypeBeg == 1 ]] ; then echo "$line" >> "$ligname"".prm" ; fi

		if [[ $atmTypeBeg == 0 || $atmTypeBeg == 2 ]] ; then echo "$line" >> "$ligname"".itp" ;	fi

	done < "$ligname"".itp.temp"
	rm "$ligname"".itp.temp"
	echo " Extract ligand parameters from the LigParGen-derived topology...DONE""$demA"

	echo "$demA"" Preparing the topology for the complex..."$'\n'

	Ligprmcount=0
	for ligprm in *.prm; do
		if [[ "$ligprm" == '' ]]; then
			echo "$demA""The current directory contains no ligand parameter (.prm) file. "\
			"Perhaps CHAPERONg could not process your Acpype-derived topology ($ligname"_GMX.itp")."
			exit 1
		else
			Ligprmcount=$(( Ligprmcount + 1))
			if [[ $Ligprmcount == 2 ]] ; then 
				echo "$demA""Multiple ligand parameter (.prm) files found in the working "\
				"directory! Please keep only the required one and try again!!"
				exit 1
			fi
			LigprmTop="$ligprm"
			echo $'\n'" Ligand parameter file found in the current directory: $LigprmTop"$'\n\n'
			sleep 1
		fi	
	done

	Ligtopcount=0
	for ligtop in *.itp; do
		if [[ "$ligtop" == '' ]]; then
			echo "$demA""The current directory contains no ligand topology (.itp) file. Perhaps you have not "\
			"copied the prodrg folder into the working directory, or the folder has not been named accordingly!"
			exit 1
		else
			if [[ "$ligtop" == 'posre.itp' || "$ligtop" == "$ligitp" || "$ligtop" == 'posre'*'.itp' ]]
			then continue
			fi
			Ligtopcount=$(( Ligtopcount + 1))
			if [[ $Ligtopcount == 2 ]] ; then 
				echo "$demA"'Multiple ligand topology (.itp) files found in the working '\
				$'directory! Please keep only the required one and try again!!'
				exit 1 ; fi
			if [[ "$ligtop" == "$ligname"".itp" ]] ; then LigitpTop="$ligtop"
				echo " Ligand topology file found in the current directory: $LigitpTop"
			fi
			sleep 1
		fi	
	done

	echo "$demA"" Now modifying topol.top file to include the ligand parameters and topology..."$'\n\n'

	catchffpar="; Include forcefield parameters"
	catchPosRes="; Include Position restraint file"
	catchWaterTOP="; Include water topology"
	#initiate counters to track target lines that determine points of insertion or replacement
	lnumffpar=7; 
	lnumposre=7; lnumMOLlist=7; lnumSysName=7; posres_pos=0; ff_pos=0
	while IFS= read -r line; do
		lnumffpar=$(( lnumffpar + 1 ))
		lnumposre=$(( lnumposre + 1 ))
		lnumMOLlist=$(( lnumMOLlist + 1 ))
		lnumSysName=$(( lnumSysName + 1 ))
		if [[ $lnumSysName == 3 ]] ; then
			echo "$coordinates""_complex" >> topol.top.temp ; continue; fi
		echo "$line" >> topol.top.temp
		if [[ $line == $catchffpar ]]; then	lnumffpar=1 ; fi
		if [[ $lnumffpar == 2 ]] ; then
			echo $'\n; Include ligand parameters\n#include "'"$LigprmTop"$'"' >> topol.top.temp ; fi
		if [[ $line == $catchPosRes ]]; then lnumposre=1 ; fi
		if [[ $lnumposre == 4 ]] ; then
			echo $'\n; Include ligand topology\n#include "'"$LigitpTop"$'"' >> topol.top.temp ; fi
		if [[ $line == "[ molecules ]" ]] ; then lnumMOLlist=1 ; fi
		if [[ $lnumMOLlist == 3 ]] ; then
			prochain="$line"
			echo "$ligname                 1" >> topol.top.temp ; fi
		if [[ "$ff_pos" == 1 ]]; then
			incl_ffWatMod_statement=$(echo "$line")
			ff_pos=$(( ff_pos + 1 )) ; fi
		if [[ $line == "$catchWaterTOP" ]]; then
			ff_pos=$(( ff_pos + 1 )) ; fi
		if [[ $line == "[ system ]" ]] ; then lnumSysName=1 ; fi

	done < topol.top

	ff_watMod_dir=$(echo "$incl_ffWatMod_statement" | awk '{print $2}')

	mv topol.top topolPreLigParmod.top && mv topol.top.temp topol.top

	echo " Modify topol.top file to include the ligand parameters and topology...DONE"

	echo $'\n\n'" Prepare the topology for the complex...DONE"
	sleep 1
	echo "$demA"$'To prepare the topology for the complex, topol.top has been modified as follows:\n\n'
	echo "$demA"$'\n    **THE LINES\n\n'\
	$'; Include ligand parameters\n#include "'"$LigprmTop"$'"\n'\
	$'\n  **HAVE BEEN WRITTEN BETWEEN THE LINES\n\n'\
	$'; Include forcefield parameters\n#include **'".ff/forcefield.itp"$'\n\n'\
	$'    **AND THE LINES\n\n'\
	$''"[ moleculetype ]"$'\n; Name            nrexcl\n'
	sleep 3

	echo $'    **TO OBTAIN:\n\n'\
	$'; Include forcefield parameters\n#include **'".ff/forcefield.itp"$'\n'\
	$'\n; Include ligand parameters\n#include "'"$LigprmTop"$'"\n'\
	$''"[ moleculetype ]"$'\n; Name            nrexcl'
	sleep 2

	echo "$demA"$'\n    **ALSO, THE LINES\n\n'\
	$'; Include ligand topology\n#include "'"$LigitpTop"$'"\n'\
	$'\n  **HAVE BEEN WRITTEN BETWEEN THE LINES\n\n'\
	$'; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'\
	$'\n  **AND THE LINES\n\n'\
	$'; Include water topology\n#include '"$ff_watMod_dir"
	sleep 3

	echo $'\n    **TO OBTAIN:\n\n'\
	$'; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'\
	$'\n; Include ligand topology\n#include "'"$LigitpTop"$'"\n'\
	$'\n; Include water topology\n#include '"$ff_watMod_dir"
	sleep 2

	echo "$demA"$'\n    **ALSO, THE LINE\n\n'"$ligname                 1"\
	$'\n\n    **HAS BEEN WRITTEN AFTER THE LINES\n\n'"[ molecules ]"\
	$'\n; Compound        #mols\n'"$prochain"$'\n'
	sleep 2

	echo $'    **TO OBTAIN:\n\n'"[ molecules ]"\
	$'\n; Compound        #mols\n'"$prochain"\
	$'\n'"$ligname                 1"
	sleep 2

	echo "$demA"$'\nThe topol.top file before the modifications has been backed up to '\
	$'topolPreLigParmod.top.\nYou may want to check the new topol.top now to verify the '\
	$'modifications.\nWHEN YOU ARE DONE, you may proceed to the next step with a "yes" below...\n\n'

	read -p 'Do you want to proceed? (yes/no): ' procds

	while [[ "$procds" != "yes" && "$procds" != "no" && "$procds" != '"yes"' && "$procds" != '"no"' ]]
	do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p 'Do you want to proceed? (yes/no): ' procds
	done
	if [[ "$procds" == "yes" || "$procds" == '"yes"' ]] ; then echo ""
	elif [[ "$procds" == "no" || "$procds" == '"no"' ]] ; then exit 0
	fi

}	

s0Charmm()
{
	echo "$demA"$'\nChoose the stage to start/resume ligand topology preparation from\n'
cat << iniGenLigTop
 Stage  Step to begin/resume ligand topology preparation 
   a    Convert CHARMM stream file (.str) to GROMACS format (cgenff_charmm2gmx)
   b    Convert ligand_ini.pdb to ligand.gro (editconf)
   c    Prepare the protein-ligand complex
   d    Prepare the topology for the complex

*ENTER A RESPONSE BELOW USING THE APPROPRIATE OPTION

iniGenLigTop

	read -p ' Initiate ligand preparation at stage: ' stepLigPrep

	while [[ "$stepLigPrep" != 'a' && "$stepLigPrep" != 'b' && \
		"$stepLigPrep" != 'c' && "$stepLigPrep" != 'd' ]] ; do
		echo $'\nYou entered: '"$stepLigPrep"
		echo $'Please enter a valid option (a, b, c or d)!!\n'
		read -p ' Initiate ligand preparation at stage (Enter a, b, c or d): ' stepLigPrep
	done

	if [[ "$stepLigPrep" == 'a' ]] ; then s0CharmmSta; s0CharmmStb; s0CharmmStc; s0CharmmStd
	elif [[ "$stepLigPrep" == 'b' ]] ; then s0CharmmStb; s0CharmmStc; s0CharmmStd
	elif [[ "$stepLigPrep" == 'c' ]] ; then	s0CharmmStc; s0CharmmStd
	elif [[ "$stepLigPrep" == 'd' ]] ; then s0CharmmStd
	fi
}

s0PrepLigTopo()
{
	echo "$demA"$'\nChoose the stage to start/resume ligand topology preparation from\n'
cat << iniGenLigTop
 Stage  Step to begin/resume ligand topology preparation 
   a    Prepare the protein-ligand complex
   b    Prepare the topology for the complex

*ENTER A RESPONSE BELOW USING THE APPROPRIATE OPTION

iniGenLigTop

	read -p ' Initiate ligand preparation at stage: ' stepLigPrep

	while [[ "$stepLigPrep" != 'a' && "$stepLigPrep" != 'b' ]] ; do
		echo $'\nYou entered: '"$stepLigPrep"
		echo $'Please enter a valid option (a or b)!!\n'
		read -p ' Initiate ligand preparation at stage (Enter a or b): ' stepLigPrep
	done

	if [[ "$stepLigPrep" == 'a' && "$pltfrm" == 2 ]] ; then s0AcpypeSta; s0AcpypeStb
	elif [[ "$stepLigPrep" == 'b' && "$pltfrm" == 2 ]] ; then s0AcpypeStb
	elif [[ "$stepLigPrep" == 'a' && "$pltfrm" == 3 ]] ; then s0ligpargenSta; s0ligpargenStb
	elif [[ "$stepLigPrep" == 'b' && "$pltfrm" == 3 ]] ; then s0ligpargenStb
	elif [[ "$stepLigPrep" == 'a' && "$pltfrm" == 4 ]] ; then s0PrdrgSta; s0PrdrgStb
	elif [[ "$stepLigPrep" == 'b' && "$pltfrm" == 4 ]] ; then s0PrdrgStb
	fi

	if [[ "$mdType" == 2 ]] ; then
		echo "$demA"$'\n Before you proceed, ensure you modify your topology to indicate the'\
		$'immobile\n reference for the pulling simulation.\n\n'\
		$'When you are done, enter "yes" below to proceed to unit cell definition...\n'

		read -p 'Do you want to proceed? (yes/no): ' procd2boxdef

		while [[ "$procd2boxdef" != "yes" && "$procd2boxdef" != "no" ]] && \
			[[ "$procd2boxdef" != '"yes"' && "$procd2boxdef" != '"no"' ]]
		do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p 'Do you want to proceed? (yes/no): ' procd2boxdef
		done
		if [[ "$procd2boxdef" == "yes" || "$procd2boxdef" == '"yes"' ]] 
			then echo ""
		elif [[ "$procd2boxdef" == "no" || "$procd2boxdef" == '"no"' ]] 
			then exit 0
		fi
	fi	
}

s0GenLigTop()
{
	echo "$demA"$'\nWhich platform have you used for your ligand parameterization?\n'
cat << ListStageGenLigTop
 Option   Platform used for ligand topology derivation
    1     CGenFF (CHARMM)
    2     ACPYPE (AMBER)
    3     LigParGen (OPLS)
    4     PRODRG2 (GROMOS)

*ENTER A RESPONSE BELOW USING THE APPROPRIATE OPTION
ListStageGenLigTop

	read -p ' Platform used for ligand topology derivation: ' pltfrm

	while [[ "$pltfrm" != 1 && "$pltfrm" != 2 && "$pltfrm" != 3 && "$pltfrm" != 4 ]] ; do
		echo $'\nYou entered: '"$pltfrm"
		echo $'Please enter a valid option (1, 2, 3 or 4)!!\n'
		read -p ' Platform used for ligand topology derivation (Enter 1, 2, 3 or 4): ' pltfrm
	done

	if [[ "$pltfrm" == 1 ]] ; then s0Charmm
	elif [[ "$pltfrm" == 2 || "$pltfrm" == 3 || "$pltfrm" == 4 ]] ; then s0PrepLigTopo
	fi

}

s1DefBox()
{
	echo "$demA"" Define box...""$demB"
	if [[ "$mdType" == 1 ]] && [[ $sysType == 1 || $sysType == 3 ]] ; then
		eval $gmx_exe_path editconf -f ${coordinates}_processed.gro -o ${coordinates}_newbox.gro -c -d ${edgeDist} -bt ${btype}
	elif [[ "$mdType" == 1 && "$sysType" == 2 ]]; then
		eval $gmx_exe_path editconf -f ${coordinates}.gro -o ${coordinates}_newbox.gro -c -d ${edgeDist} -bt ${btype}
	elif [[ "$mdType" == 2 && "$sysType" == 1 ]] ; then
		eval $gmx_exe_path editconf -f ${coordinates}_processed.gro -o ${coordinates}_newbox.gro -c -d ${edgeDist} -bt ${btype}
	elif [[ "$mdType" == 2 && $sysType == 2 ]] ; then
		eval $gmx_exe_path editconf -f ${coordinates}.gro -o ${coordinates}_newbox.gro -c -d ${edgeDist} -bt ${btype}
	fi

	if [[ "$mdType" == 2 ]] ; then
		echo "$demA"$'\n'"A placeholder unit cell of the type ${btype} has been generated."\
		$'\nInspect the unit cell and/or provide dimensions for a new box and the'\
		$'\nplacement of the center of mass of your protein.\n'

		echo "Provide dimensions or proceed to solvation?"
		echo $'\n    1.  Provide dimensions\n    2.  Proceed\n'

		read -p 'ENTER A RESPONSE HERE (1 or 2): ' respDim

		while [[ $respDim == 1 ]] ; do
			echo ""
			read -p 'Center dimensions (x y z):' centerVector
			echo ""
			read -p 'Box dimensions (x y z):' boxVector
			echo ""
			if [[ "$sysType" == 1 ]] ; then
				eval $gmx_exe_path editconf -f ${coordinates}_processed.gro -o ${coordinates}_newbox.gro -center $centerVector -box $boxVector
			elif [[ "$sysType" == 2 ]] ; then
				eval $gmx_exe_path editconf -f ${coordinates}.gro -o ${coordinates}_newbox.gro -center $centerVector -box $boxVector
			fi
			echo "Provide dimensions or proceed to solvation?"
			echo $'\n    1.  Provide dimensions\n    2.  Proceed\n'

			read -p 'ENTER A RESPONSE HERE (1 or 2): ' respDim
		done
	fi
	echo "$demA"" Define box... DONE""$demB"
sleep 2
}

s2Solvat()
{
echo "$demA"" Solvating the system...""$demB"
sleep 2
eval $gmx_exe_path solvate -cp ${coordinates}_newbox.gro -cs spc216.gro -o ${coordinates}_solv.gro -p topol.top
echo "$demA"" Solvate system... DONE""$demB"
sleep 2
}

s3AddIons1()
{
echo "$demA"" Adding ions...""$demB"
sleep 2
#1. Assemble .tpr file with grompp, using ions.mdp
eval $gmx_exe_path grompp -f ions.mdp -c ${coordinates}_solv.gro -p topol.top -o ions.tpr -maxwarn $mw
}

s4AddIons2()
{
#2. Pass .tpr file to genion to add ions
echo 'SOL' | eval $gmx_exe_path genion -s ions.tpr -o ${coordinates}_solv_ions.gro -p topol.top -neutral ${pnam_nnam}
echo "$demA"" Add ions... DONE""$demB"
sleep 2
}

s5EnMin1()
{
echo "$demA"" Now running energy minimization...""$demB"
sleep 2
#1. Assemble the binary input with grompp
for param in *m.mdp; do
	if [[ "${param}" == "minim.mdp" ]] ; then Minmdp=1
	elif [[ "${param}" == "em.mdp" ]] ; then Enmdp=1
	fi
done
if [[ "$Minmdp" == 1 ]] && [[ "$Enmdp" == 1 ]]; then
	EnMmdp="minim.mdp"
	echo "$demA"$' Parameter files "minim.mdp" and "em.mdp" found.\n *Choosing "minim.mdp" for energy minimization.'"$demB"
	sleep 2
elif [[ "$Minmdp" == 1 ]] && [[ "$Enmdp" == 0 ]]; then
	EnMmdp="minim.mdp"
	echo "$demA"" Parameter file $EnMmdp found."$'\n'" *Using $EnMmdp for energy minimization.""$demB"
	sleep 2
elif [[ "$Minmdp" == 0 ]] && [[ "$Enmdp" == 1 ]]; then
	EnMmdp="em.mdp"
	echo "$demA"" Parameter file $EnMmdp found."$'\n'" *Using $EnMmdp for energy minimization.""$demB"
	sleep 2
	
fi
sleep 1	
eval $gmx_exe_path grompp -f $EnMmdp -c ${coordinates}_solv_ions.gro -p topol.top -o em.tpr -maxwarn $mw

}

s6EnMin2()
{
	#2. Invoke mdrun to carry out the EM
	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -deffnm em
	#tail -n 8 em.log
	echo "$demA"" Run energy minimization... DONE""$demB"
	sleep 2

	if [[ "$sysType" == 3 ]] ; then
		echo "$demA""CHAPERONg will run a check to detect the group numbers of protein and DNA...""$demB"
		sleep 2
	echo q | eval $gmx_exe_path make_ndx -f em.gro -o temp.ndx  > tempNDXfile

	while IFS= read -r tNDXline; do
		if [[ "$tNDXline" == *"DNA"*":"*"atoms"* ]] ; then
			dnaNDXno=$(echo "$tNDXline" | awk '{print $1}')
			#echo $tNDXline >> check
			#echo $dnaNDXno >> check
		elif [[ "$tNDXline" == *"Protein"*":"*"atoms"* ]] ; then
			proNDXno=$(echo "$tNDXline" | awk '{print $1}')
			#echo $tNDXline >> check
			#echo $proNDXno >> check
			break
		fi
	done < tempNDXfile
	rm tempNDXfile temp.ndx
	echo "$demA"" Checking completed...""$demB"
	echo "$demA""CHAPERONg will now make a Protein-DNA index group...""$demB"
	sleep 2

eval $gmx_exe_path make_ndx -f em.gro -o index.ndx << ProvProDNAMakNdx
$proNDXno | $dnaNDXno
q
ProvProDNAMakNdx

	echo "$demA"" Make a Protein-DNA index group...DONE""$demB"

if [[ "$sysType" == 3 ]] && [[ $flw == 1 ]]; then
	echo "$demA""CHAPERONg will now check your mdp files and make corrections on the tc-grps if needed...""$demB" 
	sleep 2
	for mdparafile in nvt.mdp npt.mdp md.mdp; do
		nbkUP=1
		mdpbak="$mdparafile"".backup."
		prevmdpbak="$mdpbak""$nbkUP"
		if [[ -f "$prevmdpbak" ]]; then
			echo "$demA""$prevmdpbak"" exists!"
			nbkUP=$(( nbkUP + 1 ))
			currentmdpbak="$mdpbak""$nbkUP"
			while [[ -f "$currentmdpbak" ]]; do
				nbkUP=$(( nbkUP + 1 ))
				currentmdpbak="$mdpbak""$nbkUP"
			done
			echo $'\nCHAPERONg will back that up as'\
			"$currentmdpbak and the $mdparafile used in the last run as $mdparafile"".backup.1""$demB"
			mv "$prevmdpbak" "$currentmdpbak"
			cp "$mdparafile" "$mdparafile"".backup.1"
		elif [[ ! -f "$prevmdpbak" ]]; then
			echo "$demA"" Backing up $mdparafile used in the last run as $mdparafile"".backup.1""$demB"
			sleep 2
			cp "$mdparafile" "$mdparafile"".backup.1"
		fi
		
		filname1="$mdparafile"".1"
		cat $mdparafile > temp_mdpfile
		while IFS= read -r mdpline; do
			checktcG=$(echo "$mdpline" | awk '{print $1}')
			if [[ "$checktcG" == "tc-grps" ]] ; then
				echo "; tc-grps below has been checked and/or modified by CHAPERONg for Protein-DNA simulation" >> "$filname1"
				echo "tc-grps                = Protein_DNA Water_and_ions    ; two coupling groups - modified by CHAPERONg" \
				>> "$filname1"
			
			elif [[ "$checktcG" == ";" ]] | [[ "$checktcG" == "title" ]]; then echo $mdpline >> "$filname1"
			elif [[ "$checktcG" == "nsteps" ]] ; then
				echo $mdpline | awk '{print $1"\t\t\t"$2,$3}'  >> "$filname1"
			else
				echo $mdpline | awk '{print $1"\t\t\t"$2,$3"\t",$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' >> "$filname1"
			fi
		done < temp_mdpfile
		rm "$mdparafile" && mv "$filname1" "$mdparafile"
	done
	rm temp_mdpfile
	echo "$demA""npt.mdp, nvt.mdp and md.mdp successfully checked and/or modified accorgingly."
	sleep 2
	echo $'\nYou may want to check these files to ensure they have been modified to your satisfaction.\n'
	$'If not, you may re-modify them accordingly or restore the backups made by CHAPERONg.\n'
	
	echo $'When you are done, enter "yes" below to proceed to equilibration...'"$demB"

	read -p 'Do you want to proceed? (yes/no): ' procd2EQ

	while [[ "$procd2EQ" != "yes" ]] && [[ "$procd2EQ" != "no" ]] && \
	[[ "$procd2EQ" != '"yes"' ]] && [[ "$procd2EQ" != '"no"' ]] ; do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p 'Do you want to proceed? (yes/no): ' procd2EQ
	done
	if [[ "$procd2EQ" == "yes" || "$procd2EQ" == '"yes"' ]] ; then echo ""
	elif [[ "$procd2EQ" == "no" || "$procd2EQ" == '"no"' ]] ; then exit 0
	fi
fi	

elif [[ "$sysType" == 2 ]] ; then

	echo "$demA"$'Confirm if you want to make custom index group(s) before proceeding.'\
	$'\n\nYou should respond with a "yes" or a "no" below.'\
	$'\n\n *If you are unsure, enter "no" and CHAPERONg will make appropriate index for'\
	$'\n you while automatically updating topol.top accordingly.'"$demB"
	sleep 3
	read -p 'Do you need to make custom index group(s) before proceeding? (yes/no): ' customNDXask

	while [[ "$customNDXask" != "yes" && "$customNDXask" != "no" && "$customNDXask" != '"yes"' && "$customNDXask" != '"no"' ]]
	do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p 'Do you need to make custom index group(s) before proceeding? (yes/no): ' customNDXask
	done

	if [[ "$customNDXask" == "yes" || "$customNDXask" == $'"yes"' ]] ; then makeNDXGroup
	elif [[ "$customNDXask" == "no" || "$customNDXask" == $'"no"' ]] ; then echo ""; fi
fi
}

s6aLigRes()
{
	modTop1
	#if user did not create an index manually
	if [[ "$customNDXask" == "no" || "$customNDXask" == $'"no"' || "$customNDXask" == '' ]] ; then
		echo "$demA"" Restraining the ligand...""$demB""$demA"$'**"In the coming steps, '\
		$'CHAPERONg will create an index group for the non-hydrogen atoms of $ligname'"$demB"
		sleep 2

eval $gmx_exe_path make_ndx -f "$lignamFile".gro -o index_"$lignamFile".ndx << ProvMakNdx
0 & ! a H*
q
ProvMakNdx
	fi
	echo "$demA"" Executing the genrestr module""$demB""$demA"$'**CHAPERONg will select '\
	$'the newly created index group i.e. group 3 in the index_'"$lignamFile.ndx file""$demB"
	sleep 2
	
	if [[ "$customNDXask" == "no" || "$customNDXask" == $'"no"' || "$customNDXask" == '' ]] ; then
		echo "System_&_!H*" | eval $gmx_exe_path genrestr -f $lignamFile.gro -n index_$lignamFile.ndx -o posre_$lignamFile.itp -fc 1000 1000 1000

		if [[ -f "$ligname""_GMX.itp" ]] && [[ -d acpype ]] ; then modTop2forAcpype
		else
			modTop2
		fi
	
		echo "$demA"$' To restrain the protein and ligand simultaneously,\n'\
		$'the following lines have been written to your topol.top file:'"$demB"\
		$'\n; Ligand position restraints\n#ifdef POSRES'\
		$'\n#include '$'"'$'posre_'"$lignamFile"$'.itp'$'"'\
		$'\n#endif\n'

		sleep 2
 
		echo "$demA"$' The above were written'"$demB""$demA"$'  **BETWEEN the lines:'"$demB"\
		$'\n; Include ligand topology\n'\
		$'#include '$'"'"$lignamFile"$'.itp'$'"'
		sleep 2

		if [[ -f "$ligname""_GMX.itp" ]] && [[ -d acpype ]]; then
			echo $'\n'"$demA"$'  **AND:'"$demB"$'\n[ moleculetype ]\n; Name            nrexcl'
			sleep 2
			echo "$demA"$'  **TO OBTAIN:'"$demB"$'; Include ligand topology'
			echo "#include "$'"'$lignamFile.itp$'"'
			sleep 2
			echo $'\n; Ligand position restraints\n#ifdef POSRES\n'\
			$'#include '$'"'$'posre_'"$lignamFile"$'.itp'$'"'\
			$'\n#endif\n\n'\
			$'\n[ moleculetype ]\n; Name            nrexcl'
		else
			echo $'\n'"$demA"$'  **AND:'"$demB"$'\n; Include water topology\n#include '"$ff_watMod_dir"
			sleep 2
			echo "$demA"$'  **TO OBTAIN:'"$demB"$'; Include Position restraint file\n#ifdef POSRES\n#include '\
			$'"posre.itp"\n#endif\n\n; Include ligand topology'
			echo "#include "$'"'$lignamFile.itp$'"'
			sleep 2
			echo $'\n; Ligand position restraints\n#ifdef POSRES\n'\
			$'#include '$'"'$'posre_'"$lignamFile"$'.itp'$'"'\
			$'\n#endif\n\n'\
			$'; Include water topology\n#include '"$ff_watMod_dir"
		fi

		sleep 2
		echo "$demA"$'\n The "topol.top" file before modification was'\
		$'backed up as "topolPreResmod.top"\n\n'
		sleep 2

		echo $' You may want to check your "topol.top" file now to verify the modifications.\n'\
		$'WHEN YOU ARE DONE, you can proceed to temperature coupling with a "yes" below\n'

		if [[ "$mdType" == 1 ]] ; then
			read -p 'Do you want to proceed to temperature coupling? (yes/no): ' restr 

			while [[ "$restr" != "yes" && "$restr" != "no" && "$restr" != '"yes"' && "$restr" != '"no"' ]]
			do
				echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
				read -p ' Do you want to proceed to temperature coupling? (yes/no) ' restr 
			done
		fi

		if [[ "$restr" == "yes" || "$restr" == $'"yes"' ]] ; then continue
			
		elif [[ "$restr" == "no" || "$restr" == $'"no"' ]] ; then exit 0
		fi

	elif [[ "$customNDXask" != "no" || "$customNDXask" != $'"no"' ]] ; then
		eval $gmx_exe_path genrestr -f $lignamFile.gro -n $ndxNAME -o "posre_"$ndxNAME -fc 1000 1000 1000
	fi

}

s6bTempCoup()
{
	if [[ "$pn" == "" && "$nn" == "" ]] ; then
		echo "$demA"" Coupling the protein with $ligname...""$demB"
		sleep 2
	#Make index for temperature coupling
eval $gmx_exe_path make_ndx -f em.gro -o index.ndx << TempCoupIN1
"Protein" | 13
q
TempCoupIN1

	elif [[ "$pn" == "NA" && "$nn" == "CL" ]]
	then echo "$demA"" Coupling the protein with $ligname...""$demB"
	sleep 2
#Make index for temperature coupling
eval $gmx_exe_path make_ndx -f em.gro -o index.ndx << TempCoupIN2
"Protein" | 13
q
TempCoupIN2

	#in the case of the charmm-jul-2021 ff
	elif [[ "$pnam_nnam" != '' && "$pn" != "NA" && "$nn" != "CL" ]]
	then echo "$demA"" Coupling the protein with $ligname..."$'\n Will '\
	$'also couple water with '"$pn and\or $nn""$demB"
	sleep 3
#Make index for temperature coupling
eval $gmx_exe_path make_ndx -f em.gro -o index.ndx << TempCoupIN3
"Protein" | 13
"Water" | "$pn" | "$nn"
q
TempCoupIN3

		if [[ "$sysType" == 2 && "$flw" == 1 ]] ; then
			echo $'\n'"$demA"$'Running a check to detect the type(s) of'\
			$'ions that GMX added to your system...\n\n'
			sleep 2
			tail -n 15 topol.top > tempMolcfile
			ionAdded1=''
			ionAdded2=''
			tcgrpWt='Water'
			while IFS= read -r Molcline; do
				molecl=$(echo $Molcline | awk '{print $1}')
				if [[ "$molecl" == "$pn" ]] ; then
					tcgrpWt="$tcgrpWt""_""$molecl"
					ionAdded1="$molecl"
					sleep 2
				elif [[ "$molecl" == "$nn" ]] ; then
					tcgrpWt="$tcgrpWt""_""$molecl"
					ionAdded2="$molecl"
					sleep 2
				fi
			done < tempMolcfile
			rm tempMolcfile
			if [[ "$ionAdded1" != '' && "$ionAdded2" != '' ]] ; then
				echo $' Checking completed...\n The ions'\
				$''"$ionAdded1 and $ionAdded2 were detected in your topol.top file..."$'\n\n'
				sleep 2
			elif [[ "$ionAdded1" != '' || "$ionAdded2" != '' ]] ; then
				echo $' Checking completed...\n The ion'\
				$''"$ionAdded1""$ionAdded2 was detected in your topol.top file"$'\n'
				sleep 2
			fi	
			echo $'\nChecking your mdp files to make corrections to the tc-grps if needed...\n'
			sleep 2
			if [[ "$mdType" == 1 ]] ; then
				for mdparafile in nvt.mdp npt.mdp md.mdp; do
					nbkUP=1
					mdpbak="$mdparafile"".backup."
					prevmdpbak="$mdpbak""$nbkUP"
					if [[ -f "$prevmdpbak" ]]; then
						echo $'\n'"$prevmdpbak"" exists!"
						nbkUP=$(( nbkUP + 1 ))
						currentmdpbak="$mdpbak""$nbkUP"
						while [[ -f "$currentmdpbak" ]]; do
							nbkUP=$(( nbkUP + 1 ))
							currentmdpbak="$mdpbak""$nbkUP"
						done
						echo $'\nCHAPERONg will back that up as '"$currentmdpbak"$' and'\
						$'the '"$mdparafile used in the last run as $mdparafile"".backup.1"$'\n'
						mv "$prevmdpbak" "$currentmdpbak"
						cp "$mdparafile" "$mdparafile"".backup.1"
					elif [[ ! -f "$prevmdpbak" ]]; then
						echo $'\n'" Backing up $mdparafile used in the last run as $mdparafile"".backup.1"$'\n'
						sleep 2
						cp "$mdparafile" "$mdparafile"".backup.1"
					fi
		
					filname1="$mdparafile"".1"
					cat $mdparafile > temp_mdpfile
					while IFS= read -r mdpline; do
						checktcG=$(echo "$mdpline" | awk '{print $1}')
						if [[ "$checktcG" == "tc-grps" || "$checktcG" == "tc_grps" ]] ; then
							echo "; tc-grps below has been checked and/or modified by CHAPERONg for Protein-Ligand simulation" >> "$filname1"
							echo "tc-grps                = Protein_$ligname   $tcgrpWt    ; two coupling groups - modified by CHAPERONg" >> "$filname1"
			
						elif [[ "$checktcG" == ";" ]] | [[ "$checktcG" == "title" ]]; then echo $mdpline >> "$filname1"
						elif [[ "$checktcG" == "nsteps" ]] ; then
							echo $mdpline | awk '{print $1"\t\t\t"$2,$3}'  >> "$filname1"
						else
							echo $mdpline | awk '{print $1"\t\t\t"$2,$3"\t",$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' >> "$filname1"
						fi
					done < temp_mdpfile
					rm "$mdparafile" && mv "$filname1" "$mdparafile"
				done
				rm temp_mdpfile
				echo "$demB""$demA""npt.mdp, nvt.mdp and md.mdp successfully checked and/or modified accorgingly."
				sleep 2

			elif [[ ! -f nvt.mdp ]] && [[ "$mdType" == 2 ]] ; then
				for mdparafile in npt.mdp md_pull.mdp npt_umbrella.mdp md_umbrella.mdp
				do
					nbkUP=1
					mdpbak="$mdparafile"".backup."
					prevmdpbak="$mdpbak""$nbkUP"
					if [[ -f "$prevmdpbak" ]]; then
						echo $'\n'"$prevmdpbak"" exists!"
						nbkUP=$(( nbkUP + 1 ))
						currentmdpbak="$mdpbak""$nbkUP"
						while [[ -f "$currentmdpbak" ]]; do
							nbkUP=$(( nbkUP + 1 ))
							currentmdpbak="$mdpbak""$nbkUP"
						done
						echo $'\nCHAPERONg will back that up as '"$currentmdpbak"$' and'\
						$'the '"$mdparafile used in the last run as $mdparafile"".backup.1"$'\n'
						mv "$prevmdpbak" "$currentmdpbak"
						cp "$mdparafile" "$mdparafile"".backup.1"
					elif [[ ! -f "$prevmdpbak" ]]; then
						echo $'\n'" Backing up $mdparafile used in the last run as $mdparafile"".backup.1"$'\n'
						sleep 2
						cp "$mdparafile" "$mdparafile"".backup.1"
					fi
		
					filname1="$mdparafile"".1"
					cat $mdparafile > temp_mdpfile
					while IFS= read -r mdpline; do
						checktcG=$(echo "$mdpline" | awk '{print $1}')
						if [[ "$checktcG" == "tc-grps" || "$checktcG" == "tc_grps" ]] ; then
							echo "; tc-grps below has been checked and/or modified by CHAPERONg for Protein-Ligand simulation" >> "$filname1"
							echo "tc-grps                = Protein_$ligname   $tcgrpWt    ; two coupling groups - modified by CHAPERONg" >> "$filname1"
			
						elif [[ "$checktcG" == ";" ]] | [[ "$checktcG" == "title" ]]; then echo $mdpline >> "$filname1"
						elif [[ "$checktcG" == "nsteps" ]] ; then
							echo $mdpline | awk '{print $1"\t\t\t"$2,$3}'  >> "$filname1"
						else
							echo $mdpline | awk '{print $1"\t\t\t"$2,$3"\t",$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' >> "$filname1"
						fi
					done < temp_mdpfile
					rm "$mdparafile" && mv "$filname1" "$mdparafile"
				done
				rm temp_mdpfile
				echo "$demB""$demA""The files npt.mdp, md_pull.mdp, npt_umbrella.mdp, and md_umbrella.mdp have"\
				$'\n been checked and/or modified accorgingly.'
				sleep 2
			
			elif [[ -f nvt.mdp ]] && [[ "$mdType" == 2 ]] ; then
				for mdparafile in nvt.mdp npt.mdp md_pull.mdp npt_umbrella.mdp md_umbrella.mdp
				do
					nbkUP=1
					mdpbak="$mdparafile"".backup."
					prevmdpbak="$mdpbak""$nbkUP"
					if [[ -f "$prevmdpbak" ]]; then
						echo $'\n'"$prevmdpbak"" exists!"
						nbkUP=$(( nbkUP + 1 ))
						currentmdpbak="$mdpbak""$nbkUP"
						while [[ -f "$currentmdpbak" ]]; do
							nbkUP=$(( nbkUP + 1 ))
							currentmdpbak="$mdpbak""$nbkUP"
						done
						echo $'\nCHAPERONg will back that up as '"$currentmdpbak"$' and'\
						$'the '"$mdparafile used in the last run as $mdparafile"".backup.1"$'\n'
						mv "$prevmdpbak" "$currentmdpbak"
						cp "$mdparafile" "$mdparafile"".backup.1"
					elif [[ ! -f "$prevmdpbak" ]]; then
						echo $'\n'" Backing up $mdparafile used in the last run as $mdparafile"".backup.1"$'\n'
						sleep 2
						cp "$mdparafile" "$mdparafile"".backup.1"
					fi
		
					filname1="$mdparafile"".1"
					cat $mdparafile > temp_mdpfile
					while IFS= read -r mdpline; do
						checktcG=$(echo "$mdpline" | awk '{print $1}')
						if [[ "$checktcG" == "tc-grps" || "$checktcG" == "tc_grps" ]] ; then
							echo "; tc-grps below has been checked and/or modified by CHAPERONg for Protein-Ligand simulation" >> "$filname1"
							echo "tc-grps                = Protein_$ligname   $tcgrpWt    ; two coupling groups - modified by CHAPERONg" >> "$filname1"
			
						elif [[ "$checktcG" == ";" ]] | [[ "$checktcG" == "title" ]]; then echo $mdpline >> "$filname1"
						elif [[ "$checktcG" == "nsteps" ]] ; then
							echo $mdpline | awk '{print $1"\t\t\t"$2,$3}'  >> "$filname1"
						else
							echo $mdpline | awk '{print $1"\t\t\t"$2,$3"\t",$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' >> "$filname1"
						fi
					done < temp_mdpfile
					rm "$mdparafile" && mv "$filname1" "$mdparafile"
				done
				rm temp_mdpfile
				echo "$demB""$demA""The files nvt.mdp, npt.mdp, md_pull.mdp, npt_umbrella.mdp, and md_umbrella.mdp"\
				$'\nhave been checked and/or modified accorgingly.'
				sleep 2
			fi
			echo $'\nYou may want to verify the modifications in these files. If unsatisfied, you'\
			$'\nmay re-modify them accordingly or restore the backups made by CHAPERONg.\n'\
				
			echo $'\n *When you are done, enter "yes" below to proceed to equilibration...'"$demB"

		elif [[ "$sysType" == 2 && "$flw" == 0 ]] ; then
			echo "$demA"$'Because you are not running CHAPERONg in the "auto mode", you may\n'\
			$'need to confirm that you have set the tc-groups in your mdp files as necessary'
		fi
		
		read -p 'Do you want to proceed? (yes/no): ' procd2usEQ

		while [[ "$procd2usEQ" != "yes" && "$procd2usEQ" != "no" ]] && \
				[[ "$procd2usEQ" != '"yes"' && "$procd2usEQ" != '"no"' ]] 
		do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p 'Do you want to proceed? (yes/no): ' procd2usEQ
		done
		if [[ "$procd2usEQ" == "yes" || "$procd2usEQ" == $'"yes"' ]] ; then echo ""
		elif [[ "$procd2usEQ" == "no" || "$procd2usEQ" == $'"no"' ]] ; then exit 0
		fi
	fi	
}

s7NVTeq1()
{
	if [[ $mdType == 2 ]] ; then
		echo "$demA"$'\n'
		read -p 'Do you want to run an optional NVT equilibration? (yes/no): ' runNVTeq
		
		while [[ "$runNVTeq" != "yes" && "$runNVTeq" != "no" && \
			"$runNVTeq" != '"yes"' && "$runNVTeq" != '"no"' ]]
		do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p 'Do you want to run an optional NVT equilibration? (yes/no): ' runNVTeq
		done
	elif [[ $mdType == 1 ]] ; then runNVTeq="yes"
	fi
	echo "$demA"" Now running NVT Equilibration...""$demB"
	sleep 2

	if [[ $sysType == 1 && $runNVTeq == "yes" ]] ; then
		eval $gmx_exe_path grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn $mw
	elif [[ $sysType == 2 || $sysType == 3 ]] && [[ $runNVTeq == "yes" ]] ; then
		eval $gmx_exe_path grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn $mw
	elif [[ $runNVTeq == "no" || $runNVTeq != '"no"' ]] ; then
		echo "$demA"$' Skipping NVT equilibration, proceeding to NPT equilibration'"$demB"
	fi
}

s8NVTeq2()
{
	if [[ $runNVTeq == "yes" ]] ; then
		eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm nvt
		echo "$demA"" NVT Equilibration... DONE""$demB"
		sleep 2
	fi
}

s9NPTeq1()
{
	echo "$demA"" Now running NPT Equilibration...""$demB"
	if [[ $sysType == 1 && $runNVTeq == "yes" ]] ; then
		eval $gmx_exe_path grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn $mw
	elif [[ $sysType == 2 || $sysType == 3 ]] && [[ $runNVTeq == "yes" ]] ; then
		eval $gmx_exe_path grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn $mw
	#if NVT equilibration was skipped
	elif [[ $sysType == 1 && $runNVTeq == "no" ]] ; then
		eval $gmx_exe_path grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -o npt.tpr -maxwarn $mw
	elif [[ $sysType == 2 || $sysType == 3 ]] && [[ $runNVTeq == "no" ]] ; then
		eval $gmx_exe_path grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn $mw
	fi
}

s10NPTeq2()
{
	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm npt
	echo "$demA"" NPT Equilibration... DONE""$demB"
	sleep 2
	if [[ $sysType == 1 && $mdType == 2 ]] ; then
		read -p 'Do you need to make custom index for pulling groups? (yes/no): ' cstmndx
		
		while [[ "$cstmndx" != "yes" && "$cstmndx" != "no" && \
			"$cstmndx" != '"yes"' && "$cstmndx" != '"no"' ]]
		do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p 'Do you need to make custom index for pulling groups? (yes/no): ' cstmndx
		done
		if [[ $cstmndx == "yes" ]] ; then
			echo "$demA"" Will now make an index for the pulling groups...""$demB"
			sleep 2
			eval $gmx_exe_path make_ndx -f npt.gro
			echo "$demA"" Make an index for the pulling groups...DONE""$demB"
		fi
	fi
}

s11RelPosRe()
{
	echo "$demA"" Releasing position restraints...""$demB"
	sleep 2
	if [[ $sysType == 1 ]]; then
		eval $gmx_exe_path grompp -v -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o ${filenm}.tpr -maxwarn $mw
	elif [[ $sysType == 2 || $sysType == 3 ]] ; then
		eval $gmx_exe_path grompp -v -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ${filenm}.tpr -maxwarn $mw
	fi
	echo "$demA"" Release position restraints... DONE""$demB"
	sleep 2
}

umbre_s11_SMD1()
{
	echo "$demA"" Now executing grompp for steered MD simulation...""$demB"
	sleep 2
	eval $gmx_exe_path grompp -f md_pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr -maxwarn $mw
}

umbre_s12_SMD2()
{
	echo "$demA"" Now running steered MD simulation...""$demB"
	sleep 2
	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm pull -pf pullf.xvg -px pullx.xvg

	echo "$demA"$'Generating finished figures of key results of the pulling simulation...'"$demB"
	sleep 2
	cat pullx.xvg | grep -v "^[@#]" | awk '{ print $2 }' > displacementAxis
	cat pullf.xvg | grep -v "^[@#]" | awk '{ print $2 }' > pullforceAxis
	echo "# The Pull Force and Displacement data in this file were extracted by CHAPERONg" > displacement_pullForce.xvg
	echo "# from the pullx.xvg and pullf.xvg files generated by GROMACS..." >> displacement_pullForce.xvg
	echo "#" >> displacement_pullForce.xvg
	echo "@    title "$'"Plot of Pull force against Displacement"' >> displacement_pullForce.xvg
	echo "@    xaxis  label "$'"Displacement (nm)"' >> displacement_pullForce.xvg
	echo "@    yaxis  label "$'"Force (kJ/mol/nm)"' >> displacement_pullForce.xvg
	echo "@TYPE xy" >> displacement_pullForce.xvg
	paste -d "        " displacementAxis pullforceAxis >> displacement_pullForce.xvg
	rm displacementAxis pullforceAxis

	gracebat pullf.xvg -hdevice PNG -autoscale xy -printfile pullf.png -fixed 7500 4000 -legend load || true

	gracebat displacement_pullForce.xvg -hdevice PNG -autoscale xy -printfile displacement_pullForce.png -fixed 7500 4000 -legend load || true
	
	AnaName="steered_MD"
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
		mv pullf.png pullf.xvg pullx.xvg ./$AnaName || true
		mv displacement_pullForce.xvg displacement_pullForce.png ./$AnaName || true	
	elif [[ ! -d "$currentAnadir" ]]; then mkdir ./$AnaName
		mv pullf.png pullf.xvg pullx.xvg ./$AnaName || true
		mv displacement_pullForce.xvg displacement_pullForce.png ./$AnaName || true
	fi

	echo "$demA"$'Generate finished figures of key results of the pulling simulation... DONE'
	sleep 2
	echo "$demA"" Steered MD simulation completed...""$demB"
	sleep 2
}

s12MDrun()
{
	echo "$demA"" Now running production MD...""$demB"
	if [[ $nohp == '' ]]; then
		eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm}
	elif [[ $nohp == 1 ]]; then
		echo "$demA"$' gmx mdrun will be launched and run in the background.\n\nYou can tail the '\
		$'output file "nohup.txt" to monitor the progress of the simulation at any time...\n'"$demB"
		sleep 2
	
		nohup eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm}
	
		echo $'\n'"$demA""CHAPERONg: Simulation Completed...""$demB"
		Credit
		echo $'\n'
		read -p ' Do you want to run any post-simulation processing/analyses?(yes/no): ' psa

		while [[ "$psa" != "yes" && "$psa" != "no" && "$psa" != '"yes"' && "$psa" != '"no"' ]] ; do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p '*Do you want to run any post-simulation processing/analyses?(yes/no): ' psa 
		done

		if test "$psa" == "yes"; then Analysis ; fi

	elif [[ $nohp == '' ]]; then
	
		eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm}
	
		echo $'\n'"$demA""CHAPERONg: Simulation Completed...""$demB"
		Credit
		echo $'\n'
		read -p '*Do you want to run any post-simulation processing/analyses?(yes/no): ' psa

		while [[ "$psa" != "yes" && "$psa" != "no" && "$psa" != '"yes"' && "$psa" != '"no"' ]] ; do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p '**Do you want to run any post-simulation processing/analyses?(yes/no): ' psa 
		done

		if test "$psa" == "yes"; then Analysis; fi

	fi
	echo $'\n'"$demA""CHAPERONg: Simulation Completed...""$demB"
	Credit
	echo $'\n'
	read -p 'Do you want to run any post-simulation processing/analyses?(yes/no): ' psa

	while [[ "$psa" != "yes" ]] && [[ "$psa" != "no" ]] && [[ "$psa" != '"yes"' ]] && [[ "$psa" != '"no"' ]] ; do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p '**Do you want to run any post-simulation processing/analyses?(yes/no): ' psa 
	done

	if test "$psa" == "yes"; then Analysis; fi
}

s13MDappend()
{
	echo "$demA"" CHAPERONg will now run an extended production MD...""$demB"

if [[ $nohp == 1 ]]; then
	
	echo "$demA"$' gmx mdrun will be launched and run in the background.\n\n\
	You can tail the output file "nohup.txt" to monitor the progress of the simulation at any time...\n'"$demB"
	sleep 2
	
	nohup eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm} -s ${filenm}."tpr" -cpi ${filenm}."cpt" -append
	
	echo $'\n'"$demA"" Simulation extended/appended successfully...""$demB"
	Credit
	echo $'\n'
	read -p '**Do you want to run any post-simulation processing/analyses?(yes/no): ' psa

	while [[ "$psa" != "yes" ]] && [[ "$psa" != "no" ]] && [[ "$psa" != '"yes"' ]] && [[ "$psa" != '"no"' ]] ; do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p '**Do you want to run any post-simulation processing/analyses?(yes/no): ' psa 
	done

	if test "$psa" == "yes"; then Analysis ; fi

elif [[ $nohp == '' ]]; then
	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm} -s ${filenm}."tpr" -cpi ${filenm}."cpt" -append
		
	#echo $'\n'"$demA""CHAPERONg: Simulation Completed...""$demB"
	#Credit
	#echo $'\n'
	#read -p '**Do you want to run any post-simulation processing/analyses?(yes/no): ' psa

	#while [[ "$psa" != "yes" ]] && [[ "$psa" != "no" ]] && [[ "$psa" != '"yes"' ]] && [[ "$psa" != '"no"' ]] ; do
	#	echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
	#	read -p '*Do you want to run any post-simulation processing/analyses?(yes/no): ' psa 
	#done

	#if test "$psa" == "yes"; then Analysis
	#fi

fi
echo $'\n'"$demA""Simulation extended/appended successfully""$demB"
Credit
echo $'\n'
read -p '*Do you want to run any post-simulation processing/analyses?(yes/no): ' psa


while [[ "$psa" != "yes" ]] && [[ "$psa" != "no" ]] && [[ "$psa" != '"yes"' ]] && [[ "$psa" != '"no"' ]] ; do
	echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
	un umbrella sampling for an additional window
done

if test "$psa" == "yes"; then Analysis
fi
}

umbre_s13_xtractFrames()
{
	echo "$demA"" Now extracting frames from the steered MDS trajectory...""$demB"
	sleep 2
	echo 0 | eval $gmx_exe_path trjconv -s pull.tpr -f pull.xtc -o coordinate.gro -sep
	echo "$demA"" Extract frames from the steered MDS trajectory...DONE""$demB"
	sleep 2
}

umbre_s14_calcCOMdist()
{
	echo "$demA"" Now calculating COM distances...""$demB"
	sleep 2
	group1_name=$(cat md_pull.mdp | grep "pull_group1_name" | awk '{print $3}')
	group2_name=$(cat md_pull.mdp | grep "pull_group2_name" | awk '{print $3}')
	StructNo=0
	mkdir coordinates_SMD
	mkdir distances_SMD

	for Structure in ./"coordinate"*".gro" ; do
		#calculate distance between the groups
		eval $gmx_exe_path distance -s pull.tpr -f coordinate"$StructNo".gro -n index.ndx -select \
		"com of group $group1_name plus com of group $group2_name" -oall dist${StructNo}.xvg
		if [[ $StructNo == 50 ]] ; then sleep 2 ; fi
		#extract the distances into a summary file
		distanc=$(tail -n 1 dist${StructNo}.xvg | awk '{print $2}')
		if [[ $StructNo == 0 ]] ; then
			echo "$StructNo"$'\t'"$distanc" > distances_summary.txt
		else
			echo "$StructNo"$'\t'"$distanc" >> distances_summary.txt
		fi
		mv dist${StructNo}.xvg ./distances_SMD
		mv coordinate"$StructNo".gro ./coordinates_SMD
		StructNo=$(( StructNo + 1 ))
	done
	echo "$demA"" Calculate COM distances...DONE""$demB"
	sleep 2
}

umbre_s15_findIniConf()
{
	echo "$demA"" Now identifying initial configurations for umbrella sampling...""$demB"
	#create a new configuratns_list.txt file and backup existing one
	if [[ ! -f "configuratns_list.txt" ]]; then
		touch configuratns_list.txt
		echo $'\n'"$demA"$' An empty configuratns_list.txt file has been created\n'
	elif [[ -f "configuratns_list.txt" ]]; then
		configList="configuratns_list.txt"
		nbkUP=1
		configListbak="$configList"".backup."
		prevconfigListbak="$configListbak""$nbkUP"
		if [[ -f "$prevconfigListbak" ]]; then
			echo "$demA""$prevconfigListbak"" exists!"
			nbkUP=$(( nbkUP + 1 ))
			currentconfigListbak="$configListbak""$nbkUP"
			while [[ -f "$currentconfigListbak" ]]; do
				nbkUP=$(( nbkUP + 1 ))
				currentconfigListbak="$configListbak""$nbkUP"
			done
			echo $'\nCHAPERONg will back that up as '"$currentconfigListbak"$' and'\
			$'the '"$configList used in the last run as $configList"".backup.1""$demB"
			mv "$prevconfigListbak" "$currentconfigListbak"
			cp "$configList" "$configList"".backup.1"
		elif [[ ! -f "$prevconfigListbak" ]]; then
			echo "$demA"" Backing up $configList used in the last run as $configList"".backup.1""$demB"
			sleep 2
			cp "$configList" "$configList"".backup.1"
		fi
		touch configuratns_list.txt
	fi

	#run extract_spaced_frame_dist.py script
	echo " Running the extract_spaced_frame_dist.py script..."
	sleep 2
	python3 extract_spaced_frame_dist.py || python extract_spaced_frame_dist.py

	echo " Run extract_spaced_frame_dist.py...DONE""$demB"
	sleep 2
}

umbre_s16_USampling()
{
	echo "$demA"" Initiating umbrella sampling...""$demB"
	window=0
	while IFS= read -r line; do
		if [[ $line == *"#"* ]] ; then continue
		elif [[ $line != *"#"* ]] ; then
			us_frame=$(echo "$line" | awk '{print $1}')
			echo "$demA""Now running NPT equilibration for configuration $us_frame"
			eval $gmx_exe_path grompp -f npt_umbrella.mdp -c ./coordinates_SMD/coordinate"$us_frame".gro -p topol.top -r \
			./coordinates_SMD/coordinate"$us_frame".gro -n index.ndx -o npt_win"$window"_conf"$us_frame".tpr -maxwarn $mw

			eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm npt_win"$window"_conf"$us_frame"

			echo "$demA""Run NPT equilibration for configuration $us_frame...DONE""$demB"
		
			echo "$demA""Now running umbrella sampling for configuration $us_frame"$'\n\n'
			eval $gmx_exe_path grompp -f md_umbrella.mdp -c npt_win"$window"_conf"$us_frame".gro -t npt_win"$window"_conf"$us_frame".cpt -p \
			topol.top -r npt_win"$window"_conf"$us_frame".gro -n index.ndx -o umbrella_win"$window"_conf"$us_frame".tpr -maxwarn 1

			eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm umbrella_win"$window"_conf"$us_frame"

			echo "$demA""Run umbrella sampling for configuration $us_frame...DONE""$demB"
		fi
		if [[ $window == 0 ]] ; then
			echo "umbrella_win"$window"_conf"$us_frame".tpr" > tpr_files.dat
			echo "umbrella_win"$window"_conf"$us_frame"_pullf.xvg" > pullf_files.dat
			echo "umbrella_win"$window"_conf"$us_frame"_pullx.xvg" > pullx_files.dat
		elif [[ $window > 0 ]] ; then
			echo "umbrella_win"$window"_conf"$us_frame".tpr" >> tpr_files.dat
			echo "umbrella_win"$window"_conf"$us_frame"_pullf.xvg" >> pullf_files.dat
			echo "umbrella_win"$window"_conf"$us_frame"_pullx.xvg" >> pullx_files.dat
		fi
	
		window=$(( window + 1 ))

	done < configuratns_list.txt

	echo "$demA"$' Umbrella sampling simulation has been completed for all the windows.\n\n'\
	$'The files tpr_files.dat, pullf_files.dat and pullx_files.dat have been written.\n'\
	$'These files list, respectively, the names of the .tpr, pullf.xvg and\n'\
	$'pullx.xvg files for the windows.\n'\
	$'*The files tpr_files.dat and pullf_files.dat will be used for data analysis...'"$demB"
	#mv npt_win*_pull*.xvg umbrella_sampling/equilibration/
	#mv umbrella_win*_pullx.xvg ./umbrella_sampling/data_collection

}

umbre_s17_WHAM()
{
	echo "$demA"" Now extracting the PMF and plotting the umbrella histograms...""$demB"
	sleep 2
	eval $gmx_exe_path wham -it tpr_files.dat -if pullf_files.dat -o -hist -unit kCal

	echo "$demA"$'Generating finished figures of key results of WHAM analysis...'"$demB"
	sleep 2
	gracebat -nxy histo.xvg -hdevice PNG -autoscale xy -printfile histograms.png -fixed 7500 4000 || true
	gracebat profile.xvg -hdevice PNG -autoscale xy -printfile profile_PMF.png -fixed 7500 4000 -legend load || true

	AnaName="Data_Analysis_PMF"
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
		mv histo.xvg histograms.png ./$AnaName || true
		mv profile.xvg profile_PMF.png ./$AnaName || true
	elif [[ ! -d "$currentAnadir" ]]; then mkdir ./$AnaName
		mv histo.xvg histograms.png ./$AnaName || true
		mv profile.xvg profile_PMF.png ./$AnaName || true
	fi

	echo "$demA"$'Generate finished figures of results of WHAM analysis...DONE'
	sleep 2
	echo "$demA"" Extract the PMF and plot the umbrella histograms...DONE""$demB"
	sleep 2
}

umbre_s18_MoreWin()
{
	echo "$demA"" Running umbrella sampling for an additional window...""$demB"
	sleep 2
	read -p ' Enter the frame number of the SMD configuration to use: ' us_frame
	echo $'\n You entered: '"$us_frame"
	sleep 2

	window=0
	while [[ -f npt_win"$window"_* ]]; do
		window=$(( window + 1 ))
	done
		
	echo "$demA""Now running NPT equilibration for configuration $us_frame"
	eval $gmx_exe_path grompp -f npt_umbrella.mdp -c ./coordinates_SMD/coordinate"$us_frame".gro -p topol.top -r \
	./coordinates_SMD/coordinate"$us_frame".gro -n index.ndx -o npt_win"$window"_conf"$us_frame".tpr -maxwarn $mw

	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm npt_win"$window"_conf"$us_frame"

	echo "$demA""Run NPT equilibration for configuration $us_frame...DONE""$demB"

	echo "$demA""Now running umbrella sampling for configuration $us_frame"$'\n\n'
	eval $gmx_exe_path grompp -f md_umbrella.mdp -c npt_win"$window"_conf"$us_frame".gro -t npt_win"$window"_conf"$us_frame".cpt -p \
	topol.top -r npt_win"$window"_conf"$us_frame".gro -n index.ndx -o umbrella_win"$window"_conf"$us_frame".tpr -maxwarn 1

	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm umbrella_win"$window"_conf"$us_frame"

	echo "$demA""Run umbrella sampling for configuration $us_frame...DONE""$demB"
	echo "umbrella_win"$window"_conf"$us_frame".tpr" >> tpr_files.dat
	echo "umbrella_win"$window"_conf"$us_frame"_pullf.xvg" >> pullf_files.dat
	echo "umbrella_win"$window"_conf"$us_frame"_pullx.xvg" >> pullx_files.dat
	echo "$demA"$' Umbrella sampling simulation has been completed for the additional window.\n\n'\
	$' The names of the .tpr, pullf.xvg and pullx.xvg files have been appended to the\n'\
	$' tpr_files.dat, pullf_files.dat and pullx_files.dat, respectively.\n'\
	$'*The files tpr_files.dat and pullf_files.dat will be used in the data analysis step...'"$demB"

	read -p ' Do you want to proceed to PMF calculation? (yes/no): ' procd2pmf

	while [[ "$procd2pmf" != "yes" && "$procd2pmf" != "no" && \
		"$procd2pmf" != '"yes"' && "$procd2pmf" != '"no"' ]] ; do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p ' Do you want to proceed? (yes/no): ' procd2pmf
		done
		if [[ "$procd2pmf" == "yes" || "$procd2pmf" == '"yes"' ]] ; then echo ""
		elif [[ "$procd2pmf" == "no" || "$procd2pmf" == '"no"' ]] ; then exit 0
		fi

}