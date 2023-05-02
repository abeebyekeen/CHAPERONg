#! /bin/bash


######################################################################
#  CHAP_deffxn - The function definition module of CHAPERONg         #
#  CHAPERONg -- An automation program for GROMACS md simulation and  #
#    trajectory analysis                                             #
#  Author -- Abeeb A. Yekeen (abeeb.yekeen@hotmail.com)              #
#  Date: 2022.02.11                                                  #
######################################################################

set -e
set -o pipefail

#demA=$'\n\n'"#**********************************CHAPERONg**********************************#"$'\n'
demA=$'\n\n'"#================================= CHAPERONg =================================#"$'\n'
demB=$'\n'"#=============================================================================#"$'\n\n'
#demB=$'\n'"#*****************************************************************************#"$'\n\n'


if [[ $initiator1 != 'avail' ]] ; then
	echo "${demA}"$' Do not run modules independently!\nLaunch CHAPERONg with run_CHAPERONg-<version>!!'"${demB}"	
	exit 1
fi	

# import module with collected parameters
#. ./CHAP_colPar-4.00.sh

valid_YesNo_response=("yes" "no" '"yes"' '"no"' "y" "n")
acceptable_Yes_response=("yes" '"yes"' '"y"' "y")
acceptable_No_response=("no" '"n"' '"no"' "n")

##Defining primary functions

#defining function for topol.top modification during ligand restrain step

modTop1()
{
	catchLigTOPO="; Include ligand topology"
	ligtopSt=0
	while IFS= read -r line; do
		if [[ "$ligtopSt" == 1 ]] ; then
			ligtopSt=$(( ligtopSt + 1 ))
			# lignamFile=$(echo "$line" | awk '{print $2}' | tr -d '"' | tr -d '.itp')
			# Get name of ligand as contained in the topol.top file
			lignamFile=$(echo "$line" | awk '{print $2}' | tr -d '".' | sed 's/itp//')
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

notifyImgFail()
{
echo "${demA}"$'CHAPERONg could not generate a finished image file from the .xvg output.'\
$'\nConfirm that your xmgrace/gracebat is functional!'"${demB}"
sleep 2
}

#defining a function for make_ndx

makeNDXGroup()
{
	echo "${demA}"$' Will now make index group(s)'"${demB}"
	sleep 2

	read -p '**Provide a filename for the index group to be made: ' nameForIndex
	echo "*You entered: ${nameForIndex}"$'\n\n**Select the groups to be indexed in the coming prompts\n\n'
	sleep 2

	ndxNAME="$nameForIndex"".ndx"

	eval $gmx_exe_path make_ndx -f em.gro -o $ndxNAME

	echo -e "${demA}\033[92m Make index group ${nameForIndex}... DONE\033[m${demB}"
}


# Defining functions for different steps of the md pipeline
# Generate protein topo
s0GenTop()
{
	echo "${demA}"$' Generating protein topology...'"${demB}"
	sleep 2
	if [[ $coordinates_raw == '' ]]; then
		echo -e "${demA}\033[31m ERROR! ""No coordinates filename is given! ""Use the --input flag!""!""\033[m${demB}"
		sleep 3
		Help; sleep 3; Credit; exit 1
	fi
	if [[ "$ffUse" == "wd" || "$ffUse" == '"wd"' ]] && [[ "$wat" != "" ]]; then
		echo -e "${demA}\033[92m GROMACS will use the force-field present in the working directory!""\033[m${demB}"
		sleep 2
		echo 1 | eval $gmx_exe_path pdb2gmx -f $coordinates.pdb -o ${coordinates}_processed.gro $extr ${wmodel}
	elif [[ "$ffUse" != "" && "$ffUse" != "wd" && "$ffUse" != '"wd"' ]]; then
		echo -e "${demA}\033[92m gmx will use the specified $ffUse force-field\033[m${demB}"
		sleep 2
		eval $gmx_exe_path pdb2gmx -f $coordinates.pdb -o ${coordinates}_processed.gro $extr -ff $ffUse ${wmodel}
	elif [[ "$ffUse" != "" ]] && [[ "$ffUse" == "wd" || "$ffUse" != '"wd"' ]] && [[ "$wat" == "" ]]; then
		echo -e "${demA}\033[92m gmx will use the specified $ffUse force-field\033[m${demB}"
		sleep 2
		eval $gmx_exe_path pdb2gmx -f $coordinates.pdb -o ${coordinates}_processed.gro $extr
	elif [[ "$ffUse" == "" ]]; then
		eval $gmx_exe_path pdb2gmx -f $coordinates.pdb -o ${coordinates}_processed.gro $extr ${wmodel}
	fi
	echo -e "${demA}\033[92m Generate protein topology...DONE\033[m${demB}"
	sleep 2

	if [[ $mdType == "umbrellaSampl" ]] ; then
		echo "${demA}"$'\n Before you proceed, ensure you modify your topology to indicate the'\
		$'immobile\n reference for the pulling simulation.\n\n'\
		$'When you are done, enter "yes" below to proceed to unit cell definition...\n'

		read -p ' Do you want to proceed? (yes/no): ' procd2boxdef

		while [[ ! " ${valid_YesNo_response[@]} " =~ " ${procd2boxdef} " ]] ; do
		# while [[ "$procd2boxdef" != "yes" && "$procd2boxdef" != "no" && \
		# 	"$procd2boxdef" != '"yes"' && "$procd2boxdef" != '"no"' ]]
		# do
			echo $' Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
			read -p ' Do you want to proceed? (y/n): ' procd2boxdef
		done
		
		if [[ "${acceptable_Yes_response[@]}" =~ "$procd2boxdef" ]] ; then echo ""
		elif [[ "${acceptable_No_response[@]}" =~ "$procd2boxdef" ]] ; then exit 0
		# if [[ "$procd2boxdef" == "yes" || "$procd2boxdef" == '"yes"' ]] ; then echo ""
		# elif [[ "$procd2boxdef" == "no" || "$procd2boxdef" == '"no"' ]] ; then exit 0
		fi
	fi
}

s0CharmmSta()
{
	ligcount=0
	for ligmol2 in *.mol2; do
		if [[ "$ligmol2" == '' ]]; then
			echo $'\n'"$demA No mol2 ligand found in the current directory."
			exit 1
		else
			ligcount=$(( ligcount + 1))
			if [[ $ligcount == 2 ]] ; then 
				echo "${demA}"$'Multiple mol2 ligands found in the working directory!'\
				$'Please keep only the required ligand and try again!!'
				exit 1
			fi
			ligandmol2="$ligmol2"
			echo $'\n'"Ligand (mol2) found in the current directory: $ligandmol2"$'\n'
			sleep 3
		fi	
	done

	# Run sort_mol2_bonds.pl on ligand if found in the working directory
	if [[ -f "sort_mol2_bonds.pl" ]]; then
		echo $' sort_mol2_bonds.pl script found in the working directory.\n'\
		$' CHAPERONg will attempt to use it on '"$ligandmol2"$'\n'
		cp $ligandmol2 ${ligandmol2}.bak
		ligandmol2_basename=$(basename "$ligandmol2")
		perl sort_mol2_bonds.pl $ligandmol2 ${ligandmol2_basename}_fix.mol2
		rm $ligandmol2
		ligandmol2="${ligandmol2_basename}_fix.mol2"
	fi
	
	ffcount=0
	for ffd in charmm*.ff; do
		if [[ "$ffd" == '' ]]; then
			echo "${demA} Provide below the path/name of the charmm forcefield to be used"$'\n'
			read -p " Path (or name if in current directory) of forcefield: " forceffd
		else
			ffcount=$(( ffcount + 1))
			if [[ $ffcount == 2 ]] ; then 
				echo "${demA} Multiple forcefield folders found in the working directory! \
				Please keep only the one to be used and try again!!"
				exit 1
			fi
			forceffd="$ffd"
			echo "${demA}"$'\nForcefield found in the current directory: '"$forceffd"$'\n'
			sleep 3
		fi
		break	
	done

	charm2gmxCount=0
	for cgenffscr in cgenff_charmm2gmx*.py; do
		if [[ "$cgenffscr" == '' ]]; then
			echo "${demA} Please copy a suitable cgenff_charmm2gmx*.py script in the current \
			directory and repeat this step."$'\n'
			exit 1
		else
			charm2gmxCount=$(( charm2gmxCount + 1))
			if [[ $charm2gmxCount == 2 ]] ; then 
				echo "${demA}"$'Multiple cgenff_charmm2gmx*.py scripts found in the working directory!'\
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

	# ligcount=0
	# for ligmol2 in *.mol2; do
	# 	if [[ "$ligmol2" == '' ]]; then
	# 		echo "$demA No mol2 ligand found in the current directory."
	# 		exit 1
	# 	else
	# 		ligcount=$(( ligcount + 1))
	# 		if [[ $ligcount == 2 ]] ; then 
	# 			echo "${demA}"$'Multiple mol2 ligands found in the working directory!'\
	# 			$'Please keep only the required ligand and try again!!'
	# 			exit 1
	# 		fi
	# 		ligandmol2="$ligmol2"
	# 		echo "Ligand (mol2) found in the current directory: $ligandmol2"$'\n'
	# 		sleep 3
	# 	fi	
	# done

	strcount=0
	for strLig in *.str; do
		if [[ "$strLig" == '' ]]; then
			echo "${demA} No CHARMM stream file (.str) found in the current directory."
			exit 1
		else
		strcount=$(( strcount + 1))
			if [[ $strcount == 2 ]] ; then 
				echo "${demA}"$'Multiple stream (.str) files found in the working directory!'\
				$'Please keep only the required one and try again!!'
				exit 1
			fi
		strLigand="$strLig"
		echo "Topology stream (str) file found in the current directory: $strLigand"$'\n'
		sleep 3
		fi	
	done

	echo "Ligand name (as will be written to your topol.top file): $ligname"\
	$'\n\nThe above files will be passed to charmm2gmx...'"${demB}"
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

	echo "${demA} Will now convert CHARMM topology to GROMACS format...""${demB}"
	sleep 2

	cgenff_charmm2gmxTry1 || cgenff_charmm2gmxTry2 || cgenff_charmm2gmxTry3

	echo -e "${demA}\033[92m Convert CHARMM topology to GROMACS format...DONE\033[m${demB}"
	sleep 2

}

s0CharmmStb()
{
	LigINIcount=0
	for LigINI in *_ini.pdb; do
		if [[ "$LigINI" == '' ]]; then
			echo "${demA}"$'The current directory contains no ligand_ini.pdb file.'\
			$'\nPerhaps cgenff_charmm2gmx.py has not been run on your stream file (.str) yet.'
			exit 1
		else
		LigINIcount=$(( LigINIcount + 1))
			if [[ $LigINIcount == 2 ]] ; then 
				echo "${demA}"$'Multiple ligand_ini.pdb files found in the working directory!'\
				$'\nPlease keep only the required one and try again!!'
				exit 1
			fi
		LigINIpdb="$LigINI"
		echo "${demA} Ligand_ini.pdb file found in the current directory: $LigINIpdb""${demB}"
		sleep 3
		fi	
	done

	ligPDB=$(basename "$LigINIpdb" _ini.pdb)

	echo "${demA} Now converting ligand_ini.pdb file to ligand.gro format...""${demB}"
	sleep 2

	eval $gmx_exe_path editconf -f $LigINIpdb -o "$ligPDB"".gro"

	echo "${demA} Convert ligand_ini.pdb file to ligand.gro format...DONE"$'\n'
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

	echo "${demA} Your ""$coordinates""_processed.gro file contains $atomCountRecept atoms"$'\n'

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
	echo " Your ""$ligPDB"".gro"" file contains $atomCountLig atoms""${demB}"
	sleep 2

	TotalAtomCount=$(( atomCountRecept + atomCountLig ))

	echo "${demA} Now preparing protein-ligand complex..."$'\n'

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
	echo -e "${demA}\033[92m Prepare protein-ligand complex...DONE\033[m${demB}"
	sleep 2

	coordinates="$coordinates""-$ligname"
}

s0CharmmStd()
{
	echo "${demA} Preparing the topology for the complex..."

	Ligprmcount=0
	for ligprm in *.prm; do
		if [[ "$ligprm" == '' ]]; then
			echo "${demA}""The current directory contains no ligand parameter (.prm) file. "\
			"Perhaps cgenff_charmm2gmx.py has not been run on your stream file (.str) yet."
			exit 1
		else
			Ligprmcount=$(( Ligprmcount + 1))
			if [[ $Ligprmcount == 2 ]] ; then 
				echo "${demA}""Multiple ligand parameter (.prm) files found in the working "\
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
			echo "${demA}""The current directory contains no ligand topology (.itp) file. "\
			"Perhaps cgenff_charmm2gmx.py has not been run on your stream file (.str) yet."
			exit 1
		else
			if [[ $ligtop == 'posre.itp' || $ligtop == "posre_"*".itp" ]] ; then continue ; fi
			Ligtopcount=$(( Ligtopcount + 1))
			if [[ $Ligtopcount == 2 ]] ; then 
				echo "${demA}""Multiple ligand topology (.itp) files found in the working "\
				"directory! Please keep only the required one and try again!!"
				exit 1 ; fi
			LigitpTop="$ligtop"
			echo " Ligand topology file found in the current directory: $LigitpTop"
			sleep 1
		fi	
	done

	echo "${demA} Now modifying topol.top file to include the ligand parameters..."$'\n\n'

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

	echo -e "\n\033[92m Prepare the topology for the complex...DONE\033[m${demB}"
	sleep 2
	echo "${demA}"$'To prepare the topology for the complex, topol.top has been modified as follows:\n\n'\
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

	echo "${demA}"$'\n    **ALSO, THE LINES\n\n'\
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

	echo "${demA}"$'\n    **ALSO, THE LINE\n\n'\
	$''"$ligname                 1"\
	$'\n\n    **HAS BEEN WRITTEN AFTER THE LINES\n\n'\
	$''"[ molecules ]"$'\n; Compound        #mols\n'"$prochain"$'\n'
	sleep 2

	echo $'    **TO OBTAIN:\n\n'\
	$''"[ molecules ]"$'\n; Compound        #mols\n'"$prochain"$'\n'\
	$''"$ligname                 1"
	sleep 2

	echo "${demA}"$'\nBefore the modifications, topol.top was backed up to'\
	$'topolPreLigParmod.top.\nYou may want to check the new topol.top now to verify the'\
	$'modifications.\n\n WHEN YOU ARE DONE, you may proceed to the next step with a "yes" below...\n'

	read -p 'Do you want to proceed? (yes/no): ' procds

	while [[ ! "${valid_YesNo_response[@]}" =~ "${procds}" ]]
	# while [[ "$procds" != "yes" && "$procds" != "no" && "$procds" != '"yes"' && "$procds" != '"no"' ]]
	do
		echo $'Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
		read -p 'Do you want to proceed? (yes/no): ' procds
	done
	if [[ " ${acceptable_Yes_response[@]} " =~ "$procds" ]] ; then echo ""
	elif [[ " ${acceptable_No_response[@]} " =~ "$procds" ]] ; then exit 0
	# if [[ "$procds" == "yes" || "$procds" == '"yes"' ]] ; then echo ""
	# elif [[ "$procds" == "no" || "$procds" == '"no"' ]] ; then exit 0
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

	echo "${demA} Your ""$coordinates""_processed.gro file contains $atomCountRecept atoms"$'\n'

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

	echo " Your ""$ligname"".gro"" file contains $atomCountLig atoms""${demB}"
	sleep 2

	TotalAtomCount=$(( atomCountRecept + atomCountLig ))

	echo "${demA} Now preparing protein-ligand complex..."$'\n'

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
	echo -e "\033[92m Prepare protein-ligand complex...DONE\033[m${demB}"
	sleep 2

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
	
	echo "${demA} Preparing the topology for the complex..."$'\n'

	Ligtopcount=0
	for ligtop in *.itp; do
		if [[ "$ligtop" == '' ]]; then
			echo "${demA}""The current directory contains no ligand topology (.itp) file. Perhaps you have not "\
			"copied the prodrg folder into the working directory, or the folder has not been named accordingly!"
			exit 1
		else
			if [[ "$ligtop" == 'posre.itp' || "$ligtop" == "$ligitp" || "$ligtop" == 'posre'*'.itp' ]]
			then continue
			fi
			Ligtopcount=$(( Ligtopcount + 1))
			if [[ $Ligtopcount == 2 ]] ; then 
				echo "${demA}"$'Multiple ligand topology (.itp) files found in the working '\
				$'directory! Please keep only the required one and try again!!'
				exit 1 
			fi
			if [[ "$ligtop" == "$ligname"".itp" ]] ; then LigitpTop="$ligtop"
				echo " Ligand topology file found in the current directory: $LigitpTop"
			fi
			sleep 1
		fi	
	done

	echo "${demA} Now modifying topol.top file to include the ligand parameters and topology..."$'\n\n'

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
	echo "${demA}"$'To prepare the topology for the complex, '\
	$'topol.top has been modified as follows:\n\n'
	echo "${demA}"$'\n    **THE LINES\n\n'\
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

	echo "${demA}"$'\n    **ALSO, THE LINE\n\n'\
	$''"$ligname                 1"\
	$'\n\n    **HAS BEEN WRITTEN AFTER THE LINES\n\n'\
	$''"[ molecules ]"$'\n; Compound        #mols\n'"$prochain"$'\n'
	sleep 2

	echo $'    **TO OBTAIN:\n\n'\
	$''"[ molecules ]"$'\n; Compound        #mols\n'"$prochain"$'\n'\
	$''"$ligname                 1"
	sleep 2

	echo "${demA}"$'\nThe topol.top file before the modifications has been backed up to '\
	$'topolPreLigParmod.top.\nYou may want to check the new topol.top now to verify the '\
	$'modifications.\nWHEN YOU ARE DONE, you may proceed to the next step with a "yes" below...\n\n'

	read -p 'Do you want to proceed? (yes/no): ' procds

	while [[ ! "${valid_YesNo_response[@]}" =~ "${procds}" ]]
	# while [[ "$procds" != "yes" && "$procds" != "no" && "$procds" != '"yes"' && "$procds" != '"no"' ]]
	do
		echo $'Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
		read -p 'Do you want to proceed? (yes/no): ' procds
	done
	if [[ "${acceptable_Yes_response[@]}" =~ "$procds" ]] ; then echo ""
	elif [[ "${acceptable_No_response[@]}" =~ "$procds" ]] ; then exit 0
	# if [[ "$procds" == "yes" || "$procds" == '"yes"' ]] ; then echo ""
	# elif [[ "$procds" == "no" || "$procds" == '"no"' ]] ; then exit 0
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

	echo "${demA} Your ""$coordinates""_processed.gro file contains $atomCountRecept atoms"$'\n'

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

	echo " Your ""$ligname"".gro"" file contains $atomCountLig atoms""${demB}"
	sleep 2

	TotalAtomCount=$(( atomCountRecept + atomCountLig ))

	echo "${demA}"" Now preparing protein-ligand complex..."$'\n'

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
	echo "\033[92m Prepare protein-ligand complex...DONE\033[m${demB}"
	sleep 2

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
	
	echo "${demA}"" Extracting ligand parameters from the Acpype-derived topology..."$'\n'
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
	echo " Extract ligand parameters from the Acpype-derived topology...DONE""${demA}"

	echo "${demA}"" Preparing the topology for the complex..."$'\n'

	Ligprmcount=0
	for ligprm in *.prm; do
		if [[ "$ligprm" == '' ]]; then
			echo "${demA}""The current directory contains no ligand parameter (.prm) file. "\
			"Perhaps CHAPERONg could not process your Acpype-derived topology ($ligname"_GMX.itp")."
			exit 1
		else
			Ligprmcount=$(( Ligprmcount + 1))
			if [[ $Ligprmcount == 2 ]] ; then 
				echo "${demA}""Multiple ligand parameter (.prm) files found in the working "\
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
			echo "${demA}""The current directory contains no ligand topology (.itp) file. Perhaps you have not "\
			"copied the prodrg folder into the working directory, or the folder has not been named accordingly!"
			exit 1
		else
			if [[ "$ligtop" == 'posre.itp' || "$ligtop" == "$ligitp" || "$ligtop" == 'posre'*'.itp' ]]
			then continue
			fi
			Ligtopcount=$(( Ligtopcount + 1))
			if [[ $Ligtopcount == 2 ]] ; then 
				echo "${demA}"$'Multiple ligand topology (.itp) files found in the working '\
				$'directory! Please keep only the required one and try again!!'
				exit 1 
			fi
			if [[ "$ligtop" == "$ligname"".itp" ]] ; then LigitpTop="$ligtop"
				echo " Ligand topology file found in the current directory: $LigitpTop"
			fi
			sleep 1
		fi	
	done

	echo "${demA}"" Now modifying topol.top file to include the ligand parameters and topology..."$'\n\n'

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
	echo "${demA}"$'To prepare the topology for the complex, topol.top has been modified as follows:\n\n'
	echo "${demA}"$'\n    **THE LINES\n\n'\
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

	echo "${demA}"$'\n    **ALSO, THE LINES\n\n'\
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

	echo "${demA}"$'\n    **ALSO, THE LINE\n\n'"$ligname                 1"\
	$'\n\n    **HAS BEEN WRITTEN AFTER THE LINES\n\n'"[ molecules ]"\
	$'\n; Compound        #mols\n'"$prochain"$'\n'
	sleep 2

	echo $'    **TO OBTAIN:\n\n'"[ molecules ]"\
	$'\n; Compound        #mols\n'"$prochain"\
	$'\n'"$ligname                 1"
	sleep 2

	echo "${demA}"$'\nThe topol.top file before the modifications has been backed up to '\
	$'topolPreLigParmod.top.\nYou may want to check the new topol.top now to verify the '\
	$'modifications.\nWHEN YOU ARE DONE, you may proceed to the next step with a "yes" below...\n\n'

	read -p 'Do you want to proceed? (yes/no): ' procds

	while [[ ! "${valid_YesNo_response[@]}" =~ "${procds}" ]]
	# while [[ "$procds" != "yes" && "$procds" != "no" && "$procds" != '"yes"' && "$procds" != '"no"' ]]
	do
		echo $'Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
		read -p 'Do you want to proceed? (yes/no): ' procds
	done
	if [[ "${acceptable_Yes_response[@]}" =~ "$procds" ]] ; then echo ""
	elif [[ "${acceptable_No_response[@]}" =~ "$procds" ]] ; then exit 0
	# if [[ "$procds" == "yes" || "$procds" == '"yes"' ]] ; then echo ""
	# elif [[ "$procds" == "no" || "$procds" == '"no"' ]] ; then exit 0
	fi

}

#Process topology for LigParGen
s0ligpargenSta()
{
	#find topology file
	for ligLPGitp in ./ligpargen/* ; do
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
	for liggro in ./ligpargen/* ; do
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

	echo "${demA}"" Your ""$coordinates""_processed.gro file contains $atomCountRecept atoms"$'\n'

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

	echo " Your ""$ligname"".gro"" file contains $atomCountLig atoms""${demB}"
	sleep 2

	TotalAtomCount=$(( atomCountRecept + atomCountLig ))

	echo "${demA}"" Now preparing protein-ligand complex..."$'\n'

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
	echo "\033[92m Prepare protein-ligand complex...DONE\033[m${demB}"
	sleep 2

	coordinates="$coordinates""-$ligname"
}

#Prepare topology for LigParGen protein-lig complex
s0ligpargenStb()
{
	#find topology file
	for ligitp in ./ligpargen/* ; do
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
	
	echo "${demA}"" Extracting ligand parameters from the LigParGen-derived topology..."$'\n'
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
	echo " Extract ligand parameters from the LigParGen-derived topology...DONE""${demA}"

	echo "${demA}"" Preparing the topology for the complex..."$'\n'

	Ligprmcount=0
	for ligprm in *.prm; do
		if [[ "$ligprm" == '' ]]; then
			echo "${demA}""The current directory contains no ligand parameter (.prm) file. "\
			"Perhaps CHAPERONg could not process your Acpype-derived topology ($ligname"_GMX.itp")."
			exit 1
		else
			Ligprmcount=$(( Ligprmcount + 1))
			if [[ $Ligprmcount == 2 ]] ; then 
				echo "${demA}""Multiple ligand parameter (.prm) files found in the working "\
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
			echo "${demA}""The current directory contains no ligand topology (.itp) file. Perhaps you have not "\
			"copied the prodrg folder into the working directory, or the folder has not been named accordingly!"
			exit 1
		else
			if [[ "$ligtop" == 'posre.itp' || "$ligtop" == "$ligitp" || "$ligtop" == 'posre'*'.itp' ]]
			then continue
			fi
			Ligtopcount=$(( Ligtopcount + 1))
			if [[ $Ligtopcount == 2 ]] ; then 
				echo "${demA}"'Multiple ligand topology (.itp) files found in the working '\
				$'directory! Please keep only the required one and try again!!'
				exit 1 ; fi
			if [[ "$ligtop" == "$ligname"".itp" ]] ; then LigitpTop="$ligtop"
				echo " Ligand topology file found in the current directory: $LigitpTop"
			fi
			sleep 1
		fi	
	done

	echo "${demA}"" Now modifying topol.top file to include the ligand parameters and topology..."$'\n\n'

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
	echo "${demA}"$'To prepare the topology for the complex, topol.top has been modified as follows:\n\n'
	echo "${demA}"$'\n    **THE LINES\n\n'\
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

	echo "${demA}"$'\n    **ALSO, THE LINES\n\n'\
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

	echo "${demA}"$'\n    **ALSO, THE LINE\n\n'"$ligname                 1"\
	$'\n\n    **HAS BEEN WRITTEN AFTER THE LINES\n\n'"[ molecules ]"\
	$'\n; Compound        #mols\n'"$prochain"$'\n'
	sleep 2

	echo $'    **TO OBTAIN:\n\n'"[ molecules ]"\
	$'\n; Compound        #mols\n'"$prochain"\
	$'\n'"$ligname                 1"
	sleep 2

	echo "${demA}"$'\nThe topol.top file before the modifications has been backed up to '\
	$'topolPreLigParmod.top.\nYou may want to check the new topol.top now to verify the '\
	$'modifications.\nWHEN YOU ARE DONE, you may proceed to the next step with a "yes" below...\n\n'

	read -p 'Do you want to proceed? (yes/no): ' procds

	while [[ ! "${valid_YesNo_response[@]}" =~ "${procds}" ]]
	# while [[ "$procds" != "yes" && "$procds" != "no" && "$procds" != '"yes"' && "$procds" != '"no"' ]]
	do
		echo $'Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
		read -p 'Do you want to proceed? (y/n): ' procds
	done
	if [[ "${acceptable_Yes_response[@]}" =~ "$procds" ]] ; then echo ""
	elif [[ "${acceptable_No_response[@]}" =~ "$procds" ]] ; then exit 0
	# if [[ "$procds" == "yes" || "$procds" == '"yes"' ]] ; then echo ""
	# elif [[ "$procds" == "no" || "$procds" == '"no"' ]] ; then exit 0
	fi

}	

s0Charmm()
{
	echo "${demA}"$'\nChoose the stage to start/resume ligand topology preparation from\n'
cat << iniGenLigTop
 Stage  Step to begin/resume ligand topology preparation 
   a    Convert CHARMM stream file (.str) to GROMACS format (cgenff_charmm2gmx)
   b    Convert ligand_ini.pdb to ligand.gro (editconf)
   c    Prepare the protein-ligand complex
   d    Prepare the topology for the complex

*ENTER A RESPONSE BELOW USING THE APPROPRIATE OPTION

iniGenLigTop

	read -p ' Initiate ligand preparation at stage: ' stepLigPrep

	valid_entries=( "a" "b" "c" "d" )
	while [[ ! "${valid_entries[@]}" =~ "${stepLigPrep}" ]] ; do
	# while [[ "$stepLigPrep" != 'a' && "$stepLigPrep" != 'b' && \
	# 	"$stepLigPrep" != 'c' && "$stepLigPrep" != 'd' ]] ; do
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
	echo "${demA}"$'\nChoose the stage to start/resume ligand topology preparation from\n'
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

	if [[ $mdType == "umbrellaSampl" ]] ; then
		echo "${demA}"$'\n Before you proceed, ensure you modify your topology to indicate the'\
		$'immobile\n reference for the pulling simulation.\n\n'\
		$'When you are done, enter "yes" below to proceed to unit cell definition...\n'

		read -p 'Do you want to proceed? (yes/no): ' procd2boxdef

		while [[ ! "${valid_YesNo_response[@]}" =~ "${procd2boxdef}" ]] ; do
		# while [[ "$procd2boxdef" != "yes" && "$procd2boxdef" != "no" ]] && \
		# 	[[ "$procd2boxdef" != '"yes"' && "$procd2boxdef" != '"no"' ]] ; do
			echo $'Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
			read -p 'Do you want to proceed? (y/n): ' procd2boxdef
		done
	
		if [[ "${acceptable_Yes_response[@]}" =~ "$procd2boxdef" ]] ; then echo ""
		elif [[ "${acceptable_No_response[@]}" =~ "$procd2boxdef" ]] ; then exit 0

		# if [[ "$procd2boxdef" == "yes" || "$procd2boxdef" == '"yes"' ]]; then echo ""
		# elif [[ "$procd2boxdef" == "no" || "$procd2boxdef" == '"no"' ]]; then exit 0
		fi
	fi	
}

s0GenLigTop()
{
	echo "${demA}"$'\nWhich platform have you used for your ligand parameterization?\n'
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
	echo "${demA}"" Defining simulation box...""${demB}"
	sleep 2
	if [[ $mdType == "regularMD" ]] && [[ $sysType == "protein_only" || $sysType == "protein_dna" ]] ; then
		eval $gmx_exe_path editconf -f ${coordinates}_processed.gro -o ${coordinates}_newbox.gro -c -d ${edgeDist} -bt ${btype}
	elif [[ $mdType == "regularMD" && $sysType == "protein_lig" ]]; then
		eval $gmx_exe_path editconf -f ${coordinates}.gro -o ${coordinates}_newbox.gro -c -d ${edgeDist} -bt ${btype}
	elif [[ $mdType == "umbrellaSampl" && $sysType == "protein_only" ]] ; then
		eval $gmx_exe_path editconf -f ${coordinates}_processed.gro -o ${coordinates}_newbox.gro -c -d ${edgeDist} -bt ${btype}
	elif [[ $mdType == "umbrellaSampl" && $sysType == "protein_lig" ]] ; then
		eval $gmx_exe_path editconf -f ${coordinates}.gro -o ${coordinates}_newbox.gro -c -d ${edgeDist} -bt ${btype}
	fi

	if [[ $mdType == "umbrellaSampl" ]] ; then
		sleep 2
		echo "${demA}"$'\n'" A placeholder unit cell of the type "$'"'"${btype}"$'"'" has been generated. Inspect the"\
		$'\n unit cell (you may check '"${coordinates}_newbox.gro"') and/or provide'\
		$'\n dimensions for a new box and the placement of the COM of your protein.\n'

		echo "Provide dimensions or proceed to solvation?"
		echo $'\n    1.  Provide dimensions\n    2.  Proceed\n'

		read -p 'ENTER A RESPONSE HERE (1 or 2): ' respDim

		while [[ $respDim == 1 ]] ; do
			echo ""
			read -p 'Center dimensions (x y z):' centerVector
			echo ""
			read -p 'Box dimensions (x y z):' boxVector
			echo ""
			if [[ $sysType == "protein_only" ]] ; then
				eval $gmx_exe_path editconf -f ${coordinates}_processed.gro -o ${coordinates}_newbox.gro -center $centerVector -box $boxVector
			elif [[ $sysType == "protein_lig" ]] ; then
				eval $gmx_exe_path editconf -f ${coordinates}.gro -o ${coordinates}_newbox.gro -center $centerVector -box $boxVector
			fi
			echo "Provide dimensions or proceed to solvation?"
			echo $'\n    1.  Provide dimensions\n    2.  Proceed\n'

			read -p 'ENTER A RESPONSE HERE (1 or 2): ' respDim
		done
	fi
	echo -e "${demA}\033[92m Define box...DONE\033[m${demB}"
	sleep 2
}

s2Solvat()
{
	echo "${demA}"" Solvating the system...""${demB}"
	sleep 2
	eval $gmx_exe_path solvate -cp ${coordinates}_newbox.gro -cs spc216.gro -o ${coordinates}_solv.gro -p topol.top
	echo -e "${demA}\033[92m Solvate system...DONE\033[m${demB}"
	sleep 2
}

s3AddIons1()
{
	echo "${demA}"" Adding ions...""${demB}"
	sleep 2
	#1. Assemble .tpr file with grompp, using ions.mdp
	eval $gmx_exe_path grompp -f ions.mdp -c ${coordinates}_solv.gro -p topol.top -o ions.tpr -maxwarn $WarnMax
}

s4AddIons2()
{
	#2. Pass .tpr file to genion to add ions
	echo 'SOL' | eval $gmx_exe_path genion -s ions.tpr -o ${coordinates}_solv_ions.gro -p topol.top -neutral ${pnam_nnam}
	echo -e "${demA}\033[92m Add ions...DONE\033[m${demB}"
	sleep 2
}

s5EnMin1()
{
	echo "${demA}"" Now running energy minimization...""${demB}"
	sleep 2
	Enmdp=0; Minmdp=0
	#1. Assemble the binary input with grompp
	for param in *m.mdp; do
		if [[ "${param}" == "minim.mdp" ]] ; then Minmdp=1
		elif [[ "${param}" == "em.mdp" ]] ; then Enmdp=1
		fi
	done
	if [[ "$Minmdp" == 1 ]] && [[ "$Enmdp" == 1 ]]; then
		EnMmdp="minim.mdp"
		echo "${demA}"$' Parameter files "minim.mdp" and "em.mdp" found.\n *Choosing "minim.mdp" for energy minimization.'"${demB}"
		sleep 2
	elif [[ "$Minmdp" == 1 ]] && [[ "$Enmdp" == 0 ]]; then
		EnMmdp="minim.mdp"
		echo "${demA}"" Parameter file $EnMmdp found."$'\n'" *Using $EnMmdp for energy minimization.""${demB}"
		sleep 2
	elif [[ "$Minmdp" == 0 ]] && [[ "$Enmdp" == 1 ]]; then
		EnMmdp="em.mdp"
		echo "${demA}"" Parameter file $EnMmdp found."$'\n'" *Using $EnMmdp for energy minimization.""${demB}"
		sleep 2
	
	fi
	sleep 1	
	eval $gmx_exe_path grompp -f $EnMmdp -c ${coordinates}_solv_ions.gro -p topol.top -o em.tpr -maxwarn $WarnMax

}

s6EnMin2()
{
	#2. Invoke mdrun to carry out the EM
	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -deffnm em
	#tail -n 8 em.log
	echo -e "${demA}\033[92m Run energy minimization... DONE\033[m${demB}"
	sleep 2

	echo "${demA}"$' Now calculating post-EM thermodynamic parameters...\n\n'
	sleep 2

	echo "Potential" | eval $gmx_exe_path energy -f em.edr -o postEM_Potential.xvg
	gracebat postEM_Potential.xvg -hdevice PNG -autoscale xy -printfile postEM_Potential.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "${demA}"$' Calculate Potential energy...DONE\n\n' ; sleep 2

	currentEMthermodyndir="$(pwd)""/postEM_thermodynamics"
	nEMtherm=1
	bkupEMtherm="$(pwd)""/#postEM_thermodynamics"".""backup.""$nEMtherm"
	base_bkupEMtherm=$(basename "$bkupEMtherm")
	if [[ -d "$currentEMthermodyndir" ]]; then
		echo $'\n'"$currentEMthermodyndir"$' folder exists,\n'"backing it up as $base_bkupEMtherm"
		sleep 1
		while [[ -d "$bkupEMtherm" ]]; do
			nEMtherm=$(( nEMtherm + 1 ))
			bkupEMtherm="$(pwd)""/#postEM_thermodynamics"".""backup.""$nEMtherm"
			base_bkupEMtherm=$(basename "$bkupEMtherm")
		done
		mv "$currentEMthermodyndir" "$bkupEMtherm" && mkdir ./postEM_thermodynamics || true
		echo $'\n'"Backing up the last postEM_thermodynamics folder and its contents as $base_bkupEMtherm"
		sleep 1
	elif [[ ! -d "$currentEMthermodyndir" ]]; then mkdir postEM_thermodynamics
	fi
	mv postEM_Potential.xvg postEM_Potential.png ./postEM_thermodynamics || true

	echo -e "${demA}\033[92m Calculate post-EM thermodynamics parameters...DONE\033[m${demB}"
	sleep 2

	if [[ $sysType == "protein_dna" ]] ; then
		echo "${demA}""CHAPERONg will run a check to detect the group numbers of protein and DNA...""${demB}"
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
		echo "${demA}"" Checking completed...""${demB}"
		sleep 2
		echo "${demA}""CHAPERONg will now make a Protein-DNA index group...""${demB}"
		sleep 2

eval $gmx_exe_path make_ndx -f em.gro -o index.ndx << ProvProDNAMakNdx
$proNDXno | $dnaNDXno
q
ProvProDNAMakNdx

	echo -e "${demA}\033[92m Make a Protein-DNA index group...DONE\033[m${demB}"
	sleep 2

	if [[ $sysType == "protein_dna" && $automode == "full" ]]; then
		echo "${demA}""CHAPERONg will now check your mdp files and make corrections on the tc-grps if needed...""${demB}" 
		sleep 2
		for mdparafile in nvt.mdp npt.mdp md.mdp; do
			nbkUP=1
			mdpbak="$mdparafile"".backup."
			prevmdpbak="$mdpbak""$nbkUP"
			if [[ -f "$prevmdpbak" ]]; then
				echo "${demA}""$prevmdpbak"" exists!"
				nbkUP=$(( nbkUP + 1 ))
				currentmdpbak="$mdpbak""$nbkUP"
				while [[ -f "$currentmdpbak" ]]; do
					nbkUP=$(( nbkUP + 1 ))
					currentmdpbak="$mdpbak""$nbkUP"
				done
				echo $'\nCHAPERONg will back that up as'\
				"$currentmdpbak and the $mdparafile used in the last run as $mdparafile"".backup.1""${demB}"
				sleep 1
				mv "$prevmdpbak" "$currentmdpbak"
				cp "$mdparafile" "$mdparafile"".backup.1"
			elif [[ ! -f "$prevmdpbak" ]]; then
				echo "${demA}"" Backing up $mdparafile used in the last run as $mdparafile"".backup.1""${demB}"
				sleep 1
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
		echo "${demA}""npt.mdp, nvt.mdp and md.mdp successfully checked and/or modified accorgingly."
		sleep 2
		echo $'\nYou may want to check these files to ensure they have been modified to your satisfaction.\n'
		$'If not, you may re-modify them accordingly or restore the backups made by CHAPERONg.\n'
		
		echo $'When you are done, enter "yes" below to proceed to equilibration...'"${demB}"

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

	elif [[ $sysType == "protein_lig" ]] ; then

		echo "${demA}"$' Confirm if you want to make custom index group(s) before proceeding.'\
		$'\n\n You should respond with a "yes" or a "no" below.'\
		$'\n\n *If you are unsure, enter "no" and CHAPERONg will make appropriate index for'\
		$'\n you while automatically updating topol.top accordingly.'"${demB}"
		sleep 3
		customNDXask=''
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
		echo "${demA}"" Restraining the ligand...""${demB}""${demA}"$'**"In the coming steps, '\
		$'CHAPERONg will create an index group for the non-hydrogen atoms of '"$ligname ${demB}"
		sleep 2

eval $gmx_exe_path make_ndx -f "$lignamFile".gro -o index_"$lignamFile".ndx << ProvMakNdx
0 & ! a H*
q
ProvMakNdx
	fi
	echo "${demA}"" Executing the genrestr module""${demB}""${demA}"$'**CHAPERONg will select '\
	$'the newly created index group i.e. group 3 in the index_'"$lignamFile.ndx file""${demB}"
	sleep 2
	
	if [[ "$customNDXask" == "no" || "$customNDXask" == $'"no"' || "$customNDXask" == '' ]] ; then
		echo "System_&_!H*" | eval $gmx_exe_path genrestr -f $lignamFile.gro -n index_$lignamFile.ndx -o posre_$lignamFile.itp -fc 1000 1000 1000

		if [[ -f "$ligname""_GMX.itp" ]] && [[ -d acpype ]] ; then modTop2forAcpype
		else
			modTop2
		fi
	
		echo "${demA}"$' To restrain the protein and ligand simultaneously,\n'\
		$'the following lines have been written to your topol.top file:'"${demB}"\
		$'\n; Ligand position restraints\n#ifdef POSRES'\
		$'\n#include '$'"'$'posre_'"$lignamFile"$'.itp'$'"'\
		$'\n#endif\n'

		sleep 2
 
		echo "${demA}"$' The above were written'"${demB}""${demA}"$'  **BETWEEN the lines:'"${demB}"\
		$'\n; Include ligand topology\n'\
		$'#include '$'"'"$lignamFile"$'.itp'$'"'
		sleep 2

		if [[ -f "$ligname""_GMX.itp" ]] && [[ -d acpype ]]; then
			echo $'\n'"${demA}"$'  **AND:'"${demB}"$'\n[ moleculetype ]\n; Name            nrexcl'
			sleep 2
			echo "${demA}"$'  **TO OBTAIN:'"${demB}"$'; Include ligand topology'
			echo "#include "$'"'$lignamFile.itp$'"'
			sleep 2
			echo $'\n; Ligand position restraints\n#ifdef POSRES\n'\
			$'#include '$'"'$'posre_'"$lignamFile"$'.itp'$'"'\
			$'\n#endif\n\n'\
			$'\n[ moleculetype ]\n; Name            nrexcl'
		else
			echo $'\n'"${demA}"$'  **AND:'"${demB}"$'\n; Include water topology\n#include '"$ff_watMod_dir"
			sleep 2
			echo "${demA}"$'  **TO OBTAIN:'"${demB}"$'; Include Position restraint file\n#ifdef POSRES\n#include '\
			$'"posre.itp"\n#endif\n\n; Include ligand topology'
			echo "#include "$'"'$lignamFile.itp$'"'
			sleep 2
			echo $'\n; Ligand position restraints\n#ifdef POSRES\n'\
			$'#include '$'"'$'posre_'"$lignamFile"$'.itp'$'"'\
			$'\n#endif\n\n'\
			$'; Include water topology\n#include '"$ff_watMod_dir"
		fi

		sleep 2
		echo "${demA}"$'\n The "topol.top" file before modification was'\
		$'backed up as "topolPreResmod.top"\n\n'
		sleep 2

		echo $' You may want to check your "topol.top" file now to verify the modifications.\n'\
		$'WHEN YOU ARE DONE, you can proceed to temperature coupling with a "yes" below\n'

		if [[ $mdType == "regularMD" ]] ; then
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
		echo "${demA}"" Coupling the protein with $ligname...""${demB}"
		sleep 2
	#Make index for temperature coupling
eval $gmx_exe_path make_ndx -f em.gro -o index.ndx << TempCoupIN1
"Protein" | 13
q
TempCoupIN1

	elif [[ "$pn" == "NA" && "$nn" == "CL" ]]
	then echo "${demA}"" Coupling the protein with $ligname...""${demB}"
	sleep 2
#Make index for temperature coupling
eval $gmx_exe_path make_ndx -f em.gro -o index.ndx << TempCoupIN2
"Protein" | 13
q
TempCoupIN2

	#in the case of the charmm-jul-2021 ff
	elif [[ "$pnam_nnam" != '' && "$pn" != "NA" && "$nn" != "CL" ]]
	then echo "${demA}"" Coupling the protein with $ligname..."$'\n Will '\
	$'also couple water with '"$pn and\or $nn""${demB}"
	sleep 3
#Make index for temperature coupling
eval $gmx_exe_path make_ndx -f em.gro -o index.ndx << TempCoupIN3
"Protein" | 13
"Water" | "$pn" | "$nn"
q
TempCoupIN3

		if [[ $sysType == "protein_lig" && $automode == "full" ]] ; then
			echo $'\n'"${demA}"$' Running a check to detect the type(s) of'\
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
				echo $'  Checking completed...\n  The ions'\
				$''"$ionAdded1 and $ionAdded2 were detected in your topol.top file..."$'\n\n'
				sleep 2
			elif [[ "$ionAdded1" != '' || "$ionAdded2" != '' ]] ; then
				echo $'  Checking completed...\n  The ion'\
				$''"$ionAdded1""$ionAdded2 was detected in your topol.top file"$'\n'
				sleep 2
			fi	
			echo $'\n Checking your mdp files to make corrections to the tc-grps if needed...\n'
			sleep 2
			if [[ $mdType == "regularMD" ]] ; then
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
						echo $'\n CHAPERONg will back that up as '"$currentmdpbak"$' and'\
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
				echo "${demB}""${demA}""npt.mdp, nvt.mdp and md.mdp successfully checked and/or modified accorgingly."
				sleep 2

			elif [[ ! -f nvt.mdp ]] && [[ $mdType == "umbrellaSampl" ]] ; then
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
						echo $'\n CHAPERONg will back that up as '"$currentmdpbak"$' and'\
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
				echo "${demB}""${demA}""The files npt.mdp, md_pull.mdp, npt_umbrella.mdp, and md_umbrella.mdp have"\
				$'\n been checked and/or modified accorgingly.'
				sleep 2
			
			elif [[ -f nvt.mdp ]] && [[ $mdType == "umbrellaSampl" ]] ; then
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
						echo $'\n CHAPERONg will back that up as '"$currentmdpbak"$' and'\
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
				echo "${demA}""The files nvt.mdp, npt.mdp, md_pull.mdp, npt_umbrella.mdp, and md_umbrella.mdp"\
				$'\nhave been checked and/or modified accorgingly.'
				sleep 2
			fi
			echo $'\nYou may want to verify the modifications in these files. If unsatisfied, you'\
			$'\nmay re-modify them accordingly or restore the backups made by CHAPERONg.\n'\
				
			echo $'\n *When you are done, enter "yes" below to proceed to equilibration...'"${demB}"

		elif [[ $sysType == "protein_lig" && $automode == "semi" ]] ; then
			echo "${demA}"$'Because you are not running CHAPERONg in the "auto mode", you may\n'\
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
	if [[ $mdType == "umbrellaSampl" ]] ; then
		echo "${demA}"$'\n'
		read -p 'Do you want to run an optional NVT equilibration? (yes/no): ' runNVTeq
		
		while [[ "$runNVTeq" != "yes" && "$runNVTeq" != "no" && \
			"$runNVTeq" != '"yes"' && "$runNVTeq" != '"no"' ]]
		do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p 'Do you want to run an optional NVT equilibration? (yes/no): ' runNVTeq
		done
	elif [[ $mdType == "regularMD" ]] ; then runNVTeq="yes"
	fi

	if [[ $sysType == "protein_only" && $runNVTeq == "yes" ]] ; then
		echo "${demA}"" Now running NVT Equilibration...""${demB}"
		sleep 2
		eval $gmx_exe_path grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn $WarnMax
	elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $runNVTeq == "yes" ]] ; then
		echo "${demA}"" Now running NVT Equilibration...""${demB}"
		sleep 2
		eval $gmx_exe_path grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn $WarnMax
	elif [[ $runNVTeq == "no" || $runNVTeq != '"no"' ]] ; then
		echo "${demA}"$' Skipping NVT equilibration, proceeding to NPT equilibration'"${demB}"
		sleep 2
	fi
}

s8NVTeq2()
{
	if [[ $runNVTeq == "yes" ]] ; then
		eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm nvt
		echo -e "${demA}\033[92m NVT Equilibration... DONE\033[m${demB}"
		sleep 2
		echo "${demA}"$' Now calculating post-NVT thermodynamic parameters...\n\n'
		sleep 2
		echo "${demA}"$' Analyzing the Temperature progression of the System...\n\n'
		echo "Temperature" | eval $gmx_exe_path energy -f nvt.edr -o postNVT_Temperature.xvg
		gracebat postNVT_Temperature.xvg -hdevice PNG -autoscale xy -printfile postNVT_Temperature.png \
		-fixed 7500 4000 -legend load || notifyImgFail
		
		currentEMthermodyndir="$(pwd)""/postEM_thermodynamics"
		if [[ -d "$currentEMthermodyndir" ]]; then echo ''
		elif [[ ! -d "$currentEMthermodyndir" ]]; then mkdir postEM_thermodynamics
		fi
		mv postNVT_Temperature.xvg postNVT_Temperature.png ./postEM_thermodynamics || true
		echo "${demA}"$' Calculate Temperature progression...DONE\n\n' ; sleep 2	
	fi
}

s9NPTeq1()
{
	echo "${demA}"" Now running NPT Equilibration...""${demB}"
	sleep 2
	if [[ $sysType == "protein_only" && $runNVTeq == "yes" ]] ; then
		eval $gmx_exe_path grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn $WarnMax
	elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $runNVTeq == "yes" ]] ; then
		eval $gmx_exe_path grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn $WarnMax
	#if NVT equilibration was skipped
	elif [[ $sysType == "protein_only" && $runNVTeq == "no" ]] ; then
		eval $gmx_exe_path grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -o npt.tpr -maxwarn $WarnMax
	elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] && [[ $runNVTeq == "no" ]] ; then
		eval $gmx_exe_path grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn $WarnMax
	# if resuming simulation from NPT equil
	else
		if [[ $sysType == "protein_only" ]] ; then
			eval $gmx_exe_path grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn $WarnMax || \
			eval $gmx_exe_path grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -o npt.tpr -maxwarn $WarnMax
		elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] ; then
			eval $gmx_exe_path grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn $WarnMax || \
			eval $gmx_exe_path grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn $WarnMax
		fi
	fi
}

s10NPTeq2()
{
	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm npt
	echo -e "${demA}\033[92m NPT Equilibration... DONE\033[m${demB}"
	sleep 2
	echo "${demA}"$' Now calculating post-NPT thermodynamic parameters...\n\n'
	sleep 2

	echo "Pressure" | eval $gmx_exe_path energy -f npt.edr -o postNPT_Pressure.xvg
	gracebat postNPT_Pressure.xvg -hdevice PNG -autoscale xy -printfile postNPT_Pressure.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "${demA}"$' Calculate Pressure progression...DONE\n\n' ; sleep 2

	echo "Density" | eval $gmx_exe_path energy -f npt.edr -o postNPT_Density.xvg
	gracebat postNPT_Density.xvg -hdevice PNG -autoscale xy -printfile postNPT_Density.png \
	-fixed 7500 4000 -legend load || notifyImgFail
	echo "${demA}"$' Calculate Density progression...DONE\n\n' ; sleep 2
	
	currentEMthermodyndir="$(pwd)""/postEM_thermodynamics"
	
	if [[ ! -d "$currentEMthermodyndir" ]]; then mkdir postEM_thermodynamics
	elif [[ -d "$currentEMthermodyndir" ]]; then echo
	fi
	mv postNPT_Density.xvg postNPT_Density.png ./postEM_thermodynamics || true
	mv postNPT_Pressure.xvg postNPT_Pressure.png ./postEM_thermodynamics || true


	if [[ $sysType == "protein_only" && $mdType == "umbrellaSampl" ]] ; then
		read -p 'Do you need to make custom index for pulling groups? (yes/no): ' cstmndx
		
		while [[ "$cstmndx" != "yes" && "$cstmndx" != "no" && \
			"$cstmndx" != '"yes"' && "$cstmndx" != '"no"' ]]
		do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p 'Do you need to make custom index for pulling groups? (yes/no): ' cstmndx
		done
		if [[ $cstmndx == "yes" ]] ; then
			echo "${demA}"" Will now make an index for the pulling groups...""${demB}"
			sleep 2
			eval $gmx_exe_path make_ndx -f npt.gro
			echo -e "${demA}\033[92m Make an index for the pulling groups...DONE\033[m${demB}"
			sleep 2
		fi
	fi
}

s11RelPosRe()
{
	echo "${demA}"" Releasing position restraints...""${demB}"
	sleep 2
	if [[ $sysType == "protein_only" ]]; then
		eval $gmx_exe_path grompp -v -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o ${filenm}.tpr -maxwarn $WarnMax
	elif [[ $sysType == "protein_lig" || $sysType == "protein_dna" ]] ; then
		eval $gmx_exe_path grompp -v -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ${filenm}.tpr -maxwarn $WarnMax
	fi
	echo -e "${demA}\033[92m Release position restraints... DONE\033[m${demB}"
	sleep 2
}

ScanTRAJ_SMD()
{
if [[ ! -f "SMD_trajectDetails.log" ]]; then
	echo "${demA}"$' Checking the SMD trajectory to extract info about the number of frames and\n simulation time'"${demB}"
	sleep 2
	eval $gmx_exe_path check -f pull.xtc |& tee SMD_trajectDetails.log
	No_of_frames=$(cat SMD_trajectDetails.log | grep "Last" | awk '{print $(NF-2)}')
	simDuratnps=$(cat SMD_trajectDetails.log | grep "Last" | awk '{print $NF}')
	sim_timestep=$(cat SMD_trajectDetails.log | grep -A1 "Item" | awk '{print $NF}' | tail -n 1)
	#simDuratn_nsFloat=$(echo "${simDuratnps%\.*} / 1000" | bc -l)
	simDuratn_nsFloat=$(awk "BEGIN {print $simDuratnps / 1000}")
	simDuratnINTns=$(echo ${simDuratn_nsFloat%\.*})
	echo $simDuratnINTns > simulation_duration

	echo -e "${demA}\033[92m Extract number of frames and simulation duration from trajectory...DONE\033[m${demB}"
	sleep 2
else
	No_of_frames=$(cat SMD_trajectDetails.log | grep "Last" | awk '{print $(NF-2)}')
	simDuratnps=$(cat SMD_trajectDetails.log | grep "Last" | awk '{print $NF}')
	sim_timestep=$(cat SMD_trajectDetails.log | grep -A1 "Item" | awk '{print $NF}' | tail -n 1)
	#simDuratn_nsFloat=$(echo "${simDuratnps%\.*} / 1000" | bc -l)
	simDuratn_nsFloat=$(awk "BEGIN {print $simDuratnps / 1000}")
	simDuratnINTns=$(echo ${simDuratn_nsFloat%\.*})
	echo $simDuratnINTns > simulation_duration
fi
}

umbre_s11_SMD1()
{
	echo "${demA}"" Now executing grompp for steered MD simulation...""${demB}"
	sleep 2
	eval $gmx_exe_path grompp -f md_pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr -maxwarn $WarnMax
}

umbre_s12_SMD2()
{
	echo "${demA}"" Now running steered MD simulation...""${demB}"
	sleep 2
	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm pull -pf pullf.xvg -px pullx.xvg

	echo "${demA}"$' Generating finished figures of key results of the pulling simulation...'"${demB}"
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

	gracebat pullx.xvg -hdevice PNG -autoscale xy -printfile pullx.png -fixed 7500 4000 -legend load || true

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
			base_bkupAnadir=$(basename "$bkupAnadir")
		done
		mv "$currentAnadir" "$bkupAnadir" && mkdir ./$AnaName
		echo "Backing up the last $AnaName folder and its contents as $base_bkupAnadir"
		mv pullf.png pullx.png pullf.xvg pullx.xvg ./$AnaName || true
		mv displacement_pullForce.xvg displacement_pullForce.png ./$AnaName || true	
	elif [[ ! -d "$currentAnadir" ]]; then mkdir ./$AnaName
		mv pullf.png pullx.png pullf.xvg pullx.xvg ./$AnaName || true
		mv displacement_pullForce.xvg displacement_pullForce.png ./$AnaName || true
	fi

	echo "${demA}"$'Generate finished figures of key results of the pulling simulation... DONE'
	sleep 2
	echo -e "${demA}\033[92m Steered MD simulation completed...\033[m${demB}"
	sleep 2
}


s12MDrun()
{
	echo "${demA}"" Now running production MD...""${demB}"
	sleep 2
	if [[ $nohp == '' ]]; then
		eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm}
	elif [[ $nohp == 1 ]]; then
		echo "${demA}"$' gmx mdrun will be launched and run in the background.\n\nYou can tail the '\
		$'output file "nohup.txt" to monitor the progress of the simulation at any time...\n'"${demB}"
		sleep 2
	
		nohup eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm}
	
		echo -e "\n${demA}\033[92m Simulation Completed...\033[m${demB}"
		sleep 3
		Credit
		sleep 2
		echo $'\n'
		read -p ' Do you want to run any post-simulation processing/analyses?(yes/no): ' psa

		while [[ ! "${valid_YesNo_response[@]}" =~ "${psa}" ]] ; do
		# while [[ "$psa" != "yes" && "$psa" != "no" && "$psa" != '"yes"' && "$psa" != '"no"' ]] ; do
			echo $' Please enter the appropriate response ("yes"/"y" or a "no"/"n")!!\n'
			read -p '*Do you want to run any post-simulation processing/analyses?(y/n): ' psa 
		done

		if test "${acceptable_Yes_response[@]}" =~ "$psa" ; then Analysis ; fi

	elif [[ $nohp == '' ]]; then
	
		eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm}
	
		echo -e "\n${demA}\033[92m Simulation Completed...\033[m${demB}"
		sleep 3
		Credit
		sleep 2
		echo $'\n'
		read -p '*Do you want to run any post-simulation processing/analyses?(yes/no): ' psa

		while [[ "$psa" != "yes" && "$psa" != "no" && "$psa" != '"yes"' && "$psa" != '"no"' ]] ; do
			echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
			read -p '**Do you want to run any post-simulation processing/analyses?(yes/no): ' psa 
		done

		if test "$psa" == "yes"; then Analysis; fi

	fi
	echo -e "\n${demA}\033[92m Simulation Completed...\033[m${demB}"
	sleep 3
	Credit
	sleep 2
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
	echo "${demA}"" CHAPERONg will now run an extended production MD...""${demB}"
	sleep 2

if [[ $nohp == 1 ]]; then
	
	echo "${demA}"$' gmx mdrun will be launched and run in the background.\n\n\
	You can tail the output file "nohup.txt" to monitor the progress of the simulation at any time...\n'"${demB}"
	sleep 2
	
	nohup eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm} -s ${filenm}."tpr" -cpi ${filenm}."cpt" -append
	
	echo -e "\n${demA}\033[92m Simulation extended/appended successfully...\033[m${demB}"
	sleep 2
	Credit
	sleep 2
	echo $'\n'
	read -p '**Do you want to run any post-simulation processing/analyses?(yes/no): ' psa

	while [[ "$psa" != "yes" ]] && [[ "$psa" != "no" ]] && [[ "$psa" != '"yes"' ]] && [[ "$psa" != '"no"' ]] ; do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p '**Do you want to run any post-simulation processing/analyses?(yes/no): ' psa 
	done

	if test "$psa" == "yes"; then Analysis ; fi

elif [[ $nohp == '' ]]; then
	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm ${filenm} -s ${filenm}."tpr" -cpi ${filenm}."cpt" -append
		
	#echo $'\n'"${demA}""CHAPERONg: Simulation Completed...""${demB}"
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
	echo -e "\n${demA}\033[92m Simulation extended/appended successfully...\033[m${demB}"
	sleep 2
	Credit
	sleep 2
	echo $'\n'
	read -p '*Do you want to run any post-simulation processing/analyses?(yes/no): ' psa

	while [[ "$psa" != "yes" ]] && [[ "$psa" != "no" ]] && [[ "$psa" != '"yes"' ]] && [[ "$psa" != '"no"' ]] ; do
		echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		read -p '**Do you want to run any post-simulation processing/analyses?(yes/no): ' psa 
		# un umbrella sampling for an additional window
	done

	if test "$psa" == "yes"; then Analysis
	fi
}

ScanTRAJ()
{
if [[ ! -f "${filenm}_${wraplabel}.xtc" ]];	then
	echo -e "${demA} The trajectory hasn't been corrected for pbc yet. Running this first...${demB}"
	sleep 2
	analyser0
fi
if [[ ! -f "trajectDetails.log" ]]; then
	echo -e "${demA} Checking the trajectory to extract info about the number of frames and\n simulation time${demB}"
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

	echo -e "${demA} \033[92mExtract number of frames and simulation duration from trajectory...DONE\033[m${demB}"
	sleep 2
else
	ScanTraj_again_if_err()
	{
		echo -e "${demA} Checking the trajectory to extract info about the number of frames and\n simulation time${demB}"
		sleep 2
		eval "$gmx_exe_path" check -f "${filenm}"_${wraplabel}.xtc |& tee trajectDetails.log		
	}
	No_of_frames=$(cat trajectDetails.log | grep "Last" | awk '{print $(NF-2)}') || ScanTraj_again_if_err
	simDuratnps=$(cat trajectDetails.log | grep "Last" | awk '{print $NF}')
	simDuratnpsINT=$(echo ${simDuratnps%\.*})
	sim_timestep=$(cat trajectDetails.log | grep -A1 "Item" | awk '{print $NF}' | tail -n 1)
	#simDuratn_nsFloat=$(echo "${simDuratnps%\.*} / 1000" | bc -l)
	simDuratn_nsFloat=$(awk "BEGIN {print $simDuratnps / 1000}")
	simDuratnINTns=$(echo ${simDuratn_nsFloat%\.*})
	echo $simDuratnINTns > simulation_duration
fi
}

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

# define function for variables_for_SMD_Movie
variables_for_SMD_Movie()
{
	message_Movie="Preparing to make a summary movie of the SMD trajectory"
	trajectlog="SMD_trajectDetails.log"
	simulationcontext="steered MD"
	xtcFileMovie="pull"
	tprFileMovie="pull"
	outXTCmovie="pull_trjEvery""$skimov""skipForMovie"
	movieDIRECORY="MOVIE_SMD"
}

analyser10()
{
	[[ $mdType == "regularMD" ]] && ScanTRAJ
	echo "${demA} $message_Movie"

	if [[ $customframeNo == '' ]]; then
cat << askMovielength

  Do you want to proceed to making a movie summarized into 200-300 frames?

    1) Yes.
    2) No, I want a different number of frames for the movie.

askMovielength

	read -p '  Enter 1 or 2 here: ' movieLeng
		while [[ "$movieLeng" != 1 && "$movieLeng" != 2 ]]; do
			echo $'\n  You entered: '"$movieLeng"
			echo $'  Please enter a valid number (1 or 2)!!\n'
			read -p '  Enter 1 or 2 here: ' movieLeng
		done
	if [[ $movieLeng == 1 ]] ; then
		customframeNo_int=200
	elif [[ $movieLeng == 2 ]] ; then
		cat "$trajectlog" | grep -v "GROMACS reminds"
		sleep 2

		echo "${demA} Above is a summary of your $simulationcontext trajectory."
		sleep 1
		echo $' You may find the info useful to provide a response to the prompt below.'
		sleep 2
		echo "${demA}"$' How many frames do you want the movie to be composed of?\n'
		sleep 1
		read -p ' *Please enter a value here: ' customframeNo
		echo " You entered: $customframeNo"$'\n'
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
	echo "${demA}"" Number of frames in the trajectory: ${No_of_frames}"\
	$'\n'" Total number of frames in the trajectory is less than $customframeNo!"\
	$'\n'" Using ${No_of_frames} frames directly.""${demB}"
	customframeNo=$(echo ${No_of_frames})
	customframeNo_int=$(echo ${No_of_frames})
fi

echo "${demA}"$' Will now extract frames to be used for the movie...\n\n'
sleep 2
#if [[ $sysType == 1 ]] || [[ $sysType == 2 ]] || [[ $sysType == 3 ]] && [[ $flw == 1 ]] ; then
echo 0 | eval $gmx_exe_path trjconv -f "$xtcFileMovie".xtc -s ${tprFileMovie}.tpr -o ${outXTCmovie}.xtc -skip $skimov
#elif [[ $sysType == 1 ]] || [[ $sysType == 2 ]] && [[ $flw == 0 ]] ; then
#gmx trjconv -f "${filenm}"_"${wraplabel}".xtc -s "${filenm}".tpr -o "${outXTCmovie}.xtc" -skip $skimov
#fi
sleep 2
echo "${demA} Preparing to extract $customframeNo snapshots..."$'\n\n'
sleep 2
if [[ $automode == "full" && $sysType == "protein_only" ]]; then
	echo 1 | eval $gmx_exe_path trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -o summaryForMovie.pdb
elif [[ $automode == "semi" && $sysType == "protein_only" ]]; then
	eval $gmx_exe_path trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -o summaryForMovie.pdb
elif [[ $automode == "semi" || $automode == "full" ]] && [[ $sysType == "protein_dna" ]]; then
	echo "Protein_DNA" | eval $gmx_exe_path trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -n index.ndx -o summaryForMovie.pdb
elif [[ $automode == "semi" && $sysType == "protein_lig" ]]; then
	eval $gmx_exe_path trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -n index.ndx -o summaryForMovie.pdb
elif [[ $automode == "full" && $sysType == "protein_lig" ]]; then
	echo "Protein_$ligname" | eval $gmx_exe_path trjconv -f ${outXTCmovie}.xtc -s ${tprFileMovie}.tpr -n index.ndx -o summaryForMovie.pdb
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
$'mclear\n'"cd ${movieDIRECORY}" > prep_movie_Pyscript.pml

if [[ "$mdType" == "umbrellaSampl" ]]; then
	echo "show cell" >> prep_movie_Pyscript.pml
fi

echo $'save PyMOLsession_allSet.pse\nmpng frame_.png\nquit' >> prep_movie_Pyscript.pml
	
echo "${demA}"$' Now, PyMOL will do the job. You sit back and have a cup of tea...Cheers!'"${demB}"
sleep 2
pyM=0
pymol prep_movie_Pyscript.pml || pyM=1
echo -e "${demA}\033[92m Extract frames as images...DONE\033[m${demB}"\
"${demA}"$' Now converting images to movie...\n'
cd ./${movieDIRECORY}
mov_make=0
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.gif || mov_make=1
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.mp4 || mov_make=2

if [[ "$mov_make" == 1 ]] && [[ "$pyM" == 0 ]]; then
	echo "${demA}"$'The program "Convert/ImageMagick" could not be found.\nCHAPERONg detected'\
	$'"PyMOL" and will use it to make a movie which may,\n'\
	$'however, be of lesser quality'"${demB}"
		
	makeMoviePy1
	makeMoviePy2
		
	echo "${demA}"$'Movie (lesser quality) made with PyMOL...\n'
fi

if [[ "$mov_make" == 2 ]] && [[ "$pyM" == 0 ]]; then
	echo "${demA}"$' A mp4 movie could not be made. This may be due to the program'\
	$'\n "Convert/ImageMagick" not being found, or some library is missing.'\
	$'\n CHAPERONg detected "PyMOL" and will use it to make a mp4 movie which may,'\
	$'\n however, be of lesser quality'"${demB}"
		
	makeMoviePy1
	makeMoviePy2
		
	echo "${demA}"$'Movie (lesser quality) made with PyMOL...\n'
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
echo -e "${demA}\033[92m Convert images to movie...DONE\033[m${demB}"
cd ..

if [[ $mdType == "umbrellaSampl" ]] ; then
	read -p ' Do you want to proceed with umbrella sampling? (yes/no): ' procdWithUS

	while [[ ! " ${valid_YesNo_response[@]} " =~ " ${procdWithUS} " ]] ; do
		echo $'Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
		read -p ' Do you want to proceed with umbrella sampling? (y/n): ' procdWithUS
	done

	if [[ "${acceptable_Yes_response[@]}" =~ "$procdWithUS" ]] ; then echo ""
	elif [[ "${acceptable_No_response[@]}" =~ "$procdWithUS" ]] ; then exit 0
	fi
fi

}

analyser10update()
{

echo "${demA}"$' Preparing to make a summary movie from a preset PyMOL session\n'
sleep 2

currentMOVIEdir="$(pwd)""/${movieDIRECORY}"
if [[ ! -d "$currentMOVIEdir" ]]; then
	echo " No MOVIE directory from a previous run exists. Exiting..."
	exit 1
fi
nMOVIE=1
currentMOVIEgif="$(pwd)""/${movieDIRECORY}/dynamics_movie.gif"
bkupMOVIEgif="$(pwd)""/${movieDIRECORY}/dynamics_movie_""backup""$nMOVIE"".gif"
if [[ -f "$currentMOVIEgif" ]]; then
	base_currentMOVIEgif=$(basename "$currentMOVIEgif")
	base_bkupMOVIEgif=$(basename "$bkupMOVIEgif")
	echo $'\n'" $base_currentMOVIEgif" "exists, backing it up as $base_bkupMOVIEgif"$'\n'
	sleep 1
	while [[ -f "$bkupMOVIEgif" ]]; do
	nMOVIE=$(( nMOVIE + 1 ))
	bkupMOVIEgif="$(pwd)""/${movieDIRECORY}/dynamics_movie_""backup""$nMOVIE"".gif"
	done
	mv "$currentMOVIEgif" "$bkupMOVIEgif" || true
	echo $'\n'" Backing up the last .gif MOVIE as $base_bkupMOVIEgif"
	sleep 1
fi	

nMOVIE=1
currentMOVIEmp4="$(pwd)""/${movieDIRECORY}/dynamics_movie.mp4"
bkupMOVIEmp4="$(pwd)""/${movieDIRECORY}/dynamics_movie_""backup""$nMOVIE"".mp4"
if [[ -f "$currentMOVIEmp4" ]]; then
	base_currentMOVIEmp4=$(basename "$currentMOVIEmp4")
	base_bkupMOVIEmp4=$(basename "$bkupMOVIEmp4")
	echo $'\n'" $base_currentMOVIEmp4"" exists, backing it up as $base_bkupMOVIEmp4"$'\n'
	sleep 1
	while [[ -f "$bkupMOVIEmp4" ]]; do
	nMOVIE=$(( nMOVIE + 1 ))
	bkupMOVIEmp4="$(pwd)""/${movieDIRECORY}/dynamics_movie_""backup""$nMOVIE"".mp4"
	done
	mv "$currentMOVIEmp4" "$bkupMOVIEmp4" || true
	echo $'\n'" Backing up the last .mp4 MOVIE as $base_bkupMOVIEmp4"
	sleep 1
fi	
	
echo "cd ${movieDIRECORY}"$'\nload PyMOLsession_allSet.pse\nmpng frame_.png\nquit' > prep_movie_Pyscript.pml
	
echo "${demA}"$' Now, PyMOL will do the job. Sit back and have a cup of tea...Cheers!'"${demB}"
sleep 2
pyM=0
pymol prep_movie_Pyscript.pml || pyM=1
echo "${demA}"$' Extract frames as images...DONE\n\n Now converting images to movie...'"${demB}"
cd ./${movieDIRECORY}
mov_make=0
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.gif || mov_make=1
convert -delay 5 -loop 0 -dispose Background frame_*.png dynamics_movie.mp4 || mov_make=2

if [[ "$mov_make" == 1 ]] && [[ "$pyM" == 0 ]]; then
	echo "${demA}"$' The program ''"'"Convert/ImageMagick "'"'"could not be found. CHAPERONg detected "'"'\
	$'PyMOL''"'" and will use it to make a movie which may, however, be of lesser quality""${demB}"
		
	makeMoviePyx
	makeMoviePyy
		
	echo "${demA}"$' Movie (lesser quality) made with PyMOL...\n'
fi

if [[ "$mov_make" == 2 ]] && [[ "$pyM" == 0 ]]; then
	echo "${demA}"$' A mp4 movie could not be made. This may be due to the program'\
	$'\n "Convert/ImageMagick" not being found, or some library is missing.'\
	$'\n CHAPERONg detected "PyMOL" and will use it to make a mp4 movie which may,'\
	$'\n however, be of lesser quality'"${demB}"
		
	makeMoviePyx
	makeMoviePyy
		
	echo "${demA}"$' Movie made with PyMOL...\n'
fi

rm frame_*.png || true 
# rm ../*movie_Pyscript.pml || true
rm ./*movie_Pyscript.pml || true
echo -e "${demA}\033[92m Convert images to movie...DONE\033[m${demB}"
cd ..

if [[ $mdType == "umbrellaSampl" ]] ; then
	read -p ' Do you want to proceed with umbrella sampling? (yes/no): ' procdWithUS

	while [[ ! " ${valid_YesNo_response[@]} " =~ " ${procdWithUS} " ]] ; do
		echo $'Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
		read -p ' Do you want to proceed with umbrella sampling? (y/n): ' procdWithUS
	done

	if [[ "${acceptable_Yes_response[@]}" =~ "$procdWithUS" ]] ; then echo ""
	elif [[ "${acceptable_No_response[@]}" =~ "$procdWithUS" ]] ; then exit 0
	fi
fi

}

# umbre_s13_SMD_movie()
# {
# 	ScanTRAJ_SMD; variables_for_SMD_Movie; analyser10;

# }

umbre_s13_SMD_movie()
{
	ScanTRAJ_SMD; variables_for_SMD_Movie
	if [[ "$stage" == 13 ]] && [[ -d "$movieDIRECORY" ]]; then
cat << MovChoic
$demA
Make a new movie or adjust (e.g. the orientation of) a previously prepared one?

  a     Make a new movie
  b     Adjust a previous one

MovChoic

		read -p ' *Enter your choice here (a or b): ' moviechoic

		while [[ "$moviechoic" != "a" ]] && [[ "$moviechoic" != "b" ]] ; do
			echo $'\n You entered: '"$moviechoic"$'\n'
			echo $' Please enter a valid letter!!\n'
			read -p ' *Enter your choice here (a or b): ' moviechoic
		done
		
		if [[ "$moviechoic" == "a" ]]; then analyser10
		elif [[ "$moviechoic" == "b" ]]; then analyser10update
		fi
	elif [[ "$stage" == 13 ]] && [[ ! -d "$movieDIRECORY" ]]; then analyser10
	fi
}

umbre_s14_xtractFrames()
{
	echo "${demA}"" Now extracting frames from the steered MDS trajectory...""${demB}"
	sleep 2
	currentcoords_SMDdir="$(pwd)""/coordinates_SMD"
	ncoords_SMD=1
	bkupcoords_SMDdir="$(pwd)""/#coordinates_SMD"".""backup.""$ncoords_SMD"
	base_bkupcoords_SMDdir=$(basename "$bkupcoords_SMDdir")
	if [[ -d "$currentcoords_SMDdir" ]]; then
		echo "$currentcoords_SMDdir"$' folder exists,\n'"backing it up as $base_bkupcoords_SMDdir"
		sleep 1
		while [[ -d "$bkupcoords_SMDdir" ]]; do
			ncoords_SMD=$(( ncoords_SMD + 1 ))
			bkupcoords_SMDdir="$(pwd)""/#coordinates_SMD"".""backup.""$ncoords_SMD"
			base_bkupcoords_SMDdir=$(basename "$bkupcoords_SMDdir")
		done
		mv "$currentcoords_SMDdir" "$bkupcoords_SMDdir" && mkdir ./coordinates_SMD || true
		echo $'\nBacking up the last coordinates_SMD folder and its contents as'\
		$'\n'"$base_bkupcoords_SMDdir"$'\n\n\n'
		sleep 1
	elif [[ ! -d "$currentcoords_SMDdir" ]]; then mkdir coordinates_SMD
	fi
	# mv coordinate*.gro ./coordinates_SMD || true
	echo 0 | eval $gmx_exe_path trjconv -s pull.tpr -f pull.xtc -o ./coordinates_SMD/coordinate.gro -sep
	echo -e "${demA}\033[92m Extract frames from the steered MDS trajectory...DONE\033[m${demB}"
	sleep 2
}

umbre_s15_calcCOMdist()
{
	echo "${demA}"" Now calculating COM distances...""${demB}"
	sleep 2
	group1_name=$(cat md_pull.mdp | grep "pull_group1_name" | awk '{print $3}')
	group2_name=$(cat md_pull.mdp | grep "pull_group2_name" | awk '{print $3}')
	StructNo=0

	currentDist_SMDdir="$(pwd)""/distances_SMD"
	nDist_SMD=1
	bkupDist_SMDdir="$(pwd)""/#distances_SMD"".""backup.""$nDist_SMD"
	base_bkupDist_SMDdir=$(basename "$bkupDist_SMDdir")
	if [[ -d "$currentDist_SMDdir" ]]; then
		echo "$currentDist_SMDdir"$' folder exists,\n'"backing it up as $base_bkupDist_SMDdir"
		sleep 1
		while [[ -d "$bkupDist_SMDdir" ]]; do
			nDist_SMD=$(( nDist_SMD + 1 ))
			bkupDist_SMDdir="$(pwd)""/#distances_SMD"".""backup.""$nDist_SMD"
			base_bkupDist_SMDdir=$(basename "$bkupDist_SMDdir")
		done
		mv "$currentDist_SMDdir" "$bkupDist_SMDdir" && mkdir ./distances_SMD || true
		echo $'\nBacking up the last distances_SMD folder and its contents as'\
		$'\n'"$base_bkupDist_SMDdir"$'\n\n\n'
		sleep 1
	elif [[ ! -d "$currentDist_SMDdir" ]]; then mkdir distances_SMD
	fi

	# currentDist_SMDdir="$(pwd)""/distances_SMD"
	# if [[ -d "$currentDist_SMDdir" ]]; then
	# 	rm -r distances_SMD || true
	# 	mkdir distances_SMD || true
	# elif [[ ! -d "$currentDist_SMDdir" ]]; then mkdir distances_SMD
	# fi

	com_groups=$'"'"com of group $group1_name plus com of group $group2_name"$'"'

	currentcoords_SMDdir="$(pwd)""/coordinates_SMD"
	if [[ -d "$currentcoords_SMDdir" ]]; then
		for Structure in ./coordinates_SMD/"coordinate"*".gro" ; do
			#calculate distance between the groups
			eval $gmx_exe_path distance -s pull.tpr -f ./coordinates_SMD/coordinate"$StructNo".gro \
			-n index.ndx -select "$com_groups" -oall ./distances_SMD/dist${StructNo}.xvg
			sleep 1
			# if [[ $StructNo == 50 || $StructNo == 100 || $StructNo == 150 || \
			# 	$StructNo == 200 || $StructNo == 250 || $StructNo == 300 ]]
			# then sleep 2
			# fi

			# no_of_structure_milestone=(50 100 150 200 250 300 350 400 450 500)
			# if [[ "${no_of_structure_milestone[@]}" =~ "${StructNo}" ]]; then
			# sleep 2
			# fi

			if (( StructNo % 50 == 0 )); then sleep 2
			fi

			# extract the distances into a summary file
			distanc=$(tail -n 1 ./distances_SMD/dist${StructNo}.xvg | awk '{print $2}')
			if [[ $StructNo == 0 ]] ; then
				echo "$StructNo"$'\t'"$distanc" > distances_summary.txt
			else
				echo "$StructNo"$'\t'"$distanc" >> distances_summary.txt
			fi
			StructNo=$(( StructNo + 1 ))
		done
		rm -r distances_SMD
		echo -e "${demA}\033[92m Calculate COM distances...DONE\033[m${demB}"
		sleep 2
	elif [[ ! -d "$currentcoords_SMDdir" ]]; then
		echo " Frames have not yet been extracted from the steered MDS trajectory!"
		sleep 4
		echo " Run stage 14 first!!"
		sleep 4
		exit 0
	fi
}

umbre_s16_findIniConf()
{
	echo "${demA}"" Now identifying initial configurations for umbrella sampling...""${demB}"
	sleep 2
	# create a new configuratns_list.txt file and backup existing one
	if [[ ! -f "configuratns_list.txt" ]]; then echo ""
	elif [[ -f "configuratns_list.txt" ]]; then
		configList="configuratns_list.txt"
		nbkUP=1
		configListbak="$configList"".backup."
		prevconfigListbak="$configListbak""$nbkUP"
		if [[ -f "$prevconfigListbak" ]]; then
			echo "${demA}""$prevconfigListbak"" exists!"
			nbkUP=$(( nbkUP + 1 ))
			currentconfigListbak="$configListbak""$nbkUP"
			while [[ -f "$currentconfigListbak" ]]; do
				nbkUP=$(( nbkUP + 1 ))
				currentconfigListbak="$configListbak""$nbkUP"
			done
			echo $'\nCHAPERONg will back that up as '"$currentconfigListbak"$' and'\
			$'the '"$configList used in the last run as $configList"".backup.1"
			sleep 2
			mv "$prevconfigListbak" "$currentconfigListbak"
			cp "$configList" "$configList"".backup.1"
		elif [[ ! -f "$prevconfigListbak" ]]; then
			echo "${demA}"" Backing up $configList used in the last run as $configList"".backup.1"
			sleep 2
			cp "$configList" "$configList"".backup.1"
		fi
	fi
	touch configuratns_list.txt

	#run CHAP_set_US_starting_configs.py script
	echo $'\n'" Identifying corresponding frames at the specified window spacing..."
	sleep 2
	python3 ${CHAPERONg_PATH}/CHAP_utilities/CHAP_set_US_starting_configs.py || \
	python ${CHAPERONg_PATH}/CHAP_utilities/CHAP_set_US_starting_configs.py

	echo -e "\n\033[92m Identify initial configurations for umbrella sampling...DONE\033[m${demB}"
	sleep 2
}

US_fxn()
{
	us_frame=$(echo "$line" | awk '{print $1}')
	echo "${demA}"" Now running NPT equilibration for configuration $us_frame"
	sleep 1

	eval $gmx_exe_path grompp -f npt_umbrella.mdp -c ./coordinates_SMD/coordinate"$us_frame".gro \
	-p topol.top -r ./coordinates_SMD/coordinate"$us_frame".gro -n index.ndx -o \
	npt_win"$window"_conf"$us_frame".tpr -maxwarn $WarnMax

	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm npt_win"$window"_conf"$us_frame"

	echo -e "${demA}\033[92m Run NPT equilibration for configuration $us_frame...DONE\033[m${demB}"
	sleep 1

	echo "${demA} Now running umbrella sampling for configuration $us_frame"$'\n\n'
	sleep 1

	eval $gmx_exe_path grompp -f md_umbrella.mdp -c npt_win"$window"_conf"$us_frame".gro \
	-t npt_win"$window"_conf"$us_frame".cpt -p topol.top -r npt_win"$window"_conf"$us_frame".gro \
	-n index.ndx -o umbrella_win"$window"_conf"$us_frame".tpr -maxwarn 1

	eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm umbrella_win"$window"_conf"$us_frame"

	echo -e "${demA}\033[92m Run umbrella sampling for configuration $us_frame...DONE\033[m${demB}"
	sleep 1
	config_no=$(( config_no + 1 ))
}

umbre_s17_USampling()
{
	echo "${demA}"" Initiating umbrella sampling...""${demB}"
	sleep 2
	window=0
	config_no=0
	if [[ "$stage" == 17 ]]; then window="$resume_win"
	fi
	while IFS= read -r line; do
		if [[ $line == *"#"* ]] ; then continue
		elif [[ $line != *"#"* && $stage != 17 ]] ; then US_fxn
		elif [[ $line != *"#"* && $stage == 17 ]] ; then 
			if (( $config_no < "$resume_win" )) ; then
				config_no=$(( config_no + 1 ))
				continue
			elif (( $config_no >= "$resume_win" )) ; then
				US_fxn
			fi
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

		# Remove empty lines
		sed -i '/^$/d' tpr_files.dat
		sed -i '/^$/d' pullf_files.dat
		sed -i '/^$/d' pullx_files.dat		

	done < configuratns_list.txt

	echo "${demA}"$' Umbrella sampling simulation has been completed for all the windows.\n\n'\
	$'The files tpr_files.dat, pullf_files.dat and pullx_files.dat have been written.\n'\
	$'These files list, respectively, the names of the .tpr, pullf.xvg and\n'\
	$'pullx.xvg files for the windows.\n'\
	$'*The files tpr_files.dat and pullf_files.dat will be used for data analysis...'"${demB}"
	sleep 2
	#mv npt_win*_pull*.xvg umbrella_sampling/equilibration/
	#mv umbrella_win*_pullx.xvg ./umbrella_sampling/data_collection
}

umbre_s18_WHAM()
{
	echo "${demA} Extracting the PMF and plotting the umbrella histograms...""${demB}"
	sleep 2
	
	eval $gmx_exe_path wham -it tpr_files.dat -if pullf_files.dat -o \
	PMF_profile.xvg -hist umbrella_sampling_histograms0.xvg -unit kCal
	rm umbrella_sampling_histograms0.xvg || true
	sleep 2

	# minPMFdG=$(grep -v "^[@#]" PMF_profile.xvg | sort -gk 2,2 | head -1 | awk '{print $2}')
	# maxPMFdG=$(grep -v "^[@#]" PMF_profile.xvg | sort -gk 2,2 | tail -1 | awk '{print $2}')
	# displacentATdGmin=$(grep -v "^[@#]" PMF_profile.xvg | sort -gk 2,2 | head -1 | awk '{print $1}')

	## The above method works, but sometimes, for some reason,
	##  the script terminates after the minPMFdG line
	## So I splitted the long command into bits. Why? This is the 2nd time I experienced this with
	##  the sort command. And that was what solved the problem the 1st time (or so I believed)
	
	grep -Ev "^[@#]|^$" PMF_profile.xvg > tempPMF
	sort -gk 2,2 -o tempPMF tempPMF
	head -n 1 tempPMF > tempPMFmin
	tail -n 1 tempPMF > tempPMFmax
	minPMFdG=$(grep "[0-9]" tempPMFmin | awk '{print $2}')
	maxPMFdG=$(grep "[0-9]" tempPMFmax | awk '{print $2}')
	displacentATdGmin=$(grep "[0-9]" tempPMFmin | awk '{print $1}')
	rm tempPMF tempPMFmin tempPMFmax
	
	# gromacs .xvg outputs often contain decimals in scientific notations and 
	# sort command with the -g flag handles that.
	# LC_ALL=C is necessary for users whose locales use a comma instead of a period to indicate
	# decimals, in which case the sort -g command would fail to sort dot decimals properly

	eval $gmx_exe_path wham -it tpr_files.dat -if pullf_files.dat -o PMF_profile_YminAdjusted.xvg \
	-hist umbrella_sampling_histograms.xvg -unit kCal -zprof0 $displacentATdGmin

	# shift the PMF plot so that the displacement (x-axis) starts at 1 nm
	# write out the XVG comments and headers
	grep "^[@#]" PMF_profile_YminAdjusted.xvg > PMF_profile_XminYminAdjusted.xvg
	count_data=0
	while IFS= read -r line; do
		# skip comments and headers
		[[ "$line" =~ ^[#@] ]] && continue
		displacement=$(echo "$line" | awk '{print $1}')
		gibbsEnergy=$(echo "$line" | awk '{print $2}')
		count_data=$(( count_data + 1 ))
		if [[ "$count_data" == 1 ]] ; then startingDisplacement="$displacement" ; fi
		adjustedDisplacement=$(awk "BEGIN {print "$displacement" - "$startingDisplacement" + 1}")
		echo -e "${adjustedDisplacement}\t${gibbsEnergy}" >> PMF_profile_XminYminAdjusted.xvg
	done < PMF_profile_YminAdjusted.xvg


	echo "${demA}"$' Generating finished figures of key results of WHAM analysis...'"${demB}"
	sleep 2
	gracebat -nxy umbrella_sampling_histograms.xvg -hdevice PNG -autoscale xy -printfile \
	umbrella_sampling_histograms.png -fixed 7500 4000 || true
	gracebat PMF_profile.xvg -hdevice PNG -autoscale xy -printfile PMF_profile.png \
	-fixed 7500 4000 -legend load || true
	gracebat PMF_profile_YminAdjusted.xvg -hdevice PNG -autoscale xy -printfile \
	PMF_profile_YminAdjusted.png -fixed 7500 4000 -legend load || true
	dG_PMF=$(awk "BEGIN {print $minPMFdG - $maxPMFdG}")
	echo $'Binding Free Energy (dG) = '"$dG_PMF"$' kCal/mol' > summary_dG.dat

	gracebat PMF_profile_XminYminAdjusted.xvg -hdevice PNG -autoscale xy -printfile \
	PMF_profile_XminYminAdjusted.png -fixed 7500 4000 -legend load || true


	AnaName="Data_Analysis_PMF"
	currentAnadir="$(pwd)""/$AnaName"
	nDir=1
	bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
	if [[ -d "$currentAnadir" ]]; then
		base_currentAnadir=$(basename "$currentAnadir")
		base_bkupAnadir=$(basename "$bkupAnadir")
		echo "$base_currentAnadir" "exists, backing it up as $base_bkupAnadir"
		sleep 2
		while [[ -d "$bkupAnadir" ]]; do
			nDir=$(( nDir + 1 )); bkupAnadir="$(pwd)""/#""$AnaName"".backup.""$nDir"
			base_bkupAnadir=$(basename "$bkupAnadir")
		done
		mv "$currentAnadir" "$bkupAnadir" && mkdir ./$AnaName
		echo "Backing up the last $AnaName folder and its contents as $base_bkupAnadir"
		sleep 2
	elif [[ ! -d "$currentAnadir" ]]; then mkdir ./$AnaName
	fi
	mv umbrella_sampling_histograms.xvg PMF_profile.xvg \
	PMF_profile_YminAdjusted.xvg PMF_profile_XminYminAdjusted.xvg ./$AnaName || true
	mv umbrella_sampling_histograms.png PMF_profile.png PMF_profile_YminAdjusted.png \
	PMF_profile_XminYminAdjusted.png summary_dG.dat ./$AnaName || true
	
	echo "${demA}"$' Generate finished figures of results of WHAM analysis...DONE'
	sleep 2
	echo -e "${demA}\033[92m Extract the PMF and plot the umbrella histograms...DONE\033[m${demB}"
	sleep 2
}

umbre_s19_MoreWin()
{
	repeatUSmore="yes"

	while [[ "$repeatUSmore" == "yes" || "$repeatUSmore" == "y" || \
		"$repeatUSmore" == '"yes"' || "$repeatUSmore" == '"y"' ]]
	do
		echo "${demA} Running umbrella sampling for an additional window...""${demB}"
		sleep 2
		read -p ' Enter the frame number of the SMD configuration to use: ' us_frame
		echo $'\n You entered: '"$us_frame"
		sleep 2

		window=0
		while [[ -f npt_win"$window"_* ]]; do
			window=$(( window + 1 ))
		done
			
		echo "${demA} Now running NPT equilibration for configuration $us_frame"
		sleep 1
		eval $gmx_exe_path grompp -f npt_umbrella.mdp -c ./coordinates_SMD/coordinate"$us_frame".gro -p topol.top -r \
		./coordinates_SMD/coordinate"$us_frame".gro -n index.ndx -o npt_win"$window"_conf"$us_frame".tpr -maxwarn $WarnMax

		eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm npt_win"$window"_conf"$us_frame"

		echo -e "${demA}\033[92m Run NPT equilibration for configuration $us_frame...DONE\033[m${demB}"
		sleep 1

		echo "${demA} Now running umbrella sampling for configuration $us_frame"$'\n\n'
		sleep 1
		eval $gmx_exe_path grompp -f md_umbrella.mdp -c npt_win"$window"_conf"$us_frame".gro -t npt_win"$window"_conf"$us_frame".cpt -p \
		topol.top -r npt_win"$window"_conf"$us_frame".gro -n index.ndx -o umbrella_win"$window"_conf"$us_frame".tpr -maxwarn 1

		eval $gmx_exe_path mdrun ${threader} ${THREA} $gpidn -v -deffnm umbrella_win"$window"_conf"$us_frame"

		echo -e "${demA}\033[92m Run umbrella sampling for configuration $us_frame...DONE\033[m${demB}"
		sleep 1
		echo "umbrella_win"$window"_conf"$us_frame".tpr" >> tpr_files.dat
		echo "umbrella_win"$window"_conf"$us_frame"_pullf.xvg" >> pullf_files.dat
		echo "umbrella_win"$window"_conf"$us_frame"_pullx.xvg" >> pullx_files.dat
		
		# Remove empty lines
		sed -i '/^$/d' tpr_files.dat
		sed -i '/^$/d' pullf_files.dat
		sed -i '/^$/d' pullx_files.dat
		
		echo "${demA}"$' Umbrella sampling simulation has been completed for the additional window.\n\n'\
		$' The names of the .tpr, pullf.xvg and pullx.xvg files have been appended to the\n'\
		$' tpr_files.dat, pullf_files.dat and pullx_files.dat, respectively...'"${demB}"
		sleep 2

		read -p ' Do you want to run umbrella sampling for more windows? (yes/no): ' repeatUSmore

		# while [[ "$repeatUSmore" != "yes" && "$repeatUSmore" != "no" && \
		# 	"$repeatUSmore" != '"yes"' && "$repeatUSmore" != '"no"' ]]
		# do
		# 	echo $'Please enter the appropriate response (a "yes" or a "no")!!\n'
		# 	read -p ' Do you want to proceed? (yes/no): ' procd2pmf
		# done

		while [[ ! " ${valid_YesNo_response[@]} " =~ " ${repeatUSmore} " ]]; do
			echo $' Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
			read -p ' Do you want to run umbrella sampling for more windows? (y/n): ' repeatUSmore
		done
	done

	read -p ' Do you want to proceed to PMF calculation? (yes/no): ' procd2pmf

	while [[ ! " ${valid_YesNo_response[@]} " =~ " ${procd2pmf} " ]] ; do
		echo $'Please enter the appropriate response ("yes"/"y" or "no"/"n")!!\n'
		read -p ' Do you want to proceed to PMF calculation? (y/n): ' procd2pmf
	done

	if [[ "${acceptable_Yes_response[@]}" =~ "$procd2pmf" ]] ; then echo ""
	elif [[ "${acceptable_No_response[@]}" =~ "$procd2pmf" ]] ; then exit 0

	# if [[ "$procd2pmf" == "yes" || "$procd2pmf" == '"yes"' ]] ; then echo ""
	# elif [[ "$procd2pmf" == "no" || "$procd2pmf" == '"no"' ]] ; then exit 0
	fi
}