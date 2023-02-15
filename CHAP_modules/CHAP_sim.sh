#! /bin/bash

#CHAP_sim - The simulation module of CHAPERONg
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
#demB=$'\n'"#*****************************************************************************#"$'\n\n'


if [[ $initiator1 != 'avail' ]] ; then
	echo "$demA"$' Do not run modules independently!\n'\
	$' Launch CHAPERONg with run_CHAPERONg-<version>!!'"$demB"	
	exit 1
fi

#call module to collect parameters
. "$CHAPERONg_PATH/CHAP_modules/CHAP_colPar.sh"

#call module with defined fxns
. "$CHAPERONg_PATH/CHAP_modules/CHAP_deffxn.sh"

#call module for traj_ana
. "$CHAPERONg_PATH/CHAP_modules/CHAP_ana.sh"


md_proOnly()
{
echo "$demA"$' Enter a stage number to start/resume the pipeline'
cat << ListStage

---------------------------------------------------------------
 StageNo | Step to begin/resume workflow
---------+-----------------------------------------------------
    0    | Prepare protein topology (pdb2gmx)
    1    | Define box (editconf)
    2    | Solvate (solvate)
    3    | Add ions (grompp)
    4    | Add ions (genion)
    5    | Energy minimization (grompp)
    6    | Energy minimization (mdrun)
    7    | NVT Equilibration (grompp)
    8    | NVT Equilibration (mdrun)
    9    | NPT Equilibration (grompp)
    10   | NPT Equilibration (mdrun)
    11   | Release position restraints (for "production mdrun")
    12   | Production md (mdrun)
    13   | Extend/resume a previously run/terminated simulation
    14   | Make index group(s) (make_ndx)
    15   | Post-simulation processing & analyses
ListStage
}

md_Complex()
{
echo "$demA"$' Enter a stage number to start/resume the pipeline'
cat << ListStageComplx

---------------------------------------------------------------
 StageNo | Step to begin/resume workflow
---------+-----------------------------------------------------
    0a   | Prepare protein topology
    0b   | Prepare ligand topology
    1    | Define box (editconf)
    2    | Solvation (solvate)
    3    | Add ions (grompp)
    4    | Add ions (genion)
    5    | Energy minimization (grompp)
    6    | Energy minimization (mdrun)
    6a   | Ligand restraint
    6b   | Temperature coupling
    7    | NVT Equilibration (grompp)
    8    | NVT Equilibration (mdrun)
    9    | NPT Equilibration (grompp)
    10   | NPT Equilibration (mdrun)
    11   | Release position restraints (for "production mdrun")
    12   | Production md (mdrun)
    13   | Extend/resume a previously run/terminated simulation
    14   | Make index group(s) (make_ndx)
    15   | Post-simulation processing & analyses
ListStageComplx
}

US_mdProLig()
{
echo "$demA"$' Enter a stage number to start/resume the pipeline'
cat << ListStageUS

-----------------------------------------------------------------
 StageNo | Step to begin/resume workflow
---------+-------------------------------------------------------
    0a   | Prepare protein topology
    0b   | Prepare ligand topology
    1    | Define box (editconf)
    2    | Solvation (solvate)
    3    | Add ions (grompp)
    4    | Add ions (genion)
    5    | Energy minimization (grompp)
    6    | Energy minimization (mdrun)
    7    | NVT Equilibration (grompp) -- optional
    8    | NVT Equilibration (mdrun) -- optional
    9    | NPT Equilibration (grompp)
    10   | NPT Equilibration (mdrun)
    11   | Steered MDS (grompp)
    12   | Steered MDS (mdrun)
    13   | Make a movie of the steered MD trajectory
    14   | Extract frames from the steered MDS trajectory
    15   | Calculate COM distances
    16   | Identify initial configurations for umbrella sampling
    17   | Umbrella sampling
    18   | Calculate PMF using WHAM
    19   | Umbrella sampling for additional window(s) -- optional
ListStageUS
}

US_mdProPro()
{
echo "$demA"$' Enter a stage number to start/resume the pipeline'
cat << ListStageUS

 StageNo   Step to begin/resume workflow
---------------------------------------------------------------
 StageNo | Step to begin/resume workflow
---------+-----------------------------------------------------
    0    | Prepare protein topology (pdb2gmx)
    1    | Define box (editconf)
    2    | Solvate (solvate)
    3    | Add ions (grompp)
    4    | Add ions (genion)
    5    | Energy minimization (grompp)
    6    | Energy minimization (mdrun)
    7    | NVT Equilibration (grompp) -- optional
    8    | NVT Equilibration (mdrun) -- optional
    9    | NPT Equilibration (grompp)
    10   | NPT Equilibration (mdrun)
    11   | Steered MDS (grompp)
    12   | Steered MDS (mdrun)
    13   | Make a movie of the steered MD trajectory
    14   | Extract frames from the steered MDS trajectory
    15   | Calculate COM distances
    16   | Identify initial configurations for umbrella sampling
    17   | Umbrella sampling
    18   | Calculate PMF using WHAM
    19   | Umbrella sampling for additional window(s) -- optional
ListStageUS
}

RegMDsimulate()
{
	echo $'\n *ENTER A STAGE NUMBER BELOW\n'
	read -p ' Initiation stage: ' stage

	checkstage=""
	for i in {0..15} ; do
		if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]] && \
			[[ "$stage" == "$i" ]] ; then checkstage="yes"; break
		elif [[ $sysType == "protein_lig" ]] && [[ "$stage" == "$i" || "$stage" == "0a" || \
			"$stage" == "0b" || "$stage" == "6a" || "$stage" == "6b" ]]
			then checkstage="yes"; break
		fi
	done

	while [[ "$checkstage" != "yes" ]] ; do
		echo $'\nYou entered: '"$stage"
		echo $'Please enter a valid value!!\n'
		read -p ' Initiation stage: ' stage

		for i in {0..15} ; do
			if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]] && [[ "$stage" == "$i" ]]
				then checkstage="yes"; break
			elif [[ $sysType == "protein_lig" ]] && [[ "$stage" == "$i" || "$stage" == "0a" || \
			"$stage" == "0b" || "$stage" == "6a" || "$stage" == "6b" ]]
				then checkstage="yes"; break
			fi
		done
	done

	if [[ $sysType == "protein_lig" ]]; then
		echo "$demA"$'\n*What is the name of your ligand (same as for/in your topol.top file)?\n'
		read -p '**Ligand name (without any extension): ' ligname
		echo $'\nYou entered: '"$ligname"$'\n'
	fi

	if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]] && [[ "$stage" == 0 ]] ; then
		s0GenTop; s1DefBox; s2Solvat; s3AddIons1; s4AddIons2; s5EnMin1;
		s6EnMin2; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun

	elif [[ $sysType == "protein_lig" && "$stage" == '0a' ]] ; then s0GenTop; s0GenLigTop
		s1DefBox; s2Solvat; s3AddIons1; s4AddIons2; s5EnMin1; s6EnMin2; s6aLigRes
		s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
	
	elif [[ $sysType == "protein_lig" && "$stage" == '0b' ]] ; then s0GenLigTop
		s1DefBox; s2Solvat; s3AddIons1; s4AddIons2; s5EnMin1; s6EnMin2; s6aLigRes
		s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		
	elif [[ "$stage" == 1 ]]; then
		if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]]; then s1DefBox; s2Solvat; s3AddIons1; s4AddIons2
			s5EnMin1;s6EnMin2; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		elif [[ $sysType == "protein_lig" ]]; then s1DefBox; s2Solvat; s3AddIons1; s4AddIons2; s5EnMin1; s6EnMin2
			s6aLigRes; s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		fi
	elif [[ "$stage" == 2 ]]; then
		if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]]; then s2Solvat; s3AddIons1; s4AddIons2; s5EnMin1
			s6EnMin2; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		elif [[ $sysType == "protein_lig" ]]; then s2Solvat; s3AddIons1; s4AddIons2; s5EnMin1; s6EnMin2;
			s6aLigRes; s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		fi
	elif [[ "$stage" == 3 ]]; then
		if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]]; then s3AddIons1; s4AddIons2; 
			s5EnMin1; s6EnMin2; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		elif [[ $sysType == "protein_lig" ]]; then s3AddIons1; s4AddIons2; s5EnMin1; s6EnMin2; s6aLigRes
			s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		fi
	elif [[ "$stage" == 4 ]]; then
		if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]]; then s4AddIons2; s5EnMin1
			s6EnMin2; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		elif [[ $sysType == "protein_lig" ]]; then s4AddIons2; s5EnMin1; s6EnMin2; s6aLigRes;
			s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		fi
	elif [[ "$stage" == 5 ]]; then
		if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]]; then s5EnMin1; s6EnMin2
			s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		elif [[ $sysType == "protein_lig" ]]; then s5EnMin1; s6EnMin2; s6aLigRes; s6bTempCoup
			s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		fi
	elif [[ "$stage" == 6 ]]; then
		if [[ $sysType == "protein_only" || $sysType == "protein_dna" ]]; then s6EnMin2; s7NVTeq1
			s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		elif [[ $sysType == "protein_lig" ]]; then s6EnMin2; s6aLigRes; s6bTempCoup; s7NVTeq1
			s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
		fi
	
	elif [[ "$stage" == "6a" && "$sysType" == "2" ]]; then s6aLigRes; s6bTempCoup
		s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
	
	elif [[ "${stage}" == "6b" && $sysType == "2" ]]; then s6bTempCoup
		s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun		
	elif [[ "$stage" == 7 ]]; then s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2
		s11RelPosRe; s12MDrun
	elif [[ "$stage" == 8 ]]; then s8NVTeq2; s9NPTeq1; s10NPTeq2
		s11RelPosRe; s12MDrun
	elif [[ "$stage" == 9 ]]; then s9NPTeq1; s10NPTeq2; s11RelPosRe; s12MDrun
	elif [[ "$stage" == 10 ]]; then s10NPTeq2; s11RelPosRe; s12MDrun
	elif [[ "$stage" == 11 ]]; then s11RelPosRe; s12MDrun
	elif [[ "$stage" == 12 ]]; then s12MDrun
	elif [[ "$stage" == 13 ]]; then s13MDappend
	elif [[ "$stage" == 14 ]]; then makeNDXGroup
	elif [[ "$stage" == 15 ]]; then Analysis
	fi
}

US_simulate()
{
	echo $'\n **ENTER A STAGE NUMBER BELOW\n'
	read -p ' Initiation stage: ' stage

	# checkstage=""
	# for i in {0..18} ; do
	# 	if [[ "$stage" == "$i" || "$stage" == "0a" || "$stage" == "0b" ]]
	# 		then checkstage="yes"; break
	# 	fi
	# done

	# while [[ "$checkstage" != "yes" ]] ; do
	# 	echo $'\nYou entered: '"$stage"
	# 	echo $'Please enter a valid value!!\n'
	# 	read -p ' Initiation stage: ' stage

	# 	for i in {0..18} ; do
	# 		if [[ "$stage" == "$i" || "$stage" == "0a" || "$stage" == "0b" ]]
	# 			then checkstage="yes"; break
	# 		fi
	# 	done
	# done

	while ! [[ "$stage" =~ ^(0a|0b|[0-9]|1[0-9]{1})$ ]]
	# Condition:
		# "0a",
		# "0b",
		# any digit between 0 and 9, or
		# "1" followed by exactly one digit between 0 and 8 ([0-8]{1})
	do
		echo $'\nYou entered: '"$stage"
		echo $'Please provide a valid entry!!\n'
		read -p 'Initiation stage: ' stage
	done

	if [[ $sysType == "protein_lig" ]]; then
		while ! [[ "$stage" =~ ^(0a|0b|[0-9]|1[0-9]{1})$ ]]
		# Condition:
			# "0a",
			# "0b",
			# any digit between 0 and 9, or
			# "1" followed by exactly one digit between 0 and 8 ([0-8]{1})
		do
			echo $'\nYou entered: '"$stage"
			echo $'Please provide a valid entry!!\n'
			read -p 'Initiation stage: ' stage
		done

		echo "$demA"$'\n*What is the name of your ligand (same as for/in your topol.top file)?\n'
		read -p '**Ligand name (without any extension): ' ligname
		echo $'\nYou entered: '"$ligname"$'\n'
	elif [[ $sysType == "protein_only" ]]; then
		while ! [[ "$stage" =~ ^([0-9]|1[0-9]{1})$ ]]
		# Condition:
			# any digit between 0 and 9, or
			# "1" followed by exactly one digit between 0 and 8 ([0-8]{1})
		do
			echo $'\nYou entered: '"$stage"
			echo $'Please provide a valid entry!!\n'
			read -p 'Initiation stage: ' stage
		done
	fi

	if [[ $sysType == "protein_only" ]] && [[ "$stage" == 0 ]]; then s0GenTop; s1DefBox; s2Solvat
		s3AddIons1;	s4AddIons2; s5EnMin1; s6EnMin2; s7NVTeq1; s8NVTeq2; s9NPTeq1
		s10NPTeq2; umbre_s11_SMD1; umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames 
		umbre_s15_calcCOMdist; umbre_s16_findIniConf; umbre_s17_USampling; umbre_s18_WHAM

	elif [[ $sysType == "protein_lig" && "$stage" == '0a' ]] ; then s0GenTop; s0GenLigTop; s1DefBox
		s2Solvat; s3AddIons1; s4AddIons2; s5EnMin1; s6EnMin2; s6bTempCoup; s7NVTeq1; s8NVTeq2
		s9NPTeq1; s10NPTeq2; umbre_s11_SMD1; umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames
		umbre_s15_calcCOMdist; umbre_s16_findIniConf; umbre_s17_USampling; umbre_s18_WHAM
	
	elif [[ $sysType == "protein_lig" && "$stage" == '0b' ]] ; then s0GenLigTop; s1DefBox; s2Solvat
		s3AddIons1; s4AddIons2; s5EnMin1; s6EnMin2; s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1
		s10NPTeq2; umbre_s11_SMD1; umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames
		umbre_s15_calcCOMdist; umbre_s16_findIniConf; umbre_s17_USampling; umbre_s18_WHAM
		
	elif [[ "$stage" == 1 ]]; then s1DefBox; s2Solvat; s3AddIons1; s4AddIons2
		s5EnMin1; s6EnMin2; s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2
		umbre_s11_SMD1; umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames;	umbre_s14_calcCOMdist
		umbre_s15_findIniConf; umbre_s16_USampling; umbre_s17_WHAM

	elif [[ "$stage" == 2 ]]; then s2Solvat; s3AddIons1; s4AddIons2; s5EnMin1
		s6EnMin2; s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; umbre_s11_SMD1
		umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames;	umbre_s14_calcCOMdist
		umbre_s15_findIniConf; umbre_s16_USampling; umbre_s17_WHAM

	elif [[ "$stage" == 3 ]]; then s3AddIons1; s4AddIons2; s5EnMin1; s6EnMin2
		s6bTempCoup; s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2; umbre_s11_SMD1
		umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames;	umbre_s14_calcCOMdist
		umbre_s15_findIniConf; umbre_s16_USampling; umbre_s17_WHAM

	elif [[ "$stage" == 4 ]]; then s4AddIons2; s5EnMin1; s6EnMin2; s6bTempCoup; s7NVTeq1
		s8NVTeq2; s9NPTeq1; s10NPTeq2; umbre_s11_SMD1; umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames
		umbre_s15_calcCOMdist; umbre_s16_findIniConf; umbre_s17_USampling; umbre_s18_WHAM

	elif [[ "$stage" == 5 ]]; then s5EnMin1; s6EnMin2; s6bTempCoup; s7NVTeq1; s8NVTeq2
		s9NPTeq1; s10NPTeq2; umbre_s11_SMD1; umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames
		umbre_s15_calcCOMdist; umbre_s16_findIniConf; umbre_s17_USampling; umbre_s18_WHAM
	
	elif [[ "$stage" == 6 ]]; then s6EnMin2; s6bTempCoup; s7NVTeq1; s8NVTeq2
		s9NPTeq1; s10NPTeq2; umbre_s11_SMD1; umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames
		umbre_s15_calcCOMdist; umbre_s16_findIniConf; umbre_s17_USampling; umbre_s18_WHAM

	elif [[ "$stage" == 7 ]]; then s7NVTeq1; s8NVTeq2; s9NPTeq1; s10NPTeq2
		umbre_s11_SMD1; umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames; umbre_s14_calcCOMdist
		umbre_s15_findIniConf; umbre_s16_USampling; umbre_s17_WHAM

	elif [[ "$stage" == 8 ]]; then s8NVTeq2; s9NPTeq1; s10NPTeq2; umbre_s11_SMD1
		umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames; umbre_s14_calcCOMdist
		umbre_s15_findIniConf; umbre_s16_USampling; umbre_s17_WHAM

	elif [[ "$stage" == 9 ]]; then s9NPTeq1; s10NPTeq2; umbre_s11_SMD1
		umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames; umbre_s14_calcCOMdist
		umbre_s15_findIniConf; umbre_s16_USampling; umbre_s17_WHAM

	elif [[ "$stage" == 10 ]]; then s10NPTeq2; umbre_s11_SMD1
		umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames; umbre_s14_calcCOMdist
		umbre_s15_findIniConf; umbre_s16_USampling; umbre_s17_WHAM

	elif [[ "$stage" == 11 ]]; then umbre_s11_SMD1; umbre_s12_SMD2
		umbre_s13_SMD_movie; umbre_s14_xtractFrames; umbre_s14_calcCOMdist
		umbre_s15_findIniConf; umbre_s16_USampling; umbre_s17_WHAM

	elif [[ "$stage" == 12 ]]; then umbre_s12_SMD2; umbre_s13_SMD_movie; umbre_s14_xtractFrames
		umbre_s15_calcCOMdist; umbre_s16_findIniConf; umbre_s17_USampling; umbre_s18_WHAM

	elif [[ "$stage" == 13 ]]; then umbre_s13_SMD_movie; umbre_s14_xtractFrames
		umbre_s15_calcCOMdist; umbre_s16_findIniConf; umbre_s17_USampling; umbre_s18_WHAM

	elif [[ "$stage" == 14 ]]; then umbre_s14_xtractFrames;	umbre_s15_calcCOMdist
		umbre_s16_findIniConf; umbre_s17_USampling; umbre_s18_WHAM

	elif [[ "$stage" == 15 ]]; then umbre_s15_calcCOMdist; umbre_s16_findIniConf
		umbre_s17_USampling; umbre_s18_WHAM

	elif [[ "$stage" == 16 ]]; then umbre_s16_findIniConf; umbre_s17_USampling
		umbre_s18_WHAM

	elif [[ "$stage" == 17 ]]; then
		echo "$demA"$' Provide the window number to start/resume umbrella sampling from.'\
		$'\n To start from the 1st window (window 0), enter zero (0).\n'
		read -p ' Enter the window number here: ' resume_win
		echo $'\n You entered: '"$resume_win"$'\n'
		sleep 2
		umbre_s17_USampling; umbre_s18_WHAM
	elif [[ "$stage" == 18 ]]; then umbre_s18_WHAM
	elif [[ "$stage" == 19 ]]; then umbre_s19_MoreWin; umbre_s18_WHAM
	fi
}

system_typeReg()
{
#check if protein-only or complex mds
cat << SystemType

##*********************************CHAPERONg**********************************##

 What type of system are we simulating?
  (1)  Protein-only (including protein-protein complex)
  (2)  Protein-ligand (small molecule) complex
  (3)  Protein-DNA complex

SystemType

read -p " *Enter 1, 2 or 3 here: " sysType
echo $'\n You entered: '"$sysType"$'\n'
while [[ "$sysType" != 1 && "$sysType" != 2 && "$sysType" != 3 ]]; do
cat << SystemTypeAgain
 *Please enter either 1, 2 or 3 below!!
 
 Protein-only mds (1), Protein-ligand mds (2) or Protein-DNA mds (3)?

SystemTypeAgain
	read -p " *Enter 1, 2 or 3 here: " sysType
done

if [[ "${filenm}" == '' ]]; then filenm="md_${coordinates}"; fi
		
if [[ "$sysType" == 1 ]]; then
	sysType="protein_only"
	wraplabel=noPBC
	md_proOnly
elif [[ "$sysType" == 2 ]]; then
	sysType="protein_lig"
	md_Complex
	wraplabel=center
elif [[ "$sysType" == 3 ]]; then
	sysType="protein_dna"
	md_proOnly
	wraplabel=center
fi

}

system_typeUS()
{
#check if protein-only or complex mds
cat << SystemType

##*********************************CHAPERONg**********************************##

 What type of system are we simulating?
  (1)  Protein-protein complex
  (2)  Protein-ligand (small molecule) complex

SystemType

read -p " *Enter 1 or 2: " sysType
echo $'\n You entered: '"$sysType"$'\n'
while [[ "$sysType" != 1 && "$sysType" != 2 ]]; do
cat << SystemTypeAgain
 *Please enter either 1 or 2 below!!
 
 Protein-only mds (1) or Protein-ligand mds (2)?

SystemTypeAgain
	read -p " *Enter 1 or 2: " sysType
done

if [[ "${filenm}" == '' ]]; then filenm="md_${coordinates}"; fi

#present a list of stages to initiate/resume simulation from
if [[ "$sysType" == 1 ]]; then
	sysType="protein_only"
	US_mdProPro
elif [[ "$sysType" == 2 ]]; then
	sysType="protein_lig"
	US_mdProLig
fi

}

#check for sampling biasedness
cat << Biasedness

##*********************************CHAPERONg**********************************##

 What type of simulation do you want to run?
  (1)  Conventional MD simulation -- Unbiased/Regular
  (2)  Enhanced (umbrella) sampling simulation -- Biased

Biasedness

read -p " *ENTER YOUR CHOICE HERE (1 or 2): " mdType
echo $'\n You entered: '"$mdType"$'\n'
while [[ "$mdType" != 1 && "$mdType" != 2 ]]; do
cat << BiasednessAgain
 *Please enter either 1 or 2 below!!
 
 Conventional mds (1) or Enhanced (umbrella) sampling simulation (2)?

BiasednessAgain
	read -p " *Enter 1 or 2 here: " mdType
done

if [[ "$mdType" == 1 ]]; then
	mdType="regularMD"
	system_typeReg
	RegMDsimulate
	Credit
elif [[ "$mdType" == 2 ]]; then
	mdType="umbrellaSampl"
	system_typeUS
	US_simulate
	Credit
fi
