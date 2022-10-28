#! /bin/bash

#CHAP_colPar - The parameter collection module of CHAPERONg
#CHAPERONg - An automation program for GROMACS md simulation
#Author: Abeeb A. Yekeen
#Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date: 2022.02.11

set -e
set -o pipefail

#set version
CHAPERONg_version="v0.1"

#Initialize (default) parameters
btype='cubic' ; edgeDist="1.0"; WarnMax=0; nt=0
THREA="-nt ""${nt}"; hbthread="-nthreads ""0"
wmodel=""; wat=""; nb='' ; indexer=''
flw=0; ffUse=""; gpid=''; filenm=''
ntmpi=0; ntomp=0; skp=''; nohp=''
threader="-ntmpi ""${ntmpi}"" -ntomp ""${ntomp}"
Enmdp=0; Minmdp=0; wraplabel=''
pn=''; nn=''; ion_conc=''; pnam_nnam=''
extr="-ignh" ; customNDXask=''
PBCcorrectType='' ; trajFraction=''
mmGMX='' ; mmGMXpath='' ; coordinates_raw=''
parfilename='' ; mmpbframesNo=''
gmx_exe_path="gmx"
#gmxV=''

#Defining primary functions
Credit()
{
	echo \
  $'\n###############################################################################'\
  $'\n#--------------------------------- CHAPERONg ---------------------------------#'\
  $'\n#   An automated pipeline for GROMACS MD simulation and trajectory analyses   #'\
  $'\n#    If you use this program in your work, please cite the relevant paper:    #'\
  $'\n#                   Yekeen, A.A. et al. To be published...                    #'\
  $'\n###############################################################################'
}

Help()
{

#Shows the use of the script
cat << guide_sh

Usage:
chmod +x ./setup_CHAPERONg-<version>
./run_CHAPERONg-<version> -i inputStructure_filename

Required (int=integer; str=string):
-i, --input <str>   Input coordinate file (.pdb or .gro)
Optional (int=integer; str=string):
-h, --help          Print this help
-b, --bt <str>      Box type i.e., cubic, dodecahedron, triclinic, octahedron
                    (default: cubic)
-T, --nt <int>      Number of threads to use (default is 0: allow gmx to guess)
-g, --nb gpu        Calculate non-bonded interactions on gpu
-G, --gpu_id <str>  List ID(s) of unique GPU device(s) available for use
-p, --deffnm <str>  Define filename prefix (default is "md_filename" for outputs)
-a, --auto          Automatic mode (less prompts): Use default parameters &
                    do common analyses
--parFile <str>     Name of the CHAPERONg input parameter file
-H, --Help          Print more, advanced options
guide_sh
}

HHelp()
{
#Shows the use of the script
cat << guide_lg

Usage:
chmod +x ./setup_CHAPERONg-<version>
./run_CHAPERONg-<version> -i inputStructure_filename [-More options]

Required (int=integer; str=string):
-i, --input <str>   Input coordinate file (.pdb or .gro)
Optional (int=integer; str=string):
-h, --help          Print the shorter version of this help
-b, --bt <str>      Box type: cubic (default), dodecahedron, triclinic, octahedron
-T, --nt <int>      Number of threads to use (default is 0: allow gmx to guess)
-g, --nb gpu        Calculate non-bonded interactions on gpu
-G, --gpu_id <str>  List ID(s) of unique GPU devices available for use
-p, --deffnm <str>  Define filename prefix (default is "md_filename" for outputs)
-a, --auto          Automatic mode (less prompts): Use default parameters &
                    do common analyses
--parFile <str>     Name of the CHAPERONg input parameter file
-H, --Help          Print more, advanced options
-s, --water <str>   Water model to use i.e. tip3p, spc, etc. (ff-dependent)
-f, --ff <str>      Force-field e.g. charmm27, amber94, etc.
                    Enter "wd" if ff in working directory
-R, --ntmpi <int>   Number of thread-MPI ranks (default is 0: gmx guesses)
-m, --ntomp <int>   Number of OpenMP threads per MPI rank (default is 0: guess)
-P, --pname <str>   Name of the positive ion (default is NA)
-N, --nname <str>   Name of the negative ion (default is CL)
-c, --conc <int>    Set salt concentration (mol/L) for the system
-W, --maxwarn <int> Number of allowed warnings (default is 0)
-M, --mmgpath <str> Absolute path to gmx binary to use for g_mmpbsa
-E, --gmx_exe <str> Path to gmx to use for all gmx runs except g_mmpbsa
                    (Default is to use the gmx set in the environment)
-v, --version       Print the installed version of CHAPERONg
--mmFrame <int>     Number of frames to be extracted for g_mmpbsa calculations
--trFrac <int>      Fraction of trajectory to use for g_mmpbsa calculations
                    (e.g. enter 2 for 2nd half, 3 for last 3rd, etc.)
--dist <float>      Solute-box distance i.e. distance to box edge (default is 1.0)
--bg                Run production mdrun with nohup ("no hang up")
--inputtraj <str>   Corrected trajectory to use for analyses (noPBC, nojump, fit)
--ter <prompt>      Interactively choose the N- & C-termini protonation states
                    (default: ionized with NH3+ & COO-)
guide_lg
}

#demA=$'\n\n'"#**********************************CHAPERONg**********************************#"$'\n'
demA=$'\n\n'"#================================= CHAPERONg =================================#"$'\n'
demB=$'\n'"#=============================================================================#"$'\n\n'
#demB=$'\n'"#*****************************************************************************#"$'\n\n'


#function to check if the parFile flag is used and then read the provided CHAPERONg parFile
read_parFile()
{
	if [[ "$parfilename" != '' ]] ; then 
		while IFS= read -r line; do
			par=$(echo "$line" | awk '{print $1}')
			if [[ "$par" == "mmgpath" ]] ; then
				mmGMXpath=$(echo "$line" | awk '{print $3}')
			elif [[ "$par" == "gmx_exe" ]] ; then
				gmx_exe_path=$(echo "$line" | awk '{print $3}')
			fi
		done < "$parfilename"
	fi
}

#then check other flags
#flags provided on the terminal overwrite parameters in parFile in case of conflicts
while [ "$1" != "" ]; do	
	case "$1" in
	--parFile) shift; parfilename="$1"; read_parFile;;	
	-a | --auto) flw=1;;
	-b | --bt) shift; btype="$1";;
	-c | --conc) shift; ion_conc="$1";;
	--dist) shift; edgeDist="$1";;
	--bg) nohp=1;;
	-p | --deffnm) shift; filenm="$1";;
	-E | --gmx_exe) shift; gmx_exe_path="$1";;
	-f | --ff) shift; ffUse="$1";;
	-F | --mmFrame) shift; mmpbframesNo="$1";;
	-g | --nb) nb=1;;
	-G | --gpu_id) shift; gpid="$1";;
	-h | --help) Help; Credit; exit 0;;
	-H | --Help) HHelp; Credit; exit 0;;
	-i | --input) shift; coordinates_raw="$1";;
	--inputtraj) shift; PBCcorrectType="$1";;
	-m | --ntomp) shift; ntomp="$1" ;;
	-M | --mmgpath) shift; mmGMXpath="$1"; mmGMX="1";;
	-N | --nname) shift; nn="$1";;
	-P | --pname) shift; pn="$1";;
	--parFile) shift; parfilename="$1";;	
	-R | --ntmpi) shift; ntmpi="$1";;
	-s | --water) shift; wat="$1";;
	--trFrac) shift; trajFraction="$1";;	
	-T | --nt) shift; nt="$1";;
	--ter) termini=1;;
	-v | --version) echo "$demA"$' CHAPERON version: '"$CHAPERONg_version"; Credit; echo $''; exit 0 ;;
	-W | --maxwarn) shift; WarnMax="$1";;
	*) echo "Invalid option: $1"; Help; echo $''; exit 1;;
	esac
	shift
done


if [[ "$edgeDist" != *"0."* && "$edgeDist" != *"1."* && "$edgeDist" != *"2."* && "$edgeDist" != *"3."* && \
		"$edgeDist" != *"4."* && "$edgeDist" != *"5."* && "$edgeDist" != *"6."* && "$edgeDist" != *"7."* && \
		"$edgeDist" != *"8."* && "$edgeDist" != *"9."* && "$edgeDist" != *"."*"0"* && "$edgeDist" != *".0"* ]];then
	echo $'\n'"$demA"$'\nThe solute-box distance i.e. minimum distance to the edge of the box you entered: '"$edgeDist"
	echo $'Please enter a valid number !!\n'
	exit 1
fi


if [[ "$#" == 1 ]] || [[ "$#" == 2 ]] && [[ "$flag" != "h" ]] && [[ "$flag" != "H" ]]; then
	echo "$demA"" No arguments are given. Default parameters will be used...""$demB"
	sleep 1
fi

if [[ $coordinates_raw == *".pdb" ]]; then
	coordinates=$(basename "$coordinates_raw" .pdb)
	echo "$demA"$' Your input coordinate filename has the extension ".pdb"\n The corrected filename is "'"$coordinates"'"'"$demB"
	sleep 2
elif [[ $coordinates_raw == *".gro" ]]; then
	coordinates=$(basename "$coordinates_raw" .gro)
	echo "$demA"$' Your input coordinate filename has the extension ".gro"\n The corrected filename is "'"$coordinates"'"'"$demB"
	sleep 2
else coordinates=$coordinates_raw
fi

if test "$wat" != ""; then wmodel="-water ""${wat}"
# elif test "$wat" == ""; then
# 	echo "$demA"" No water model is provided. You maybe be prompted to choose later.$demB"
# 	sleep 1
fi


if [[ "${filenm}" == '' ]]; then filenm="md_${coordinates}"; fi

if [[ "$gpid" != '' ]] && [[ "$nb" == '' ]] ; then gpidn="-gpu_id $gpid"
elif [[ "$gpid" != '' ]] && [[ "$nb" == 1 ]] ; then gpidn="-gpu_id $gpid -nb gpu"
elif [[ "$gpid" == '' ]] && [[ "$nb" == 1 ]] ; then gpidn="-nb gpu"
elif [[ "$gpid" == '' ]] && [[ "$nb" == '' ]] ; then gpidn=''
fi

if [[ "$termini" == 1 ]]; then extr="-ignh -ter"; fi

if [[ "$pn" != '' ]] || [[ "$nn" != '' ]] && [[ "$ion_conc" == '' ]]; then
	pnam_nnam="-pname $pn -nname $nn"
elif [[ "$pn" != '' ]] || [[ "$nn" != '' ]] && [[ "$ion_conc" != '' ]]; then
	pnam_nnam="-pname $pn -nname $nn -conc ${ion_conc}"
elif [[ "$pn" == '' ]] || [[ "$nn" == '' ]] && [[ "$ion_conc" == '' ]]; then
	pnam_nnam=''
elif [[ "$pn" == '' ]] || [[ "$nn" == '' ]] && [[ "$ion_conc" != '' ]]; then
	pnam_nnam="-conc ${ion_conc}"
fi

if [[ "$nt" == 0 ]] && [[ "$ntmpi" == 0 ]] && [[ "$ntomp" == 0 ]]; then threader='' && THREA=''
elif [[ "$nt" != 0 ]] && [[ "$ntmpi" != 0 ]] && [[ "$ntomp" == 0 ]]; then
	threader="-ntmpi ""${ntmpi}" && THREA="-nt ""${nt}" && hbthread="-nthreads ""${nt}"
elif [[ "$nt" == 0 ]] && [[ "$ntmpi" != 0 ]] && [[ "$ntomp" != 0 ]]; then
	threader="-ntmpi ""${ntmpi}"" -ntomp ""${ntomp}" && THREA='' && hbthread="-nthreads ""${ntomp}"
elif [[ "$nt" != 0 ]] && [[ "$ntmpi" != 0 ]] && [[ "$ntomp" != 0 ]]; then
	threader="-ntmpi ""${ntmpi}"" -ntomp ""${ntomp}" && THREA=''	&& hbthread="-nthreads ""${ntomp}"
elif [[ "$nt" != 0 ]] && [[ "$ntmpi" == 0 ]] && [[ "$ntomp" == 0 ]]; then
	threader='' && THREA="-nt ""${nt}" && hbthread="-nthreads ""${nt}"
fi


echo \
  $'\n###############################################################################'\
  $'\n#--------------------------------- CHAPERONg ---------------------------------#'\
  $'\n#   An automated pipeline for GROMACS MD simulation and trajectory analyses   #'\
  $'\n#    If you use this program in your work, please cite the relevant paper:    #'\
  $'\n#                    Yekeen, A.A. et al. To be published...                   #'\
  $'\n###############################################################################'

sleep 2
# cat << usageSt

# #-----------------------------------------------------------------------------#
# ######## ======================== BASIC USAGE ======================== ########
# #-----------------------------------------------------------------------------#

# chmod +x ./run_CHAPERONg-<version>
# ./run_CHAPERONg-<version> -i inputStructure_filename [-More options]

# #-----------------------------------------------------------------------------#
# ######## ========================= IMPORTANT ========================= ########
# #-----------------------------------------------------------------------------#
#  MAKE SURE the following are in the current working directory:
#   (1) Input structure (.pdb or .gro)
#   (2) mdp files (named as minim.mdp/em.mdp, nvt.mdp, npt.md, ions.mdp, md.mdp)
#     *PLEASE READ THE HIGHLIGHTED NOTES PRINTED ON THE TERMINAL DURING RUNS!!
#        *THIS WAY, YOU WON'T MISS ANY INFO YOU MAY FIND IMPORTANT... ENJOY!
# #-----------------------------------------------------------------------------#

# usageSt

# sleep 2

if [[ "$PBCcorrectType" != '' && "$PBCcorrectType" == 'noPBC' ]] ; then touch pbcmol
elif [[ "$PBCcorrectType" != '' && "$PBCcorrectType" == 'nojump' ]] ; then touch pbcjump
elif [[ "$PBCcorrectType" != '' && "$PBCcorrectType" == 'fit' ]] ; then touch pbcfit
elif [[ "$PBCcorrectType" != '' && "$PBCcorrectType" == 'center' ]] ; then touch pbccenter
fi