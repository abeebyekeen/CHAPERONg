#! /bin/bash

######################################################################
#  CHAP_colPar -- The parameter collection module of CHAPERONg       #
#  CHAPERONg -- An automation program for GROMACS md simulation and  #
#    trajectory analysis                                             #
#  Author -- Abeeb A. Yekeen                                         #
#  Contact -- abeeb.yekeen@hotmail.com                               #
#  Date -- 2022.02.11                                                #
######################################################################

set -e
set -o pipefail

#set version
CHAPERONg_version="v0.1"


#Defining primary functions
# Credit()
# {
# 	echo \
#   $'\n###############################################################################'\
#   $'\n#--------------------------------- CHAPERONg ---------------------------------#'\
#   $'\n#   An automated pipeline for GROMACS MD simulation and trajectory analyses   #'\
#   $'\n#    If you use this program in your work, please cite the relevant paper:    #'\
#   $'\n#                   Yekeen, A.A. et al. To be published...                    #'\
#   $'\n###############################################################################'
# }

Credit()
{
	echo -e \
  '\033[92m'\
  '\n###############################################################################'\
  '\n#\033[5m--------------------------------- CHAPERONg ---------------------------------\033[25m#'\
  '\n#   An automated pipeline for GROMACS MD simulation and trajectory analyses   #'\
  '\n#   \033[92;7m'\
  'If you use this program in your work, please cite the relevant paper: \033[m '\
  '\033[92m #'\
  '\n#                  \033[92;7m Yekeen, A.A. et al. To be published... \033[m '\
  '\033[92m                 #'\
  '\n###############################################################################'\
  '\033[m'
}

Help()
{

#Shows the use of the script
cat << guide_sh

Usage:
chmod +x ./setup_CHAPERONg-<version>
./run_CHAPERONg-<version> -i inputStructure_filename

Required (int=integer; str=string):
-i, --input <str>    Input coordinate file (.pdb or .gro)
Optional (int=integer; str=string):
-h, --help           Print this help
-b, --bt <str>       Box type: cubic (default), dodecahedron, triclinic, etc.
-T, --nt <int>       Number of threads to use [default: 0 (gmx guesses)]
-g, --nb gpu         Calculate non-bonded interactions on gpu
-G, --gpu_id <str>   List ID(s) of unique GPU device(s) available for use
-p, --deffnm <str>   Set filename prefix (default for outputs: "md_filename")
-a, --auto_mode      Automation mode [options: full, semi(default)]. full: Use
                     default parameters & do common analyses (less prompts)
--paraFile <str>     Name of the CHAPERONg input parameter file
-H, --Help           Print more, advanced options
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
-i, --input <str>    Input coordinate file (.pdb or .gro)
Optional (int=integer; str=string):
-h, --help           Print the shorter version of this help
-b, --bt <str>       Box type: cubic (default), dodecahedron, triclinic, etc.
-T, --nt <int>       Number of threads to use [default: 0 (gmx guesses)]
-g, --nb gpu         Calculate non-bonded interactions on gpu
-G, --gpu_id <str>   List ID(s) of unique GPU devices available for use
-p, --deffnm <str>   Set filename prefix (default for outputs: "md_filename")
-a, --auto_mode      Automation mode [options: full, semi(default)]. full: Use
                     default parameters & do common analyses (less prompts)
-H, --Help           Print more, advanced options
-s, --water <str>    Water model: tip3p, spc, spce, etc. (ff-dependent)
-f, --ff <str>       Force-field: charmm27, amber94, oplsaa, gromos54a7, etc.
                     (Enter "wd" if forcefield is in working directory)
-P, --posname <str>  Name of the positive ion (default: NA)
-N, --negname <str>  Name of the negative ion (default: CL)
-c, --conc <int>     Set salt concentration (mol/L) for the system
-W, --maxwarn <int>  Number of allowed warnings (default is 0)
-M, --mmgpath <str>  Absolute path to gmx binary to use for g_mmpbsa
-E, --gmx_exe <str>  Path to gmx to use for all gmx runs except g_mmpbsa
                     (Default is to use the gmx set in the environment)
-v, --version        Print the installed version of CHAPERONg
-t, --temp <int>     Simulation temperature in kelvin
--ntmpi <int>        Number of thread-MPI ranks [default: 0 (gmx guesses)]
--ntomp <int>        Number of OpenMP threads per MPI rank; default: 0 (guess)
--paraFile <str>     Name of the CHAPERONg input parameter file
--inputtraj <str>    Corrected trajectory to generate and use for analyses
                     (options: noPBC, nojump, center, fit, combo)
--clustr_cut <float> RMSD cut-off (nm) for cluster membership (default: 1.0)
--clustr_methd <str> Method for cluster determination: gromos (default),
                     linkage, jarvis-patrick, monte-carlo, diagonalization
--frame_beginT <int> Time (ps) of first frame to read from trajectory
--frame_endT <int>   Time (ps) of last frame to read from trajectory
--dt <int>           Interval (ps) at which frames are taken from trajectory
--mmFrame <int>      Number of frames to be extracted for g_mmpbsa
--movieFrame <int>   Number of frames to extract and use for movie
--trFrac <int>       Fraction of trajectory to use for g_mmpbsa
                     (enter 1 for all, 2 for 2nd half, 3 for last 3rd, etc.)
--kde_opt <int>      Range (above and below the estimate) to test for the  
                     optimization of histogram number of bins for KDE
--path_av_plot<str>  Path to input files for average of replica plots
--dist <float>       Solute-box distance (distance to box edge; default: 1.0)
--bg                 Run production mdrun in the background with "nohup"
--ter <prompt>       Interactively choose the N- & C-termini protonation 
                     states (default: ionized with NH3+ & COO-)
guide_lg
}

demA=$'\n\n'"#================================= CHAPERONg =================================#"$'\n'
demB=$'\n'"#=============================================================================#"$'\n\n'


# Initialize (default) parameters
btype='cubic' ; edgeDist="1.0"; WarnMax=0; nt=0
wat=""; nb='' ; termini=0 ; gmx_exe_path="gmx"
automode="semi"; ffUse=""; gpid=''; filenm=''; Temp=300
ntmpi=0; ntomp=0; skp=''; nohp='' ; mmpbframesNo=''
pn=''; nn=''; ion_conc='' ; customframeNo=''
PBCcorrectType='' ; trajFraction='' ; dt=1
mmGMX='' ; mmGMXpath='' ; coordinates_raw=''
parfilename='' ; frame_b=0 ; frame_e=0
method_clust='gromos' ; cut_cl='0.1'
bin_number_range='' ; customNDXask=''
mmpb_begin='' ; path_av=''
#gmxV=''

# check if the paraFile flag is used and then read the provided parameter file
read_paraFile()
{
	if [[ "$parfilename" != '' ]] ; then 
		while IFS= read -r line; do
			par=$(echo "$line" | awk '{print $1}')
			par_input=$(echo "$line" | awk '{print $3}')
			if [[ "$par" == "input" ]]; then coordinates_raw="$par_input"
			elif [[ "$par" == "bt" ]]; then btype="$par_input"
			elif [[ "$par" == "nt" ]]; then nt="$par_input"
			elif [[ "$par" == "nb" && "$par_input" == "gpu" ]]; then nb=1
			elif [[ "$par" == "gpu_id" ]]; then gpid="$par_input"
			elif [[ "$par" == "deffnm" ]]; then filenm="$par_input"
			elif [[ "$par" == "water" ]]; then wat="$par_input"
			elif [[ "$par" == "ff" ]]; then ffUse="$par_input"
			elif [[ "$par" == "ntmpi" ]]; then ntmpi="$par_input"
			elif [[ "$par" == "ntomp" ]]; then ntomp="$par_input"
			elif [[ "$par" == "mmgpath" ]]; then mmGMX="1"; mmGMXpath="$par_input"
			elif [[ "$par" == "frame_beginT" ]]; then frame_b="$par_input"
			elif [[ "$par" == "frame_endT" ]]; then frame_e="$par_input"
			elif [[ "$par" == "movieFrame" ]]; then customframeNo="$par_input"
			elif [[ "$par" == "posname" ]]; then pn="$par_input"
			elif [[ "$par" == "negname" ]]; then nn="$par_input"
			elif [[ "$par" == "conc" ]]; then ion_conc="$par_input"
			elif [[ "$par" == "temp" ]]; then Temp="$par_input"
			elif [[ "$par" == "maxwarn" ]]; then WarnMax="$par_input"
			elif [[ "$par" == "dist" ]]; then edgeDist="$par_input"
			elif [[ "$par" == "inputtraj" ]]; then PBCcorrectType="$par_input"
			elif [[ "$par" == "trFrac" ]]; then trajFraction="$par_input"
			elif [[ "$par" == "mmFrame" ]]; then mmpbframesNo="$par_input"
			elif [[ "$par" == "mmBegin" ]]; then mmpb_begin="$par_input"
			elif [[ "$par" == "gmx_exe" ]]; then gmx_exe_path="$par_input"
			elif [[ "$par" == "clustr_methd" ]]; then method_clust="$par_input"
			elif [[ "$par" == "clustr_cut" ]]; then cut_cl="$par_input"
			elif [[ "$par" == "dt" ]]; then dt="$par_input"
			elif [[ "$par" == "auto_mode" ]]; then automode="$par_input"
			elif [[ "$par" == "path_av_plot"]]; then path_av="$part_input"
			fi
		done < "$parfilename"
	fi
}

# then check other flags
# flags provided on the terminal overwrite parameters in paraFile in case of conflicts
while [ "$1" != "" ]; do	
	case "$1" in
	--paraFile) shift; parfilename="$1"; read_paraFile;;	
	-a | --auto_mode) shift; automode="$1";;
	-b | --bt) shift; btype="$1";;
	-c | --conc) shift; ion_conc="$1";;
	--clustr_cut) shift; cut_cl="$1";;
	--clustr_methd) shift; method_clust="$1";;
	--dist) shift; edgeDist="$1";;
	--dt) shift; dt="$1";;
	--bg) nohp=1;;
	-p | --deffnm) shift; filenm="$1";;
	-E | --gmx_exe) shift; gmx_exe_path="$1";;
	-f | --ff) shift; ffUse="$1";;
	-F | --mmFrame) shift; mmpbframesNo="$1";;
	--mmBegin) shift; mmpb_begin="$1";;
	--frame_beginT) shift; frame_b="$1";;
	--frame_endT) shift; frame_e="$1";;
	-g | --nb) nb=1;;
	-G | --gpu_id) shift; gpid="$1";;
	-h | --help) Help; Credit; exit 0;;
	-H | --Help) HHelp; Credit; exit 0;;
	-i | --input) shift; coordinates_raw="$1";;
	--inputtraj) shift; PBCcorrectType="$1";;
	--kde_opt) shift; bin_number_range="$1";;
	--ntomp) shift; ntomp="$1" ;;
	--movieFrame) shift; customframeNo="$1" ;;
	-M | --mmgpath) shift; mmGMXpath="$1"; mmGMX="1";;
	-N | --negname) shift; nn="$1";;
	-P | --posname) shift; pn="$1";;
	--path_av_plot) shift; path_av="$1";;
	--paraFile) shift; parfilename="$1";;	
	--ntmpi) shift; ntmpi="$1";;
	-s | --water) shift; wat="$1";;
	-T | --nt) shift; nt="$1";;
	-t | --temp) shift; Temp="$1";;
	--ter) termini=1;;
	--trFrac) shift; trajFraction="$1";;
	-v | --version) echo "$demA"$' CHAPERON version: '"$CHAPERONg_version"; Credit; echo $''; exit 0 ;;
	-W | --maxwarn) shift; WarnMax="$1";;
	*) echo "Invalid option: $1"; Help; echo $''; exit 1;;
	esac
	shift
done

pattern="^[0-9]+(\.[0-9]+)?$"

if ! [[ "$edgeDist" =~ $pattern ]]; then
  echo -e "\n$demA\nThe solute-box distance i.e. minimum distance to the edge of the box you entered: $edgeDist"
  echo "Please enter a valid number !!\n"
  exit 1
fi

# if [[ "$edgeDist" != *"0."* && "$edgeDist" != *"1."* && "$edgeDist" != *"2."* && "$edgeDist" != *"3."* && \
# 		"$edgeDist" != *"4."* && "$edgeDist" != *"5."* && "$edgeDist" != *"6."* && "$edgeDist" != *"7."* && \
# 		"$edgeDist" != *"8."* && "$edgeDist" != *"9."* && "$edgeDist" != *"."*"0"* && "$edgeDist" != *".0"* ]];then
# 	echo $'\n'"$demA"$'\nThe solute-box distance i.e. minimum distance to the edge of the box you entered: '"$edgeDist"
# 	echo $'Please enter a valid number !!\n'
# 	exit 1
# fi

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
elif test "$wat" == ""; then wmodel=""
# elif test "$wat" == ""; then
# 	echo "$demA"" No water model is provided. You maybe be prompted to choose later.$demB"
# 	sleep 1
fi

if [[ "${filenm}" == '' ]]; then filenm="md_${coordinates}"; fi

if [[ "$gpid" != '' && "$nb" == '' ]] ; then gpidn="-gpu_id $gpid"
elif [[ "$gpid" != '' && "$nb" == 1 ]] ; then gpidn="-gpu_id $gpid -nb gpu"
elif [[ "$gpid" == '' && "$nb" == 1 ]] ; then gpidn="-nb gpu"
elif [[ "$gpid" == '' && "$nb" == '' ]] ; then gpidn=''
fi

if [[ $termini == 1 ]]; then extr="-ignh -ter"
elif [[ $termini == 0 ]]; then extr="-ignh"
fi

pnam_nnam=''

if [[ "$pn" != '' || "$nn" != '' ]] && [[ "$ion_conc" == '' ]]; then
	pnam_nnam="-pname $pn -nname $nn"
elif [[ "$pn" != '' || "$nn" != '' ]] && [[ "$ion_conc" != '' ]]; then
	pnam_nnam="-pname $pn -nname $nn -conc ${ion_conc}"
elif [[ "$pn" == '' || "$nn" == '' ]] && [[ "$ion_conc" == '' ]]; then
	pnam_nnam=''
elif [[ "$pn" == '' || "$nn" == '' ]] && [[ "$ion_conc" != '' ]]; then
	pnam_nnam="-conc ${ion_conc}"
fi


THREA="-nt ""${nt}"; hbthread="-nthreads ""0"
threader="-ntmpi ""${ntmpi}"" -ntomp ""${ntomp}"

if [[ "$nt" == 0 && "$ntmpi" == 0 && "$ntomp" == 0 ]]; then threader='' && THREA=''
elif [[ "$nt" != 0 ]] && [[ "$ntmpi" != 0 ]] && [[ "$ntomp" == 0 ]]; then
	threader="-ntmpi ""${ntmpi}" && THREA="-nt ""${nt}" && hbthread="-nthreads ""${nt}"
elif [[ "$nt" == 0 ]] && [[ "$ntmpi" != 0 ]] && [[ "$ntomp" != 0 ]]; then
	threader="-ntmpi ""${ntmpi}"" -ntomp ""${ntomp}" && THREA='' && hbthread="-nthreads ""${ntomp}"
elif [[ "$nt" != 0 ]] && [[ "$ntmpi" != 0 ]] && [[ "$ntomp" != 0 ]]; then
	threader="-ntmpi ""${ntmpi}"" -ntomp ""${ntomp}" && THREA=''	&& hbthread="-nthreads ""${ntomp}"
elif [[ "$nt" != 0 ]] && [[ "$ntmpi" == 0 ]] && [[ "$ntomp" == 0 ]]; then
	threader='' && THREA="-nt ""${nt}" && hbthread="-nthreads ""${nt}"
fi


	# '\033[92m'\
	# '\n###############################################################################'\
	# '\n#--------------------------------- CHAPERONg ---------------------------------#'\
	# '\n#   An automated pipeline for GROMACS MD simulation and trajectory analyses   #'\
	# '\n#    If you use this program in your work, please cite the relevant paper:    #'\
	# '\n#                   Yekeen, A.A. et al. To be published...                    #'\
	# '\n###############################################################################'\
	# '\033[m'

echo -e \
'\033[92m'\
'\n###############################################################################'\
'\n#\033[5m--------------------------------- CHAPERONg ---------------------------------\033[25m#'\
'\n#   An automated pipeline for GROMACS MD simulation and trajectory analyses   #'\
'\n#   \033[92;7m'\
' If you use this program in your work, please cite the relevant paper: \033[m  '\
'\033[92m #'\
'\n#                  \033[92;7m Yekeen, A.A. et al. To be published... \033[m '\
'\033[92m                  #'\
'\n###############################################################################'\
'\033[m'

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
elif [[ "$PBCcorrectType" != '' && "$PBCcorrectType" == 'combo' ]] ; then touch pbccombo
fi