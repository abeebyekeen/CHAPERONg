#! /bin/bash


###############################################################################
# install_CHAPERONg.sh -- The installation & preparatory script for CHAPERONg #
# CHAPERONg -- An automation program for GROMACS md simulation and trajectory #
#    analysis                                                                 #
# Author -- Abeeb A. Yekeen                                                   #
# Contact -- contact@abeebyekeen.com                                          #
# Date -- 2022.02.11                                                          #
###############################################################################

set -e
set -o pipefail

# get user's consent for setting the environment

cat << seekConsent

 This action will allow CHAPERONg to safely modify your environment for installation.
 This is a one-time process. For details, visit:
 https://abeebyekeen.com/chaperong-online-documentation/#installation

seekConsent

read -p ' Do you agree to proceed with the installation? (y/n) ' prepCHAP
while [[ ${prepCHAP} != "y" && ${prepCHAP} != "yes" && ${prepCHAP} != "no" && ${prepCHAP} != "n" ]]; do
    echo $'\n You entered: '"$prepCHAP"
    echo $'\n Please enter a valid response (a "yes"/"y" or a "no"/"n" )!!\n'
    read -p ' Enter a response here: ' prepCHAP
done
if [[ "$prepCHAP" == "no" || "$prepCHAP" == "n" ]]
    then echo $'User quitted the installation!\nThanks' ; sleep 1 ; exit 0
fi

# add CHAPERONg installation path to the ~/.bashrc
CHAPERONg_PATH=$(pwd)
echo $'\n'"#=========== This section has been added by CHAPERONg upon installation ==============#"\
$'\n'"#================================ Delete to uninstall ================================#"\
$'\n'"export CHAPERONg_PATH=$(pwd)"$'\n'"export PATH=$(echo $CHAPERONg_PATH):"$'$PATH'\
$'\n'"#=====================================================================================#" >> ~/.bashrc

demA=$'\n\n'"#================================= CHAPERONg =================================#"$'\n'
demB=$'\n'"#=============================================================================#"$'\n\n'

echo "$demA"$' Installation completed'"$demB"	
sleep 2

# set version
# CHAPERONg_version="beta3"

cp "$CHAPERONg_PATH/CHAP_modules/run_CHAPERONg.sh" "$CHAPERONg_PATH"

chmod a+x "$CHAPERONg_PATH/CHAP_modules/CHAP_colPar.sh" "$CHAPERONg_PATH/CHAP_modules/CHAP_deffxn.sh" \
"$CHAPERONg_PATH/CHAP_modules/CHAP_ana.sh" "$CHAPERONg_PATH/CHAP_modules/CHAP_sim.sh" \
"$CHAPERONg_PATH/conda_env_setup.sh" "$CHAPERONg_PATH/conda_dependencies.yml" \
"$CHAPERONg_PATH/CHAP_modules/run_CHAPERONg.sh" "$CHAPERONg_PATH/CHAP_utilities/g_mmpbsa_pkg/g_mmpbsa" \
"$CHAPERONg_PATH/CHAP_utilities/g_mmpbsa_pkg/energy2bfac" "$CHAPERONg_PATH/CHAP_utilities/dssp-x64" \
"$CHAPERONg_PATH/run_CHAPERONg.sh"

  
echo -e \
  '\033[92m'\
  '\n###############################################################################'\
  '\n#\033[5m--------------------------------- CHAPERONg ---------------------------------\033[25m#'\
  '\n#              If you use this program in your work, please cite:             #'\
  '\n# \033[92;7m'\
  'Yekeen A.A., Durojaye O.A., Idris M.O., Muritala H.F., Arise R.O. (2023). \033[m'\
  '\033[92m#'\
  '\n#     \033[92;7m CHAPERONg: A tool for automated GROMACS-based molecular dynamics \033[m'\
  '\033[92m     #'\
  '\n#  \033[92;7m simulations  and  trajectory  analyses,  Computational  and  Structural \033[m'\
  '\033[92m #'\
  '\n#                   \033[92;7m Biotechnology Journal 21: 4849-4858. \033[m'\
  '\033[92m                   #'\
  '\n###############################################################################'\
  '\033[m'



sleep 4
cat << usageSt

######## ======================== BASIC USAGE ======================== ########   
#-----------------------------------------------------------------------------#
  run_CHAPERONg.sh -i inputStructure_filename [-More options]
#-----------------------------------------------------------------------------#

######## ========================= IMPORTANT ========================= ########
#-----------------------------------------------------------------------------#
     *DO READ THE HIGHLIGHTED NOTES PRINTED ON THE TERMINAL DURING RUNS!!          
     **THIS WAY, YOU WON'T MISS ANY INFO YOU MAY FIND IMPORTANT... ENJOY!
#-----------------------------------------------------------------------------#

usageSt
sleep 3
. ~/.bashrc