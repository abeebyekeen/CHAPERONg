#! /bin/bash

#install_CHAPERONg - The automated installation and preparatory script for CHAPERONg
#CHAPERONg - An automation program for GROMACS md simulation
#Author: Abeeb A. Yekeen
#Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date: 2022.10.11  # Modified: 2022.10.26

set -e
set -o pipefail

# get user's consent for setting the environment

cat << seekConsent

 This action will allow CHAPERONg to safely modify your environment for installation.
 This is a one-time process. For details, visit:
 https://abeebyekeen.com/resources/chaperong/installation

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
CHAPERONg_version="beta3"

chmod a+x "$CHAPERONg_PATH/CHAP_modules/CHAP_colPar.sh" "$CHAPERONg_PATH/CHAP_modules/CHAP_deffxn.sh" \
"$CHAPERONg_PATH/CHAP_modules/CHAP_ana.sh" "$CHAPERONg_PATH/CHAP_modules/CHAP_sim.sh" \
"$CHAPERONg_PATH/conda_env_setup.sh" "$CHAPERONg_PATH/conda_dependencies.yml" \
"$CHAPERONg_PATH/CHAP_modules/run_CHAPERONg.sh" "$CHAPERONg_PATH/CHAP_utilities/g_mmpbsa_pkg/g_mmpbsa" \
"$CHAPERONg_PATH/CHAP_utilities/g_mmpbsa_pkg/energy2bfac" "$CHAPERONg_PATH/CHAP_utilities/dssp-x64"

cp "$CHAPERONg_PATH/CHAP_modules/run_CHAPERONg.sh" "$CHAPERONg_PATH"

echo \
  $'\n###############################################################################'\
  $'\n#--------------------------------- CHAPERONg ---------------------------------#'\
  $'\n#   An automated pipeline for GROMACS MD simulation and trajectory analyses   #'\
  $'\n#    If you use this program in your work, please cite the relevant paper:    #'\
  $'\n#                    Yekeen, A.A. et al. To be published...                   #'\
  $'\n###############################################################################'

sleep 4
cat << usageSt

######## ======================== BASIC USAGE ======================== ########   
#-----------------------------------------------------------------------------#
run_CHAPERONg -i inputStructure_filename [-More options]
#-----------------------------------------------------------------------------#

######## ========================= IMPORTANT ========================= ########
#-----------------------------------------------------------------------------#
     *DO READ THE HIGHLIGHTED NOTES PRINTED ON THE TERMINAL DURING RUNS!!          
     **THIS WAY, YOU WON'T MISS ANY INFO YOU MAY FIND IMPORTANT... ENJOY!
#-----------------------------------------------------------------------------#

usageSt
sleep 3
. ~/.bashrc