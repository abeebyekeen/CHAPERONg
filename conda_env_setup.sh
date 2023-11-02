

#install_CHAPERONg.sh - The installation and preparatory script for CHAPERONg
#CHAPERONg - An automation program for GROMACS md simulation
#Author: Abeeb A. Yekeen
#Contact: contact@abeebyekeen.com
#Date: 2022.10.11  # Modified: 2022.10.26

###############################################################################
# conda_env_setup.sh -- The conda environment setup script for CHAPERONg      #
# CHAPERONg -- An automation program for GROMACS md simulation and trajectory #
#    analysis                                                                 #
# Author -- Abeeb A. Yekeen                                                   #
# Contact -- contact@abeebyekeen.com                                          #
# Date -- 2022.02.11                                                          #
###############################################################################

set -e
set -o pipefail


# conda config --add channels conda-forge

# detect the path to the anaconda base environment
baseENV=$(conda env list | grep "base" | awk '{print $3}')
conda config --append envs_dirs "${baseENV}/envs"

# create the chaperong conda environment
conda env create --prefix "${baseENV}/envs/chaperong" --file ./conda_dependencies.yml