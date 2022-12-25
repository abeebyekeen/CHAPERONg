

#install_CHAPERONg - The automated installation and preparatory script for CHAPERONg
#CHAPERONg - An automation program for GROMACS md simulation
#Author: Abeeb A. Yekeen
#Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date: 2022.10.11  # Modified: 2022.10.26

set -e
set -o pipefail


# conda config --add channels conda-forge

# detect the path to the anaconda base environment
baseENV=$(conda env list | grep "base" | awk '{print $3}')
conda config --append envs_dirs "${baseENV}/envs"

# create the chaperong conda environment
conda env create --prefix "${baseENV}/envs/CHAPERONg" --file ./conda_dependencies.yml