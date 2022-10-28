

#install_CHAPERONg - The automated installation and preparatory script for CHAPERONg
#CHAPERONg - An automation program for GROMACS md simulation
#Author: Abeeb A. Yekeen
#Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date: 2022.10.11  # Modified: 2022.10.26

set -e
set -o pipefail

# conda deactivate
# conda deactivate
# conda activate base
# conda config --add channels conda-forge

baseENV=$(conda env list | grep "base" | awk '{print $3}')

conda env create --prefix "${baseENV}/CHAPERONg" --file ./conda_dependencies.yml