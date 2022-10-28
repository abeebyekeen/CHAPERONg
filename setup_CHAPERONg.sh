#! /bin/bash

#setup_CHAPERONg - The preparatory script for setting up CHAPERONg
#CHAPERONg - An automation program for GROMACS md simulation and trajectory analyses
#Author: Abeeb A. Yekeen
#Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date: 2022.02.11

#set version
# CHAPERONg_version="beta3"

chmod a+x "$CHAPERONg_PATH/CHAP_modules/CHAP_colPar.sh" "$CHAPERONg_PATH/CHAP_modules/CHAP_deffxn.sh" \
"$CHAPERONg_PATH/CHAP_modules/CHAP_ana.sh" "$CHAPERONg_PATH/CHAP_modules/CHAP_sim.sh" \
"$CHAPERONg_PATH/CHAP_modules/run_CHAPERONg.sh" "$CHAPERONg_PATH/CHAP_utilities/g_mmpbsa_pkg/g_mmpbsa" \
"$CHAPERONg_PATH/CHAP_utilities/g_mmpbsa_pkg/energy2bfac" "$CHAPERONg_PATH/CHAP_utilities/dssp-x64"

cp "$CHAPERONg_PATH/CHAP_modules/run_CHAPERONg.sh" "$CHAPERONg_PATH/"