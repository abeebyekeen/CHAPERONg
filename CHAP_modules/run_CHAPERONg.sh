#! /bin/bash

#CHAPERONg - An automation program for GROMACS md simulation
#Author: Abeeb A. Yekeen
#Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
#Date: 2022.02.11

set -e
set -o pipefail

#call module to collect parameters
#. "$CHAPERONg_PATH/CHAP_modules/CHAP_colPar-4.01.sh"
initiator1=avail
#set CHAPERONg_version
# CHAPERONg_version="beta3"

#call module with defined fxns
. "$CHAPERONg_PATH/CHAP_modules/CHAP_deffxn.sh"
initiator2=avail

#call module for traj_ana
. "$CHAPERONg_PATH/CHAP_modules/CHAP_ana.sh"
initiator3=avail

#call module for mdsim
. "$CHAPERONg_PATH/CHAP_modules/CHAP_sim.sh"



