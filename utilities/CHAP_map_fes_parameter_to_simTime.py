
# CHAP_map_fes_parameter_to_simTime.py -- A python script to map fes data point and
#   free energy to simulation time and order parameter in the md trajectory
# Part of the CHAPERONg suite of scripts
# Input parameters are generated by other scripts in CHAPERONg and
#   are read by this script
# CHAPERONg -- An automation program for GROMACS md simulation
# Author -- Abeeb A. Yekeen
# Contact -- yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
# Date: 2022.10.17

import os

with open("OrderParameters1_2_dG_nogap-sorted.dat", "r") as paramet:
    Energymin = paramet.readline().strip("\n")
    filecount=1
    spltEnMinData = Energymin.split("\t")
    # print(spltEnMinData)
    Par1 = spltEnMinData[0]
    Par1float = float(Par1)
    Par1float_decimal = float("{:.7f}".format(Par1float))
    Par2 = spltEnMinData[1]
    Par2float = float(Par2)
    Par2float_decimal = float("{:.7f}".format(Par2float))
    apprxPar1 = float(Par1[:4])
    # print(apprxPar1)
    apprxPar2 = float(Par2[:4])
    # print(apprxPar2)
    with open("SimTime_OrderParameters1_2.dat", "r") as timedparamet:
        timedPairLines = timedparamet.readlines()
        linecount = 0
        for timeddatapair in timedPairLines:
            # print(timeddatapair)
            spltimeddatapair = timeddatapair.rstrip("\n").split("\t")
            if (float(spltimeddatapair[1]) <= (apprxPar1+0.009)) and (float(spltimeddatapair[1]) >= (apprxPar1-0.009)) and \
                (float(spltimeddatapair[2]) <= (apprxPar2+0.005)) and (float(spltimeddatapair[2]) >= (apprxPar2-0.005)):
                os.chdir('collect_mappings')
                readinFile = str(filecount)+".txt"
                if linecount == 0:
                    with open(readinFile , "w") as mapping:
                        mapping.write("Time"+"\t"+"OrderPar1"+"\t"+"OrderPar2"+"\t"+"Energy"+"\n")
                        fesdat = str("----")+"\t"+str(Par1float_decimal)+"\t"+str(Par2float_decimal)+"\t"+str(spltEnMinData[2])+"\n"
                        mapping.write(fesdat)
                        mapping.write(str((timeddatapair).rstrip())+"\t"+str(spltEnMinData[2])+"\n")
                        linecount+=1
                        os.chdir('../')
                    # continue
                elif linecount > 0:
                    with open(readinFile , "a") as mapping:
                        mapping.write(str((timeddatapair).rstrip())+"\t"+str(spltEnMinData[2])+"\n")
                        linecount+=1
                        os.chdir('../')
        filecount+=1
