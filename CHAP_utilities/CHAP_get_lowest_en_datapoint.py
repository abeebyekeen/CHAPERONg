
##########################################################################
#  CHAP_get_lowest_en_datapoint.py -- A python script to get the approx. #
#    datapoint for the lowest energy structure                           #
#  CHAP_get_lowest_en_datapoint.py is part of the CHAPERONg package      #
#  Input parameters are generated by other scripts in CHAPERONg and are  #
#    read by this script                                                 #
#  CHAPERONg -- An automation program for GROMACS MD simulations and     #
#    trajectory analyses                                                 #
##########################################################################

__author__  = 'Abeeb A. Yekeen'
__email__   = 'contact@abeebyekeen.com'
__date__    = '2022.10.17'
__version__ = '1.0'
__status__  = 'Production'

import os

os.chdir('collect_mappings')
with open("EnergyMinim.txt") as paramet:
    Energyminim = paramet.readline().strip("\n")
    Energyminimsplt = Energyminim.split("\t")
with open("sorted_1.txt", "r") as dataPts:
    prevdatapt = "nil"
    fesMin = 0
    lineno = 0
    for dataP in dataPts.readlines():
        lineno+=1
        # print(lineno)
        # print(dataP)
        dataPstrip = dataP.strip("\n")
        dataPstripsplt = dataPstrip.split("\t")
        if fesMin == 0:
            if str(dataPstrip) == str(Energyminim):
                fesMin+=1
                # print("Found minimum: "+str(dataPstrip)+"\n")
                # print(lineno)
                continue
            elif str(dataPstrip) != str(Energyminim):
                prevdatapt = dataPstrip
        elif fesMin > 0:
            if prevdatapt != "nil":
                # print(lineno)
                prevdataptsplt = prevdatapt.split("\t")
                # print(prevdataptsplt)
                with open("lowest_energy_datapoints_timed.dat", "w") as lowEn:
                    if (abs(float(prevdataptsplt[2]) - float(Energyminimsplt[2])) < 
                        abs(float(dataPstripsplt[2]) - float(Energyminimsplt[2]))):
                        lowEn.write("Time"+"\t"+"OrderPar1"+"\t"+"OrderPar2"+"\t"+"Energy"+"\n")
                        lowEn.write(prevdatapt)
                    elif (abs(float(prevdataptsplt[2]) - float(Energyminimsplt[2])) > 
                        abs(float(dataPstripsplt[2]) - float(Energyminimsplt[2]))):
                        lowEn.write("Time"+"\t"+"OrderPar1"+"\t"+"OrderPar2"+"\t"+"Energy"+"\n")
                        lowEn.write(dataPstrip)
                os.chdir('../')
                break
            elif prevdatapt == "nil":
                with open("lowest_energy_datapoints_timed.dat", "w") as lowEn:
                    lowEn.write("Time"+"\t"+"OrderPar1"+"\t"+"OrderPar2"+"\t"+"Energy"+"\n")
                    lowEn.write(dataPstrip)
                os.chdir('../')
                break
