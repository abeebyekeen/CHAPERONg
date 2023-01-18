
# CHAP_construct_free_en_surface.py -- A python script to construct and
#   plot a free energy surface from a pair of ordered parameters
# Part of the CHAPERONg suite of scripts
# Input parameters are generated by other scripts in CHAPERONg and
#   are read by this script
# CHAPERONg -- An automation program for GROMACS md simulation
# Author -- Abeeb A. Yekeen
# Contact -- yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
# Date: 2022.10.16

import matplotlib
# Configure the non-interactive backend
matplotlib.use('Cairo')
from matplotlib import pyplot as plt
import math
import time
import numpy as np
import pandas
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Read in parameters for FES calculations
print (" Reading in parameters for FES calculations"+"\n")
with open("CHAP_fes_Par.in") as in_par:
	for parameter in in_par.readlines():
		if "minPar1" in parameter:
			para_data = str(parameter).split(",")
			para1_min = float(para_data[1])
		elif "maxPar1" in parameter:
			para_data = str(parameter).split(",")
			para1_max = float(para_data[1])
		elif "minPar2" in parameter:
			para_data = str(parameter).split(",")
			para2_min = float(para_data[1])
		elif "maxPar2" in parameter:
			para_data = str(parameter).split(",")
			para2_max = float(para_data[1])
		elif "XaxisL" in parameter:
			para_data = parameter.rstrip('\n').split(",")
			xaxis_label = str(para_data[1])
		elif "YaxisL" in parameter:
			para_data = parameter.rstrip('\n').split(",")
			yaxis_label = str(para_data[1])
		elif "no_of_frames" in parameter:
			para_data = str(parameter).split(",")
			no_of_frames = int(para_data[1])
		elif "Temp" in parameter:
			para_data = str(parameter).split(",")
			Temp = float(para_data[1])
		elif "outFilename" in parameter:
			para_data = parameter.rstrip('\n').split(",")
			out_file = str(para_data[1])+str(".png")
		elif "plotTitle" in parameter:
			para_data = parameter.rstrip('\n').split(",")
			plotTitle = str(para_data[1])+str(" Free Energy Surface")
		elif "x_bin_count" in parameter:
			para_data = parameter.rstrip('\n').split(",")
			xbin_custom = int(para_data[1])
		elif "y_bin_count" in parameter:
			para_data = parameter.rstrip('\n').split(",")
			ybin_custom = int(para_data[1])
time.sleep(2)

print (" Reading in data of order parameters"+"\n")
# Initialize the order parameter lists
order_p1 = []
order_p2 = []

# Read in data of order parameters
with open("OrderParameterPair.dat") as alldata:
	alldata_lines = alldata.readlines()
	for line in alldata_lines:
		data_point = str(line).split(",")
		order_p1.append(float(data_point[0]))
		order_p2.append(float(data_point[1]))
time.sleep(2)

if "PCA-derived" in plotTitle:
	print (" Binning and generating a 2D histogram"+"\n")
	time.sleep(2)
	# Range of data from order parameters
	para1_range = para1_max - para1_min
	para2_range = para2_max - para2_min

	# Preset number of bins for PCA data
	xbin = xbin_custom
	ybin = ybin_custom
	
		
	# Create a 2D histogram using the numpy histogram2d function
	hist, x_edges, y_edges = np.histogram2d(order_p1, order_p2, bins=(xbin, ybin), \
		range=[[para1_min, para1_max], [para2_min, para2_max]])
	time.sleep(2)

else:
	print (" Estimating the optimal number of bins"+"\n")
	time.sleep(2)
	# Determine the number of bins using the Freedman-Diaconis (1981) method
	dist_p1 = pandas.Series(order_p1)
	para1_max = dist_p1.max()
	para1_min = dist_p1.min()
	para1_range = para1_max - para1_min
	qt1 = dist_p1.quantile(0.25)
	qt3 = dist_p1.quantile(0.75)
	iqr1 = qt3 - qt1
	bin_width_p1 = (2 * iqr1) / (len(dist_p1) ** (1 / 3))
	bin_count1_p1 = int(np.ceil((para1_range) / bin_width_p1))
	xbin = bin_count1_p1
	print(f'  Number of bins deduced using the Freedman-Diaconis (1981) rule')
	time.sleep(2)
	print(f'\n    bin_count_x = {xbin}')

	# Scott (1979) method
	stdev1 = dist_p1.std()
	bin_width_p1_scott = (3.5 * stdev1) / (len(dist_p1) ** (1 / 3))
	bin_count1_p1_scott = int(np.ceil((para1_range) / bin_width_p1_scott))
	time.sleep(1)

	dist_p2 = pandas.Series(order_p2)
	para2_max = dist_p2.max()
	para2_min = dist_p2.min()
	para2_range = para2_max - para2_min
	qt1 = dist_p2.quantile(0.25)
	qt3 = dist_p2.quantile(0.75)
	iqr2 = qt3 - qt1
	bin_width_p2 = (2 * iqr2) / (len(dist_p2) ** (1 / 3))
	bin_count2_p2 = int(np.ceil((para2_range) / bin_width_p2))
	ybin = bin_count2_p2
	print(f'    bin_count_y = {ybin}')
	# Scott (1979) method
	stdev2 = dist_p2.std()
	bin_width_p2_scott = (3.5 * stdev2) / (len(dist_p2) ** (1 / 3))
	bin_count2_p2_scott = int(np.ceil((para2_range) / bin_width_p2_scott))
	time.sleep(2)

	with open("CHAP_fes_Par.in", "a") as in_par:
		in_par.write(f'x_bin_count,{xbin}\ny_bin_count,{ybin}')

	# Determine the number of bins using the sqrt method
	num_of_bins_sqrt = int(np.ceil(math.sqrt(len(dist_p1))))
	num_of_bins_rice = int(np.ceil( 2 * (len(dist_p1) ** (1 / 3))))

	print(f"\n Optimal binning parameters have been estimated.\
		\n These parameters have been written to file (CHAP_fes_Par.in).\
		\n\n Do you want to proceed?\n  (1) Yes\n  (2) No\n")

	prmpt = " Enter a response here (1 or 2): "
	response = int(input(prmpt))

	while response != 1 and response != 2:
		print("\n ENTRY REJECTED!\n **Please enter the appropriate option (1 or 2)\n")
		response = int(input(prmpt))

	with open("binning_summary.dat", "w") as bin_summary:
		bin_summary.write(f"Binning Method\t  | Number of bins\n")
		bin_summary.write(f"------------------|---------------\n")
		bin_summary.write(f"Freedman-Diaconis | x={xbin}, y={ybin}\t(*)\n")
		bin_summary.write(f"Square root\t\t  | x={num_of_bins_sqrt}, y={num_of_bins_sqrt}\n")
		bin_summary.write(f"Rice\t\t\t  | x={num_of_bins_rice}, y={num_of_bins_rice}\n")
		bin_summary.write(f"Scott\t\t\t  | x={bin_count1_p1_scott}, y={bin_count2_p2_scott}")

	if response == 2:
		sys.exit(0)
	elif response == 1:
		print (f"\n Updating input parameters for FES calculations\n")
		with open("CHAP_fes_Par.in") as in_par:
			for parameter in in_par.readlines():
				if "XaxisL" in parameter:
					para_data = parameter.rstrip('\n').split(",")
					xaxis_label = str(para_data[1])
				elif "YaxisL" in parameter:
					para_data = parameter.rstrip('\n').split(",")
					yaxis_label = str(para_data[1])
				elif "Temp" in parameter:
					para_data = str(parameter).split(",")
					Temp = float(para_data[1])
				elif "outFilename" in parameter:
					para_data = parameter.rstrip('\n').split(",")
					out_file = str(para_data[1])+str(".png")
				elif "plotTitle" in parameter:
					para_data = parameter.rstrip('\n').split(",")
					plotTitle = str(para_data[1])+str(" Free Energy Surface")
				elif "x_bin_count" in parameter:
					para_data = parameter.rstrip('\n').split(",")
					xbin_custom = int(para_data[1])
				elif "y_bin_count" in parameter:
					para_data = parameter.rstrip('\n').split(",")
					ybin_custom = int(para_data[1])
		xbin = xbin_custom
		ybin = ybin_custom
		print (f" Generating a 2D histogram\n")		
		# Create a 2D histogram using the numpy histogram2d function
		hist, x_edges, y_edges = np.histogram2d(order_p1, order_p2, bins=(xbin, ybin), \
			range=[[para1_min, para1_max], [para2_min, para2_max]])
		time.sleep(2)

print (" Identifying the highest probability bin"+"\n")
# Flatten the 2D histogram into a 1D array
Prob = hist.flatten()

# Identify the most populated bin
max_bin = np.max(Prob)
time.sleep(2)

# Constant -> product of kilocal conversion factor, Avogadro's number, Boltzmann constant & temperature
RT = -0.001 * 6.02214E23 * 3.29763E-24 * Temp

print (" Calculating delta_G values by Boltzmann inversion of the histogram"+"\n")

# Initialize the dG array
dG = np.zeros((xbin,ybin))

# Estimate dG values using Boltzmann inversion
with open("OrderParameters1_2_dG.dat", "w") as dGoutFile:
	for x in range(xbin):
		for y in range(ybin):
			if hist[x][y] == 0:
				dG[x][y] = 10
				dGoutFile.write((f"{(2*para1_min+(2*x+1)*para1_range/xbin)/2}\t{(2*para2_min+(2*y+1)*para2_range/ybin)/2}\t{dG[x][y]}\n"))
				continue
			else:
				dG[x][y] = RT*(np.log(hist[x][y]) - np.log(max_bin))
				dGoutFile.write((f"{(2*para1_min+(2*x+1)*para1_range/xbin)/2}\t{(2*para2_min+(2*y+1)*para2_range/ybin)/2}\t{dG[x][y]}\n"))
				dGoutFile.write("\n")
time.sleep(2)

print (" Generating and saving FES plot")
# Plot figure
plt.figure()
plt.xlabel(xaxis_label)
plt.ylabel(yaxis_label)
plt.title(plotTitle)
ext = [para1_min,para1_max,para2_min,para2_max]
im = plt.gca().imshow(dG.T, origin='lower', aspect='auto', cmap="gnuplot", extent=ext)
c_ax = make_axes_locatable(plt.gca()).append_axes("right", size="2.5%", pad=0.1)
plt.colorbar(im, cax=c_ax).set_label(r'$\Delta G$'+' (kcal/mol)', size=12)
plt.savefig(out_file,dpi=600)