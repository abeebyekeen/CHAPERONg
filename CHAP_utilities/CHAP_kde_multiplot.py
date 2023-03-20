
# CHAP_kde_multiplot.py -- A python script to estimate probability density
#   function using the kernel density estimator with the Gaussian kernel.
#   The estimation is made for multiple comparable data and plotted together.
# Part of the CHAPERONg suite of scripts
# Input parameters are generated by other scripts in CHAPERONg and are read
#   by this script
# CHAPERONg -- An automation program for GROMACS md simulation
# Author -- Abeeb A. Yekeen
# Contact -- abeeb.yekeen@hotmail.com
# Date: 2023.02.14


import math
import os
import shutil
import sys
import time

# def check_and_import_lib
missingLib = []
try:
	import matplotlib
except ModuleNotFoundError:
	print("The matplotlib library has not been installed!\n")
	missingLib.append("matplotlib")
else:
	# Configure the non-interactive backend
	matplotlib.use('AGG')
	import matplotlib.pyplot as plt
try:
	import numpy as np
except ModuleNotFoundError:
	print(" The numpy library has not been installed!\n")
	missingLib.append("numpy")
try:
	import pandas as pd
except ModuleNotFoundError:
	print(" The pandas library has not been installed!\n")
	missingLib.append("pandas")
try:
	import scipy.stats as st
except ModuleNotFoundError:
	print(" The scipy library has not been installed!\n")
	missingLib.append("scipy")

if len(missingLib) >= 1 :
	print('\n\n#================================= CHAPERONg =================================#\n')
	print(" One or more required libraries have not been installed\n")
	# if len(missingLib) == 1 : print(" Missing library:\n")
	# elif len(missingLib) > 1 : print(" Missing libraries:\n")
	# Ternary Operator
	# print(" Missing library:\n") if len(missingLib) == 1 else print(" Missing libraries:\n")
	print(" Missing library:\n" if len(missingLib) == 1 else " Missing libraries:\n")
	suggestmsg = "pip install "
	for lib in missingLib:
		print(f"  {lib}\n")
		if lib == "numpy" : lib = "numpy==1.23.4"
		suggestmsg = suggestmsg + f'{lib} '
	print(
		f' {suggestmsg}\n\n or\n\n'
		' See https://www.abeebyekeen.com/chaperong-online-documentation/'
		)
	sys.exit(0)

def make_dir_for_KDE(motherDir='Kernel_Density_Estimation_multi_plot'):
	if os.path.exists(motherDir):
		backup_count = 1
		backupDir=f'#{motherDir}.backup.{backup_count}'
		while os.path.exists(backupDir):
			backup_count += 1
			backupDir = f'#{motherDir}.backup.{backup_count}'			
		shutil.move(motherDir, backupDir)
	os.mkdir(motherDir)

def store_data_label_name():
	# Read in the data and label
	print (" Reading in parameters for density estimation\n")
	with open("CHAP_kde_dataset_list.dat") as in_par:
		alldatasets = in_par.readlines()
		# Create a dictionary of label-data pair
		input_data_dict = {}
		
		for lineNo, line in enumerate(alldatasets):
			if int(lineNo) == 0 and "auto mode" in line:
				global auto_mode
				auto_mode_raw = str(line).rstrip("\n").split(",")
				auto_mode = auto_mode_raw[1]
			elif int(lineNo) == 2:
				# Get the type of data from the header
				dataName = str(line).rstrip("\n")
				continue
			elif int(lineNo) >= 3:
				data_info = str(line).rstrip("\n").split(",")
				data_label = str(data_info[0])
				data = str(data_info[1])
				
				# Add label and data to the dictionary
				input_data_dict[data_label] = data
	return dataName, input_data_dict

output_and_para_files = []

# def estimate_PDF_with_KDE_multiplot():
def plot_multidata_hist(dataName, input_data_dict):
	data_count = 1
	plt.figure() # Create a new figure
	for key, value in input_data_dict.items():
		dataLabel = key
		extracted_data = value
		data_in = []
		with open(extracted_data) as alldata:
			alldata_lines = alldata.readlines()
			for line in alldata_lines:
				data_point = str(line).rstrip("\n")
				data_in.append(float(data_point))

		# Determine the number of bins automatically
		print (
			"#=============================================================================#\n"
			f"\n Estimating the optimal number of histogram bins for {dataLabel}\n"
			)
		time.sleep(2)
		# Determine the number of bins using the Freedman-Diaconis (1981) method
		dist = pd.Series(data_in)
		data_max, data_min = dist.max(), dist.min()
		data_range = data_max - data_min
		qt1 = dist.quantile(0.25)
		qt3 = dist.quantile(0.75)
		iqr = qt3 - qt1
		bin_width = (2 * iqr) / (len(dist) ** (1 / 3))
		bin_count = int(np.ceil((data_range) / bin_width))
		print(f'  Number of bins deduced using the Freedman-Diaconis (1981) rule')
		time.sleep(2)
		print(f'\n    bin_count = {bin_count}')

		def writeOut_parameters():
			in_par.write(
				f'{dataLabel}\n'
				f'bin_count,{bin_count}\n'
				"bandwidth_method,silverman\n\n"
				)			

		if data_count == 1 :
			with open("CHAP_kde_Par.in", "w") as in_par:
				in_par.write(f'=>{dataName}\n')
				writeOut_parameters()

		elif data_count > 1 :
			with open("CHAP_kde_Par.in", "a") as in_par:
				writeOut_parameters()				

		# output_and_para_files.append(f'CHAP_kde_Par.in')

		# Scott (1979) method
		stdev = dist.std()
		bin_width_scott = (3.5 * stdev) / (len(dist) ** (1 / 3))
		bin_count_scott = int(np.ceil((data_range) / bin_width_scott))
		time.sleep(1)

		# Determine the number of bins using the sqrt method
		num_of_bins_sqrt = int(np.ceil(math.sqrt(len(dist))))
		num_of_bins_rice = int(np.ceil( 2 * (len(dist) ** (1 / 3))))

		def write_binning_parameters():
			bin_summary.write(
				f'{dataLabel}\n'
				"----------------------------------\n"
				"Binning Method\t  | Number of bins\n"
				"------------------+---------------\n"
				f'Freedman-Diaconis | {bin_count}\t(*)\n'
				f'Square root\t\t  | {num_of_bins_sqrt}\n'
				f'Rice\t\t\t  | {num_of_bins_rice}\n'
				f'Scott\t\t\t  | {bin_count_scott}\n'
				"----------------------------------\n\n"
				)
		
		if data_count == 1:
			with open("kde_bins_estimated_summary.dat", "w") as bin_summary:
				bin_summary.write(f'=> {dataName}\n')
				write_binning_parameters()

		elif data_count > 1:
			with open("kde_bins_estimated_summary.dat", "a") as bin_summary:
				write_binning_parameters()

		if auto_mode == 'semi':
			print(
				"\n  Optimal binning parameters have been estimated."
				'\n  Parameters have been written to the file "CHAP_kde_Par.in".'
				'\n\n  You can modify the parameters if required.'
				'\n   Enter "Yes" below when you are ready.'
				"\n\n   Do you want to proceed?\n    (1) Yes\n    (2) No\n"
				)

			prmpt = "  Enter a response here (1 or 2): "
			response = int(input(prmpt))

			while response != 1 and response != 2:
				print(
					"\n ENTRY REJECTED!"
					"\n **Please enter the appropriate option (1 or 2)\n"
					)
				response = int(input(prmpt))

		elif auto_mode == 'full':
			response = 1
			print(
				"\n   CHAPERONg is running in full-auto mode"
				"\n   The estimated number of bins above will be used"
				"\n   To use a different number or estimator,"
				"\n   run CHAPERONg in the semi-auto mode. For details, see"
				"\n   https://www.abeebyekeen.com/post-sim-analysis-1/"
			)
			time.sleep(2)

		if response == 2:
			sys.exit(0)
		elif response == 1:
			print ("\n Updating input parameters for density estimation\n")
			time.sleep(2)
			with open("CHAP_kde_Par.in" , 'r') as in_par:
				for parameter in in_par.readlines():
					if "bin_count" in parameter:
						para_data = parameter.rstrip('\n').split(",")
						bin_custom = int(para_data[1].strip())
			bin_set = bin_custom

			print (f" Generating and plotting the histogram of the {dataLabel}\n")
			time.sleep(2)
		
			if "RMSD" in dataName: 
				XaxisLabelXVG = r'RMSD (\cE\C)' 
				XaxisLabelPNG = 'RMSD' + r' ($\AA$)' # Using Latex in matplotlib
			elif "Rg" in dataName:
				XaxisLabelXVG = r'Radius of gyration (\cE\C)'
				XaxisLabelPNG = 'Radius of gyration' + r' ($\AA$)'
			elif "Hbond" in dataName:
				XaxisLabelXVG = "Number of hydrogen bonds"
				XaxisLabelPNG = XaxisLabelXVG
			elif "SASA" in dataName:
				XaxisLabelXVG = r'SASA (nm\S2\N)'
				XaxisLabelPNG = 'SASA' + r' ($nm^{2}$)'

			# Generate and plot the histogram of the data
			histo = plt.hist(data_in, bins=bin_set, label=dataLabel, alpha=0.8)

			# The first elements are the ys, the second are the xs.
			# ys = histo[0]; xs = histo[1]

			# The x-axis values are boundaries, starting with the lower bound 
			# of the first and ending with the upper bound of the last.
			# So, length of the x-axis > length of the y-axis by 1.

			out_hist = dataLabel + "_histogram.xvg"
			with open (out_hist, 'w') as out_his_file:
				out_his_file.write(
					f'# This file contains the histogram values of the {dataLabel}'
					'\n# data calculated by CHAPERONg from the output of GROMACS\n#\n'
					f'@    title Histogram of {dataName}\n'
					f'@    xaxis  label "{XaxisLabelXVG}"\n'
					'@    yaxis  label "Count"\n'
					'@TYPE bar\n'
					f'@ s0 legend "{dataLabel}"\n'
					'@    s0 symbol size 0.200000\n'
					'@    s0 line type 0\n'
					)
			
			# Get the upper bounds and write out the histogram			
			pd.DataFrame({'x_upper':histo[1][1:], 'y': histo[0]}).to_csv(
				out_hist, header=False, index=False, sep="\t", mode='a'
				)
			
			output_and_para_files.append(out_hist)

			# Increase counter for additional data
			data_count += 1

	print (f" Generating the combined histogram plots of the {dataName}\n")
	time.sleep(2)
	
	plt.xlabel(XaxisLabelPNG) # using Latex expression in matplotlib
	plt.ylabel('Count')
	plt.legend()
	plt.title("Histogram of the " + dataName)
	figname = dataName + "_histogram_multi_plot.png"
	plt.savefig(figname, dpi=600)
		
	output_and_para_files.append(figname)
	return dataName, XaxisLabelXVG, XaxisLabelPNG

def estimate_PDF_with_KDE_multiplot(dataName, XaxisLabelXVG, XaxisLabelPNG):		
	data_count = 1
	bins_number_dict = {}
	bandwidth_dict = {}
	bins_number_count = 1
	print (f" Extracting pre-calculated number of bins for plotting histogram\n")
	time.sleep(2)	
	with open('CHAP_kde_Par.in', 'r') as in_par:
		for line in in_par.readlines():
			if "bin_count" in line:
				bin_w = line.rstrip("\n").split(",")
				bins_number_dict[bins_number_count] = bin_w[1]
			elif "bandwidth_method" in line:
				band_w = line.rstrip("\n").split(",")
				bandwidth_dict[bins_number_count] = band_w[1]
				bins_number_count += 1

	# Reset bin count
	bins_number_count2 = 1

	plt.figure() # Create a new figure for KDE

	for key, value in input_data_dict.items():
		dataLabel = key
		extracted_data = value
		data_in = []
		with open(extracted_data, 'r') as alldata:
			alldata_lines = alldata.readlines()
			for line in alldata_lines:
				data_point = str(line).rstrip("\n")
				data_in.append(float(data_point))

		bin_set = int(bins_number_dict[bins_number_count2])
		bandwidth = bandwidth_dict[bins_number_count2]
		bins_number_count2 += 1
		print(f'  Number of bins deduced using the Freedman-Diaconis (1981) rule')
		time.sleep(2)
		print(f'\n    bin_count = {bin_set}\n    bandwidth = {bandwidth}')

		# Generate and plot the histogram of the data
		plt.hist(data_in, density=True, bins=bin_set, label=dataLabel, alpha=0.6)

		print (f"\n Estimating the probability density function for {dataLabel}\n")
		time.sleep(2)
		kde_xs = np.linspace(min(data_in), max(data_in), 300)
		kde = st.gaussian_kde(data_in, bw_method=bandwidth)
		kde_ys = kde.pdf(kde_xs)
		kdeLabel = dataLabel + "_PDF"
		plt.plot(kde_xs, kde.pdf(kde_xs), label=kdeLabel)
		plt.legend()

		out_kde = dataLabel + "_KDEdata.xvg"
		with open (out_kde, 'w') as out_kde_file:
			out_kde_file.write(
				f'# This file contains the KDE-estimated PDF values of the {dataLabel}'
				'\n# data calculated by CHAPERONg from the output of GROMACS\n#\n'
				f'@    title "KDE-estimated Probability Density of {dataName}"\n'
				f'@    xaxis  label "{XaxisLabelXVG}"\n'
				'@    yaxis  label "Density"\n'
				'@TYPE xy\n'
				f'@ s0 legend "{dataLabel}_PDF"\n'
				)
			
		pd.DataFrame({'x':kde_xs, 'y': kde_ys}).to_csv(
					out_kde, header=False, index=False, sep="\t", mode='a'
					)
		
		output_and_para_files.append(out_kde)
		# Increase counter for additional data
		data_count += 1
		plt.ylabel("Density")
		plt.xlabel(XaxisLabelPNG)

		output_and_para_files.append(value)

	plt.title(f'Kernel Density Estimation Plots of {dataName} Data')
	figname = dataName + "_KDE_multi_plot.png"
	plt.savefig(figname, dpi=600)
	
	output_and_para_files.append(figname)

	print (f"\033[92m Estimate probability density function for {dataName}...DONE\033[00m\n"
			"#=============================================================================#\n")

	for file in output_and_para_files:
		try: shutil.move(file, './Kernel_Density_Estimation_multi_plot')
		except: pass

	para_summary = [
		'CHAP_kde_Par.in', 
		'kde_bins_estimated_summary.dat', 
		'CHAP_kde_dataset_list.dat',
		]

	for file in para_summary:
		try: shutil.move(file, './Kernel_Density_Estimation_multi_plot')
		except FileNotFoundError: pass				

make_dir_for_KDE()
dataName, input_data_dict = store_data_label_name()
dataName, XaxisLabelXVG, XaxisLabelPNG = plot_multidata_hist(dataName, input_data_dict)
estimate_PDF_with_KDE_multiplot(dataName, XaxisLabelXVG, XaxisLabelPNG)
