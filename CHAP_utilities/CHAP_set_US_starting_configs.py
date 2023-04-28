'''
A python script to identify initial configurations for umbrella sampling
Part of the CHAPERONg suite of scripts
CHAPERONg -- An automation program for GROMACS md simulation
Author -- Abeeb A. Yekeen
Contact -- yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com
Date: 2022.02.11
'''

##lines with ## are actual comments
def get_spaced_frame_dist(phrase='us_window_spacing'):
	with open("paraFile.par") as par:
		for parameter in par.readlines():
			if phrase in parameter:
				paraData=[]
				paraData_temp = str(parameter).strip().split(" ")
				for paraTemp in paraData_temp:
					if paraTemp != ' ' and paraTemp != '' and paraTemp != '\t' and paraTemp != '\n':
						paraData.append(paraTemp)
				spacing = float(paraData[2])
	with open("distances_summary.txt") as alldata:
		alldataLines = alldata.readlines()
		for lineNo, line in enumerate(alldataLines):
			##using enumerate to map each line of the file to
			##it's line_number starting line number from zero
			datalist = str(line).split("\t")
			frame = int(datalist[0])
			dist = float(datalist[1])
			#bkup_dst = dist 
			if int(lineNo) == 0:
				##register and write out the distance for frame zero
				#regFrame = frame
				regDist = dist
				with open("configuratns_list.txt", "w") as config:
					capture1 = "frame#"+"\t"+"dist"+"\t"+"d_dist"+"\n"
					capture2 = str(frame)+"\t"+str(regDist)+"\t"+"nil"+"\n"
					config.write(capture1)
					config.write(capture2)
					#print ("lineNo "+str(frame)+" "+str(dist))	
				continue

			elif int(lineNo) > 0:
				##this condition applies to other frames after frame zero
				diff = float("{:.3f}".format(dist - regDist))
				##setting a range for configurations to save
				spacing_upper = float(spacing) + float(spacing/16)
				spacing_lower = float(spacing) - float(spacing/16)
				#print ("lineNo "+str(frame)+" "+str(dist))			
				if diff > spacing_upper:
					##if the distance for the current frame jumps too higher, higher
					##than the interval, use the distance for the frame just before it
					diff = float("{:.3f}".format(bkup_dst - regDist))
					with open("configuratns_list.txt", "a") as config:
						capture = str(bkup_frame)+"\t"+str(bkup_dst)+"\t"+str(diff)+"\n"
						config.write(capture)
						##now use the bkup as the new regDist
						regDist = bkup_dst
					diff = float("{:.3f}".format(dist - regDist))

				if diff <= spacing_upper and diff >= spacing_lower:
					with open("configuratns_list.txt", "a") as config:
						capture = str(frame)+"\t"+str(dist)+"\t"+str(diff)+"\n"
						config.write(capture)
						regDist = dist
						#print ("lineNo "+str(frame)+" "+str(dist)+" "+str(diff))
						
				if diff < spacing_upper and int(lineNo)+1 == len(alldataLines):
					##write the last frame and the difference with the last regDist
					with open("configuratns_list.txt", "a") as config:
						config.write(str(frame)+"\t"+str(dist)+"\t"+str(diff)+"\n")
			bkup_dst = dist
			bkup_frame = frame
					
get_spaced_frame_dist()

