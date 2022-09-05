#############################################################################
#-------------------------------- CHAPERONg --------------------------------#
#  An automated pipeline for GROMACS MD simulation and end-point analysis   #
#      If you find this program useful please cite the relevant paper:      #
#                   Yekeen, A.A. et al. To be published...                  #
#############################################################################

#---------------------------------------------------------------------------#
#### ================== SETTING UP & RUNNING CHAPERONg ================= ####
#---------------------------------------------------------------------------#

1) After extracting the compressed file CHAPERONg-<version>, you will find a
   folder CHAPERONg-<version> and this README.txt file.

2) Go into the CHAPERONg-<version> folder and you will find the binary
   "setup_CHAPERONg-<version>", a folder "CHAP_modules" and another
   copy of this README.txt file.

3) Copy "setup_CHAPERONg-<version>" and the folder "CHAP_modules" into
   your md simulation working directory (i.e the directory where you have
   your input files ready for mds).
   
4) In your linux terminal, and while in the mds working directory, execute
   the following in order:
   
	chmod +x setup_CHAPERONg-<version>
	./setup_CHAPERONg-<version>

5) The script above will make all the CHAPERONg binaries executable and will
   also copy the binary "run_CHAPERONg-<version>" into the working 
   directory.
   
6) And that's it! You can then launch CHAPERONg by executing:

   	./run_CHAPERONg-<version>


#---------------------------------------------------------------------------#
#  Author: Abeeb A. Yekeen                                                  #
#  Contact: yekeenaa@mail.ustc.edu.cn, abeeb.yekeen@hotmail.com             #
#---------------------------------------------------------------------------#


