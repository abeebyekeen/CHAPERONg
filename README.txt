###############################################################################
#--------------------------------- CHAPERONg ---------------------------------#
#   An automated pipeline for GROMACS MD simulation and trajectory analyses   #
#    If you use this program in your work, please cite the relevant paper:    #
#                   Yekeen, A.A. et al. To be published...                    #
###############################################################################

#-----------------------------------------------------------------------------#
#### ========================= INSTALLATION GUIDE ======================== ####
#-----------------------------------------------------------------------------#

1) Downloaded the CHAPERONg package (if you haven't done so)
   wget ""

2) Extract the downloaded CHAPERONg package
   tar xfz CHAPERONg-v0.1.tar.gz
   OR
   gunzip CHAPERONg-v0.1.zip

3) Change into the CHAPERONg-v1.0 and run the installation script
   cd CHAPERONg-v1.0
   chmod +x install_CHAPERONg.sh
   ./install_CHAPERONg.sh

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


