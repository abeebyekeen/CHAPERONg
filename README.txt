CHAPERONg: A Tool for Automated GROMACS-based MD Simulations and Trajectory Analyses

## 1. Prerequisites/Requirements

The only essential requirement of CHAPERONg is, of course, a functional installation of the GROMACS software run on Linux (Ubuntu) or Windows (via WSL).
You can download and install a version by visiting the GROMACS documentation page (https://manual.gromacs.org/),
or you could find the latest version of GROMACS at https://manual.gromacs.org/current/index.html.

CHAPERONg has been tested with GROMACS versions 2020.x, 2021.x, and 2022.x, but should generally work with GROMACS 5.x versions and later.


## 2. Installation

Step 1 -- Copy the CHAPERONg package to your preferred installation directory.

Step 2 -- Extract the package with:

      tar xfz CHAPERONg-v0.1.tar.gz

or with:

      unzip CHAPERONg-v0.1.zip

Step 3 -- Enter the extracted CHAPERONg folder:

      cd CHAPERONg-v0.1

Step 4 -- Run the installation script:

      chmod +x install_CHAPERONg.sh && ./install_CHAPERONg.sh

   Enter `y` or `yes` when prompted for permission.
   During installation, CHAPERONg will append the installation path to your `$PATH` via the `~/.bashrc` file.
   After installation, you should see a completion message.


## 3. Uninstall

At any time, you can uninstall CHAPERONg by simply, and CAREFULLY, deleting the lines that have been added your `~/.bashrc` file.


## 4. Optional Dependencies

CHAPERONg will run the entire GROMACS simulation pipelines without additional dependencies installed.
However, some other dependencies are required to run some specific analyses.
The more of these additional tools you install, the more you are able to maximize the usage and functionalities of CHAPERONg.

Existing dependencies do not need to be re-installed, i.e., you do not need to re-install any active packages.

There are two options by which you can install the dependencies:
- installation in conda environment (recommended), or
- global (system-wide) installation.

### 4.1. Installing dependencies with conda (recommended)

Installation in conda environment is recommended for most of the dependencies.

Step 1 -- Change to your `base` anaconda (default) environment:

      conda activate base

Step 2 -- Run the conda environment setup script `conda_env_setup.sh` while in the CHAPERONg folder:

      chmod +x conda_env_setup.sh && ./conda_env_setup.sh

   This step requires an internet connection and could take a while to complete depending on the stength of your connection.
   So, relax and give it time.

You should see a message indicating that the `chaperong` environment has been set up successfully.

You can confirm the newly created environment with:

      conda env list

You should see `chaperong` listed among your conda environments.

You can then change into the `chaperong` environment with:

      conda activate chaperong

If, after setting up the conda environment as described above, some python modules/libraries are, for some reasons,
not accessible, you can further run the commands described in the system-wide option 1 below.

### 4.2. Installing dependencies system-wide (optional)

Step 1 -- Install all of the dependencies at once:

      sudo apt install grace pymol libpng-dev libjpeg-dev ghostscript-x imagemagick-6.q16 ffmpeg

Step 2 -- Install/update python libraries as well as a few other packages:

      pip3 install numpy=1.23.4 scipy matplotlib pandas networkx=2.3

If you don't have pip installed, you can first install it with:

      sudo apt install python3-pip

And then, run the command:

      pip3 install numpy=1.23.4 scipy matplotlib pandas networkx=2.3


## 5. Citation

If you use this code in your work, kindly cite it as:

Yekeen A. A., Durojaye O. A., Idris M. O., Muritala H. F., and Arise R. O. (2023). CHAPERONg: A tool for automated GROMACS-based molecular dynamics simulations and trajectory analyses, bioRxiv 2023.2007.2001.546945. DOI: https://doi.org/10.1101/2023.07.01.546945

For BibTeX, you can insert the following:

@article {yekeenbiorxiv2023chaperong,
	author = {Abeeb Abiodun Yekeen and Olanrewaju Ayodeji Durojaye and Mukhtar Oluwaseun Idris and Hamdalat Folake Muritala and Rotimi Olusanya Arise},
	title = {{CHAPERONg}: A tool for automated GROMACS-based molecular dynamics simulations and trajectory analyses},
	elocation-id = {2023.07.01.546945},
	year = {2023},
	doi = {10.1101/2023.07.01.546945},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2023/07/02/2023.07.01.546945},
	eprint = {https://www.biorxiv.org/content/early/2023/07/02/2023.07.01.546945.full.pdf},
	journal = {bioRxiv}
}

