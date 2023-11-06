## CHAPERON*g*: Automated GROMACS-based MD Simulations and Trajectory Analyses

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/abeebyekeen/CHAPERONg?style=flat-square)](https://github.com/abeebyekeen/CHAPERONg/releases)
[![DOI](https://zenodo.org/badge/doi/10.1101/2023.07.01.546945.svg?style=svg)](https://doi.org/10.1101/2023.07.01.546945)

<p align="center">
  <a href="https://github.com/abeebyekeen/CHAPERONg">
	<img alt="CHAPERONg" src="CHAP_utilities/chaperong-banner-dir-lowRes.png">
  </a>
</p>

### For detailed documentation and tutorials, see the [CHAPERON*g* homepage](https://abeebyekeen.com/chaperong-online/).

## 0. Contents

0. [Contents](contents)
1. [Prerequisites/Requirements](#1-prerequisitesrequirements)
2. [Installation](#2-installation)
3. [Uninstall](#3-uninstall)
4. [Optional Dependencies](#4-optional-dependencies)
5. [Usage](#5-usage)
6. [Updating the Code](#6-updating-the-code)
7. [Citation](#7-citation)

## 1. Prerequisites/Requirements

<font><p align="justify">The only essential requirement of **CHAPERON*****g*** is, of course, a functional installation of the <a href=https://www.gromacs.org/ target="_blank">GROMACS software</a> run on Linux (Ubuntu) or Windows (via WSL). You can download and install a version by visiting the [GROMACS documentation page](https://manual.gromacs.org/), or you could find the latest version of GROMACS [here](https://manual.gromacs.org/current/index.html).

> **CHAPERON*****g*** has been tested with GROMACS versions 2020.x, 2021.x, and 2022.x, but should generally work with GROMACS 5.x versions and later.
</p></font>

## 2. Installation

> ***In the instructions that follow, replace "xxx" with "main" or "v0.1" as you have it in the downloaded package***

**Step 1:** Copy the **CHAPERON*****g*** package to your preferred installation directory.

**Step 2:** Extract the package with:

```bash
tar xzvf CHAPERONg-xxx.tar.gz
```

or with:

```bash
unzip CHAPERONg-xxx.zip
```

**Step 3:** Enter the extracted **CHAPERON*****g*** folder:

```bash
cd CHAPERONg-xxx
```

**Step 4:** Run the installation executable:

```bash
chmod +x install_CHAPERONg.sh && ./install_CHAPERONg.sh
```

Enter `y` or `yes` when prompted for permission. During installation, **CHAPERON*****g*** will append the installation path to your `$PATH` via the `~/.bashrc` file. After installation, you should see a completion message.
<br>

## 3. Uninstall

<font><p align="justify">At any time, you can uninstall **CHAPERON*****g*** by simply, and **CAREFULLY**, deleting the lines that have been added your `~/.bashrc` file.
</p></font>

## 4. Optional Dependencies

<font><p align="justify">**CHAPERON*****g*** will run the entire GROMACS simulation pipelines without additional dependencies installed. However, some other dependencies are required to run some specific analyses. The more of these additional tools you install, the more you are able to maximize the usage and functionalities of **CHAPERON*****g***.</p></font>

> **Existing dependencies do not need to be re-installed**, i.e., you do not need to re-install any active packages.

**There are two options by which you can install the dependencies:**
- ***installation in conda environment (recommended), or***
- ***global (system-wide) installation.***

### 4.1. Installing dependencies with conda (recommended)

> Installation in conda environment is recommended for most of the dependencies.

**Step 1**: Change to your `base` anaconda (default) environment:

```bash
conda activate base
```

**Step 2**: Run the conda environment setup script `conda_env_setup.sh` while in the **CHAPERON*****g*** folder:

```bash
chmod +x conda_env_setup.sh && ./conda_env_setup.sh
```

> This step requires an internet connection and could take a while to complete depending on the stength of your connection. So, relax and give it time.

> In case of errors, see the [installation page](https://abeebyekeen.com/chaperong-online-documentation/#installation) of the [CHAPERON*g* documentation](https://abeebyekeen.com/chaperong-online-documentation/) for guides on troubleshooting.

You should see a message indicating that the `chaperong` environment has been set up successfully.

You can confirm the newly created environment with:
```bash
conda env list
```

You should see `chaperong` listed among your conda environments.

You can then change into the `chaperong` environment with:
```bash
conda activate chaperong
```

> If, after setting up the conda environment as described above, some python modules/libraries are, for some reasons, not accessible, you can further run the commands described in the *system-wide option 1* below.

### 4.2. Installing dependencies system-wide
#### 4.2.1. Option 1: Batch installation

**Step 1**: Install all of the dependencies at once:
```bash
sudo apt install grace pymol libpng-dev libjpeg-dev ghostscript-x imagemagick-6.q16 ffmpeg
```

**Step 2**: Install/update python libraries as well as a few other packages:
```bash
pip3 install numpy=1.23.4 scipy matplotlib pandas networkx=2.3
```

> If you don't have pip installed, you can first install it with:

> `sudo apt install python3-pip`

> And then, run the command:

> `pip3 install numpy=1.23.4 scipy matplotlib pandas networkx=2.3`

#### 4.2.2. Option 2: One-by-one installation

See the [relevant section](https://abeebyekeen.com/chaperong-online-documentation/#install-system-wide-option2) of the **CHAPERON*g*** [documentation](https://abeebyekeen.com/chaperong-online-documentation/) and [installation page](https://abeebyekeen.com/chaperong-online-documentation/#installation).

## 5. Usage

The basic usage of **CHAPERON*****g*** is:

```bash
run_CHAPERONg.sh -i inputStructure_filename [more options]
```

For more details, see the **CHAPERON*****g*** [documentation](https://abeebyekeen.com/chaperong-online-documentation/) and [tutorials](https://abeebyekeen.com/chaperong-online-tutorials/)

## 6. Updating the Code

With the updates of **CHAPERON*****g*** pushed since Nov. 5 2023, it is now possible to automatically update **CHAPERON*****g*** for users who installed **CHAPERON*****g*** from this date onwards. To update the code, simply navigate to the **CHAPERON*****g*** installation directory:

```bash
cd $CHAPERONg_PATH
```

Run the `install_CHAPERONg.sh` script in the installation directory:

```bash
./install_CHAPERONg.sh
```

**CHAPERON*****g*** will detect the existing installation, and you will be prompted to confirm if you would like to proceed with an update of the code.

## 7. Citation

If you use this code in your work, kindly cite it as:

**Yekeen A. A., Durojaye O. A., Idris M. O., Muritala H. F., and Arise R. O. (2023). CHAPERON<em>g</em>: A tool for automated GROMACS-based molecular dynamics simulations and trajectory analyses, <em>Computational and Structural Biotechnology Journal</em> 21: 4849-4858. DOI: https://doi.org/10.1016/j.csbj.2023.09.024**

For BibTeX, you can insert the following:

```latex
@article {yekeenCSBJ2023chaperong,
	author = {Abeeb Abiodun Yekeen and Olanrewaju Ayodeji Durojaye and Mukhtar Oluwaseun Idris and Hamdalat Folake Muritala and Rotimi Olusanya Arise},
	title = {{CHAPERONg}: A tool for automated {GROMACS}-based molecular dynamics simulations and trajectory analyses},
	journal = {{Computational and Structural Biotechnology Journal}},
	volume = {21},
	pages = {4849-4858},
	year = {2023},
	issn = {2001-0370},
	doi = {https://doi.org/10.1016/j.csbj.2023.09.024},
	publisher = {Elsevier},
	url = {https://www.sciencedirect.com/science/article/pii/S2001037023003367}
}

```
