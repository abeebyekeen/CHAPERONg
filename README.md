## CHAPERON*g*: Automated GROMACS-based MD Simulations and Trajectory Analyses

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/abeebyekeen/CHAPERONg?style=flat-square)](https://github.com/abeebyekeen/CHAPERONg/releases)
[![DOI](https://zenodo.org/badge/doi/10.1101/2023.07.01.546945.svg?style=svg)](https://doi.org/10.1101/2023.07.01.546945)


### For detailed documentation and tutorials, see the [CHAPERON*g* homepage](https://abeebyekeen.com/chaperong-online/).

## 1. Prerequisites/Requirements

<font><p align="justify">The only essential requirement of **CHAPERON*****g*** is, of course, a functional installation of the <a href=https://www.gromacs.org/ target="_blank">GROMACS software</a> run on Linux (Ubuntu) or Windows (via WSL). You can download and install a version by visiting the [GROMACS documentation page](https://manual.gromacs.org/), or you could find the latest version of GROMACS [here](https://manual.gromacs.org/current/index.html).

> **CHAPERON*****g*** has been tested with GROMACS versions 2020.x, 2021.x, and 2022.x, but should generally work with GROMACS 5.x versions and later.
</p></font>

## 2. Installation

**Step 1:** Copy the **CHAPERON*****g*** package to your preferred installation directory.

**Step 2:** Extract the package with:

```bash
tar xfz CHAPERONg-v0.1.tar.gz
```

or with:

```bash
unzip CHAPERONg-v0.1.zip
```

**Step 3:** Enter the extracted **CHAPERON*****g*** folder:

```bash
cd CHAPERONg-v0.1
```

**Step 4:** Run the installation script:
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

> For potential error troubleshooting, see the [installation page](https://abeebyekeen.com/chaperong-online-documentation/#installation) of the [CHAPERON*g* documentation](https://abeebyekeen.com/chaperong-online-documentation/).

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

### 4.2. Installing dependencies system-wide (optional)
#### 4.2.1. Option 1: Batch installation {#install-system-wide-option1}

**Step 1**: To install all of the dependencies at once, enter the following command:
```bash
sudo apt install grace pymol libpng-dev libjpeg-dev ghostscript-x imagemagick-6.q16 ffmpeg
```
<br>

**Step 2**: Install/update python libraries as well as a few other packages:
```bash
pip3 install numpy=1.23.4 scipy matplotlib pandas networkx=2.3
```

{{< admonition note "If the pip package is not available" >}}
If you don't have pip installed, you can first install it with:

`sudo apt install python3-pip`

And then, run the command:

`pip3 install numpy=1.23.4 scipy matplotlib pandas networkx=2.3`

{{< /admonition >}}

<br>

#### 4.2.2. Option 2: One-by-one installation {#install-system-wide-option1}

**&nbsp; 1\. Python3**
{{< admonition note "Note" >}}
The analyses run with **CHAPERON*****g*** have been tested with Python versions 3.8.x and 3.9.x.
{{< /admonition >}}

<br>

**&nbsp; 2\. [Grace](https://plasma-gate.weizmann.ac.il/Grace/)**
<br>
```bash
sudo apt install grace
```
{{< admonition note "Note" >}}
The system-wide installation of xmgrace as above is probably more convenient and is recommended when possible.
{{< /admonition >}}

<br>

**&nbsp; 3\. [PyMOL](https://pymol.org/)**
<br>
You can install PyMOL using apt (as shown below) or visit the [PyMOL page](https://pymol.org/) to download a distribution of your choice.
```bash
sudo apt install pymol
```

{{< admonition tip "Installation in conda environment is recommended for most of the tools" /admonition >}}

<br>

**&nbsp; 4\. [ImageMagick](https://imagemagick.org/index.php)**
<br>
{{< admonition note "Note" >}}
The `convert` program needed in the ImageMagick package is often already installed in most Linux installations. You can run the `convert --version` command first to check before the next command below. Either way, it shouldn't do any harm.
{{< /admonition >}}

<br>

**Step 1**: Install the `ghostscript` interpreter and some libraries that might be missing
```bash
sudo apt install -y libpng-dev libjpeg-dev ghostscript-x
```

<br>

**Step 2**: Install ImageMagick
```bash
sudo apt install imagemagick-6.q16
```
{{< admonition tip "Tip" >}}
Confirm that the program `convert` works
```bash
convert
```
or with:
```bash
convert --version
```
{{< /admonition >}}

<br>

**&nbsp; 5\. [Networkx2.3](https://networkx.org/)**

```bash
pip3 install numpy scipy matplotlib networkx=2.3
```

{{< admonition note "Note" >}}
If you don't have pip installed, you can first install it with:
`sudo apt install python3-pip`
{{< /admonition >}}

<br>

**&nbsp; 6\. [ffmpeg](https://networkx.org/)**

```bash
sudo apt install ffmpeg
```

<br><br><br>
