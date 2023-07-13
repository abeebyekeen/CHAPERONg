# CHAPERON*g*: Automated GROMACS-based MD Simulations and Trajectory Analyses

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/abeebyekeen/CHAPERONg?style=flat-square)](https://github.com/abeebyekeen/CHAPERONg/releases)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.7049711.svg?style=svg)](https://zenodo.org/record/7049711#.YxWvrHZBzi0)


## For detailed documentations and tutorials, see the [CHAPERON*g* homepage](https://abeebyekeen.com/chaperong-online/)

## 1. Prerequisites/Requirements

<font><p align="justify">The only essential requirement of **CHAPERON*****g*** is, of course, a functional installation of the <a href=https://www.gromacs.org/ target="_blank">GROMACS software</a> run on Linux (Ubuntu) or Windows (via WSL). You can download and install a version by visiting the [GROMACS documentation page](https://manual.gromacs.org/), or you could find the latest version of GROMACS [here](https://manual.gromacs.org/current/index.html).

> **CHAPERON*****g*** has been tested with GROMACS versions 2020.x, 2021.x, and 2022.x, but should generally work with GROMACS 5.x versions and later.
</p></font>

## 2. Installation

**Step 1:** Copy the downloaded **CHAPERON*****g*** &#xf1c6; &#xf1c6; package to your preferred installation directory.

<br>

**Step 2:** Extract the package with:

```bash
tar xfz CHAPERONg-v0.1.tar.gz
```

or with:

```bash
unzip CHAPERONg-v0.1.zip
```
<br>

**Step 3:** Enter the extracted **CHAPERON*****g*** folder:

```bash
cd CHAPERONg-v0.1
```
<br>

**Step 4:** Run the installation script:
```bash
chmod +x install_CHAPERONg.sh && ./install_CHAPERONg.sh
```

<br>

Enter `y` or `yes` when prompted for permission. During installation, **CHAPERON*****g*** will append the installation path to your `$PATH` via the `~/.bashrc` file. After installation, you should see a completion message as shown below:<br><br>

{{< image src="chaperong-installation2.png" caption="CHAPERONg installation" width="800">}}
<br>
At any time, you can uninstall **CHAPERON*****g*** by simply, and CAREFULLY, deleting the lines that have been added your `~/.bashrc` file.
<br><br>

## 4. Optional Dependencies {#optional-dependencies}
<font><p align="justify">**CHAPERON*****g*** will run the entire GROMACS simulation pipelines without additional dependencies installed. However, some other dependencies are required to run some specific analyses facilitated by **CHAPERON*****g***. The more of these additional tools you install, the more you are able to maximize the usage and functionalities that **CHAPERON*****g*** offers.</p></font>

{{< admonition note "Existing dependencies do not need to be re-installed" >}}
If you have been using GROMACS before, you probably have some or most of these dependencies installed already. You do not need to re-install any active packages.
{{< /admonition >}}

<br>

***There are two options by which you can install the dependencies:***

- ***installation in conda environment (recommended), or***
- ***global (system-wide) installation.***

<br>

### 4.1. Installing dependencies with conda (recommended) {#installing-dependencies-with-conda}

{{< admonition tip "Tip" >}}
Installation in conda environment is recommended for most of the dependencies.
{{< /admonition >}}

<font><p align="justify">This approach simplifies the installation of **all** the optional dependencies (including `python-3.8`, `pymol-2.5.0`, `imagemagick`, `numpy`, `matplotlib`, `networkx-2.3`, `md-davis`, `pandas`, `scipy`, etc.), hence maximizing the usage of **CHAPERON*****g***.</p></font><br>

**Step 1**: Change to your `base` anaconda (default) environment:
```bash
conda activate base
```

<br>

**Step 2**: Run the conda environment setup script `conda_env_setup.sh` while in the **CHAPERON*****g*** folder:
```bash
chmod +x conda_env_setup.sh && ./conda_env_setup.sh
```
Note that this step requires an internet connection and could take a while to complete depending on the size of the available `RAM` and the stength of your connection. So, relax and give it time.

{{< admonition info "Troubleshooting errors" >}}
In case the above command terminates with the following error:

`Collecting package metadata (current_repodata.json): failed`<br>
`CondaHTTPError: HTTP 000 CONNECTION FAILED for url <https://conda.anaconda.org/conda-forge/linux-64/current_repodata.json>`<br>
`Elapsed: -`<br>
`An HTTP error occurred when trying to retrieve this URL.`,<br>
try the following sequence of commands instead, entering them in order, line-by-line:
`conda config --set ssl_verify false`<br>
`conda_env_setup.sh`<br>
`conda config --set ssl_verify true`<br>
{{< /admonition >}}


{{< admonition tip "Tip" >}}
The input parameters for the `conda_env_setup.sh` script are specified in the `conda_dependencies.yml` configuration file. Although the default parameters specified in the <em>`yaml`</em> file contained in the package as you have downloaded it are the most suited for **CHAPERON*****g***, you can however alter the specifications based on your other needs. (But do this only if you know what you are doing. &#x1F923; )
{{< /admonition >}}

<br>

You should see a message like the one shown below, indicating that the `chaperong` environment has been set up successfully.

<br>

{{< image src="chaperong-env.png" caption="CHAPERONg conda environment" width="800">}}
<br>

{{< admonition tip "Tip" >}}
You can confirm the newly created environment named `chaperong` with:
```bash
conda env list
```
You should see `chaperong` listed among your conda environments.

You can then change into the `chaperong` environment with:
```bash
conda activate chaperong
```
{{< /admonition >}}

{{< admonition note "Note" >}}
If after setting up the conda environment as described above, some python modules/libraries are, for some reasons, not accessible, you can further run the commands described in system-wide [Option 1](#install-system-wide-option1) below.
{{< /admonition >}}

<br>

### 4.2. Installing dependencies system-wide (optional) {#install-dependencies-system-wide}
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
