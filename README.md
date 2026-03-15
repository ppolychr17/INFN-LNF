# $B^0 \to D^+ \mu^- \bar{\nu}_\mu$ Exclusive Semileptonic Analysis 

This repository contains the analysis scripts for the exclusive semileptonic decay $B^0 \to D^+ \mu^- \bar{\nu}_\mu$. The analysis explicitly accounts for feed-down contributions from excited charm resonances, such as D* and D** states, that decay into the reconstructed D+ final state.

The primary goal of this project is to fit current MC simulations to real data to accurately describe the data from Run 3 2024. To achieve this, the codebase is designed to test and compare several different methodologies and statistical tools.


## References

The statistical tools and analysis methodologies implemented in this repository rely heavily on the following foundational works:

* **Background Extraction (sWeights):** M. Pivk and F. R. Le Diberder, *"sPlot: A Statistical tool to unfold data distributions"* [arXiv:physics/0402083](https://arxiv.org/abs/physics/0402083)

* **Template Fitting with Finite MC (Beeston-Barlow):** R. Barlow and C. Beeston, *"Fitting using finite Monte Carlo samples"* [DOI:10.1016/0010-4655(93)90005-W](https://doi.org/10.1016/0010-4655(93)90005-W)

* **Approximate Maximum-Likelihood Method:** H. P. Dembinski and A. Abdelmotteleb, *"A fast and stable approximate maximum-likelihood method for template fits"* [arXiv:2206.12346](https://arxiv.org/abs/2206.12346)

* **Python Fitting Library (iminuit / Scikit-HEP):** E. Rodrigues, H. Dembinski, et al., *"The Scikit HEP Project -- overview and prospects"* [arXiv:2007.03577](https://arxiv.org/abs/2007.03577)

* **Fitting Framework (RooFit):** W. Verkerke and D. Kirkby, *"The RooFit toolkit for data modeling"* [arXiv:physics/0306116](https://arxiv.org/abs/physics/0306116)


## Analysis Features & Methodologies

This codebase allows for the comparison of:
* **Fit Types:** Binned vs. unbinned fits
* **Background Extraction:** Sideband extrapolation vs. sWeights extraction
* **Fitting Tools:** ROOT's RooFit framework (C++) vs. the `iminuit` library (Python)

**Additional checks included:**
* Kinematic cut checks
* Estimation of the MC statistics required to reach specific target accuracies

## Repository Structure

To ensure the scripts run correctly and output files are saved in the appropriate locations, your working directory should follow this structure:

    INFN/
    ├── scripts/                  # main_MC.cpp, roofit-fit.cpp, iminuit-fit.ipynb
    ├── root_files/          
    │   ├── MC/                   # Contains input MC root files
    │   │   └── templates_MC/     # Used for storing different MC signal histograms (final signal shapes)
    │   └── Data/                 # Contains input Data root files
    └── plots/                    # Output directory for generated fit plots

*Note: You can change the structure as you please, changing the file names on top of each file. Detailed explanations of specific code logic, cut and variable definitions are documented directly in the script files.*

## Prerequisites

### 1. Software Dependencies
* **C++ / ROOT:** Required for compiling macros and running RooFit
* **Python:** Required for Jupyter notebooks and the `iminuit` fitting library

### 2. Required Input Files

Before running the analysis, you should have these files:
* `MC_eff.root`: This file contains the MC data with the PID efficiencies already included as branches 
* `file_XXXXXXX_XXXXX`: These files contain the raw data without any cuts as acquired online

---

## Workflow Instructions

Follow these steps in order to reproduce the analysis:

### Step 1: Process MC and Apply Cuts
* Run `main_MC.cpp` (located in the `scripts/` directory)
* This script will read `MC_eff.root` and apply the necessary kinematic and PID selection
* The output of this script will be the final MC `.root` file that is ready for fitting

### Step 2: Merge Data and Extract sWeights
* Next, you need to merge the raw data (`file_XXXXXXX_XXXXX`), using the `roofit-fit.cpp` script where preselection and cuts are also applied
* After this, you can use the merged data with cuts applied for the rest of the analysis
* In order to use sWeights, you need to extract them using the `roofit-fit.cpp` script, which is found in the end of the script. Make sure to have correct file names
* **It is highly recommended to include the sWeights in the final .root file you use before performing any fits!**

### Step 3: Perform the Fits
Once the background extraction/sWeights step is complete, you can perform the final fits. In the `scripts/` folder, you can choose between two fitting frameworks:
* `iminuit-fit.ipynb` (Python / Jupyter Notebook)
* `roofit-fit.cpp` (ROOT / C++ Macro)

*Note: If you need to change the names of the input/output files or the variables being fitted, you can easily modify them in the configuration blocks at the very beginning of the script files.*

