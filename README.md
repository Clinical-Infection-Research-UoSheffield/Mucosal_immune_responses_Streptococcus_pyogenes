# Mucosal immune responses to *Streptococcus pyogenes* across the life course in a high burden setting: a focus on candidate vaccine antigens

Code to reproduce analyses from the SpyCATS longitudinal study of mucosal immunity to *Streptococcus pyogenes* in The Gambia.

---

## **Description**

This repository contains the code to reproduce analyses from the manuscript entitled:  
**"Mucosal immune responses to *Streptococcus pyogenes* across the life course in a high burden setting: a focus on candidate vaccine antigens"**

1. **Data Anonymization:**  
    Data including antibody measurements, clinical events, and demographics have been anonymized to protect participant confidentiality.

2. **Analytical Methods:**  
    Analyses include mixed-effects logistic regression to model event probabilities by IgA levels, correlation matrices, Principal Component Analysis (PCA), with results formatted into comprehensive tables and plots.

The analyses contained within extend the analyses of systemic blood IgG responses which can be found at 
[[https://github.com/Clinical-Infection-Research-UoSheffield/Development_natural_protective_immunity_Streptococcus_pyogenes/tree/main](https://github.com/Clinical-Infection-Research-UoSheffield/Development_natural_protective_immunity_Streptococcus_pyogenes/tree/main)  


---

## **Data Anonymisation**

- IDs are anonymized by mapping them to randomly generated codes.
- All dates are uniformly offset by a constant to preserve time intervals while concealing actual dates.
- Exact ages are replaced with pseudo-ages randomly generated within defined age groups.

Upon publication, data to reproduce analyses will be made publicly available, hosted on Zenodo

---

## **Repository Structure**

```
.
├── README.md        # Project documentation  
├── data             # Contains processed data files (empty initially)  
├── R_output         # Output from R scripts (figures, summaries)  
└── scripts          # R scripts for data analysis and visualization  
```
- **data:** Contains processed data files.  
- **R_output:** Contains output from R scripts, such as figures and summaries.  
- **scripts:** Contains R scripts used for data analysis and visualization.  

---

## **Prerequisites**

Ensure you have R installed on your system. You will also need the librarian package version (1.8.1) to manage dependencies.

- **R Version:** Requires **R version 4.4.0** or later.  
- **RStudio:** Version 2024.04.2+764 (2024.04.2+764) of later   
- **Operating System:** The code has been tested on:  
  - macOS 14.5  

Typical installation and run time is several minutes. 

---

## **Installing Dependencies**

The following R packages are required for the analyses. The `librarian` package will ensure all dependencies are installed and loaded without needing to install each package individually:

- `broom (1.0.6)`  
- `broom.mixed (0.2.9.5)`  
- `CorrMixed (1.1)`  
- `corrplot (0.92)`  
- `cowplot (1.1.3)` 
- `devtools (2.4.5)`
- `dplyr (1.1.4)`  
- `dunn.test (1.3.6)`  
- `flextable (0.9.6)`  
- `forcats (1.0.0)`  
- `FSA (0.9.5)`  
- `ggdist (3.3.2)`  
- `ggplot2 (3.5.1)`  
- `ggpubr (0.6.0)`  
- `grid (4.4.0)`  
- `gridExtra (2.3)`  
- `gtsummary (2.0.0)` 
- `inborutils (0.4.0)`
- `librarian (1.8.1)`  
- `lme4 (1.1-35.5)`  
- `lmerTest (3.1-3)`  
- `mfp (1.5.4.1)`  
- `officer (0.6.6)`  
- `patchwork (1.2.0)`  
- `pheatmap (1.0.12)`  
- `psych (2.4.6.26)`  
- `splines (4.4.0)`  
- `survival (3.7-0)`  
- `tidyr (1.3.1)`  
- `tidyverse (2.0.0)`  
- `UpSetR (1.4.0)`  
- `wesanderson (0.3.7)`  

---

## **Instructions for Use**

### **Import Data:**

There are two ways to import the data into this repository:

#### 🔁 **Option 1: Run the automated script**

Run the `load_data.R` script to automatically download and extract the anonymized dataset from the Zenodo depository.

`source("load_data.R")`

> The data has been made publicly available upon publication via DOI: <insert at publication>

#### 🛠 **Option 2: Manual download (fewer dependencies required)**

1. Go to the dataset's Zenodo page: <insert at publication>
2. Click the **“Download all”** button to download the full archive (e.g., `14887949.zip`).  
3. Unzip the contents.  
4. Move the extracted folder or files into a new directory in the project directory called `data/`.

Your directory structure should now look like:

```
repo/
├── data/
│   ├── file1.RDS
│   ├── file2.RDS
│   └── ...
├── scripts/
└── ...etc
```



### **Reproduce Analyses:**

Work through the scripts sequentially to reproduce the analyses presented in the manuscript. The recommended order is as follows:

```
- `01_oral_fluid_baseline.R`  
- `02_etc

```

Analyses should be viewed directly within **RStudio**. The output will be displayed through plots, tables, and console outputs within the RStudio environment, or outputted to a directory called "R_output"

---

### **Review Results:**

To review results, load the scripts in **RStudio** and run them interactively. The figures, tables, and statistical summaries will be displayed directly in the **Plots** and **Viewer** panels or within the **Console**. 

---

## **License**

This repository is licensed under the **Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)** license.

### **You are free to:**

- **Share:** Copy and redistribute the material in any medium or format.

### **Under the following terms:**

- **Attribution:** You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.  
- **NonCommercial:** You may not use the material for commercial purposes.  
- **NoDerivatives:** If you remix, transform, or build upon the material, you may not distribute the modified version.  

For full license details, see:  
**[Creative Commons License](https://creativecommons.org/licenses/by-nc-nd/4.0/)**  
