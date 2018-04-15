# RoCA


## Table of Contents
*  [Overview](#overview)
*  [Details](#details)
*  [Requirements](#requirements)
*  [Usage](#usage)
*  [Troubleshooting](#troubleshooting)

## Overview
Software for identifying co-evolutionary sectors in proteins using "Robust Co-evolutional Analysis (RoCA)"

## Details
#### Title of paper
Co-evolution networks of HIV/HCV are modular with direct association to structure and function
#### Authors
Ahmed A. Quadeer, David Morales-Jimenez, and Matthew R. McKay

## Requirements
1. A PC with MATLAB (preferrably v2017a or later) installed on it with the following additional toolboxes:
    * Bioinformatics Toolbox
    * Statistics and Machine Learning Toolbox
    
2. For running codes related to statistical coupling analysis (SCA), register and download the SCA software from https://ais.swmed.edu/rrlabs/register.htm
 
3. For mapping predicted sector residues on crystal structures, download Pymol available at https://pymol.org/ 

## Usage
* Inferring co-evolutionary networks for a protein using RoCA
   * Open MATLAB
   * Run the script main_RoCA.m and provide the MSA matrix as an input

* Reproducing results in the paper for HIV and HCV viral proteins
   * Run the following scripts to generate RoCA (and PCA [Quadeer et al. 2014]) results
      * main_gag.m for HIV Gag
      * main_nef.m for HIV Nef
      * main_ns34a.m for HCV NS3-4A
      * main_ns4b.m for HCV NS4B
      
   * Run the following scripts to generate SCA results
      * main_gag_sca.m for HIV Gag
      * main_nef_sca.m for HIV Nef
      * main_ns34a_sca.m for HCV NS3-4A
      * main_ns4b_sca.m for HCV NS4B
      
   * Run the following script (in the GT folder) to compare the performance of RoCA and PCA using binary synthetic data
      * main_GT.m
      
* To visualizing the step-by-step procedure and the corresponding output
   * Download the html folder
   * Open the "main.html" file in your browser

---
[Quadeer et al. 2014] Quadeer AA, Louie RHY, Shekhar K, Chakraborty AK, Hsing I-M, McKay MR. 2014. Statistical linkage analysis of substitutions in patient-derived sequences of genotype 1a hepatitis C virus non-structural protein 3 exposes targets for immunogen design. J. Virol. 88:7628–44. doi:10.1128/JVI.03812-13.

## Troubleshooting
For any questions or comments, please email at ahmedaq@gmail.com. 
