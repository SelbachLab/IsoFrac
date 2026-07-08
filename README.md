# IsoFrac

**IsoFrac** is a collection of tools for protein fractionation-based isoform proteomics. The repository provides an implementation of **Peptide Correlation Profiling (PCP)**, a strategy for extracting protein- and isoform-level information from peptide-centric (bottom-up) proteomics data.

## Repository Structure

### Proteomics

Code and workflows for the analysis of bottom-up proteomics data.

#### Analysis_Pipeline

This directory contains cleaned and documented versions of the original analysis scripts used in the Kedan *et al.* manuscript:

> **An Integrated Landscape of mRNA and Protein Isoforms**

The scripts are provided for transparency and reproducibility. They are **not intended to run out of the box**, as users will need to adapt file paths, working directories, and references to external datasets according to their own computing environment.

**Input:** MaxQuant `evidence.txt` files.

#### IsoFrac_Wrapper_Function

A simplified wrapper that streamlines the processing of peptide correlation profiling data, providing an accessible starting point for IsoFrac analyses.

#### Additional_Data

Contains preprocessed reference files and auxiliary datasets required by several analysis scripts.

#### LOPDIA

Scripts for preprocessing and preparing the LOPDIA dataset. These scripts require input files from the `Additional_Data` directory.

### Transcriptomics

Code for the processing and analysis of long-read transcriptomic data generated using **PacBio** sequencing.

## Disclaimer

The scripts in this repository are intended for **research use only**. They are provided as a reference implementation and starting point for peptide correlation profiling workflows. Depending on the experimental design, instrumentation, software versions, and biological application, users are expected to adapt the code to their specific requirements.

The repository represents an actively evolving implementation of the IsoFrac workflow and may differ from the exact analysis pipeline used in the associated publication. Users are responsible for validating all results generated with these scripts.