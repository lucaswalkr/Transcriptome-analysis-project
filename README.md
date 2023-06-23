-----------------------------------------------------------------------------------------
University of Bristol - MSc Bioinformatics - Group Project
by Lucy Sibbring, Yilong Yang and Lucas Walker (Group 3) 
--------------------------------------------------------------

## Overview 

Here you will find all the necessary files and guidance to follow a bioinformatic pipeline to assist with analysing transcriptomic data from the pathogen Trypanosoma congolense.

This pipeline assumes you have access to or are considering the use of expressed sequence tag (EST) data for transcriptome analysis and may want to:
1. Extract protein sequences from the EST data and convert sequence names into a human-readable format.
2. Retrieve specific protein sequences based on a unique ID.
3. Understand the output files produced.
4. Be able to repeat analyses on future EST datasets.

The EST data used for this project was obtained from the [sanger institute](https://www.sanger.ac.uk/resources/downloads/protozoa/trypanosoma-congolense.html).

## List of files 

- [Transcriptome Analysis of Trypanosoma Congolense.ipynb] : This document acts as a detailed manual for those looking to replicate the analyses, 
containing all the relevant information on the approach, a visual and step-by-step guide of the pipeline and some biological findings.

- [Images] : This folder contains all images included in the notebook document.

- [Raw input data] : This folder contains 7 fasta files with the raw EST data for each developmental stage of T. congolense, downloaded from the sanger institute website.

- [Adjusted input data] : This folder contains the fasta files with the EST data for each developmental stage of T. congolense, with all "-" characters replaced with "N".

- [Protein coding sequences] : This folder contains all the protein coding sequences (.pep fasta files), generated by TransDecoder, for each developmental stage of T. congolense.

- [Output data] : This folder contains all the renamed fasta files of the protein coding sequences, as well as an example of the output from the retrieve_seq function.

- [Tcongo_functions.py] : This python file contains all the code required for the analysis of the EST data. It is used as a module with two key functions:
    - `rename_seq()` - renames all protein sequence names within a file into a human-readable format. 
    - `retrieve_seq()` - retrieves a specific protein sequence based on a sequence ID value.
  
- [Pipeline.pdf] : This pdf contains a graphical overview of the pipeline followed for the bioinformatic analysis of the Trypanosoma congolense.

- [README.txt] : This document.

## Getting started

As the analysis will need to be performed on your local system, there are some requirements before getting started.

1. (For windows users) Install Ubuntu (v22.04) to enable work in a Linux environment.
	- To install Ubuntu you need Windows Subsystem for Linux (WSL) and Ubuntu itself.
	- To install both WSL and Ubuntu at the same time, search for and run 'Windows PowerShell' as an administrator on your system, 
          then copy and paste this code into PowerShell: `wsl --install -d ubuntu`.
	- Ubuntu is a distribution running on the Linux operating system (OS), 
          for more information and alternative methods of installing see: https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10#1-overview

2. Restart your system, run Ubuntu and create a username and password, then log in.

3. Update Ubuntu by entering `sudo apt update` into the command prompt.

4. Install Python (v3.10.6) within Ubuntu/Linux by entering `sudo apt install python3` on the command line.
	- to enable python to run on the command line you may also need to enter `sudo apt install python3-pip`.

5. Install biopython and pytest by entering these lines in Ubuntu (one at a time): `sudo pip3 install pytest` and `sudo pip3 install biopython`.

6. Install TransDecoder (v5.6.0) by entering `sudo apt install transdecoder` in Ubuntu.

7. Install Anaconda by following this link and scrolling to the bottom of the page: https://www.anaconda.com/products/distribution#Downloads
	- Once installed, you can run Anaconda and open "JupyterLab" to look at the code and tweak it if necessary.
	- (This step is not necessary if you already have an IDE such as pycharm or visual studio installed to read and write python code).

8. You're all set! Open the document named "Transcriptome Analysis of Trypanosoma Congolense" in this archive and follow on with our analysis.
