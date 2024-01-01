# ProteoMIC
Proteomic Modification Intensity Calculator

This Python script provides a simple user interface to perform proteome analysis based on the output file from Proteome Discoverer/Metamorpheus. It allows users to select a Proteome Discoverer/Metamorpheus output file, enter a UniProt ID, and analyze amino acid modifications.

# Requirements
Python 3.x
Required Python libraries: tkinter, pandas, numpy, requests


# How to use the script
download from Release tab (https://github.com/benfass/ProteoMIC/releases/tag/v2.0) the latest EXE file run the sorce code using python IDE/cmd propmt

# GUI instructions
1) Select Proteome Discoverer/Metamorpheus output file to analyse
2) Enter UniProt ID and press "Get FASTA file" - The button will turn green if the FASTA file is successfully loaded.
3) (optional step) - select folder for output files(by default output files will be created in input file folder)
4) Select Amino Acid (AA) to Analyze for Modifications
5) Click Run analysis


# Output:

The script will generate an output CSV file with detailed analysis results, including information about peptide modifications, amino acid coverage, and intensity percentages.

Note: Ensure that the required Python libraries are installed before running the script. You can install the necessary libraries using the following command:

bash
pip install pandas numpy requests

Feel free to customize and adapt the script to suit your specific needs or integrate it into your workflow.
