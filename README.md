# ProteoMIC
Proteomic Modification Intensity Calculator

This Python script provides a simple user interface to perform proteome analysis based on the output file from Proteome Discoverer. It allows users to select a Proteome Discoverer output file, enter a UniProt ID, and analyze amino acid modifications within the specified range.

# Directions to Usage
download from Release tab (https://github.com/benfass/ProteoMIC/releases/tag/v1.0) the latest EXE file or follow instructions bellow to use as python script

# Requirements
Python 3.x
Required Python libraries: tkinter, pandas, numpy, requests
How to Use
Run the Script: Execute the script using Python. The graphical user interface (GUI) will appear.

bash
python proteome_analysis_tool.py
Select Proteome Discoverer Output File:

Click the "Select File" button to choose the Proteome Discoverer output file (either in Excel or CSV format).
The selected file will be displayed, and the button will turn green.
Enter UniProt ID:

Enter the UniProt ID in the provided entry field.
Click the "Get FASTA file" button to fetch the corresponding FASTA file from UniProt.
The button will turn green if the FASTA file is successfully loaded.
Select Amino Acid (AA) to Analyze for Modifications:

Choose the amino acid to be analyzed for modifications from the dropdown menu.
Click the "Run Analysis" button to perform the analysis.

# Output:

The script will generate an output CSV file with detailed analysis results, including information about peptide modifications, amino acid coverage, and intensity percentages.
The output file will be saved in the same directory as the input Proteome Discoverer file.
View Results:

A message box will appear indicating that the output file has been saved successfully.
Open the generated CSV file to view the detailed analysis results.
Note: Ensure that the required Python libraries are installed before running the script. You can install the necessary libraries using the following command:

bash
pip install pandas numpy requests

Feel free to customize and adapt the script to suit your specific needs or integrate it into your workflow.
