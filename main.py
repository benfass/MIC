import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import re

import pandas as pd


AMINO_ACID_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


class UI:
    def __init__(self):
        self.root = tk.Tk()
        self.menubar = tk.Menu(self.root)
        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Import Input file", command=self.select_file)
        self.filemenu.add_command(label="add FASTA", command=self.select_fasta_file)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.root.config(menu=self.menubar)

        self.instruction_frame = tk.Frame(self.root)
        self.instruction_frame.grid(row=1,column=0)
        ttk.Label(self.instruction_frame,text="Please select AA to be checked for modifications:").grid(row=0,column=0)
        self.aa_selection = tk.StringVar()
        tk.OptionMenu(self.instruction_frame, self.aa_selection , *AMINO_ACID_LIST).grid(row=0,column=1)
        ttk.Button(self.instruction_frame,text="Run analysis", command = self.run_analysis).grid(row=1,column=0)

    
    def select_file(self):
        file_path = filedialog.askopenfilename()
        if file_path == "":
            return
        if ".xlsx" in file_path:
            self.input_df = pd.read_excel(file_path)
        elif ".csv" in file_path:
            self.input_df = pd.read_csv(file_path)
        self.input_df.dropna(subset = ['Sequence'], inplace=True)

    def select_fasta_file(self):
        file_path = filedialog.askopenfilename()
        if file_path == "":
            return
        with open(file_path) as fasta:
            self.fasta = fasta.read()
        self.fasta = self.fasta.splitlines()[1]


    def run_analysis(self):
        # if not self.input_df:
        #     return
        
        self.input_df['Sequence'] =  self.input_df['Sequence'].str.upper()
        self.filtered_df = self.input_df[self.input_df["Sequence"].str.contains(self.aa_selection.get())]
        self.sum_of_inteseties = self.filtered_df["Intensity"].sum()
        self.filtered_df["Modifications"].unique().tolist()
        mod_dict = dict()
        pattern = r'\((.*?)\)'
        for value in self.filtered_df["Modifications"].unique().tolist():
            try:
                if value.strip().startswith(self.aa_selection.get()):
                    mod_dict[re.findall(pattern, value)[0]] = 1
            except:
                pass
        modification_df = self.filtered_df
        modification_df.dropna(subset = ['Modifications'], inplace = True)
        for key in mod_dict:
            mod_dict[key] = modification_df[modification_df["Modifications"].str.contains(key)]["Intensity"].sum()*100/self.sum_of_inteseties
        print(mod_dict)

        


def main():
    tk_ui = UI()
    # tk_ui.select_file()
    #print(tk_ui.file_path)
    tk_ui.root.mainloop()
    #UI.root.mainloop()

if __name__ == "__main__":
    main()