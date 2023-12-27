import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
from tkinter  import messagebox
from tkinter import ttk
import os
import re
import pandas as pd
import numpy as np
import requests

# TODO:
# in case FASTA fails to load from uniprot - add option to attach FASTA manually


AMINO_ACID_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def fetch_fasta(uniprot_id):
    base_url = "https://www.uniprot.org/uniprot/"
    query_url = f"{base_url}{uniprot_id}.fasta"

    response = requests.get(query_url)

    if response.ok:
        fasta_data = response.text
        return parse_fasta(fasta_data)
    else:
        print(f"Failed to retrieve FASTA for {uniprot_id}. Status code: {response.status_code}")
        return None
    
def parse_fasta(fasta_data):
    metadata = None
    sequence = ""

    lines = fasta_data.splitlines()
    
    if lines and lines[0].startswith(">"):
        metadata = lines[0][1:]  # Exclude the ">" symbol from the header
        sequence = "".join(lines[1:])  # Concatenate sequence lines

    return metadata, sequence


def numbers_between(lst, lower_limit, upper_limit):
    result = [num for num in lst if lower_limit < num < upper_limit]
    return result


class UI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.minsize(600,200)

        self.fasta = ""
        self.input_df = pd.DataFrame()

        self.instruction_frame = tk.Frame(self.root)
        self.instruction_frame.grid(sticky = "W",row=1,column=0, padx=10)
        self.input_prog_var = tk.StringVar()
        self.input_prog_frame = tk.Frame(self.instruction_frame)
        self.input_prog_frame.grid(row=2, column=0, columnspan= 3, sticky = "W")

        ttk.Label(self.input_prog_frame, text= "Please select program used to create input file: ").grid(row=0, column=0, padx=2)
        ttk.Radiobutton(self.input_prog_frame,text="Proteom discoverer", variable=self.input_prog_var, value= "Proteom discoverer", command=self.chenge_input_prog).grid(sticky = "W",row=0,column=1, padx=2)
        ttk.Radiobutton(self.input_prog_frame,text="Metamorpheus", variable=self.input_prog_var, value= "Metamorpheus", command=self.chenge_input_prog).grid(sticky = "W",row=0,column=2, padx=2)


        ttk.Label(self.instruction_frame,text="Please select Proteome Discoverer output file:").grid(sticky = "W",row=3,column=0)
        self.select_file_button_text = tk.StringVar()
        self.select_file_button_text.set("please select file")
        self.file_select_button = tk.Button(self.instruction_frame,textvariable=self.select_file_button_text,command=self.select_file)
        self.file_select_button.grid(sticky = "W",row=3,column=2)

        ttk.Label(self.instruction_frame,text="Please enter UNIPROT ID:").grid(sticky = "W",row=4,column=0)
        self.select_FASTA_result_text = tk.StringVar()
        self.select_FASTA_button_text = tk.StringVar()
        self.select_FASTA_button_text.set("Get FASTA file")

        tk.Entry(self.instruction_frame, textvariable=self.select_FASTA_result_text).grid(sticky = "W",row=4,column=1)
        self.get_fasta_button = tk.Button(self.instruction_frame,textvariable=self.select_FASTA_button_text ,command=self.select_fasta_file)
        self.get_fasta_button.grid(sticky = "W",row=4,column=2)

        ttk.Label(self.instruction_frame,text="Please select AA to be checked for modifications:").grid(sticky = "W",row=6,column=0)
        self.aa_selection = tk.StringVar()
        tk.OptionMenu(self.instruction_frame, self.aa_selection , *AMINO_ACID_LIST).grid(sticky = "W",row=6,column=2)
        tk.Button(self.instruction_frame,text="Run analysis", command = self.run_analysis).grid(sticky = "W",row=7,column=0)

    def chenge_input_prog(self):
        selection = self.input_prog_var.get()
        if selection == "Proteom discoverer":
            self.input_prog_selection = "Proteom discoverer"


    def select_file(self):
        self.file_path_input = filedialog.askopenfilename()
        if self.file_path_input == "":
            return
        if ".xlsx" in self.file_path_input:
            self.input_df = pd.read_excel(self.file_path_input)
        elif ".csv" in self.file_path_input:
            self.input_df = pd.read_csv(self.file_path_input)
        self.input_df.dropna(subset = ['Sequence'], inplace=True)
        self.select_file_button_text.set(self.file_path_input)
        self.file_select_button.configure(bg="green")

    def select_fasta_file(self):
        self.uniprot_id = self.select_FASTA_result_text.get()
        try:
            self.fasta = fetch_fasta(self.select_FASTA_result_text.get().strip())[1]
        except:
            self.select_FASTA_button_text.set("Failed to load FASTA!")
            self.get_fasta_button.configure(bg="red")
        if self.fasta:
            self.select_FASTA_button_text.set("FASTA file loaded sucssesfully")
            self.get_fasta_button.configure(bg="green")
        else:
            self.select_FASTA_button_text.set("Failed to load FASTA!")
            self.get_fasta_button.configure(bg="red")


    def run_analysis(self):
        if not self.fasta or self.input_df.empty or self.aa_selection.get() == "":
            messagebox.showerror(title="error running analisys", message="please fill all parameters and try again")
            return
        self.input_df = self.input_df[self.input_df["Protein Group Accessions"] == self.uniprot_id]
        self.aa_to_analyise = self.aa_selection.get()
        self.input_df['Sequence_upper'] =  self.input_df['Sequence'].str.upper()
        self.filtered_df = self.input_df[self.input_df["Sequence_upper"].str.contains(self.aa_to_analyise)]
        self.sum_of_inteseties = self.filtered_df["Intensity"].sum()
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
            mod_dict[key] = modification_df[modification_df["Modifications"].str.contains(key)]["Intensity"].sum()*(100/self.sum_of_inteseties)
        
        # total number of peptides found
        self.num_peptides_total = len(self.input_df['Sequence'].unique().tolist())
        # total number of peptides containing selected AA
        self.num_peptides_aa_selection = len(self.filtered_df['Sequence'].unique().tolist())
        
        # AA analysis
        selected_aa_dict = dict()
        self.selected_aa_locations = [m.start() for m in re.finditer(self.aa_to_analyise, self.fasta)]
        for peptide in self.filtered_df['Sequence_upper'].unique().tolist():
            start = [re.search(peptide,self.fasta)][0].span()[0]
            end = [re.search(peptide,self.fasta)][0].span()[1]
            match_list = numbers_between(self.selected_aa_locations, start, end)
            for match in match_list:
                if self.aa_to_analyise + "_" +str(match) in selected_aa_dict:
                    selected_aa_dict[self.aa_to_analyise + "_" +str(match)]["unique_peptides"][peptide] = (start,end)
                else:
                    selected_aa_dict[self.aa_to_analyise + "_" +str(match)] = {"unique_peptides": {peptide:(start,end)}}
        
        df_cols=["AA number","location_in_FASTA", "total_peptides","total_unique_peptides","total_intensity","no_modification"]
        for mod in mod_dict:
            df_cols.append("peptides_with_" + mod)
            df_cols.append("percentage_" + mod + "%")
            df_cols.append("total_intensity_" + mod)
            df_cols.append("intensity_pecentage_" + mod)

        df_to_save = pd.DataFrame(columns = df_cols)

        for aa in selected_aa_dict:
            selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"] = len(self.input_df[self.input_df["Sequence_upper"].isin(selected_aa_dict[aa]["unique_peptides"].keys())])
            selected_aa_dict[aa]["total_unique_peptides"] = len(selected_aa_dict[aa]["unique_peptides"].keys())
            selected_aa_dict[aa]["total_intensity"] = sum(self.input_df[self.input_df["Sequence_upper"].isin(selected_aa_dict[aa]["unique_peptides"].keys())]["Intensity"].tolist())
            
            
            selected_aa_dict[aa]["no_modification"] = 0
            total_mod_amount = 0
            for mod in mod_dict:
                mod_intensity = 0
                selected_aa_dict[aa]["total_modification_" + mod] = 0
                for peptide_iter in selected_aa_dict[aa]["unique_peptides"].keys():
                    phrase_to_search = aa[0:1] + str( 1+ int(aa[2:]) - selected_aa_dict[aa]["unique_peptides"][peptide_iter][0]) + "(" + mod + ")"
                    selected_aa_dict[aa]["total_modification_" + mod] += len(self.input_df[(self.input_df["Sequence_upper"] == peptide_iter) & (self.input_df["Modifications"].str.contains(phrase_to_search,na=False,regex=False))])
                    mod_intensity += sum(self.input_df[(self.input_df["Sequence_upper"] == peptide_iter) & (self.input_df["Modifications"].str.contains(phrase_to_search,na=False,regex=False))]["Intensity"])
                total_mod_amount += selected_aa_dict[aa]["total_modification_" + mod]
                selected_aa_dict[aa]["total_intensity_" + mod] = mod_intensity
            selected_aa_dict[aa]["no_modification"] = selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"] - total_mod_amount



            dict_to_append = {"location_in_FASTA": [aa], 
                "total_peptides": [selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"]],
                "total_unique_peptides": [selected_aa_dict[aa]["total_unique_peptides"]],
                "total_intensity":[selected_aa_dict[aa]["total_intensity"]],
                "no_modification":[selected_aa_dict[aa]["no_modification"]]
                }  
            for mod in mod_dict:
                dict_to_append["AA number"] = int(aa[2:])
                selected_aa_dict[aa]["percentage_" + mod] = 100*selected_aa_dict[aa]["total_modification_" + mod]/selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"]
                dict_to_append["peptides_with_" + mod] = selected_aa_dict[aa]["total_modification_" + mod]
                dict_to_append["percentage_" + mod + "%"] = selected_aa_dict[aa]["percentage_" + mod]
                dict_to_append["total_intensity_" + mod] = selected_aa_dict[aa]["total_intensity_" + mod]
                dict_to_append["intensity_pecentage_" + mod] = 100*(dict_to_append["total_intensity_" + mod]/selected_aa_dict[aa]["total_intensity"])
            # order data in dict to save as a df
            
            
            entry = pd.DataFrame.from_dict(dict_to_append)
            df_to_save = pd.concat([df_to_save, entry], ignore_index=True)
        file_name, file_extension = os.path.splitext(self.file_path_input)
        df_to_save.sort_values(by=["AA number"],inplace=True)
        df_to_save.to_csv(file_name+ "_output.csv", index = False)
        messagebox.showinfo(title="file saved!", message="output file saved sucsessfully at:\n{0}".format(file_name+ "_output.csv"))




        


def main():
    tk_ui = UI()
    # tk_ui.select_file()
    #print(tk_ui.file_path)
    tk_ui.root.mainloop()
    #UI.root.mainloop()

if __name__ == "__main__":
    main()