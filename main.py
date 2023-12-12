import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
from tkinter import ttk
import os
import re
import pandas as pd


AMINO_ACID_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def numbers_between(lst, lower_limit, upper_limit):
    result = [num for num in lst if lower_limit < num < upper_limit]
    return result


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
        self.file_path_input = filedialog.askopenfilename()
        if self.file_path_input == "":
            return
        if ".xlsx" in self.file_path_input:
            self.input_df = pd.read_excel(self.file_path_input)
        elif ".csv" in self.file_path_input:
            self.input_df = pd.read_csv(self.file_path_input)
        self.input_df.dropna(subset = ['Sequence'], inplace=True)

    def select_fasta_file(self):
        file_path = filedialog.askopenfilename()
        if file_path == "":
            return
        with open(file_path) as fasta:
            self.fasta = fasta.read()
        self.fasta = self.fasta.splitlines()[1]
        self.uniprot_id = simpledialog.askstring("Input", "please enter UNIPROT ID",parent=self.root)


    def run_analysis(self):
        # if not self.input_df:
        #     return
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
        
        df_to_save_dict = dict()
        df_cols=["location_in_FASTA", "total_peptides_covarage_w_repetitions","total_unique_peptides","total_intensity","no_modification"]
        for mod in mod_dict:
            df_cols.append("total_modification_" + mod)
            df_cols.append("percentage_" + mod)

        df_to_save = pd.DataFrame(columns = df_cols)

        for aa in selected_aa_dict:
            selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"] = len(self.input_df[self.input_df["Sequence_upper"].isin(selected_aa_dict[aa]["unique_peptides"].keys())])
            selected_aa_dict[aa]["total_unique_peptides"] = len(selected_aa_dict[aa]["unique_peptides"].keys())
            selected_aa_dict[aa]["total_intensity"] = sum(self.input_df[self.input_df["Sequence_upper"].isin(selected_aa_dict[aa]["unique_peptides"].keys())]["Intensity"].tolist())
            
            
            selected_aa_dict[aa]["no_modification"] = 0
            for mod in mod_dict:
                selected_aa_dict[aa]["total_modification_" + mod] = 0

                
                for peptide_iter in selected_aa_dict[aa]["unique_peptides"].keys():
                    # C# - start(mod)
                    phrase_to_search = aa[0:1] + str( 1+ int(aa[2:]) - selected_aa_dict[aa]["unique_peptides"][peptide_iter][0]) + "(" + mod + ")"
                    selected_aa_dict[aa]["total_modification_" + mod] += len(self.input_df[(self.input_df["Sequence_upper"] == peptide_iter) & (self.input_df["Modifications"].str.contains(phrase_to_search,na=False,regex=False))])
                    selected_aa_dict[aa]["no_modification"] += len(self.input_df[(self.input_df["Sequence_upper"] == peptide_iter) & (self.input_df["Modifications"].isna())])

            dict_to_append = {"location_in_FASTA": [aa], 
                "total_peptides_covarage_w_repetitions": [selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"]],
                "total_unique_peptides": [selected_aa_dict[aa]["total_unique_peptides"]],
                "total_intensity":[selected_aa_dict[aa]["total_intensity"]],
                "no_modification":[selected_aa_dict[aa]["no_modification"]]
                }  
            for mod in mod_dict:
                selected_aa_dict[aa]["percentage_" + mod] = 100*selected_aa_dict[aa]["total_modification_" + mod]/selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"]
                dict_to_append["total_modification_" + mod] = selected_aa_dict[aa]["total_modification_" + mod]
                dict_to_append["percentage_" + mod] = selected_aa_dict[aa]["percentage_" + mod]
            # order data in dict to save as a df
            
            
            entry = pd.DataFrame.from_dict(dict_to_append)
            df_to_save = pd.concat([df_to_save, entry], ignore_index=True)
            file_name, file_extension = os.path.splitext(self.file_path_input)
            df_to_save.sort_values(by=["location_in_FASTA"])
            df_to_save.to_csv(file_name+ "_output.csv", index = False)
            simpledialog.info




        


def main():
    tk_ui = UI()
    # tk_ui.select_file()
    #print(tk_ui.file_path)
    tk_ui.root.mainloop()
    #UI.root.mainloop()

if __name__ == "__main__":
    main()