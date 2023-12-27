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

def handle_devide_by_zero(number, devider):
    try:
        return round(number/devider*100,2)
    except ZeroDivisionError:
        return "NA"

def numbers_between(lst, lower_limit, upper_limit):
    """
    Filters a list of numbers and returns only those within a specified range.

    Parameters:
    - lst (list of numeric): The input list of numbers.
    - lower_limit (numeric): The lower limit of the range (exclusive).
    - upper_limit (numeric): The upper limit of the range (exclusive).

    Returns:
    list of numeric: A new list containing only the numbers that fall within the specified range.

    Example:
    >>> numbers_between([1, 5, 10, 15, 20], 5, 15)
    [10]
    """
    result = [num for num in lst if lower_limit < num < upper_limit]
    return result


class UI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.minsize(600,200)

        self.fasta = ""
        self.file_output_path = ""
        self.input_df = pd.DataFrame()

        self.instruction_frame = tk.Frame(self.root)
        self.instruction_frame.grid(sticky = "W",row=1,column=0, padx=10)
        self.input_prog_var = tk.StringVar()
        self.input_prog_var.set("Proteom discoverer")
        self.input_prog_selection = "Proteom discoverer"
        self.input_prog_frame = tk.Frame(self.instruction_frame)
        self.input_prog_frame.grid(row=2, column=0, columnspan= 3, sticky = "W")

        ttk.Label(self.input_prog_frame, text= "Select program used to create input file: ").grid(row=0, column=0,pady=5)
        ttk.Radiobutton(self.input_prog_frame,text="Proteom discoverer", variable=self.input_prog_var, value= "Proteom discoverer", command=self.chenge_input_prog).grid(sticky = "W",row=0,column=1, padx=2)
        ttk.Radiobutton(self.input_prog_frame,text="Metamorpheus", variable=self.input_prog_var, value= "Metamorpheus", command=self.chenge_input_prog).grid(sticky = "W",row=0,column=2, padx=2)


        ttk.Label(self.instruction_frame,text="Select Proteome Discoverer/Metamorpheus output file:").grid(sticky = "W",row=3,column=0,pady=5)
        self.select_file_button_text = tk.StringVar()
        self.select_file_button_text.set("Select file")
        self.file_select_button = tk.Button(self.instruction_frame,textvariable=self.select_file_button_text,command=self.select_file)
        self.file_select_button.grid(sticky = "W",row=3,column=2,padx=5)

        ttk.Label(self.instruction_frame,text="Enter UNIPROT ID:").grid(sticky = "W",row=4,column=0,pady=5)
        self.select_FASTA_result_text = tk.StringVar()
        self.select_FASTA_button_text = tk.StringVar()
        self.select_FASTA_button_text.set("Get FASTA file")

        tk.Entry(self.instruction_frame, textvariable=self.select_FASTA_result_text).grid(sticky = "W",row=4,column=1)
        self.get_fasta_button = tk.Button(self.instruction_frame,textvariable=self.select_FASTA_button_text ,command=self.select_fasta_file)
        self.get_fasta_button.grid(sticky = "W",row=4,column=2, padx=5)

        ttk.Label(self.instruction_frame,text="Select AA to be checked for modifications:").grid(sticky = "W",row=6,column=0,pady=5)
        self.aa_selection = tk.StringVar()
        tk.OptionMenu(self.instruction_frame, self.aa_selection , *AMINO_ACID_LIST).grid(sticky = "W",row=6,column=2,padx=3)

        self.select_output_file_button_text = tk.StringVar()
        self.select_output_file_button_text.set("Select file")
        ttk.Label(self.instruction_frame,text="Select location for output file(optional)").grid(sticky = "W",row=7,column=0,pady=5)
        self.select_output_path_button = tk.Button(self.instruction_frame,textvariable=self.select_output_file_button_text ,command=self.select_output_path)
        self.select_output_path_button.grid(sticky = "W",row=7,column=2,padx=5)

        tk.Button(self.instruction_frame,text="Run analysis", command = self.run_analysis_per_program).grid(sticky = "W",row=8,column=0)

    def select_output_path(self):
        self.file_output_path = filedialog.askdirectory()
        self.select_output_file_button_text.set(self.file_output_path)
        self.select_output_path_button.configure(bg="green")
        pass

    def chenge_input_prog(self):
        selection = self.input_prog_var.get()
        if selection == "Proteom discoverer":
            self.input_prog_selection = "Proteom discoverer"
        elif selection == "Metamorpheus":
            self.input_prog_selection = "Metamorpheus"


    def select_file(self):
        self.file_path_input = filedialog.askopenfilename()
        if self.file_path_input == "":
            return
        if ".xlsx" in self.file_path_input:
            self.input_df = pd.read_excel(self.file_path_input)
        elif ".csv" in self.file_path_input:
            self.input_df = pd.read_csv(self.file_path_input)
        elif ".tsv" in self.file_path_input:
            self.input_df = pd.read_csv(self.file_path_input, sep='\t', header=0)
        if self.input_prog_selection == "Metamorpheus":
            try:
                self.metamorpheus_input_df_prework()
            except KeyError:
                messagebox.showerror(title="bad input file", message="please check that file output matches the selected program \nand that no changes made to file post creation")
                return
        try:
            self.input_df.dropna(subset = ['Sequence'], inplace=True)
        except KeyError:
            messagebox.showerror(title="bad input file", message="please check that file output matches the selected program \nand that no changes made to file post creation")
            return
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

    def extract_sequence_and_modifications(self, value):
        # Use regular expressions to find matches
        # sequence_match = re.search(r'([A-Za-z]+)', value)
        

        # Join the sequence and modifications
        # sequence = sequence_match.group(1) if sequence_match else ""
        sequence = value
        modifications = []
        legth_to_substrunct = 0
        modification_matches = re.finditer(r'\[([^\]]+)\]', value)
        if modification_matches:
            for modification_info in modification_matches:
                sequence = sequence.replace(modification_info.group(),"")
                location = modification_info.start()
                mod = modification_info.group()[modification_info.group().find(r':')+1:modification_info.group().find(r' on ')]
                AA = modification_info.group()[modification_info.group().find(r' on ')+4:modification_info.group().find(r']')]
                modifications.append(AA+str(location - legth_to_substrunct)+"("+mod+")")
                legth_to_substrunct = legth_to_substrunct + modification_info.span()[1]-modification_info.span()[0]
                continue
                

        modifications_str = '; '.join(modifications) if modifications else ""

        return sequence, modifications_str

    def metamorpheus_input_df_prework(self):
        if self.input_prog_selection != "Metamorpheus":
            return False
        self.input_df.rename(columns={"Protein Group" : "Protein Group Accessions"},inplace=True)
        self.input_df.rename(columns={"Peak intensity" : "Intensity"},inplace=True)
        
        self.input_df[['Sequence', 'Modifications']] = self.input_df['Full Sequence'].apply(self.extract_sequence_and_modifications).apply(pd.Series)
        pass

    def run_analysis_per_program(self):
        if self.input_prog_selection == "Proteom discoverer":
            output_path = self.run_analysis(self.input_df)
        elif self.input_prog_selection == "Metamorpheus":
            for sample in self.input_df["File Name"].unique().tolist():
                output_path = self.run_analysis(df_to_analyze = self.input_df[self.input_df["File Name"]==sample],sample_name = "_"+sample)
        messagebox.showinfo(title="file saved!", message="output file saved sucsessfully at:\n{0}".format(output_path))
            

    def run_analysis(self , df_to_analyze, sample_name = ""):
        if not self.fasta or df_to_analyze.empty or self.aa_selection.get() == "":
            messagebox.showerror(title="error running analisys", message="please fill all parameters and try again")
            return
        df_to_analyze = df_to_analyze[df_to_analyze["Protein Group Accessions"].str.contains(self.uniprot_id,na=False)]
        self.aa_to_analyise = self.aa_selection.get()
        df_to_analyze['Sequence_upper'] =  df_to_analyze['Sequence'].str.upper()
        filtered_df = df_to_analyze[df_to_analyze["Sequence_upper"].str.contains(self.aa_to_analyise)]
        sum_of_inteseties = filtered_df["Intensity"].sum()
        mod_dict = dict()
        pattern = r'\((.*?)\)'
        for value in filtered_df["Modifications"].unique().tolist():
            try:
                if value.strip().startswith(self.aa_selection.get()):
                    mod_dict[re.findall(pattern, value)[0]] = 1
            except:
                pass
        modification_df = filtered_df
        modification_df.dropna(subset = ['Modifications'], inplace = True)
        for key in mod_dict:
            mod_dict[key] = modification_df[modification_df["Modifications"].str.contains(key)]["Intensity"].sum()*(100/sum_of_inteseties)
        
        # # total number of peptides found
        # num_peptides_total = len(df_to_analyze['Sequence'].unique().tolist())
        # total number of peptides containing selected AA
        # num_peptides_aa_selection = len(filtered_df['Sequence'].unique().tolist())
        
        # AA analysis
        selected_aa_dict = dict()
        selected_aa_locations = [m.start() for m in re.finditer(self.aa_to_analyise, self.fasta)]
        for peptide in filtered_df['Sequence_upper'].unique().tolist():
            start = [re.search(peptide,self.fasta)][0].span()[0]
            end = [re.search(peptide,self.fasta)][0].span()[1]
            match_list = numbers_between(selected_aa_locations, start, end)
            for match in match_list:
                if self.aa_to_analyise + "_" +str(match) in selected_aa_dict:
                    selected_aa_dict[self.aa_to_analyise + "_" +str(match)]["unique_peptides"][peptide] = (start,end)
                else:
                    selected_aa_dict[self.aa_to_analyise + "_" +str(match)] = {"unique_peptides": {peptide:(start,end)}}
        
        df_cols=["AA number","location_in_FASTA", "total_peptides","total_unique_peptides","total_intensity","peptides_with_no_modification"]
        for mod in mod_dict:
            df_cols.append("peptides_with_" + mod)
            df_cols.append("percentage_" + mod + "%")
            df_cols.append("total_intensity_" + mod)
            df_cols.append("intensity_pecentage_" + mod)

        df_to_save = pd.DataFrame(columns = df_cols)

        for aa in selected_aa_dict:
            selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"] = len(df_to_analyze[df_to_analyze["Sequence_upper"].isin(selected_aa_dict[aa]["unique_peptides"].keys())])
            selected_aa_dict[aa]["total_unique_peptides"] = len(selected_aa_dict[aa]["unique_peptides"].keys())
            selected_aa_dict[aa]["total_intensity"] = sum(df_to_analyze[df_to_analyze["Sequence_upper"].isin(selected_aa_dict[aa]["unique_peptides"].keys())]["Intensity"].tolist())
            
            
            selected_aa_dict[aa]["no_modification"] = 0
            total_mod_amount = 0
            total_mod_intensity = 0
            for mod in mod_dict:
                mod_intensity = 0
                selected_aa_dict[aa]["total_modification_" + mod] = 0
                for peptide_iter in selected_aa_dict[aa]["unique_peptides"].keys():
                    phrase_to_search = aa[0:1] + str( 1+ int(aa[2:]) - selected_aa_dict[aa]["unique_peptides"][peptide_iter][0]) + "(" + mod + ")"
                    selected_aa_dict[aa]["total_modification_" + mod] += len(df_to_analyze[(df_to_analyze["Sequence_upper"] == peptide_iter) & (df_to_analyze["Modifications"].str.contains(phrase_to_search,na=False,regex=False))])
                    mod_intensity += sum(df_to_analyze[(df_to_analyze["Sequence_upper"] == peptide_iter) & (df_to_analyze["Modifications"].str.contains(phrase_to_search,na=False,regex=False))]["Intensity"])
                total_mod_amount += selected_aa_dict[aa]["total_modification_" + mod]
                selected_aa_dict[aa]["total_intensity_" + mod] = mod_intensity
                total_mod_intensity += mod_intensity
            selected_aa_dict[aa]["no_modification"] = selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"] - total_mod_amount


            dict_to_append = {"location_in_FASTA": [aa], 
                "total_peptides": [selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"]],
                "total_unique_peptides": [selected_aa_dict[aa]["total_unique_peptides"]],
                "total_intensity":[selected_aa_dict[aa]["total_intensity"]],
                "peptides_with_no_modification":[selected_aa_dict[aa]["no_modification"]],
                "no_modification_intensity" : selected_aa_dict[aa]["total_intensity"] - total_mod_intensity
                }  
            for mod in mod_dict:
                dict_to_append["AA number"] = int(aa[2:])
                selected_aa_dict[aa]["percentage_" + mod] = 100*selected_aa_dict[aa]["total_modification_" + mod]/selected_aa_dict[aa]["total_peptides_covarage_w_repetitions"]
                dict_to_append["peptides_with_" + mod] = selected_aa_dict[aa]["total_modification_" + mod]
                dict_to_append["percentage_of_" + mod + "modification[%]"] = round(selected_aa_dict[aa]["percentage_" + mod],2)
                dict_to_append["intensity_" + mod] = round(selected_aa_dict[aa]["total_intensity_" + mod],2)
                if selected_aa_dict[aa]["total_intensity"]!= 0:
                    dict_to_append["intensity_pecentage_" + mod + "[%]"] = round(100*(dict_to_append["intensity_" + mod]/selected_aa_dict[aa]["total_intensity"]),2)
                else:
                    dict_to_append["intensity_pecentage_" + mod + "[%]"] = "NA"
            
            
            entry = pd.DataFrame.from_dict(dict_to_append)
            df_to_save = pd.concat([df_to_save, entry], ignore_index=True)

        df_to_save["percentage_of_no_modifications[%]"] = df_to_save.apply(lambda x: round(x["peptides_with_no_modification"]/x["total_peptides"]*100,2), axis=1)
        try:
            df_to_save["intensity_pecentage_no_modifications[%]"] = df_to_save.apply(lambda x: handle_devide_by_zero(x["no_modification_intensity"],x["total_intensity"]), axis=1)
        except ZeroDivisionError:
            df_to_save["intensity_pecentage_no_modifications[%]"] = "NA"
        ordered_columns = ["AA number", "location_in_FASTA", "total_unique_peptides", "total_peptides", "peptides_with_no_modification"]
        for mod in mod_dict:
            ordered_columns.append("peptides_with_" + mod)
        
        ordered_columns.append("percentage_of_no_modifications[%]")
        for mod in mod_dict:
            ordered_columns.append("percentage_of_" + mod + "modification[%]")
        
        ordered_columns.append("total_intensity")
        ordered_columns.append("no_modification_intensity")
        for mod in mod_dict:
            ordered_columns.append("intensity_" + mod)
        
        ordered_columns.append("intensity_pecentage_no_modifications[%]")
        for mod in mod_dict:
            ordered_columns.append("intensity_pecentage_" + mod + "[%]")

        df_to_save = df_to_save.reindex(columns=ordered_columns)
        df_to_save.sort_values(by=["AA number"],inplace=True)

        if self.file_output_path == "":
            file_name, file_extension = os.path.splitext(self.file_path_input)
            df_to_save.to_csv(file_name + sample_name + "_" +self.uniprot_id + "_output.csv", index = False)
            return os.path.dirname(file_name + sample_name + "_" +self.uniprot_id + "_output.csv")
        else:
            file_name = os.path.basename(self.file_path_input)
            df_to_save.to_csv(self.file_output_path +"/" + os.path.splitext(file_name)[0] + sample_name + "_" +self.uniprot_id + "_output.csv", index = False)
            return self.file_output_path
        




        


def main():
    tk_ui = UI()
    # tk_ui.select_file()
    #print(tk_ui.file_path)
    tk_ui.root.mainloop()
    #UI.root.mainloop()

if __name__ == "__main__":
    main()