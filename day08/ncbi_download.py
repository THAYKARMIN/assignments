
import os
from datetime import datetime
import tkinter as tk
from tkinter import ttk, messagebox
from Bio import Entrez
import re

Entrez.email = "your.email@example.com"


def get_input():
    term = term_entry.get()
    number = number_entry.get()
    db = db_selector.get()
    file_format = format_selector.get()

    if not term or not number or not db or not file_format:
        messagebox.showerror("Input Error", "All fields are required.")
        return None, None, None, None
    
    return term, number, db, file_format


def search_files(term, number, db):
    handle = Entrez.esearch(db=db, term=term, idtype="acc", retmax=number)
    record = Entrez.read(handle)
    
    handle.close()
    
    if len(record["IdList"]) == 0:
        messagebox.showerror("Not Found", "No records found.")
        return None
    
    return record


def download_files(term, number, db, file_format):
    record = search_files(term, number, db)
    
    if record is None:
        return None, None

    data = []
    entry_names = []

    for id in record["IdList"]:
        handle = Entrez.efetch(db=db, id=id, rettype=file_format.lower(), retmode="text")
        item = handle.read()
        handle.close()
        data.append(item)

        if file_format.lower() == 'fasta':
            entry_name = extract_name_from_fasta(item)
        elif file_format.lower() == 'gb':
            entry_name = extract_name_from_genbank(item)
        elif file_format.lower() == 'txt':
            entry_name = extract_name_from_txt(item)
        else:
            entry_name = "unknown"

        entry_names.append(entry_name)


    return data, entry_names


def extract_name_from_fasta(fasta_data):
    lines = fasta_data.split('\n')
    header = lines[0]
    if header.startswith('>'):
        parts = header.split(None, 1)  
        name = re.split(r'\[|\(', parts[1])[0].strip()
        return name
    return "unknown"



def extract_name_from_genbank(genbank_data):
    lines = genbank_data.split('\n')
    for line in lines:
        if line.startswith("DEFINITION"):
            parts = line.split(None, 1)  
            name = re.split(r'\[|\(', parts[1])[0].strip()

            return name
    return "unknown"

def extract_name_from_txt(txt_data):

    lines = txt_data.split('\n')
    
    for line in lines:
        if line.startswith("1."):
            name = line.split("1.")[1].strip()
            return name
    return "unknown"


def save_files(term, data, entry_names, file_format):
    output_dir = f"output_{term}"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    filenames = []
    for entry_name, entry in zip(entry_names, data): 
        if file_format == "gb":
            filename = f"{term}_{entry_name}.gb"
        elif file_format == "fasta":
            filename = f"{term}_{entry_name}.fasta"
        else:
            filename = f"{term}_{entry_name}.txt"

        file_path = os.path.join(output_dir, filename)
        with open(file_path, "w") as fh:
            fh.write(entry)
            filenames.append(filename)
    return filenames


def csv_save(term, number, db):
    record = search_files(term, number, db)
    if record is None:
        return

    current_time = datetime.now()
    output_dir = f"output_{term}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
    filename = f"{term}.csv"
    file_path = os.path.join(output_dir, filename)
    with open(file_path, "a") as fh:
        fh.write("date,term,max,total\n")
        fh.write(f"{formatted_time},{term},{number},{record['Count']}\n")

    return filename

def main():
    term, number, db, file_format = get_input()
    if not term:
        return

    data, entry_names = download_files(term, number, db, file_format)
    if data is None:
        return

    filenames = save_files(term, data, entry_names, file_format)
    csv_filename = csv_save(term, number, db)

    messagebox.showinfo("Download Complete", f"Files saved:\n{', '.join(filenames)}\nCSV saved as: {csv_filename}")
    
    root.destroy()

def update_file_format_options(event):
    selected_db = db_selector.get()
    if selected_db in ["Taxonomy", "Gene"]:
        format_selector['values'] = ["txt"]
        format_selector.set("txt")
    else:
        format_selector['values'] = ["gb", "fasta"]
        format_selector.set("gb")

# GUI setup
root = tk.Tk()
root.title("NCBI Data Downloader")

# Search term
ttk.Label(root, text="Search Term:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
term_entry = ttk.Entry(root, width=50)
term_entry.grid(row=0, column=1, padx=5, pady=5)

# Number of records
ttk.Label(root, text="Number of Records:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
number_entry = ttk.Entry(root, width=50)
number_entry.grid(row=1, column=1, padx=5, pady=5)

# Database selector
ttk.Label(root, text="Database:").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
db_selector = ttk.Combobox(root, values=["Taxonomy", "Protein", "Gene", "Nucleotide"], width=47)
db_selector.grid(row=2, column=1, padx=5, pady=5)
db_selector.current(0)  # Set default selection
db_selector.bind("<<ComboboxSelected>>", update_file_format_options)  # Bind event handler

# File format selector
ttk.Label(root, text="File Format:").grid(row=3, column=0, padx=5, pady=5, sticky=tk.W)
format_selector = ttk.Combobox(root, values=["txt"], width=47)
format_selector.grid(row=3, column=1, padx=5, pady=5)
format_selector.current(0)  # Set default selection

# Download button
download_button = ttk.Button(root, text="Download Data", command=main)
download_button.grid(row=4, column=0, columnspan=2, pady=10)

# Initialize file format options based on default database selection
update_file_format_options(None)

root.mainloop()


