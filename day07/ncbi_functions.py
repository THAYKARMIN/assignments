
import sys
from Bio import Entrez
from datetime import datetime
import os


def get_input():

    if len(sys.argv) != 3:
        exit("Please enter in the command line: ncbi_download.py TERM NUMBER")

    term = sys.argv[1]
    number = sys.argv[2]
    return term, number


def search_files():

    term, number = get_input()

    handle = Entrez.esearch(db="Taxonomy", term=term, idtype="acc", retmax=number)
    record = Entrez.read(handle)
    if len(record["IdList"]) == 0:
        exit("Taxonomy not found")

    handle.close()
    return record


def download_files():

    record = search_files()
    data = []
    species_names = []

    for id in record["IdList"]:
        handle = Entrez.efetch(db="Taxonomy", id=id, rettype="gb", retmode="text")
        item = handle.read()
        handle.close()
        data.append(item)

        lines = item.split("\n")
        for line in lines:
            if line.startswith("1."):
                species_name = line.split("(", 1)[0].strip()[3:]
                species_names.append(species_name)

    return data, species_names


def save_files():

    term, _ = get_input()
    data, species_names = download_files()
    output_dir = f"output_{term}"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    filenames = []
    for species_name, entry in zip(species_names, data):
        filename = f"{term}_{species_name}"
        file_path = os.path.join(output_dir, filename)
        with open(file_path, "w") as fh:
            fh.write(entry)
            filenames.append(filename)

    return filenames


def csv_save():

    term, number = get_input()
    record = search_files()
    current_time = datetime.now()
    output_dir = f"output_{term}"
    formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
    filename = f"{term}.csv"
    file_path = os.path.join(output_dir, filename)
    with open(file_path, "a") as fh:
        fh.write("date,term,max,total\n")
        fh.write(f"{formatted_time},{term},{number},{record["Count"]}\n")

    return filename


def printing():

    filenames = save_files()
    print("Files saved:")
    for filename in filenames:
        print(filename)

    record = search_files()
    current_time = datetime.now()
    formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
    print("date,term,max,total")
    print(formatted_time, sys.argv[1], sys.argv[2], record["Count"])
