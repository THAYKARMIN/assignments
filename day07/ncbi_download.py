import sys
from Bio import Entrez
from datetime import datetime
import ncbi_functions

Entrez.email = "Thay.karmin@weizmann.ac.il"


def main():

    ncbi_functions.save_files()
    ncbi_functions.csv_save()
    ncbi_functions.printing()


main()
