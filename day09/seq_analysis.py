import argparse
import seq_functions

def main():
    parser = argparse.ArgumentParser(description='Find the longest repeated substring or potential start codon indexes in a DNA sequence from a FASTA file.')
    parser.add_argument('--fasta_path', type=str, required=True, help='Path to the FASTA file')
    parser.add_argument('--analysis', nargs='+', choices=['longest_repetition', 'start_codon'], required=True, help='Choose analysis options: longest_repetition, start_codon, or both')
    
    args = parser.parse_args()
    fasta_path = args.fasta_path
    
    sequence = seq_functions.fasta_input(fasta_path)

    if 'longest_repetition' in args.analysis:
        result, len_seq = seq_functions.longest_repetition(sequence)
        print(f"Longest repeated substring: {result}")
        print(f"length: {len_seq}")

    if 'start_codon' in args.analysis:
        codon_i = seq_functions.potential_start_codon(sequence)
        print(f"List of potential start codon indexes: {codon_i}")
    

if __name__ == '__main__':
    main()
