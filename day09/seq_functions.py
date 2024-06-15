import re

def fasta_input(file_path):
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
          
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def longest_repetition(sequence):

    length = 1
    result = ''
    while True:
        regex = r'([GATC]{' + str(length) + r'}).*\1'
        m = re.search(regex, sequence)
        if m:
            result = m.group(1)
            length += 1
        else:
            break
    len_seq = len(result)
    return result, len_seq

def potential_start_codon(sequence):

    pattern = r'ATG'
    matches = re.finditer(pattern, sequence)
    indexes = [match.start() for match in matches]
    return indexes