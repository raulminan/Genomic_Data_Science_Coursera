""" 
@author: raulminan

This script contains the code necessary to complete the final project  for
the "Python for Genomic Data Science" course in Coursera, which is part of 
the Genomic Data Science Specialization from the Johns Hopkins University.
"""

from json import tool
import regex as re

DATA = "C:/Users/raulm/Documents/repos/genomic_data_science/python_for_genomic_DS/data/dna2.fasta"
    
class DNAToolKit():
    """This class contains the functions necessary to complete the course
    
    Methods
    -------
    record_count():
        Counts the number of records
    seq_length():
        Checks the length of each record
    orf_finder():
        Finds all ORFs in a sequence
    repeat_finder():
        Finds all repeats in a sequence
    """
    def __init__(self, file) -> dict:
        """Opens a FASTA file and saves it in a dictionary, where the keys are
        the sequence identifiers and the values are the DNA sequencies.

        Args:
            file (str): filepath to the fasta file with the sequences
        """
        self.file = file
        self.dict = {} # init empty dict to store the data
        
        with open(self.file) as f:
            for line in f:
                line = line.rstrip()
                if line[0] == ">":
                    identifier = line
                    self.dict[identifier] = ""
                else:
                    self.dict[identifier] += line
    
    def record_count(self):
        """Counts the number of records in a FASTA file
        """
        n = len(self.dict)
        print(f"There are {n} records in the file")
        
    def seq_length(self):
        """Prints the sequences with the max and min lengths, and their respective
        lengths
        """
        max_len = 0
        min_len = 1e100
        for identifier, sequence in self.dict.items():
            if len(sequence) > max_len:
                max_len = len(sequence)
                max_len_id = identifier
            if len(sequence) < min_len:
                min_len = len(sequence)
                min_len_id = identifier
                
        print(f"Max length: {max_len}")
        print(f"Min length: {min_len}")
        print(f"Identifier with max_len: {max_len_id}")
        print(f"Identifier with min_len: {min_len_id}")
        
    def orf_finder(self, sequence, frame):
        """Finds all open reading frames (ORFs) present in each sequence of the 
        FASTA file. 

        Args:
            sequence (str): sequence to find the ORFs in
            frame (int): Reading frame used to find ORFs. Must be 1, 2 or 3.

        Returns:
            Longest ORF, its length and its position
        """
        if frame not in [1, 2, 3]:
            raise ValueError("Frame must be 1, 2 or 3")
        
        start = "ATG"
        stop = ["TAA", "TAG", "TGA"]
        
        max_orf_len = 0
        max_orf_start = 0
        max_orf = "" 
        temp_max_orf = "" # to store the temporary ORFs
        sequence = sequence[frame - 1:]
        
        for i in range(0, len(sequence)+3, 3):
            if sequence[i:i+3] == start:
                for j in range(i, len(sequence)+3, 3):
                    codon = sequence[j:j+3]
                    temp_max_orf += codon #add current codon to the ORF sequence
                    orf_len = len(sequence[i:j+3]) 
                    if codon in stop:
                        if orf_len > max_orf_len:
                            max_orf_len = orf_len
                            max_orf_start = i + frame # add the bases removed in the beginning
                                                      # (sequence = sequence[frame - 1:])
                            max_orf = temp_max_orf
                            temp_max_orf = ""
                        break
        
        if max_orf == "":
            #print("No ORF found")
            pass                 
        else:
            #print(f"longest orf is {max_orf}, it has a length of {max_orf_len} and starts at position {max_orf_start}")
            return max_orf, max_orf_len, max_orf_start
                   
    def repeats(self, sequence, n):
        """Given a length n, this function identifies all repeats of length n in
        all sequences of a FASTA f

        Args:
            n (int): length of the repeats to find
            sequence(str): sequnce to find the repeats in 
        """
        repeats = {}
        
        for i in range(0, len(sequence)):
            repeat = sequence[i:i+n]
            matches = re.findall(repeat, sequence, overlapped=True)
            if len(repeat) < n:
                break
            else:
                repeats[repeat] = len(matches)
            
        return repeats    
        
        
if __name__ == "__main__":
    toolkit = DNAToolKit(DATA)
    
    # Q1
    # How many records are in the multi FASTA file?
    toolkit.record_count()
    print("==================================================================\n")
    
    # Q2
    # What is the length of the longest sequence in the file?
    # Q3
    # What is the length of the shortest sequence in the file?
    toolkit.seq_length()
    print("==================================================================\n")
    # Q4
    # What is the longest ORF appearing in reading frame 2 of any of the sequences
    sequences = list(toolkit.dict.values())
    max_len = 0
    
    for sequence in sequences:
        orfs = toolkit.orf_finder(sequence=sequence, frame=2)
        try:
            if orfs[1] > max_len:
                max_len = orfs[1]
        except TypeError: # raises when orf_finder() can't find an ORF
            continue
        
    print(f"The longest ORF in any sequence is {max_len}")
    print("==================================================================\n")
    # Q5
    # What is the starting position of the longest ORF in reading frame 3 in any 
    # of the sequences?
    
    sequences = list(toolkit.dict.values())
    max_len = 0
    start = 0
    
    for sequence in sequences:
        orfs = toolkit.orf_finder(sequence=sequence, frame=3)
        try:
            if orfs[1] > max_len:
                max_len = orfs[1]
                start = orfs[2]
        except TypeError: # raises when orf_finder() can't find an ORF
            continue
    
    print(f"The starting position of the longest ORF in reading frame 3 is {start}")
    print("==================================================================\n") 
    
    # Q6
    # What is the length of the longest ORF appearing in any sequence and in any 
    # reading frame?
    
    sequences = sequences = list(toolkit.dict.values())
    frames = [1,2,3]
    max_len = 0
    
    for frame in frames:
        for sequence in sequences:
            orfs = toolkit.orf_finder(sequence=sequence, frame=frame)
            try:
                if orfs[1] > max_len:
                    max_len = orfs[1]
            except TypeError: # raises when orf_finder() can't find an ORF
                continue
            
    print(f"The longest ORF in any sequence, in any reading frame is {max_len}")
    print("==================================================================\n") 
    
    # Q7
    # What is the length of the longest forward ORF that appears in the seq with 
    # identifier "gi|142022655|gb|EQ086233.1|16"

    
 
    # find the key with that identifier
    id = "gi|142022655|gb|EQ086233.1|16"
    sequences = [seq for key, seq in toolkit.dict.items() if id in key]
    
    # find longest ORF
    frames = [1,2,3]
    max_len = 0
    
    for frame in frames:
        for sequence in sequences:
            orfs = toolkit.orf_finder(sequence=sequence, frame=frame)
            try:
                if orfs[1] > max_len:
                    max_len = orfs[1]
            except TypeError:
                continue
        
    print(f"The longest ORF in that sequence, in any reading frame is {max_len}")
    print("==================================================================\n") 
    
    # Q8
    # Find the most frequently occuring repeat of length 6 in all sequences. How
    # many times does it occur in all?

    max_repeat_global = 0
    most_frequent_repeat = ""

    # find the repeat that appears the most times in a sequence
    for id, seq in toolkit.dict.items():
        repeats = toolkit.repeats(seq, 6) # dict with all repeats of length 6 in sequence seq
        max_repeat = max(repeats.values()) # get number of times the most frequent repeat appears in that seq
        max_repeat_seq = max(repeats, key = repeats.get)
        if max_repeat > max_repeat_global:
            max_repeat_global = max_repeat
            most_frequent_repeat = max_repeat_seq
    
    # find how many times that repeat appears in all other sequences
    appearences = []
    for id, seq in toolkit.dict.items():
        appearences.append(re.findall(most_frequent_repeat, seq, overlapped=True))
    
    flat = [item for sublist in appearences for item in sublist]
    
    print(f"The most frequent repeat of length 6 is {most_frequent_repeat} repeats itself {len(flat)} times")
    print("==================================================================\n") 
    
    # Q9
    # Find all repeats of length 12 in the input file. Let's use Max to specify
    # the number of copies of the most frequent repeat of length 12. How many
    # different 12-base sequences occur Max times?
    
    Max = 0
    all_repeats = {}
    
    # find all repeats of length 12:
    for id, seq in toolkit.dict.items():
        repeats = toolkit.repeats(seq, 12)
        for repeat, n in repeats.items():
            try:
                all_repeats[repeat] += n
            except KeyError:
                all_repeats[repeat] = n # create key if it's the 1st time the repeat appears
    Max = max(all_repeats.values())
    
    # how many different repeats occur Max times?
    counter = 0
    repeats = []
    for repeat, n in all_repeats.items():
        if n == Max:
            counter += 1
            repeats.append(repeat)
        
    print(
        f"There are {counter} different repeats of length 12 that repeat themselves\
        {Max} times and they are: {repeats}"
    )
    print("==================================================================\n")         
    
    # Q10
    # Which one of the following repeats of length 7 has a maximum number of occurences
    
    CANDIDATES = {
        "AATGGCA": 0,
        "TGCGCGC": 0,
        "CATCGCC": 0,
        "CGCGCCG": 0
    }
    
    for id, seq in toolkit.dict.items():
        repeats = toolkit.repeats(seq, 7)
        for repeat, n in repeats.items():
            if repeat in CANDIDATES:
                CANDIDATES[repeat] += n
                
    print(f"The repeat that appears the most times is {max(CANDIDATES, key=CANDIDATES.get)}")
    print("==================================================================\n") 