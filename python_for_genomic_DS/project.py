""" 
@author: raulminan

This script contains the code necessary to complete the final project  for
the "Python for Genomic Data Science" course in Coursera, which is part of 
the Genomic Data Science Specialization from the Johns Hopkins University.
"""

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
        
        for i in range(0, len(sequence), 3):
            if sequence[i:i+3] == start:
                for j in range(i, len(sequence), 3):
                    codon = sequence[j:j+3]
                    temp_max_orf += codon #add current codon to the ORF sequence
                    orf_len = len(sequence[i:j+3]) 
                    if codon in stop:
                        if orf_len > max_orf_len:
                            max_orf_len = orf_len
                            max_orf_start = i
                            max_orf = temp_max_orf
                            temp_max_orf = ""
                        break
        
        if max_orf == "":
            print("No ORF found")                 
        else:
            print(f"longest orf is {max_orf}, it has a length of {max_orf_len} and starts at position {max_orf_start}")
            return max_orf, max_orf_len, max_orf_start
        
        

                    
      
        

    