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
        for k, v in self.dict.items():
            if len(v) > max_len:
                max_len = len(v)
                max_len_id = k
            if len(v) < min_len:
                min_len = len(v)
                min_len_id = k
                
        print(f"Max length: {max_len}")
        print(f"Min length: {min_len}")
        print(f"Identifier with max_len: {max_len_id}")
        print(f"Identifier with min_len: {min_len_id}")
        
    def orf_finder(self, frame):
        """Finds all ORFs present in each sequence of the FASTA file

        Args:
            frame (int): Reading frame used. Must be 1, 2 or 3.
        """
                            
        