import matplotlib.pyplot as plt
import numpy as np

class hwWeek1():
    """Class to solve the homework for week 1 of Algorithms for DNA sequencing 
    in Coursera"""
    
    def __init__(self, filename):
        self.filename = filename
    
    def read_fastq(self):
        sequences = []
        qualities = []
        with open(self.filename) as fh:
            while True:
                fh.readline()  # skip name line
                seq = fh.readline().rstrip()  # read base sequence
                fh.readline()  # skip placeholder line
                qual = fh.readline().rstrip() # base quality line
                if len(seq) == 0:
                    break
                sequences.append(seq)
                qualities.append(qual)
        return sequences, qualities

    def read_genome(self):
        genome = ''
        with open(self.filename, 'r') as f:
            for line in f:
                # ignore header line with genome information
                if not line[0] == '>':
                    genome += line.rstrip()
        return genome
    
    def reverse_complement(self, s):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        t = ''
        for base in s:
            t = complement[base] + t
        return t
    
    
    def naive_with_rc(self, p, t):
        occurrences = []
        for i in range(len(t) - len(p) + 1):  # loop over alignments
            for seq in (p, hwWeek1.reverse_complement(self, p)):
                match = True
                for j in range(len(seq)):  # loop over characters
                    if t[i+j] != seq[j] :  # compare characters
                        match = False
                        break
                if match:
                    occurrences.append(i)
                    break    
    
        return occurrences
    
    def naive_2mm(self, p, t):
        occurrences = []
        for i in range(len(t) - len(p) + 1):
            mismatches = 0 # loop over alignments
            for j in range(len(p)):  # loop over characters
                if t[i+j] != p[j]:  # compare characters
                    mismatches += 1
                    if mismatches > 2:
                        break
            if mismatches <= 2:
                occurrences.append(i)  # all chars matched; record
        return occurrences

    def find_GC_by_pos(self, reads):
        gc = [0] * 100
        totals = [0] * 100
        
        for read in reads:
            for i in range(len(read)):
                if read[i] == "C" or read[i] == "G":
                    gc[i] += 1
                totals[i] += 1
        
        for i in range(len(gc)):
            if len(totals) != 0:
                gc[i] /= float(totals[i])
                
        return gc
        
    
if __name__ == "__main__":
    FILENAME = R"Algorithms_for_DNA_sequencing\lambda_virus.fa"
    
    wk1 = hwWeek1(FILENAME)
    GENOME = wk1.read_genome()
    
    # Q1
    # How many times does AGGT or its revers complement occur in the lambda virus
    # genome?
    PATTERN_1 = 'AGGT'
    occurences_1 = wk1.naive_with_rc(PATTERN_1, GENOME)
    
    print(len(occurences_1))
    
    # Q2
    # How many times does AGGT or its revers complement occur in the lambda virus
    # genome?

    PATTERN_2 = 'TTAA'
    occurences_2 = wk1.naive_with_rc(PATTERN_2, GENOME)
    
    print(len(occurences_2))
    
    # Q3
    # What is the leftmost occurrence of ACTAAGT or its reverse complement in the 
    # Lambda virus genome
    
    PATTERN_3 = 'ACTAAGT'
    occurences_3 = wk1.naive_with_rc(PATTERN_3, GENOME)
    print(occurences_3[0])

    
    # Q4
    # What is the leftmost occurrence of AGTCGA or its reverse complement in the 
    # Lambda virus genome
    
    PATTERN_4 = 'AGTCGA'
    occurences_4 = wk1.naive_with_rc(PATTERN_4, GENOME)
    print(occurences_4[0])
    
    
    # Q5
    # How many times does TTCAAGCC occur in the Lambda virus genome when 
    # allowing up to 2 mismatches? 
    
    PATTERN_5 = 'TTCAAGCC'
    occurences_5 = wk1.naive_2mm(PATTERN_5, GENOME)
    print(len(occurences_5))
    
    # Q6
    # What is the offset of the leftmost occurrence of AGGAGGTT
    # in the Lambda virus genome when allowing up to 2 mismatches?
    
    PATTERN_6 = 'AGGAGGTT'
    occurences_6 = wk1.naive_2mm(PATTERN_6, GENOME)
    print(occurences_6[0])
    
    # Q7
    # This dataset has something wrong with it; one of the sequencing cycles is 
    # poor quality. Report which sequencing cycle has the problem.  Remember that a 
    # sequencing cycle corresponds to a particular offset in all the reads. For 
    # example, if the leftmost read position seems to have a problem consistently 
    # across reads, report 0. If the fourth position from the left has the problem, 
    # report 3. Do whatever analysis you think is needed to identify the bad cycle. 
    # It might help to review the "Analyzing reads by position" video.
    
    FILENAME_2 = R"Algorithms_for_DNA_sequencing\ERR037900_1.first1000.fastq"
    q7 = hwWeek1(FILENAME_2)

    GENOME, _ = q7.read_fastq()
    gc = q7.find_GC_by_pos(GENOME)

    plt.plot(range(len(gc)), gc)
    plt.show()
    
    error = [i for i in range(len(gc)) if gc[i] < 0.3]
    print(error)