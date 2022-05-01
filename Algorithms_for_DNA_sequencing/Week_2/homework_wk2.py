#!/usr/bin/env python

from bm_preproc import BoyerMoore
from kmer_index import Index, SubseqIndex


__author__ = "raulminan"

""" 
This script contains the code necessary to complete the homework for week 1
of "Algorithms for DNA sequencing" course in Coursera, which is part of 
the Genomic Data Science Specialization from the Johns Hopkins University.

Instructions:
implement versions of the naive exact matching and Boyer-Moore algorithms 
that additionally count and return (a) the number of character comparisons 
performed and (b) the number of alignments tried. Roughly speaking, these 
measure how much work the two different algorithms are doing.
"""

def read_genome(filename):
        genome = ''
        with open(filename, 'r') as f:
            for line in f:
                # ignore header line with genome information
                if not line[0] == '>':
                    genome += line.rstrip()
        return genome

def boyer_moore_with_counts(p, p_bm, t):
    """Do Boyer-Moore matching

    Args:
        p (str): pattern to match
        p_bm (bm_preproc.BoyerMoore): BoyerMoore object
        t (str): text to match

    Returns:
        list: list of occurences
        int: number of alingments
        int: number of character comparisons
    """
    i = 0
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        num_alignments += 1
        for j in range(len(p)-1, -1, -1):
            num_character_comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, num_alignments, num_character_comparisons


def naive_with_counts(p, t):
    """Do naive exact matching matching

    Args:
        p (str): pattern to match
        t (str): text to match

    Returns:
        list: list of occurences
        int: number of alingments
        int: number of character comparisons
    """
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        num_alignments += 1
        for j in range(len(p)):  # loop over characters
            num_character_comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, num_alignments, num_character_comparisons

def naive_2mm(p, t):
    """Do the naive exact match algorithms allowing up to 2 mismatches

    Args:
        p (str): pattern to match
        t (str): text to match the pattern to

    Returns:
        list: list with the indeces of occurences of the pattern p in text t
    """
    occurrences = set() 
    
    for i in range(len(t) - len(p) + 1):
        mismatches = 0 # loop over alignments
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.add(i)  # all chars matched; record
    return list(occurrences)

def boyer_approximate__matching(p, t, n):
    """Do approximate marching using the Boyer-Moore algorith and the pigeon
    principle

    Args:
        p (str): pattern to match
        t (str): text to match the pattern to
        n (int): maximum amount of mismatches allowed
    
    """
    segment_length = round(len(p) // (n+1))
    all_matches = set()
 
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = BoyerMoore(p[start:end], alphabet="ACGT")
        matches, _, _ = boyer_moore_with_counts(p[start:end], p_bm, t)
        
        for m in matches:
            if m - start < 0 or m - start + len(p) > len(t):
                continue
            
            mismatches = 0
            # verify to the left
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break  
                    
            # verify to the right
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
        
            if mismatches <= n:
                all_matches.add(m - start)
    
    return list(all_matches), len(matches)
    
def index_approximate_matching(p, t, n, k):
    """Do approximate matching using the pigeonhole principle

    Args:
        p (str): pattern to match
        t (str): text to match the pattern to
        n (int): maximum amount of mismatches allowed
        k (int): k-mer to use for indexing
    """
    
    segment_length = round(len(p) // (n+1))
    all_matches = set()
    idx = Index(t, k)
    n_hits = 0
 
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        hits = idx.query(p[start:end])
        n_hits += len(hits)
        
        for hit in hits:
            if hit - start < 0 or hit - start + len(p) > len(t):
                continue
            
            mismatches = 0
            # verify to the left
            for j in range(0, start):
                if not p[j] == t[hit-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break  
                    
            # verify to the right
            for j in range(end, len(p)):
                if not p[j] == t[hit-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
        
            if mismatches <= n:
                all_matches.add(hit - start)
    
    return list(all_matches), n_hits

def subseq_approximate_matching(p, t, n, k, ival):
    """Do approximate matching using subsequences

    Args:
        p (str): pattern to match
        t (str): text to match the pattern to
        n (int): maximum amount of mismatches allowed
        k (int): k-mer to use for indexing
        ival (int): interval to use for create sub sequences
    """
    all_matches = set()
    idx = SubseqIndex(t, k, ival)
    n_hits = 0
    
    for i in range(n+1):
        start = i 
        hits = idx.query(p[start:])
        n_hits += len(hits)
        for hit in hits:
            if hit - start < 0 or hit - start + len(p) > len(t):
                continue
            
            mismatches = 0
            t
            for j in range(0, len(p)):
                if not p[j] == t[hit-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break  
        
            if mismatches <= n:
                all_matches.add(hit - start)
    
    return list(all_matches), n_hits
    
if __name__ == "__main__":
    
    FILENAME = R"Algorithms_for_DNA_sequencing\Week_2\chr1.GRCh38.excerpt.fasta"
    TEXT = read_genome(FILENAME)

    # Q1, Q2
    # How many alignments and character comparisons does the naive exact matching 
    # algorithm try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG
    # (derived from human Alu sequences) to the excerpt of human chromosome 1?  
    # (Don't consider reverse complements.)
    
    PATTERN_1 = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
    occurrences, num_aligments, num_character_comparisons =\
        naive_with_counts(PATTERN_1, TEXT)
        
    print(f"The naive exact matching algorithm tries {num_aligments} alignments and "
        f"{num_character_comparisons} character comparisons")
    
    # Q3
    # How many alignments does Boyer-Moore try when matching the string 
    # GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu 
    # sequences) to the excerpt of human chromosome 1?  
    # (Don't consider reverse complements.)
    p_bm = BoyerMoore(PATTERN_1)
    occurrences_bm, num_aligments_bm, num_character_comparisons_bm =\
        boyer_moore_with_counts(PATTERN_1, p_bm, TEXT)
    
    print(f"The naive exact matching algorithm tries {num_aligments_bm } alignments")
    
    # Q4
    # How many times does the string "GGCGCGGTGGCTCACGCCTGTAAT",  which is 
    # derived from a human Alu sequence, occur with up to 2 substitutions in the 
    # excerpt of human chromosome 1?  (Don't consider reverse complements here.)

    # Q5
    # How many total index hits are there when searching for occurences with up
    # to 2 substitutions in the excerpt of human chromosome 1?
    
    
    PATTERN_2 = "GGCGCGGTGGCTCACGCCTGTAAT"
    n = 2
    k = 8
    
    matches, hits = index_approximate_matching(PATTERN_2, TEXT, 2, 8)
    print(f"The {len(matches)} matches found where {matches} "
          f"for a total of {hits} hits")
    
    # Q6
    # Write a function that, given a length-24 pattern P and given a SubseqIndex
    # object built with k = 8 and ival = 3, finds all approximate occurrences of 
    # P within T with up to 2 mismatches.
    # When using this function, how many total index hits are there when searching 
    # for GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of 
    # human chromosome 1?  (Again, don't consider reverse complements.)
    
    n = 2
    k = 8
    ival = 3
    
    occurrences, n_hits = subseq_approximate_matching(PATTERN_2, TEXT, n, k, ival)
    
    print(f"The number of total hits is {n_hits}")

    
    
    