import unittest

def reverse_complement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def naive_2mm(p, t):
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

def naive_with_rc(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        for seq in (p, reverse_complement(p)):
            match = True
            for j in range(len(seq)):  # loop over characters
                if t[i+j] != seq[j] :  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)
                break    

    return occurrences

class TestHomeWorkWk1(unittest.TestCase):
    """Tests for the functions used for Week 1's assignment"""
    def test_naive_2mm(self):
        p = 'CTGT'
        ten_as = 'AAAAAAAAAA'
        t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
        occurences = naive_2mm(p,t)
        
        self.assertEqual(occurences, [10, 24, 38])

    def test_naive_with_rc(self):
        p = 'CCC'
        ten_as = 'AAAAAAAAAA'
        t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
        occurrences = naive_with_rc(p, t)
        
        self.assertEqual(occurrences, [10, 23])
    
if __name__ == "__main__":
    unittest.main()
    