import unittest
import homework_wk2
import bm_preproc


class TestHomeworkWk2(unittest.TestCase):
    """Tests for the functions used for Week 2's assigment"""
    def test_naive_with_counts(self):
        p = 'word'
        t = 'there would have been a time for such a word'
        occurrences, num_alignments, num_character_comparisons =\
            homework_wk2.naive_with_counts(p, t)
            
        self.assertEqual(occurrences, [40])
        self.assertEqual(num_alignments, 41)
        self.assertEqual(num_character_comparisons, 46)
        
    def test_naive_with_counts_2(self):
        p = 'needle'
        t = 'needle need noodle needle'
        occurrences, num_alignments, num_character_comparisons =\
            homework_wk2.naive_with_counts(p, t)
            
        self.assertEqual(occurrences, [0, 19])
        self.assertEqual(num_alignments, 20)
        self.assertEqual(num_character_comparisons, 35)
        
    def test_boyer_moore_with_counts(self):
        p = 'word'
        t = 'there would have been a time for such a word'
        lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
        p_bm = bm_preproc.BoyerMoore(p, lowercase_alphabet)
        occurrences, num_alignments, num_character_comparisons =\
            homework_wk2.boyer_moore_with_counts(p, p_bm, t)
            
        self.assertEqual(occurrences, [40])
        self.assertEqual(num_alignments, 12)
        self.assertEqual(num_character_comparisons, 15)

    def test_boyer_moore_with_counts_2(self):
        p = 'needle'
        t = 'needle need noodle needle'
        lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
        p_bm = bm_preproc.BoyerMoore(p, lowercase_alphabet)
        occurrences, num_alignments, num_character_comparisons =\
            homework_wk2.boyer_moore_with_counts(p, p_bm, t)
            
        self.assertEqual(occurrences, [0, 19])
        self.assertEqual(num_alignments, 5)
        self.assertEqual(num_character_comparisons, 18)
    
    def test_index_approximate_matching(self):
        FILENAME = R"Algorithms_for_DNA_sequencing\Week_2\chr1.GRCh38.excerpt.fasta"
        p = "GGCGCGGTGGCTCACGCCTGTAAT"
        t = homework_wk2.read_genome(FILENAME)
        n = 2
        k = 8
        
        occurences_naive = homework_wk2.naive_2mm(p, t)
        occurences_index, total_matches_index =\
            homework_wk2.index_approximate_matching(p, t, n, k)
        occurences_boyer, total_matches_boyer =\
            homework_wk2.boyer_approximate__matching(p, t, n)
             
        self.assertEqual(occurences_naive, occurences_index)
        self.assertEqual(occurences_naive, occurences_boyer)
        self.assertEqual(total_matches_boyer, total_matches_boyer)

    def test_subseq_approximate_matching(self):
        t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
        p = 'to-morrow and to-morrow '
        n = 2 
        k = 8
        ival = 3 
        occurrences, num_index_hits  = \
            homework_wk2.subseq_approximate_matching(p, t, n, k, ival)
            
        self.assertEqual(occurrences, [0, 14])
        self.assertEqual(num_index_hits, 6)
          

if __name__ == "__main__":
    unittest.main()