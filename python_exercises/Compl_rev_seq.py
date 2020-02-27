# -*- coding: utf-8 -*-
"""
Functions to get the complementary, reverse and complementary reverse of a given sequence
"""

def compl(seq):
    """Creates the complementary of a given DNA sequence
      str -> str"""
    base_compl = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
    compl_seq = ''
    for base in seq:
        compl_seq += base_compl[base]
    return compl_seq

def reverse(seq):
    """Creates the reverse of a given DNA sequence
      str -> str"""
    return seq[::-1]

reverse('ATTCCAAGC')


def rev_compl(seq):
    """Creates the reverse from the complementary of a given sequence
    str -> str"""
    return compl(reverse(seq))
