# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from superfocus_app.superfocus import is_wanted_file, get_denominators, aggregate_level, is_valid_number

import pytest


def test_is_valid_number():
    assert is_valid_number("1") == True
    assert is_valid_number("1.0") == True
    assert is_valid_number("0.00001") == True
    assert is_valid_number("-0.00001") == False
    assert is_valid_number("-1") == False
    assert is_valid_number("-1.0") == False
    assert is_valid_number("aa") == False
    assert is_valid_number("1a") == False
    assert is_valid_number("-1a") == False


def test_is_wanted_file():
    assert is_wanted_file(["a.fasta", "b.fastq", "x.FASTq", "y.FASTA", "n.fna"]) == ['a.fasta', 'b.fastq', 'n.fna',
                                                                                     'x.FASTq', 'y.FASTA']
    assert is_wanted_file(["f.png", "a.fasta", "b.fastq", "x.FASTq", "y.FASTA", "n.fna"]) == ['a.fasta', 'b.fastq',
                                                                                              'n.fna', 'x.FASTq',
                                                                                              'y.FASTA']
    assert is_wanted_file(["f.png", "queries/"]) == []


def test_get_denominators():
    test_1 = {
        "a": [1, 2, 3],
        "b": [4, 8, 10]
    }
    assert list(get_denominators(test_1)) == [5, 10, 13]

    test_2 = {
        "a": [10, 5, 30],
        "b": [2, 7, 11],
        "c": [4, 0, 1]
    }
    assert list(get_denominators(test_2)) == [16, 12, 42]

    test_3 = {
        "a": [0, 0, 0],
        "b": [0, 0, 0]
    }
    assert list(get_denominators(test_3)) == [0, 0, 0]


def test_aggregate_level():
    test_1 = {"a1\ta2\ta3\tas": [1, 2, 3], "b1\tb2\tb3\tbs": [4, 8, 10], "c1\tc2\tc3\tcs": [1, 2, 3]}
    normalizer_1 = [6, 12, 16]
    position = 0

    assert aggregate_level(test_1, position, normalizer_1) == {
        'a1': [1, 2, 3, 16.666666666666664, 16.666666666666664, 18.75],
        'b1': [4, 8, 10, 66.666666666666657, 66.666666666666657, 62.5],
        'c1': [1, 2, 3, 16.666666666666664, 16.666666666666664, 18.75]
    }
