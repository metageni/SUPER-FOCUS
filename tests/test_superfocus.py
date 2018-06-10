# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import random

from superfocus_app.superfocus import is_wanted_file

import pytest


def test_is_wanted_file():
    assert is_wanted_file(["a.fasta", "b.fastq", "x.FASTq", "y.FASTA", "n.fna"]) == ['a.fasta', 'b.fastq', 'n.fna', 'x.FASTq', 'y.FASTA']
    assert is_wanted_file(["f.png", "a.fasta", "b.fastq", "x.FASTq", "y.FASTA", "n.fna"]) == ['a.fasta', 'b.fastq', 'n.fna', 'x.FASTq', 'y.FASTA']
    assert is_wanted_file(["f.png" , "queries/"]) == []
