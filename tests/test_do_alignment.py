# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from superfocus_app.do_alignment import normalise_counts, update_results, parse_alignments


def test_parse_alignment():
    pass


def test_normalise_counts():
    # example 1: empty input
    assert normalise_counts({}) == {}

    # example 2: one function
    assert normalise_counts({"function_1": 1}) == {"function_1": 1}

    # example 3: more than one function
    assert normalise_counts({"function_1": 1, "function_2": 1}) == {"function_1": 0.5, "function_2": 0.5}


def test_update_results():
    sample_index = 0
    normalise = 0
    number_samples = 1

    ## One function - not normalised

    # example 1
    results = {}
    data = {"function_1": 1}

    true_answer = {"function_1": [1]}
    assert update_results(results, sample_index, data, normalise, number_samples) == true_answer

    # example 2
    results = {"function_1": [1]}
    data = {"function_1": 1}

    true_answer = {"function_1": [2]}
    assert update_results(results, sample_index, data, normalise, number_samples) == true_answer

    # example 3
    results = {"function_1": [1]}
    data = {"function_2": 1}

    true_answer = {"function_1": [1], "function_2": [1]}
    assert update_results(results, sample_index, data, normalise, number_samples) == true_answer

    ## Two functions - not normalised

    # example 1
    results = {}
    data = {"function_1": 1, "function_2": 1}

    true_answer = {"function_1": [1], "function_2": [1]}
    assert update_results(results, sample_index, data, normalise, number_samples) == true_answer

    # example 2
    results = {"function_1": [1]}
    data = {"function_1": 1, "function_2": 1}

    true_answer = {"function_1": [2], "function_2": [1]}
    assert update_results(results, sample_index, data, normalise, number_samples) == true_answer

    # normalised example
    normalise = 1

    results = {}
    data = {"function_1": 1, "function_2": 1}

    true_answer = {"function_1": [0.5], "function_2": [0.5]}
    assert update_results(results, sample_index, data, normalise, number_samples) == true_answer
