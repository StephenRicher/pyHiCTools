#!/usr/bin/env python3

import pytest, os

from hictools.hictools_truncate import *

arguments = (
    [('^GATC',  ('GATCGATC', 'GATC')), 
    ( 'A^GATCT',   ('AGATCGATCT', 'AGATCT')), 
    pytest.param('GATC', None, 
        marks = pytest.mark.xfail(raises = AssertionError)),
    pytest.param(None, None, 
        marks = pytest.mark.xfail(raises = AssertionError)),
    pytest.param(1, None, 
        marks = pytest.mark.xfail(raises = AssertionError))])
    
@pytest.mark.parametrize('test_input,   expected', arguments)
def test_process_restriction(test_input, expected):
    assert(process_restriction(test_input) == expected)

# Test entire pipeline from input to output
@pytest.mark.parametrize(
     'test_in,     exp_out,        exp_err,        re', [
    ('test.txt',  'test_out.txt', 'test_err.txt', '^GATC')])
def test_truncate(capsys, test_in, exp_out, exp_err, re):
    
    # Tool to test
    truncate(infile = test_in, output = None, 
        read_gzip = False, write_gzip = False, 
        sample = None, restriction = re)
    
    # Capture stdout and stderr
    out, err = capsys.readouterr()
    
    # Compare stdour and stderr against expected
    assert out == open(exp_out).read()
    assert err == open(exp_err).read()
