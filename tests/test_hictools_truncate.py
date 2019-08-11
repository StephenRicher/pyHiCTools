#!/usr/bin/env python3

import pytest 
from hictools.hictools_truncate import *

@pytest.mark.parametrize(
    'test_input, expected', [
    ('^GATC', ('GATCGATC', 'GATC')),
    ('A^GATCT', ('AGATCGATCT', 'AGATCT')),
    pytest.param('GATC', None, marks=pytest.mark.xfail(raises = AssertionError)),
    pytest.param(None, None, marks=pytest.mark.xfail(raises = AssertionError)),
    pytest.param(1, None, marks=pytest.mark.xfail(raises = AssertionError))])
def test_process_restriction(test_input, expected):
    assert(process_restriction(test_input) == expected)

