#!/usr/bin/env python3

""" Extract HiC read pair information encodiing by hictools process.
    Useful for assessing contamination and determining parameters for
    hictools filter.
"""

import sys
import fileinput
import pyHiCTools as hic
import pyCommonTools as pct
from contextlib import ExitStack

def extract(infile, out, sample):

    log = pct.create_logger()

    if not sample:
        if infile == []:
            sample = 'stdin'
        else:
            sample = infile

    with ExitStack() as stack:
        input = stack.enter_context(
            fileinput.input(files=infile))
        output = stack.enter_context(
            pct.open_smart(out, mode='w'))
        log.info(f'Writing output to {out}.')

        output.write(
            'sample\torientation\tinteraction_type\tditag_length\t'
            'insert_size\tfragment_seperation\n')

        for i, line in enumerate(input):
            if line.startswith('@'):
                continue
            else:
                try:
                    read1 = pct.Sam(line)
                    read2 = pct.Sam(next(input))
                except StopIteration:
                    log.exception('Odd number of alignments in file')
                    sys.exit(1)
                output.write(
                    f'{sample}\t{read1.optional["or:Z"]}\t'
                    f'{read1.optional["it:Z"]}\t{read1.optional["dt:i"]}\t'
                    f'{read1.optional["is:i"]}\t{read1.optional["fs:i"]}\n')
