#!/usr/bin/env python3


''' Truncate proximity ligated restriction digest fragments at the
    ligation junction.
'''

import sys
import pyCommonTools as pct


def process_restriction(restriction):
    assert(isinstance(restriction, str))
    assert(restriction.count('^') == 1)
    restriction_seq = restriction.upper().replace('^', '')
    cut_site1 = restriction.index('^')
    ligation2 = restriction_seq[cut_site1:]
    cut_site2 = len(restriction) - cut_site1 - 1
    ligation1 = restriction_seq[0:cut_site2]
    ligation_seq = ligation1 + ligation2
    return ligation_seq, restriction_seq


def truncate(infile, output, read_gzip, write_gzip, sample, restriction):

    ''' Run main loop. '''

    log = pct.create_logger()

    ligation_seq, restriction_seq = process_restriction(restriction)
    total = 0
    truncated = 0
    truncated_length = 0

    if not sample:
        sample = infile

    with pct.open_gzip(output, 'wt', write_gzip) as out_obj, \
            pct.open_gzip(infile, 'rt', read_gzip) as in_obj:

        is_truncated = False
        for index, line in enumerate(in_obj):
            line = line.rstrip('\n')
            # Sequence line
            if index % 4 == 1:
                total += 1
                line = line.upper()
                if ligation_seq in line:
                    line = line[0: line.index(ligation_seq)] + restriction_seq
                    is_truncated = True
                seq_length = len(line)
            # Quality line
            elif index % 4 == 3:
                line = line[0:seq_length]
                if is_truncated:
                    truncated += 1
                    truncated_length += seq_length
                    is_truncated = False
            out_obj.write(f'{line}\n')
        try:
            mean_truncated_length = truncated_length/truncated
        except ZeroDivisionError:
            mean_truncated_length = 'na'
        sys.stderr.write(
            f'{sample}\tTotal\t{total}\n'
            f'{sample}\tTruncated\t{truncated}\n'
            f'{sample}\tNot truncated\t{total-truncated}\n'
            f'{sample}\tMean truncated length\t{mean_truncated_length}\n')
