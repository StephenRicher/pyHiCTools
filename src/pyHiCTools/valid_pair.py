#!/usr/bin/env python3

import pyCommonTools as pct


def is_valid(read1, read2):

    log = pct.create_logger()

    if read1.qname != read2.qname:
        log.error(
            f'Qname mismatch: {read1.qname} {read2.qname}. '
            'Is file name sorted?')
    elif not read1.is_paired and not read2.is_paired:
        log.error(f'{read1.qname} is not paired')
    elif read1.is_read1 == read2.is_read1:
        log.error(f'R1 and R2 flags in {read1.qname} not correctly set')
    elif read1.pnext != read2.left_pos or read2.pnext != read1.left_pos:
        log.error(f'Mate position mismatch in {read1.qname}.')
    else:
        return True
    return False
