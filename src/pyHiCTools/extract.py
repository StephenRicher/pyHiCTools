#!/usr/bin/env python3

""" Extract HiC read pair information encodiing by hictools process.
    Useful for assessing contamination and determining parameters for
    hictools filter.
"""

import sys
import pyCommonTools as pct
import pyHiCTools as hic


def extract(infile, output, samtools, sample, write_gzip):

    log = pct.create_logger()

    if not sample:
        sample = infile

    with pct.open_sam(infile, samtools=samtools) as in_obj, \
            pct.open_gzip(output, 'wt', write_gzip) as out_obj:
        log.info(f'Writing output to {output}.')
        out_obj.write(
            'sample\torientation\tinteraction_type\tditag_length\t'
            'insert_size\tfragment_seperation\n')
        for i, line in enumerate(in_obj):
            if line.startswith('@'):
                continue
            else:
                try:
                    read1 = pct.Sam(line)
                    read2 = pct.Sam(next(in_obj))
                except StopIteration:
                    log.exception('Odd number of alignments in file')
                    sys.exit(1)
                if not hic.valid_pair.is_valid(read1, read2):
                    log.error(f'Invalid format in {read1.qname}.')
                out_obj.write(
                    f'{sample}\t{read1.optional["or:Z"]}\t'
                    f'{read1.optional["it:Z"]}\t{read1.optional["dt:i"]}\t'
                    f'{read1.optional["is:i"]}\t{read1.optional["fs:i"]}\n')
