#!/usr/bin/env python3

""" Filter HiC read pairs to remove potential sources of contamination.
    Input should be in SAM/BAM format and have been processed by
    hictools process.
    """

import sys
import pyCommonTools as pct
import pyHiCTools as hic


def filter(infile, output, sample, samtools, sam_out,
           min_inward, min_outward, min_ditag, max_ditag):

    ''' Iterate through each infile. '''

    log = pct.create_logger()
    inputs = [min_inward, min_outward, max_ditag, min_ditag]
    if all(i is None for i in inputs):
        log.error('No filter settings defined.')
        sys.exit(1)

    if not sample:
        sample = infile

    mode = 'wt' if sam_out else 'wb'

    with pct.open_sam(output, mode, samtools=samtools) as out_obj, \
            pct.open_sam(infile, samtools=samtools) as in_obj:
        log.info(f'Writing output to {output}.')
        total = 0
        retained = 0
        invalid = 0
        above_ditag = 0
        below_ditag = 0
        same_fragment = 0
        below_min_inward = 0
        below_min_outward = 0

        for line in in_obj:
            if line.startswith("@"):
                out_obj.write(line)
            else:
                try:
                    read1 = pct.Sam(line)
                    read2 = pct.Sam(next(in_obj))
                    total += 1
                except StopIteration:
                    log.exception('Odd number of alignments in file.')
                    sys.exit(1)
                if not hic.valid_pair.is_valid(read1, read2):
                    log.error(f'Invalid format in {read1.qname}.')
                    invalid += 1
                    continue
                if max_ditag is not None:
                    if read1.optional['dt:i'] > max_ditag:
                        above_ditag += 1
                        continue
                if min_ditag is not None:
                    if read1.optional['dt:i'] < min_ditag:
                        below_ditag += 1
                        continue
                if read1.optional['it:Z'] == "cis":
                    if read1.optional['fs:i'] == 0:
                        same_fragment += 1
                        continue
                    if read1.optional['or:Z'] == 'Inward':
                        if min_inward is not None:
                            if read1.optional['is:i'] < min_inward:
                                below_min_inward += 1
                                continue
                    elif read1.optional['or:Z'] == 'Outward':
                        if min_outward is not None:
                            if read1.optional['is:i'] < min_outward:
                                below_min_outward += 1
                                continue
                retained += 1
                out_obj.write(read1.get_record())
                out_obj.write(read2.get_record())

        sys.stderr.write(
            f'{sample}\tTotal\t{total}\n'
            f'{sample}\tRetained\t{retained}\n'
            f'{sample}\tFiltered\t{total - retained}\n'
            f'{sample}\tInvalid\t{invalid}\n'
            f'{sample}\tDitag < {min_ditag}bp\t{above_ditag}\n'
            f'{sample}\tDitag > {max_ditag}bp\t{below_ditag}\n'
            f'{sample}\tSame fragment\t{same_fragment}\n'
            f'{sample}\tInward insert < {min_inward}bp\t{below_min_inward}\n'
            f'{sample}\tOutward insert < {min_outward}bp\t{below_min_outward}\n')
