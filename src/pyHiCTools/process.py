#!/usr/bin/env python3

""" Process named-sorted SAM/BAM alignment files to identify fragment
    mapping, insert size, ditag size and relative orientation of pairs.
"""

import argparse, collections, bisect, time, re, sys, math, logging

import pyCommonTools as pct
import pyHiCTools as hic


def process(
    infile, output, read_gzip, samtools, sam_out, digest):

    log = pct.create_logger()

    mode = 'wt' if sam_out else 'wb'

    with pct.open_sam(output, mode, samtools = samtools) as out_obj, \
            pct.open_sam(infile, samtools = samtools) as in_obj, \
            pct.open_gzip(digest, 'rt', read_gzip) as digest:
        d = process_digest(digest)
        for line in in_obj:
            if line.startswith("@"):
                out_obj.write(line)
                continue
            else:
                try:
                    read1 = pct.Sam(line)
                    read2 = pct.Sam(next(in_obj))
                except StopIteration:
                    log.exception("Odd number of alignments in file")
                orientation, ditag_length, insert_size, interaction, fragment_seperation, read1_fragment, read2_fragment = run_filter(read1, read2, d)
                read1.optional['or:Z'] = orientation
                read2.optional['or:Z'] = orientation
                read1.optional['it:Z'] = interaction
                read2.optional['it:Z'] = interaction
                read1.optional['dt:i'] = ditag_length
                read2.optional['dt:i'] = ditag_length
                read1.optional['is:i'] = insert_size
                read2.optional['is:i'] = insert_size
                read1.optional['fs:i'] = fragment_seperation
                read2.optional['fs:i'] = fragment_seperation
                read1.optional['fn:i'] = read1_fragment
                read2.optional['fn:i'] = read2_fragment
                out_obj.write(read1.get_record())
                out_obj.write(read2.get_record())


def run_filter(read1, read2, digest):
    if not hic.valid_pair.is_valid(read1, read2):
       log.error(f'Invalid format in {read1.qname}.')
    read1, read2 = reorder_read_pair(read1, read2)
    orientation = get_orientation(read1, read2)
    read1_fragment = get_fragment(read1, digest)
    read2_fragment = get_fragment(read2, digest)
    ditag_length = tag_length(read1, read1_fragment) + tag_length(read2, read2_fragment)
    insert_size = read2.right_pos - read1.left_pos + 1
    interaction = interaction_type(read1, read2)
    fragment_seperation = abs(read2_fragment.number - read1_fragment.number)
    # ADD this to a dictionary??
    return orientation, ditag_length, insert_size, interaction, fragment_seperation, read1_fragment.number, read2_fragment.number


def get_orientation(read1, read2):
    """
    Return relative orientation of read pairs.
    Assumes read pairs have been ordered such that read 1 is five prime of read 2.
    """

    if read1.is_reverse:
        if read2.is_reverse:
            orientation = "Same-reverse"
        else:
            orientation = "Outward"
    else:
        if read2.is_reverse:
            orientation = "Inward"
        else:
            orientation = "Same-forward"
    return orientation

def reorder_read_pair(read1, read2):
    """
    Return a pair of reads such that read1 is left of read 2.
    Read pairs aligning to different chromosomes are returned unchanged.
    """

    if interaction_type(read1, read2) == "cis" and read1.left_pos > read2.left_pos:
        r1_reorder = read2
        r2_reorder = read1
    else:
        r1_reorder = read1
        r2_reorder = read2
    return r1_reorder, r2_reorder

def get_fragment(read, digest):
    rf_num = bisect.bisect_left(digest[read.rname], read.middle_pos)
    rf_start = 1 if rf_num == 0 else digest[read.rname][rf_num - 1] + 1
    rf_end = digest[read.rname][rf_num]
    return fragment(rf_num, rf_start, rf_end)

def interaction_type(read1, read2):
    if read1.rname != read2.rname:
        interaction = "trans"
    else:
        interaction = "cis"
    return interaction

def tag_length(read, fragment):
    if read.is_reverse:
        return read.five_prime_pos - fragment.start + 1
    else:
        return fragment.end - read.five_prime_pos + 1



class fragment:

    number = pct.IntRange(0, math.inf)
    start = pct.IntRange(0, math.inf)
    end = pct.IntRange(0, math.inf)

    def __init__(self, number, start, end):
        self.number = number
        self.start = start
        self.end = end

def process_digest(digest):
    d = {}
    for fragment in digest:
        [ref, start, end, number] = fragment.split()
        if not (int(start) > 0 and int(end) > 0):
            log.error(f'Negative fragment start/end positions on ref {ref}.')
        if ref not in d.keys():
            if not (int(start) == 1 and int(number) == 1):
                log.error(f'Invalid first fragment in ref {ref}.')
            d[ref] = []
        d[ref].append(int(end))
    return(d)

def pysam_test():
    import pysam
    samfile = pysam.AlignmentFile("/media/stephen/Data/hic_analysis/data/test.sam", "rb")
    mysamfile = open("/media/stephen/Data/hic_analysis/data/test2.sam")
    for read_pysam, read_mysam in zip(samfile.fetch(), mysamfile):
        myread = pct.Sam(read_mysam.split())
        assert read_pysam.get_reference_positions()[-1] + 1 == myread.right_pos
        assert read_pysam.get_reference_positions()[0] + 1 == myread.left_pos
        assert read_pysam.reference_length == myread.reference_length
        assert read_pysam.is_paired == myread.is_paired
        assert read_pysam.is_reverse == myread.is_reverse
        if myread.is_reverse:
            assert read_pysam.get_reference_positions()[-1] + 1 == myread.five_prime_pos
            assert read_pysam.get_reference_positions()[0] + 1 == myread.three_prime_pos
        else:
            assert read_pysam.get_reference_positions()[-1] + 1 == myread.three_prime_pos
            assert read_pysam.get_reference_positions()[0] + 1 == myread.five_prime_pos
    mysamfile.close()

