#!/usr/bin/env python3


''' HiCTools is a set of tools for analysing HiC data. '''

import re
import time
import argparse
import fileinput
import pyCommonTools as pct
import pyHiCTools as hic

def main():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser = argparse.ArgumentParser(
        prog='pyHiCTools',
        description=__doc__,
        formatter_class=formatter_class,
        epilog=epilog)

    subparsers, base_parser = pct.set_subparser(parser)

    # Parent parser options for commands with multi-threading.
    parallel_parser = argparse.ArgumentParser(
        formatter_class=formatter_class,
        add_help=False)
    parallel_parser.add_argument(
        '-@', '--threads', default=1,
        type=pct.positive_int,
        help='Threads for parallel processing.')

    # Parent parser for gzipping output.
    gzip_parser = argparse.ArgumentParser(
        formatter_class=formatter_class,
        add_help=False)
    gzip_parser.add_argument(
        '-z', '--gzip',
        action='store_true', dest='write_gzip',
        help='Compress output using gzip')

    # Parent parser for gun-zipping input.
    gunzip_parser = argparse.ArgumentParser(
        formatter_class=formatter_class,
        add_help=False)
    gunzip_parser.add_argument(
        '-u', '--gunzip',
        action='store_true', dest='read_gzip',
        help='Read gzip compressed input.')

    # Parent parser options for SAM processing.
    sam_parser = argparse.ArgumentParser(
        formatter_class=formatter_class,
        add_help=False)
    sam_parser.add_argument(
        '--samtools', default='samtools',
        help='Set path to samtools installation.')

    sam_input_parser2 = argparse.ArgumentParser(
        formatter_class=formatter_class,
        add_help=False)
    sam_input_parser2.add_argument(
        'infile', nargs='?', default=[],
        help='Input file in SAM format.')

    # Parent parser options for SAM input.
    sam_input_parser = argparse.ArgumentParser(
        formatter_class=formatter_class,
        add_help=False)
    sam_input_parser.add_argument(
        'infile', nargs='?', default='-',
        help='Input file in SAM/BAM format.')

    qc_parser = argparse.ArgumentParser(
        formatter_class=formatter_class,
        add_help=False)
    qc_parser.add_argument(
        '--qc', nargs='?',
        help='Output file for QC statistics.')

    # Parent parser options for SAM output.
    sam_output_parser = argparse.ArgumentParser(
        formatter_class=formatter_class,
        add_help=False)
    sam_output_parser.add_argument(
        '-S', '--sam', dest='sam_out',
        action='store_true',
        help='Output alignments in SAM format.')

    # Parent parser options for bowtie2 processing.
    bowtie2_parser = argparse.ArgumentParser(
        formatter_class=formatter_class,
        add_help=False)
    bowtie2_parser.add_argument(
        '--bowtie2', default='bowtie2',
        help='Set path to bowtie2 installation.')

    # Digest sub-parser
    digest_parser = subparsers.add_parser(
        'digest',
        description=hic.digest.__doc__,
        help='Generate in silico restriction digest of reference FASTA.',
        parents=[base_parser, gzip_parser, gunzip_parser],
        formatter_class=formatter_class,
        epilog=epilog)
    digest_parser.add_argument(
        'infile', nargs='?', default='-',
        help='Input reference in FASTA format.')
    digest_parser.add_argument(
        '-o', '--output', nargs='?', default='-',
        help='Reference sequence digest file.')
    requiredNamed_digest = digest_parser.add_argument_group(
        'required named arguments')
    requiredNamed_digest.add_argument(
        '-r', '--restriction', required=True,
        type=restriction_seq,
        help='''Restriction cut sequence with "^" to indicate cut site.
                  e.g. Mbol = ^GATC''')
    digest_parser.set_defaults(function=hic.digest.digest)

    # Truncate sub-parser
    truncate_parser = subparsers.add_parser(
        'truncate',
        description=hic.truncate.__doc__,
        help='Truncate FASTQ sequences at restriction enzyme ligation site.',
        parents=[base_parser, gzip_parser, gunzip_parser],
        formatter_class=formatter_class,
        epilog=epilog)
    truncate_parser.add_argument(
        'infile', nargs='?', default='-',
        help='Input file in FASTQ format.')
    truncate_parser.add_argument(
        '-o', '--output', nargs='?', default='-',
        help='Truncated FASTQ file.')
    truncate_parser.add_argument(
        '-n', '--sample', default=None,
        help='Sample name in case infile name cannot be detected.')
    requiredNamed_truncate = truncate_parser.add_argument_group(
        'required named arguments')
    requiredNamed_truncate.add_argument(
        '-r', '--restriction', required=True,
        type=restriction_seq,
        help=('Restriction cut sequence with "^" to indicate cut site.'
              'e.g. Mbol = ^GATC'))
    truncate_parser.set_defaults(function=hic.truncate.truncate)

    # Map sub-parser
    map_parser = subparsers.add_parser(
        'map',
        description=hic.map.__doc__,
        help='Map R1 and R2 of HiC paired-end reads.',
        parents=[base_parser, parallel_parser, bowtie2_parser,
                 sam_parser, sam_output_parser],
        formatter_class=formatter_class,
        epilog=epilog)
    map_parser.add_argument(
        'infiles', nargs=2,
        help='Input R1 and R2 FASTQ files.')
    map_parser.add_argument(
        '-o', '--output', nargs='?', default='-',
        help='R1 and R2 aligned sequences in SAM/BAM format.')
    map_parser.add_argument(
        '-i', '--intermediate', default=None,
        help='Path to write intermediate BAM prior to filtering.')
    map_parser.add_argument(
        '-n', '--sample',
        default=f'sample_{time.strftime("%Y%m%d-%H%M%S")}',
        help='Sample name to prefix R1 and R2 BAMs.')
    sensitivity_options = ['very-fast', 'fast', 'sensitive', 'very-sensitive']
    map_parser.add_argument(
        '--sensitivity',
        default='very-sensitive',
        choices=sensitivity_options,
        help='Set bowtie2 alignment sensitivity.')
    requiredNamed_map = map_parser.add_argument_group(
        'required named arguments')
    requiredNamed_map.add_argument(
        '-x', '--index', required=True,
        help='Bowtie2 index of reference sequence.')
    map_parser.set_defaults(function=hic.map.map)

    # Deduplicate sub-parser
    deduplicate_parser = subparsers.add_parser(
        'deduplicate',
        description=hic.deduplicate.__doc__,
        help='Deduplicate aligned HiC sequences processed by pyHiCTools map.',
        parents=[base_parser, parallel_parser, sam_parser,
                 sam_input_parser, sam_output_parser],
        formatter_class=formatter_class,
        epilog=epilog)
    deduplicate_parser.add_argument(
        '-o', '--output', nargs='?', default='-',
        help='Deduplicated sequences in SAM/BAM format.')
    deduplicate_parser.set_defaults(function=hic.deduplicate.deduplicate)

    # Process sub-parser
    process_parser = subparsers.add_parser(
        'process',
        description=hic.process.__doc__,
        help='Determine HiC fragment mappings from named-sorted SAM/BAM file.',
        parents=[base_parser, sam_input_parser2],
        formatter_class=formatter_class,
        epilog=epilog)
    process_parser.add_argument(
        '-u', '--gunzip',
        action='store_true', dest='read_gzip',
        help='Read gzip compressed digest file.')
    requiredNamed_process = process_parser.add_argument_group(
        'required named arguments')
    requiredNamed_process.add_argument(
        '-d', '--digest', required=True,
        help='Output of pyHiCTools digest using same '
             'reference genome as used to map reads.')
    process_parser.set_defaults(function=hic.process.process)

    # Extract sub-parser
    extract_parser = subparsers.add_parser(
        'extract',
        description=hic.extract.__doc__,
        help='Extract HiC information encoded by hic process from SAM/BAM.',
        parents=[base_parser, gzip_parser, sam_parser, sam_input_parser],
        formatter_class=formatter_class,
        epilog=epilog)
    extract_parser.add_argument(
        '-o', '--output', nargs='?', default='-',
        help='HiC read pair information in long data format.')
    extract_parser.add_argument(
        '-n', '--sample', default=None,
        help='Sample name for input.')
    extract_parser.set_defaults(function=hic.extract.extract)

    # Filter sub-parser
    filter_parser = subparsers.add_parser(
        'filter',
        description=hic.filter.__doc__,
        help='Filter SAM/BAM file processed with pyHiCTools process.',
        parents=[base_parser, qc_parser, sam_input_parser2],
        formatter_class=formatter_class,
        epilog=epilog)
    filter_parser.add_argument(
        '--min_inward', default=None,
        type=pct.positive_int,
        help='Specify mininum insert size for inward facing read pairs.')
    filter_parser.add_argument(
        '--min_outward', default=None,
        type=pct.positive_int,
        help='Specify mininum insert size for outward facing read pairs.')
    filter_parser.add_argument(
        '--min_ditag', default=None,
        type=pct.positive_int,
        help='Specify minimum ditag size for read pairs.')
    filter_parser.add_argument(
        '--max_ditag', default=None,
        type=pct.positive_int,
        help='Specify maximum ditag size for read pairs.')
    filter_parser.add_argument(
        '-n', '--sample', default=None,
        help='Sample name in case infile name cannot be detected.')
    filter_parser.set_defaults(function=hic.filter.filter)

    return (pct.execute(parser))


def restriction_seq(value):

    ''' Custom argument type for restriction enzyme argument. '''

    if value.count('^') != 1:
        raise argparse.ArgumentTypeError(
            f'Restriction site {value} must contain one "^" at cut site.')
    elif re.search('[^ATCG^]', value, re.IGNORECASE):
        raise argparse.ArgumentTypeError(
            f'Restriction site {value} must only contain "ATCG^".')
    else:
        return value.upper()
