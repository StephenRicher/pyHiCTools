#!/usr/bin/env python3


''' HiCTools is a set of tools for analysing HiC data. '''

import re
import argparse
import pyCommonTools as pct
import pyHiCTools as hic
from version import __version__

def main():

    parser = pct.make_parser(prog='pyHiCTools', version=__version__)

    base_args = pct.get_base_args()
    subparser = pct.make_subparser(parser)

    # Parent parser options for commands with multi-threading.
    parallel_parser = argparse.ArgumentParser(add_help=False)
    parallel_parser.add_argument(
        '-@', '--threads', default=1,
        type=pct.positive_int,
        help='Threads for parallel processing.')

    fastq_input_arg = pct.get_in_arg(in_type = 'FASTQ')
    sam_input_arg = pct.get_in_arg(in_type = 'SAM')

    qc_arg = argparse.ArgumentParser(add_help=False)
    qc_arg.add_argument(
        '--qc', metavar='FILE', help='Output file for QC statistics.')

    # Digest sub-parser
    digest_parser = subparser.add_parser(
        'digest',
        description=hic.digest.__doc__,
        help='Generate in silico restriction digest of reference FASTA.',
        parents=[base_args, fastq_input_arg],
        epilog=parser.epilog)
    requiredNamed_digest = digest_parser.add_argument_group(
        'required named arguments')
    requiredNamed_digest.add_argument(
        '-r', '--restriction', required=True,
        type=restriction_seq,
        help='''Restriction cut sequence with "^" to indicate cut site.
                  e.g. Mbol = ^GATC''')
    digest_parser.set_defaults(function=hic.digest.digest)

    # Truncate sub-parser
    truncate_parser = subparser.add_parser(
        'truncate',
        description=hic.truncate.__doc__,
        help='Truncate FASTQ sequences at restriction enzyme ligation site.',
        parents=[base_args, qc_arg, fastq_input_arg],
        epilog=parser.epilog)
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

    # Process sub-parser
    process_parser = subparser.add_parser(
        'process',
        description=hic.process.__doc__,
        help='Determine HiC fragment mappings from named-sorted SAM/BAM file.',
        parents=[base_args, sam_input_arg],
        epilog=parser.epilog)
    requiredNamed_process = process_parser.add_argument_group(
        'required named arguments')
    requiredNamed_process.add_argument(
        '-d', '--digest', required=True,
        help='Output of pyHiCTools digest using same '
             'reference genome as used to map reads.')
    process_parser.set_defaults(function=hic.process.process)

    # Extract sub-parser
    extract_parser = subparser.add_parser(
        'extract',
        description=hic.extract.__doc__,
        help='Extract HiC information encoded by hic process from SAM/BAM.',
        parents=[base_args, sam_input_arg],
        epilog=parser.epilog)
    extract_parser.add_argument(
        '-n', '--sample', default=None,
        help='Sample name for input.')
    extract_parser.set_defaults(function=hic.extract.extract)

    # Filter sub-parser
    filter_parser = subparser.add_parser(
        'filter',
        description=hic.filter.__doc__,
        help='Filter SAM/BAM file processed with pyHiCTools process.',
        parents=[base_args, qc_arg, sam_input_arg],
        epilog=parser.epilog)
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
