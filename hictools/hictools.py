#!/usr/bin/env python3


''' HiCTools is a set of tools for analysing HiC data. '''

import sys, argparse, logging, select, re, time

from common_tools.exception_logger import *
from common_tools.gzip_opener import *

import  hictools_digest, hictools_truncate, \
        hictools_map, hictools_filter, \
        hictools_extract, hictools_process, \
        hictools_deduplicate

def main():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser = argparse.ArgumentParser(
        prog = 'hictools',
        description = __doc__,
        formatter_class = formatter_class,
        epilog = epilog)
        
    # Parent parser for options common to all sub commands.
    base_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)
    base_parser.add_argument(
        '-v', '--verbose', 
        action = 'store_true',
        help = 'Verbose logging for debugging.')
    base_parser.add_argument(
        '-l', '--log', nargs = '?',
        help = 'Log output file.')
        
    # Parent parser options for commands with multi-threading.
    parallel_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)   
    parallel_parser.add_argument(
        '-@', '--threads', default = 1,
        type = positive_int,
        help = 'Threads for parallel processing.')

    # Parent parser for gzipping output.
    gzip_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)
    gzip_parser.add_argument(
        '-z', '--gzip', 
        action = 'store_true', dest = 'write_gzip',
        help = 'Compress output using gzip')
    
    # Parent parser for gun-zipping input.
    gunzip_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)
    gunzip_parser.add_argument(
        '-u', '--gunzip', 
        action = 'store_true', dest = 'read_gzip',
        help = 'Read gzip compressed input.')
    
    # Parent parser options for SAM processing.
    sam_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)
    sam_parser.add_argument(
        '--samtools', default = 'samtools',
        help = 'Set path to samtools installation.')
    
    # Parent parser options for SAM input.
    sam_input_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)   
    sam_input_parser.add_argument(
        'infile', nargs = '?', default = '-',
        help = 'Input file in SAM/BAM format.')
        
    # Parent parser options for SAM output.
    sam_output_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)
    sam_output_parser.add_argument(
        '-S', '--sam', dest = 'sam_out',
        action = 'store_true',
        help = 'Output alignments in SAM format.')

    # Parent parser options for bowtie2 processing.
    bowtie2_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)   
    bowtie2_parser.add_argument(
        '--bowtie2', default = 'bowtie2',
        help = 'Set path to bowtie2 installation.')

    # Define subparser
    subparsers = parser.add_subparsers(
        title = 'required commands',
        description = '',
        dest = 'command',
        metavar = 'Commands',
        help = 'Description:')
    
    # Digest sub-parser
    digest_parser = subparsers.add_parser('digest',
        description = hictools_digest.description(),
        help = 'Generate in silico restriction digest of reference FASTA.', 
        parents = [base_parser, gzip_parser, gunzip_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    digest_parser.add_argument(
        'infile', nargs = '?', default = '-',
        help = 'Input reference in FASTA format.')
    digest_parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'Reference sequence digest file.')
    requiredNamed_digest = digest_parser.add_argument_group(
        'required named arguments')
    requiredNamed_digest.add_argument(
        '-r', '--restriction', required = True,
        type = restriction_seq, 
        help = '''Restriction cut sequence with "^" to indicate cut site.
                  e.g. Mbol = ^GATC''')
    digest_parser.set_defaults(function = hictools_digest.digest)
    
    # Truncate sub-parser
    truncate_parser = subparsers.add_parser('truncate',
        description = hictools_truncate.description(),
        help = 'Truncate FASTQ sequences at restriction enzyme ligation site.', 
        parents = [base_parser, gzip_parser, gunzip_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    truncate_parser.add_argument(
        'infile', nargs = '?', default = '-',
        help = 'Input file in FASTQ format.')
    truncate_parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'Truncated FASTQ file.')
    truncate_parser.add_argument(
        '-n', '--sample', default = None,
        help = 'Sample name in case infile name cannot be detected.')
    requiredNamed_truncate = truncate_parser.add_argument_group(
        'required named arguments')
    requiredNamed_truncate.add_argument(
        '-r', '--restriction', required = True,
        type = restriction_seq, 
        help = '''Restriction cut sequence with "^" to indicate cut site.
                  e.g. Mbol = ^GATC''')
    truncate_parser.set_defaults(function = hictools_truncate.truncate)
    
     # Map sub-parser
    map_command = 'map'
    map_parser = subparsers.add_parser('map',
        description = hictools_map.description(),
        help = 'Map R1 and R2 of HiC paired-end reads.',
        parents = [base_parser, parallel_parser, bowtie2_parser, 
            sam_parser, sam_output_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    map_parser.add_argument(
        'infiles', nargs = 2,
        help = 'Input R1 and R2 FASTQ files.')
    map_parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'R1 and R2 aligned sequences in SAM/BAM format.')
    map_parser.add_argument(
        '-i', '--intermediate', default = None,
        help = 'Path to write intermediate BAM prior to filtering.')
    map_parser.add_argument(
        '-n', '--sample', 
        default = f'sample_{time.strftime("%Y%m%d-%H%M%S")}',
        help = 'Sample name to prefix R1 and R2 BAMs.')
    requiredNamed_map = map_parser.add_argument_group(
        'required named arguments')
    requiredNamed_map.add_argument(
        '-x', '--index', required = True,
        help = 'Bowtie2 index of reference sequence.')
    map_parser.set_defaults(function = hictools_map.map)
    
    # Map sub-parser
    deduplicate_parser = subparsers.add_parser('deduplicate',
        description = hictools_deduplicate.description(),
        help = 'Deduplicate aligned HiC sequences processed by hictools map.',
        parents = [base_parser, parallel_parser, sam_parser, 
            sam_input_parser, sam_output_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    deduplicate_parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'Deduplicated sequences in SAM/BAM format.')
    deduplicate_parser.set_defaults(function = hictools_deduplicate.deduplicate)
    
    # Process sub-parser
    process_parser = subparsers.add_parser('process',
        description = hictools_process.description(),
        help = 'Determine HiC fragment mappings from '
               'named-sorted SAM/BAM file.', 
        parents = [base_parser, sam_parser, 
            sam_input_parser, sam_output_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    process_parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'Processed SAM/BAM file.')
    process_parser.add_argument(
        '-u', '--gunzip', 
        action = 'store_true', dest = 'read_gzip',
        help = 'Read gzip compressed digest file.')
    requiredNamed_process = process_parser.add_argument_group(
        'required named arguments')
    requiredNamed_process.add_argument(
        '-d', '--digest', required = True,
        help = 'Output of hictools digest using same '
               'reference genome as used to map reads.')
    process_parser.set_defaults(function = hictools_process.process)
                  
    # Extract sub-parser
    extract_parser = subparsers.add_parser('extract',
        description = hictools_extract.description(),
        help = 'Extract HiC information encoded by hic process from SAM/BAM.', 
        parents = [base_parser, gzip_parser, sam_parser, sam_input_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    extract_parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'HiC read pair information in long data format.')
    extract_parser.add_argument(
        '-n', '--sample', default = None,
        help = 'Sample name for input.')
    extract_parser.set_defaults(function = hictools_extract.extract)
    
    # Filter sub-parser
    filter_parser = subparsers.add_parser('filter',
        description = hictools_filter.description(),
        help = 'Filter named-sorted SAM/BAM file processed with hictools process.', 
        parents = [base_parser, sam_parser, 
            sam_input_parser, sam_output_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    filter_parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'Filtered SAM/BAM alignment file.')
    filter_parser.add_argument(
        '--min_inward', default = None, 
        type = positive_int, 
        help = 'Specify mininum insert size for inward facing read pairs.')
    filter_parser.add_argument(
        '--min_outward', default = None, 
        type = positive_int, 
        help = 'Specify mininum insert size for outward facing read pairs.')
    filter_parser.add_argument(
        '--max_ditag', default = None, 
        type = positive_int, 
        help = 'Specify maximum ditag size for read pairs.')
    filter_parser.set_defaults(function = hictools_filter.filter)
                  
    args = parser.parse_args()

    try:
        func = args.function
    except AttributeError:
        parser.print_help()
        sys.exit()
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    log_level = logging.DEBUG if args.verbose else None
    (logging.basicConfig(
        filename = args.log,
        format = log_format, 
        level = log_level))
    sys.excepthook = handle_exception

    args_dict = vars(args)
    [args_dict.pop(key) for key in ['command', 'function', 'verbose', 'log']]
    return func(**vars(args))

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

def positive_int(value):
    ''' Custom argparse type for positive integer. '''
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(
            f'{value} is not a positive integer.')
    return ivalue

if __name__ == '__main__':
    try:
        sys.exit(main())
    finally:
        logging.shutdown()

