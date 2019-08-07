#!/usr/bin/env python3


''' Deduplicate aligned HiC sequences processed by hictools map.
'''

import os, sys, time, logging, tempfile
from subprocess import Popen, PIPE
from contextlib import ExitStack
from common_tools.gzip_opener import *

def description():
    
    ''' Returns top-level docstring. Useful for providing
        descriptions to sub-parsers after importing.
        '''
        
    return __doc__

def deduplicate(infile, output, threads, samtools, sam_out):

    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    out_format = 'SAM' if sam_out else 'BAM'
    stdin = sys.stdin if infile == '-' else None

    cmd1 = [f'{samtools}', 'sort', '-l', '0', '-m', '1G', '-@', f'{threads}', 
        f'{infile}']
    cmd2 = [f'{samtools}', 'markdup', '-sr', '-@', f'{threads}', '-', '-']
    cmd3 = [f'{samtools}', 'sort', '-n', '-m', '1G', '-@', f'{threads}', 
        '-O', f'{out_format}', '-o', f'{output}', '-'] 
        
    with ExitStack() as stack:
        try:
            tmp = stack.enter_context(tempfile.TemporaryFile())
            p1 = stack.enter_context(
                Popen(cmd1, stdin = stdin, stdout = PIPE, stderr = tmp))
            p2 = stack.enter_context(
                Popen(cmd2, stdin = p1.stdout, stdout = PIPE, 
                    stderr = sys.stderr))
            p1.stdout.close()
            p3 = stack.enter_context(
                Popen(cmd3, stdin = p2.stdout, stderr = tmp))
            p2.stdout.close()
                    
            exit_codes = [p.wait() for p in [p1, p2, p3]]
            log.debug(f'Exit_codes for p1, p2, p3: {exit_codes}.')
            if not all(ec is 0 for ec in exit_codes):
                log.error('A sub-process returned a non-zero exit code.')
    
        # Ensure tmp file is always written to log.
        finally:
            tmp.seek(0)
            for line in tmp:
                log.error(line.decode().rstrip())
