#!/usr/bin/env python3


''' Map R1 and R2 of HiC paired-end reads to reference genome seperately.
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

def map(infiles, output, index, threads, sample,
        samtools, bowtie2, bowtie2_log, sam_out):
    
    ''' Map R1 and R2 of HiC paired-end reads seperately. '''
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    out_format = 'S' if sam_out else 'b'

    with tempfile.TemporaryFile() as tmp:
        try:
            for i, fastq in enumerate(infiles):
        
                flag = '0x80' if i else '0x40'
                read = 'R2' if i else 'R1'
                
                cmd1 = [f'{bowtie2}', '-x', f'{index}', '-U', f'{fastq}', 
                    '-p', f'{threads}', '--very-sensitive']
                cmd2 = ['awk', '-v', 'OFS=\t', 
                    f'!/^ *@/ {{$2 = $2+{flag}}} {{print}}']
                cmd3 = [f'{samtools}', 'sort', '-n', '-O', 'bam', '-m', '100M', 
                    '-@', f'{threads}', '-o', f'{sample}-{read}.sorted.bam']
        
                with ExitStack() as stack:
                    bowtie_log = stack.enter_context(
                        open(f'{bowtie2_log}.{read}.txt', 'w'))
                    p1 = stack.enter_context(
                        Popen(cmd1, stdout = PIPE, stderr = bowtie_log))
                    p2 = stack.enter_context(
                        Popen(cmd2, stdin = p1.stdout, 
                            stdout = PIPE, stderr = tmp))
                    p1.stdout.close()
                    p3 = stack.enter_context(
                        Popen(cmd3, stdin = p2.stdout, stderr = tmp))
                    p2.stdout.close()
                    
                    exit_codes = [p.wait() for p in [p1, p2, p3]]
                    log.debug(f'Exit_codes for p1, p2, p3: {exit_codes}.')
                    if not all(ec is 0 for ec in exit_codes):
                        log.error(
                            'A sub-process returned a non-zero exit code.')
                        
            cmd4 = [f'{samtools}', 'merge', '-n', '-@', f'{threads}', 
                '-', f'{sample}-R1.sorted.bam', f'{sample}-R2.sorted.bam']
            cmd5 = [f'{samtools}', 'fixmate', '-pr', '-@', f'{threads}', 
                '-', '-'] 
            cmd6 = [f'{samtools}', 'view', '-u', '-F', '8', '-q', '15']
            cmd7 = [f'{samtools}', 'fixmate', '-pm', '-@', f'{threads}', 
                '-', '-']
            cmd8 = [f'{samtools}', 'view', f'-{out_format}h', '-f', '1', 
                '-o', f'{output}']

            with ExitStack() as stack:
                p4 = stack.enter_context(
                    Popen(cmd4, stdout = PIPE, stderr = tmp))
                p5 = stack.enter_context(
                    Popen(cmd5, stdin = p4.stdout, stdout = PIPE, stderr = tmp))
                p4.stdout.close()
                p6 = stack.enter_context(
                    Popen(cmd6, stdin = p5.stdout, stdout = PIPE, stderr = tmp))
                p5.stdout.close()
                p7 = stack.enter_context(
                    Popen(cmd7, stdin = p6.stdout, stdout = PIPE, stderr = tmp))
                p6.stdout.close()
                p8 = stack.enter_context(
                    Popen(cmd8, stdin = p7.stdout, stderr = tmp))
                p7.stdout.close()

                exit_codes = [p.wait() for p in [p4, p5, p6, p7]]
                
                log.debug(f'Exit_codes for p4, p5, p6, p7: {exit_codes}.')
                if not all(ec is 0 for ec in exit_codes):
                    log.error('A sub-process returned a non-zero exit code.')

        # Ensure tmp file is always written to log.
        finally:
            tmp.seek(0)
            for line in tmp:
                log.error(line.decode().rstrip())

