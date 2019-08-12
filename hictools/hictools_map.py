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
        intermediate, samtools, bowtie2, sam_out):
    
    ''' Map R1 and R2 of HiC paired-end reads seperately. '''
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    out_format = 'SAM' if sam_out else 'BAM'


    if not intermediate:
        remove_intermediate = True
        intermediate = f'{sample}.fixmate.tmp.bam'
    else:
        remove_intermediate = False

    with tempfile.TemporaryFile() as tmp:
        try:
            for i, fastq in enumerate(infiles):
        
                flag = '0x80' if i else '0x40'
                read = 'R2' if i else 'R1'
                
                cmd1 = [f'{bowtie2}', '-x', f'{index}', '-U', f'{fastq}', 
                    '-p', f'{threads}', '--very-sensitive']
                cmd2 = ['awk', '-v', 'OFS=\t', 
                    f'!/^ *@/ {{$2 = $2+{flag}}} {{print}}']
                cmd3 = [f'{samtools}', 'sort', '-n', '-O', 'SAM', '-m', '1G', 
                    '-@', f'{threads}', '-o', f'{sample}-{read}.sorted.tmp.sam']
                
                sys.stderr.write(f'Mapping {fastq}.\n')
                sys.stderr.flush()
                with ExitStack() as stack:
                    
                    p1 = stack.enter_context(
                        Popen(cmd1, stdout = PIPE))
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
                        
            cmd4 = [f'{samtools}', 'merge', '-un', '-@', f'{threads}', '-', 
                f'{sample}-R1.sorted.tmp.sam', f'{sample}-R2.sorted.tmp.sam']
            cmd5 = [f'{samtools}', 'fixmate', '-pm', '-O', f'{out_format}',
                '-@', f'{threads}', '-', f'{intermediate}'] 

            with ExitStack() as stack:
                p4 = stack.enter_context(
                    Popen(cmd4, stdout = PIPE, stderr = tmp))
                p5 = stack.enter_context(
                    Popen(cmd5, stdin = p4.stdout, stderr = tmp))
                p4.stdout.close()
 
                exit_codes = [p.wait() for p in [p4, p5]]
                
                log.debug(f'Exit_codes for p4, p5: {exit_codes}.')
                if not all(ec is 0 for ec in exit_codes):
                    log.error('A sub-process returned a non-zero exit code.')
             
            os.remove(f'{sample}-R1.sorted.tmp.sam')
            os.remove(f'{sample}-R2.sorted.tmp.sam')       
            if remove_intermediate:
                os.remove(intermediate)
    
            cmd6 = [f'{samtools}', 'sort', '-l', '0', '-m', '2G', 
                '-@', f'{threads}', f'{intermediate}']
            cmd7 = [f'{samtools}', 'markdup', '-sr', '-O', 'SAM',
                '-@', f'{threads}', '-', '-']
            cmd8 = [f'{samtools}', 'sort', '-l', '0', '-n', '-m', '2G', 
                '-@', f'{threads}']
            cmd9 = [f'{samtools}', 'view', '-u', '-F', '12', '-q', '15', 
                '-@', f'{threads}']
            cmd10 = [f'{samtools}', 'fixmate', '-pr', '-O', 'SAM',
                '-@', f'{threads}', '-', '-']
            cmd11 = [f'{samtools}', 'view', '-O', f'{out_format}', '-f', '1', 
                '-@', f'{threads}', '-o', f'{output}' ]

            
            with ExitStack() as stack:
                p6 = stack.enter_context(
                    Popen(cmd6, stdout = PIPE, stderr = tmp))
                p7 = stack.enter_context(
                    Popen(cmd7, stdin = p6.stdout, stdout = PIPE))
                p6.stdout.close()
                p8 = stack.enter_context(
                    Popen(cmd8, stdin = p7.stdout, stdout = PIPE, stderr = tmp))
                p7.stdout.close()
                p9 = stack.enter_context(
                    Popen(cmd9, stdin = p8.stdout, stdout = PIPE, stderr = tmp))
                p8.stdout.close()
                p10 = stack.enter_context(
                    Popen(cmd10, stdin = p9.stdout, stdout = PIPE, stderr = tmp))
                p9.stdout.close()
                p11 = stack.enter_context(
                    Popen(cmd11, stdin = p10.stdout, stderr = tmp))
                p10.stdout.close()
                
                exit_codes = [p.wait() for p in [p6, p7, p8, p9, p10]]
                
                log.debug(f'Exit_codes for p6, p7, p8, p9, p10: {exit_codes}.')
                if not all(ec is 0 for ec in exit_codes):
                    log.error('A sub-process returned a non-zero exit code.')

        # Ensure tmp file is always written to log.
        finally:
            tmp.seek(0)
            for line in tmp:
                log.error(line.decode().rstrip())

