#!/usr/bin/env python3


''' Map R1 and R2 of HiC paired-end reads to reference genome seperately.
'''

import os, sys, time, logging, subprocess

def description():
    
    ''' Returns top-level docstring. Useful for providing
        descriptions to sub-parsers after importing.
        '''
        
    return __doc__

def map(infiles, output, sample, index, threads, samtools, bowtie2, sam_out):
    
    ''' Map R1 and R2 of HiC paired-end reads seperately. '''
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    # Get file directory of this script where bash script also kept
    path = os.path.dirname(os.path.realpath(__file__)) + '/'
    bash_script = 'hictools_map_bash.sh'

    cmd = ([path + bash_script, sample, output, 
            index, str(threads), samtools, bowtie2] 
           + infiles)
    
    with subprocess.Popen(
        cmd, stderr = subprocess.PIPE) as p:
        # Capture stderr from bash script and redirect 
        for err in p.stderr:
            log.error(err.decode().rstrip())
    
    while p.poll() == None:
        time.sleep(0.5)
    
    if p.returncode:
        log.error(f'{bash_script} returned non-zero exit status')
        sys.exit(p.returncode)
        
