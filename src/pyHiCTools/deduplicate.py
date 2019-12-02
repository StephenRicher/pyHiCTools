#!/usr/bin/env python3


''' Deduplicate aligned HiC sequences processed by hictools map.
'''

import sys
import tempfile
from subprocess import Popen, PIPE
from contextlib import ExitStack
import pyCommonTools as pct


def deduplicate(infile, output, threads, samtools, sam_out):

    log = pct.create_logger()

    out_format = 'SAM' if sam_out else 'BAM'
    stdin = sys.stdin if infile == '-' else None
    cmd1 = [f'{samtools}', 'sort', '-O', 'SAM', '-m', '1G',
            '-@', f'{threads}', infile]
    cmd2 = [f'{samtools}', 'markdup', '-sr', '-O', 'SAM',
            '-@', f'{threads}', '-', '-']
    cmd3 = [f'{samtools}', 'sort', '-l', '0', '-n', '-m', '1G',
            '-@', f'{threads}']
    cmd4 = [f'{samtools}', 'fixmate', '-p', '-O', 'SAM',
            '-@', f'{threads}', '-', '-']
    cmd5 = [f'{samtools}', 'view', '-O', f'{out_format}', '-f', '1',
            '-@', f'{threads}', '-o', f'{output}']

    with ExitStack() as stack:
        try:
            tmp = stack.enter_context(tempfile.TemporaryFile())
            p1 = stack.enter_context(
                Popen(cmd1, stdin=stdin, stdout=PIPE, stderr=tmp))
            p2 = stack.enter_context(
                Popen(cmd2, stdin=p1.stdout, stdout=PIPE))
            p1.stdout.close()
            p3 = stack.enter_context(
                Popen(cmd3, stdin=p2.stdout, stdout=PIPE, stderr=tmp))
            p2.stdout.close()
            p4 = stack.enter_context(
                Popen(cmd4, stdin=p3.stdout, stdout=PIPE, stderr=tmp))
            p3.stdout.close()
            p5 = stack.enter_context(
                Popen(cmd5, stdin=p4.stdout, stderr=tmp))
            p4.stdout.close()
            exit_codes = [p.wait() for p in [p1, p2, p3, p4, p5]]
            log.debug(f'Exit_codes for p1, p2, p3, p4, p5: {exit_codes}.')
            if not all(ec is 0 for ec in exit_codes):
                log.error(exit_codes)
                log.error('A sub-process returned a non-zero exit code.')

        # Ensure tmp file is always written to log.
        finally:
            tmp.seek(0)
            for line in tmp:
                log.error(line.decode().rstrip())
