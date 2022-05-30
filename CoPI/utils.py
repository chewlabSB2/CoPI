#!/usr/bin/env python3

import math
import string
import os
import numpy as np
import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from collections import namedtuple

table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

BASIC_AA = set(list(table.values()))
BASIC_AA.remove('_')
BASIC_AA = list(BASIC_AA)

ENCODING_MAP = {'A': 0, 'G': 1, 'C': 2, 'T': 3}
DECODING_LST = ['A', 'G', 'C', 'T']
LIMIT = 512
MISMATCH = 1

def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10**(q / -10) for q in range(n+1)]

def ave_qual(quals, qround=False, tab=errs_tab(128)):
    """
    Calculate average basecall quality of a read.
    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if quals:
        mq = -10 * math.log(sum([tab[ord(q) - 33] for q in quals]) / len(quals), 10)
        if qround:
            return round(mq)
        else:
            return mq
    else:
        return None

def update_qual(phred, add, c):
    phred = [(ord(a)-33)*c for a in phred]
    phred = [chr(math.floor(((ord(add[i]) - 33) + phred[i])/(c+1))+33) for i in range(len(phred))]
    return phred

qual_set = namedtuple('qual_set', 'q c')
def update_qual2(qual):
    q_updated = []
    for i, qu in enumerate(qual):
        qu = [(ord(q) - 33) * qu.c for q in qu.q]
        q_updated.append(qu)

    qual_sum = sum([q.c for q in qual])
    qual = [sum(x) for x in zip(*q_updated)]
    qual = [chr(math.floor(qual[i]/qual_sum) + 33) for i in range(len(qual))]

    return qual

def encode_DNA(seq):    
    code = 0
    for ch in seq:
        code <<= 2
        code |= ENCODING_MAP[ch]
    return code

def decode_DNA(code, seqlen):
    ret = ''
    for _ in range(seqlen):
        index = code & 3
        code >>= 2
        ret = DECODING_LST[index] + ret
    return ret

def reverseC(seq):
    RT = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
    reverseComplement = ''
    for i in seq:
        nt = RT.get(i)
        reverseComplement += nt
    return reverseComplement[::-1]

def readfq(fp, gzipped = True): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                l = l.decode() if gzipped else l
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            l = l.decode() if gzipped else l
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield (name, ''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                l = l.decode() if gzipped else l
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield (name, seq, ''.join(seqs)) # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield (name, seq, None) # yield a fasta record instead
                break

def is_tool(name):
    ## Check whether `name` is on PATH and marked as executable.
    from shutil import which
    return which(name) is not None

def print_error(message, **kwargs):
    """
    Formats *message* with *kwargs* using :meth:`str.format` and
    :func:`textwrap.dedent` and uses it to print an error message to
    ``sys.stderr``.
    """
    print("\nERROR: " + dedent(message.format(**kwargs)).lstrip("\n")+"\n", file = sys.stderr) 

def run_shell_command(cmd, raise_errors = False, extra_env = None):
    """
    Run the given command string via Bash with error checking.
    Returns True if the command exits normally.  Returns False if the command
    exits with failure and "raise_errors" is False (the default).  When
    "raise_errors" is True, exceptions are rethrown.

    If an *extra_env* mapping is passed, the provided keys and values are
    overlayed onto the default subprocess environment.
    """
    env = os.environ.copy()

    if extra_env:
        env.update(extra_env)

    shargs = ['-c', "set -euo pipefail; " + cmd]

    if os.name == 'posix':
        shellexec = ['/bin/bash']
    else:
        # We try best effort on other systems. For now that means nt/java.
        shellexec = ['env', 'bash']

    try:
        # Use check_call() instead of run() since the latter was added only in Python 3.5.
        subprocess.check_output(
            shellexec + shargs,
            shell = False,
            stderr = subprocess.STDOUT,
            env = env)

    except subprocess.CalledProcessError as error:
        print_error(
            "{out}\nshell exited {rc} when running: {cmd}{extra}",
            out = error.output,
            rc  = error.returncode,
            cmd = cmd,
            extra = "\nAre you sure this program is installed?" if error.returncode==127 else "",
        )
        if raise_errors:
            raise
        else:
            return False

    except FileNotFoundError as error:
        print_error(
            """
            Unable to run shell commands using {shell}!
            Module requires {shell} to be installed.
            """,
            shell = ' and '.join(shellexec)
        )
        if raise_errors:
            raise
        else:
            return False

    else:
        return True  