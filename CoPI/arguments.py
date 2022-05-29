#!/usr/bin/env python3

import argparse
from pathlib import Path
import warnings
import logging
import logging.handlers
import sys

THRESHOLD = 50
THREADS = 2
INDEL_RANGE = 16
BASIC_NT = ['A','C','G','T']

from CoPI.utils import table, reverseC

class Query():
    def __init__(self, query):
        self.seq = query
        self.bitap = self._preprocess(query)
        self.len = len(query)
        self.aa_len = self.len // 3
        self.aa = [table[query[i:i+3]] for i in range(0, len(query), 3)]

    def _preprocess(self, query):
        
        for_bitap_dict = {}
        for letter in BASIC_NT:
            letterPositionInQuery = 0
            for symbol in query:
                letterPositionInQuery = letterPositionInQuery << 1
                letterPositionInQuery |= int(letter != symbol)
            for_bitap_dict[letter] = letterPositionInQuery
        return for_bitap_dict

class Parameters():
    ALPHABETS = ['A', 'C', 'G', 'T']

    def __init__(self, 
                 reference, 
                 anchor1, 
                 anchor2, 
                 read1=None, 
                 read2=None, 
                 merge=True, 
                 paired=True, 
                 unpaired=[], 
                 threads = 2, 
                 collapse = True,
                 threshold = THRESHOLD,
                 limit=INDEL_RANGE, 
                 plot=False, 
                 prefix='CoPI', 
                 label1 = "T589", 
                 label2 = "D590",
                 min_indel_len = 7,):
        
        self.ref_seq = reference.upper()
        self.seq_start = Query(anchor1.upper())
        self.seq_end = Query(anchor2.upper())
        self.wt_distance = -1
        self.anchor_dist = -1
        self.anchor_diff_seq = None
        self.seq_start_rev = Query(reverseC(self.seq_start.seq))
        self.seq_end_rev = Query(reverseC(self.seq_end.seq))
        self.collapse = collapse
        self.threshold = threshold
        self.threads = threads
        self.tolerate = min([self.seq_start.len, self.seq_end.len]) // 9
        self.indel_range = limit
        self.plot = plot
        self.prefix = prefix
        self.min_indel_len = min_indel_len * 3
        
        self.read1 = read1
        self.read2 = read2
        self.paired = paired
        self.merge = merge
        self.readM = None
        self.readU = unpaired
        self.fastqC = 0 #self._calculate_number_fastq()
        self._check_paired()

        self.label1 = label1    
        self.label2 = label2

    def _check_paired(self):
        if self.paired: 
            return

        self.readU.append(self.read1)
        self.readU.append(self.read2)

        self.read1 = None
        self.read2 = None

    def _calculate_number_fastq(self):
        count = 0
        if self.read1 and self.read2:
            if self.paired: 
                count += 1
            else: 
                count += 2

        if self.merge or self.readM:
            count += 1

        if self.readU:
            for u in self.readU:
                count += 1

        self.fastqC = count

required_args = ['reference',
                 'anchor1',
                 'anchor2',
]

read_args = [ 'read1',
              'read2']

bool_args = {'collapse' : True,
             'plot' : True,
             'paired' : True,
             'merge' : True,
}

int_args = ['threads',
            'threshold',
            'limit',
]

unpaired_args = ['unpaired']

optional_args = {'prefix' : "CoPI",
                 'label1' : "label1",
                 'label2' : "label2",
}

ARGUMENT = required_args + read_args + int_args + unpaired_args

class softwareError(Exception):
    pass

def print_args(args, logger=None, exclude=None):
    """Print and format all arguments from the command line"""
    if exclude is None:
        exclude = []
    
    dirs = dir(args)
    m = max([len(a) for a in dirs if a[0] != "_"])
    for attr in dirs:
        if attr[0] != "_" and attr not in exclude and attr.lower() == attr:
            rec = getattr(args, attr)
            if isinstance(rec, Query):
                record = "{a}={b}".format(a=attr, m=m, b=rec.seq)
            else:
                record = "{a}={b}".format(a=attr, m=m, b=rec)

            if logger is not None:
                logger.info(record)
            else:
                print(record)

def get_toml(args, validate=True, logger=None):
    
    def check_args(toml_dict, arg):
        if not toml_dict:
            raise softwareError(f'{arg} is required for CoPI. Please check the argument file!')

    p = Path(args.config)
    if not p.is_file():
        raise FileNotFoundError("Argument File not found at '{}'".format(args.config))

    toml_dict = {k:None for k in ARGUMENT}
    toml_dict['unpaired'] = set()
    toml_dict.update(optional_args)
    toml_dict.update(bool_args)

    with open(args.config, 'r') as f:
        for line in f:
            line = line.strip('\n')
            col = line.split('=')
            try:
                arg = col[0].lower()
                if arg == 'unpaired':
                    toml_dict['unpaired'].add(col[1])
                else:
                    toml_dict[arg] = col[1]
                
            except:
                warnings.warn(f"{col[0]} not part of the necessary arguments. Did you misspell?")

    [check_args(toml_dict[a], a) for a in required_args]

    if not toml_dict['unpaired']:
        if not toml_dict['read2'] or not toml_dict['read1']:
            raise softwareError(f'A Fastq File (Read1 + Read2 or Unpaired) is required for CoPI. Please check the argument file!')    

    # Check reads path
    for r in read_args:
        if not Path(toml_dict[r]): continue
        reference_path = Path(toml_dict[r])
        if not reference_path.is_file():
            raise FileNotFoundError("{} not found at '{}'".format(r, reference_path))

    if toml_dict['unpaired']:
        for c, r in enumerate(list(toml_dict['unpaired'])):
            reference_path = Path(r)
            if not reference_path.is_file():
                raise FileNotFoundError("{} not found at '{}'".format(f"{c}th unpaired", reference_path))
        toml_dict['unpaired'] = list(toml_dict['unpaired'])

    for r in bool_args:
        if toml_dict[r].isdigit():
            toml_dict[r] = int(toml_dict[r])
        else:
            toml_dict[r] = int(toml_dict[r].lower() == 'true')

    for r in int_args:
        if not toml_dict[r]:
            if r == 'threads':
                toml_dict[r] = THREADS
            elif r == 'threshold':
                toml_dict[r] = THRESHOLD
            elif r == 'limit':
                toml_dict[r] = INDEL_RANGE

        toml_dict[r] = int(toml_dict[r])

    para = Parameters(**toml_dict)
    
    switch, anchor_dist = check_anchors(para)
    if switch:
        tmp = para.seq_start
        para.seq_start = para.seq_end
        para.seq_end = tmp

    para.anchor_dist = anchor_dist 
    para.anchor_diff_seq = para.ref_seq[para.ref_seq.find(para.seq_start.seq) + para.seq_start.len: para.ref_seq.find(para.seq_end.seq)]
    para.wt_distance = para.ref_seq.find(para.seq_end.seq) + para.seq_end.len - para.ref_seq.find(para.seq_start.seq) 

    print_args(para, logger=logger)
    return para

def pipeline(func, parameters, software):
    cmdline = func(parameters)
    success = run_shell_command(cmdline)
    if not success:
        raise softwareError(f"Error while running: {software}")

def fastp_process(para):
    if para.merge:
        cmd = f"fastp -w {para.threads} -i {para.read1} -I {para.read2} --merge --merged_out {para.prefix}.merged.fastq --out1 {para.prefix}.unmerged-1.fastq --out2 {para.prefix}.unmerged-2.fastq"
    else:
        cmd = f"fastp -w {para.threads} -i {para.read1} -I {para.read2}"
    return cmd

def gzip_process(file):
    cmd = f"gzip {file}"
    return cmd

def execute_fastp(para):
    if sys.platform == 'win32':
        warnings.warn("fastp not installed due to Windows OS! Skipping Merging")
        return 0, 0

    if not is_tool('fastp'):
        warnings.warn("fastp not installed! Check again! Skipping Merging")
        return 0, 0
    
    gzipped = 1
    if not is_tool('gzip'):
        warnings.warn("gzip not installed! Skipping gzip")
        gzipped = 0

    pipeline(fastp_process, para, "fastp")
    if para.merge and gzipped:
        for f in [f"{para.prefix}.merged.fastq", f"{para.prefix}.unmerged-1.fastq", f"{para.prefix}.unmerged-2.fastq"]:
            pipeline(gzip_process, f, "gzip")
        
    return 1, gzipped

def check_anchors(para):

    def check_nt(nt):
        for n in nt: 
            if n not in BASIC_NT:
                raise softwareError(f"{n} nucleotide not accepted. Please use generic Nucleotides (A C G T)!")

    for_anchor = para.seq_start
    rev_anchor = para.seq_end
    reference = para.ref_seq

    if for_anchor.len <= 6:
        raise softwareError("Length of Anchor1 should be at least 9")

    if rev_anchor.len <= 6:
        raise softwareError("Length of Anchor2 should be at least 9")

    check_nt(list(set(reference)))
    check_nt(list(set(for_anchor.seq)))
    check_nt(list(set(rev_anchor.seq)))

    for_anchor_pos = reference.find(for_anchor.seq)
    rev_anchor_pos = reference.find(rev_anchor.seq)
    
    for_anchor_end = for_anchor_pos + for_anchor.len
    if for_anchor.seq == rev_anchor.seq:
        warnings.warn("Anchor1 and Anchor2 are the same!")
        rev_anchor_pos = reference[for_anchor_end:]

    if for_anchor_pos < 0:
        raise softwareError("Anchor1 not found")
    
    if rev_anchor_pos < 0:
        raise softwareError("Anchor2 not found")

    switch = 0
    if rev_anchor_pos < for_anchor_pos:
        tmp = for_anchor_pos    
        for_anchor_pos = rev_anchor_pos 
        rev_anchor_pos = tmp

        tmp = for_anchor    
        for_anchor = rev_anchor
        rev_anchor = tmp 
        for_anchor_end = for_anchor_pos + for_anchor.len
        switch = 1

    overlap = 0 if for_anchor_end < rev_anchor_pos else 1
    if overlap:
        raise softwareError("Anchors overlap, adjust Anchors in reference!")

    anchor_dist = rev_anchor_pos - for_anchor_end 
    if anchor_dist % 3 != 0 or anchor_dist < 0:
        raise softwareError("Distance between Anchors is not divisble by 3, adjust Anchors in reference!")

    if anchor_dist < 6:
        warnings.warn("Distance between Anchors should preferably be 6. Will still proceed!")

    return switch, anchor_dist

def print_usage():
    print ("Will Update in Future Version!")

PROG = 'Counting Peptide Inserts (CoPI)'
DESCRIPTION  = f'''
{PROG} quantify Peptide Variants for Protein Engineering. 
Options include to Collapse variants with low counts.
'''

def get_args():
    parser = argparse.ArgumentParser(
        prog = PROG,
        description = DESCRIPTION
    )
    parser.add_argument('-p', '--parameters', dest = 'config',
                        help="Parameters of Analysis")
    parser.add_argument('-d', '--debug',
                        help='Print lots of debugging statements',
                        action="store_const", dest="loglevel",const=logging.DEBUG,
                        default=logging.INFO)
    parser.add_argument('--usage', dest = 'usage', action="store_true", #required=True,
                        help="Print Usage")

    args = parser.parse_args()
    return args