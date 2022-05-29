#!/usr/bin/env python3

import gzip
from itertools import zip_longest
from CoPI.arguments import *
from CoPI.collapse import *
from CoPI.utils import *
from recordclass import recordclass
from collections import defaultdict

import time
logger = logging.getLogger(__name__)
FORMAT = '%(asctime)s - %(module)-16s - %(levelname)s - %(message)s'

class Variant():
    def __init__(self, aa):
        self.aa_seq = aa
        self.aa_len = len(aa)
        self.tot_c = 0
        self.seq_s = set()
        self.seq_C = dict()
        self.seq_Q = dict()

    def __add__(self, v):
        self.tot_c += v.tot_c
        self.seq_s.update(v.seq_s)
        for k1, v1 in v.seq_C.items():
            ## k1: nt_seq_decoded v1: count
            if k1 in self.seq_C.keys():
                tmp_count = self.seq_C[k1] 
                tmp_qual = self.seq_Q[k1] 
                self.seq_Q[k1] = self._updateQual(tmp_qual, v.seq_Q[k1], tmp_count, v1)
                self.seq_C[k1] += v1
            else:
                self.seq_Q[k1] = v.seq_Q[k1]
                self.seq_C[k1] = v1

        tmpVariant = Variant(v.aa_seq) 
        tmpVariant.tot_c = self.tot_c
        tmpVariant.seq_s = self.seq_s
        tmpVariant.seq_C = self.seq_C
        tmpVariant.seq_Q = self.seq_Q

        return tmpVariant


    def _updateQual(self, qual1, qual2, cnt1, cnt2, tab=errs_tab(128)):

        qual1 = [(ord(q) - 33) * cnt1 for q in qual1]
        qual2 = [(ord(q) - 33) * cnt2 for q in qual2]
        total = cnt1 + cnt2
        #qual = [chr(round(-10 * math.log((qual1[i] + qual2[i])/(),10))) for i in range(len(qual1))]
        qual = [chr(math.floor((qual1[i] + qual2[i])/(total)) + 33) for i in range(len(qual1))]
        #print (qual, ave_qual(qual))
        return qual


def __naive_backtrace(B_matrix):
    i, j = len(B_matrix)-1, len(B_matrix[0])-1
    dC = 0
    iC = 0

    while (i, j) != (0, 0):
        if B_matrix[i][j] & 2:
            i, j = i-1, j-1
        elif B_matrix[i][j] & 4:
            i, j = i-1, j
            if i>0 and j>0:
                dC += 1
        elif B_matrix[i][j] & 1:
            i, j = i, j-1
            if i>0 and j>0:
                iC += 1
    
    return iC, dC

def __edit_distance(string_x, string_y, mm=1, get_backtrace=True):
    matrix_a = np.zeros([2, len(string_y) + 1], dtype = int)
    matrix_a[0, :] = np.arange(0,len(string_y)+1,1)
    col = np.arange(0,len(string_x)+1,1)

    if get_backtrace:
            B = [[0 if x > 0 else 4 for x in range(len(string_y) + 1)] for y in range(len(string_x) + 1)]
            B[0] = [0 if x == 0 else 1 for x in range(len(string_y) + 1)]

    for i in range(1,len(col)):
        matrix_a[1, 0] = col[i]
        for j in range(1,len(matrix_a[1,:])):
            if abs(i-j) > mm: 
                matrix_a[1, j] = mm+1
                continue
                
            deletion = matrix_a[0, j] + 1
            insertion = matrix_a[1, j - 1] + 1
            if (string_x[i - 1] == string_y[j - 1]):
                substitution = matrix_a[0, j - 1]
            else:
                substitution = matrix_a[0, j - 1] + 1
                
            score = min(([deletion, insertion, substitution]))  
            if get_backtrace:
                temp = 0
                if score == deletion:
                    temp |= 4
                if score == substitution:
                    temp |= 2
                if score == insertion:
                    temp |= 1

                B[i][j] = temp
            matrix_a[1, j] = score
        matrix_a[0, :] = matrix_a[1, :]

    return __naive_backtrace(B)

def _calc_indel(reference, query):
    '''
    if deletion: returns -offset
    if insertion: returns +offset
    '''
    iC, dC = __edit_distance(reference, query)
    return iC - dC

_hamming_distance = lambda x, y: sum([(x1 != y1) for x1, y1 in zip(x, y)])
def _bitap_search(query_bitap_dict, query, reference, mm = 1, indels = False):
    referenceLen = len(reference)
    alphabet = query_bitap_dict

    queryLen = len(query)
    emptyColumn = (2 << (queryLen - 1)) - 1
    matrix = [[]]

    for k in range(0, mm + 1):
        matrix.append([emptyColumn])
        for columnNum in range(1, referenceLen + 1):
            prevColumn = (matrix[1][columnNum - 1]) >> 1
            letterPattern = alphabet[reference[columnNum - 1]]
            curColumn = prevColumn | letterPattern

            if k >= 1:
                replaceColumn = curColumn & (matrix[0][columnNum - 1] >> 1)
                if indels:
                    insertColumn = curColumn & (matrix[0][columnNum - 1])
                    deleteColumn = curColumn & (matrix[0][columnNum] >> 1)
                    curColumn = insertColumn & deleteColumn & replaceColumn
                else:
                    curColumn = replaceColumn
                    
            matrix[1].append(curColumn)

            if (curColumn & 0x1) == 0:
                offset = 0
                startPos = max(0, columnNum - queryLen) 
                if indels:
                    temp = reference[startPos: columnNum]
                    if _hamming_distance(temp, query) > mm:
                        offset = _calc_indel(temp, query)
                
                return startPos, offset

        matrix = [matrix[1]]

    return -1, 0


def write_to_file_and_print(message, filename):
    filename.write(f"{message}\n")

def is_gzipped(path):
    with open(path, "rb") as f:
        return f.read(2) == b'\x1f\x8b'

Counter = recordclass('Counter', 'readC startC endC missingC indelC editC')
Anchor = namedtuple('Anchor', 'start end')

def check_reverse(para, path, gzipped = True):
    limit = LIMIT
    mismatch = para.tolerate
    threshold = para.threshold

    for_anchor = Anchor(para.seq_start, para.seq_end)
    rev_anchor = Anchor(para.seq_end_rev, para.seq_start_rev)

    for_count = 0
    rev_count = 0

    open_fn = gzip.open if gzipped else open
    with open_fn(path) as fin:
        for i, (name, seq, quality) in enumerate(readfq(fin, gzipped)):
            if i >= limit: 
                fc = for_count 
                rc = rev_count
                if fc > threshold*rc or rc > threshold*fc: break
                limit <<= 1

            seq = seq.upper()
            for_read_hit = (seq.find(for_anchor.start.seq), seq.find(for_anchor.end.seq))
            rev_read_hit = (seq.find(rev_anchor.start.seq), seq.find(rev_anchor.end.seq))

            for j, (read_start, read_end) in enumerate([for_read_hit, rev_read_hit]):
                if read_start >= 0 and read_end >= 0:
                    if read_start < read_end:
                        if j == 0:
                            for_count += 1
                        else:
                            rev_count += 1
                        
                elif read_start >= 0 and read_end < 0:
                    if mismatch: 
                        if j == 0:
                            temp_seq = seq[read_start + rev_anchor.start.len:]
                            read_end, _ = _bitap_search(rev_anchor.end.bitap, rev_anchor.end.seq, temp_seq, mm = mismatch)
                        else:
                            temp_seq = seq[read_start + for_anchor.end.len:]                            
                            read_end, _ = _bitap_search(for_anchor.end.bitap, for_anchor.end.seq, temp_seq, mm = mismatch)
                        
                        if read_end > 0: 
                            if j == 0:
                                for_count += 1
                            else:
                                rev_count += 1

                elif read_start < 0 and read_end >=0:
                    if mismatch: 
                        temp_seq = seq[:read_end]
                        if j == 0:
                            read_start, _ = _bitap_search(rev_anchor.start.bitap, rev_anchor.start.seq, temp_seq, mm = mismatch)
                        else:
                            read_start, _ = _bitap_search(for_anchor.start.bitap, for_anchor.start.seq, temp_seq, mm = mismatch)
                        
                        if read_start > 0: 
                            if j == 0:
                                for_count += 1
                            else:
                                rev_count += 1

    reverse = 0 if for_count > rev_count else 1
    return reverse

def searchUnpaired(para, file, aa_dict, indel_dict, reverse = 0, check = True, number_of_fastq = 2, curr_gz = 0, QC = False):
    gzipped = is_gzipped(file) 
    open_fn = gzip.open if gzipped else open
    mm = para.tolerate
    
    if check: reverse = check_reverse(para, file, gzipped)

    logger.info(f'{file} is reverse: {True if reverse else False}')

    anchor = Anchor(para.seq_end_rev, para.seq_start_rev) if reverse else Anchor(para.seq_start, para.seq_end)
    counter = Counter(0,0,0,0,0,0)
    count = 0

    with open_fn(file) as fin:
        for i, (name, seq, quality) in enumerate(readfq(fin, gzipped)):
            count += 1
            seq = seq.upper()
            read_start, read_end = seq.find(anchor.start.seq), seq.find(anchor.end.seq)
            for_offset = 0
            rev_offset = 0
            if read_start >= 0 and read_end >= 0:
                if read_start > read_end:
                    counter.missingC += 1
                    continue
                counter.readC += 1
                read_end += anchor.end.len
                    
            elif read_start >= 0 and read_end < 0:
                if MISMATCH < 1: 
                    counter.startC += 1 
                    continue
                
                temp_seq = seq[read_start + anchor.start.len:]
                read_end, rev_offset = _bitap_search(anchor.end.bitap, anchor.end.seq, temp_seq, mm=mm)
                if read_end < 0: 
                    counter.startC += 1
                    continue
                counter.readC += 1
                counter.editC += 1
                read_end = read_end + anchor.end.len + len(temp_seq)

            elif read_start < 0 and read_end > 0:
                if MISMATCH < 1: 
                    counter.endC += 1 
                    continue
                
                temp_seq = seq[:read_end]
                read_end += anchor.end.len
                read_start, for_offset = _bitap_search(anchor.start.bitap, anchor.start.seq, temp_seq, mm=mm)
                if read_start < 0: 
                    counter.endC += 1 
                    continue

                counter.readC += 1
                counter.editC += 1
                    
            else:
                counter.missingC += 1 
                continue

            sequence = seq[read_start:read_end]
            mm = read_end - read_start - para.wt_distance
            if mm < 0: 
                counter.missingC += 1
                continue

            if mm > para.indel_range*3: 
                counter.missingC += 1
                continue 

            '''
            #Update Read End -> Assuming Read End may lie in the actual variant
            if mm < para.min_indel_len:
                temp_seq = seq[read_end:]
                tmp_read_end, rev_offset = _bitap_search(anchor.end.bitap, anchor.end.seq, temp_seq, mm=mm)
                if tmp_read_end > read_start:
                    counter.editC += 1
                    read_end = tmp_read_end + anchor.end.len + len(temp_seq)
                    sequence = seq[read_start:read_end]
                    mm = read_end - read_start - para.wt_distance
            '''

            indel_dict[mm] += 1

            if mm > 30: continue
            if mm % 3 != 0: 
                counter.indelC += 1 
                continue

            quality = quality[read_start:read_end]
            if reverse: 
                sequence = reverseC(sequence)
                quality = quality[::-1]

            sequence = sequence[para.seq_start.len: -para.seq_end.len] 
            quality = quality[para.seq_start.len: -para.seq_end.len] 

            mm_aa = mm/3
            aa_sequence = [table[sequence[i:i+3]] for i in range(0, len(sequence), 3)]
            aa_sequence = ''.join(aa_sequence)
            if aa_sequence not in aa_dict[mm_aa].keys():
                aa_dict[mm_aa][aa_sequence] = [Variant(aa_sequence) for x in range(number_of_fastq)]
            
            aa_dict[mm_aa][aa_sequence][curr_gz].tot_c += 1
            nt_seq_decoded = encode_DNA(sequence)
            if nt_seq_decoded in aa_dict[mm_aa][aa_sequence][curr_gz].seq_s:
                tmp_qual = aa_dict[mm_aa][aa_sequence][curr_gz].seq_Q[nt_seq_decoded]
                tmp_count = aa_dict[mm_aa][aa_sequence][curr_gz].seq_C[nt_seq_decoded]
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_Q[nt_seq_decoded] = update_qual(tmp_qual, quality, tmp_count)
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_C[nt_seq_decoded] += 1
            else:
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_s.add(nt_seq_decoded)
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_C[nt_seq_decoded] = 1
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_Q[nt_seq_decoded] = quality


    logger.info(f"{count} reads were processed from {file}")
    filename = open(f'{para.prefix}.summary.log', 'a')
    write_to_file_and_print(f'*****************************************************', filename)
    write_to_file_and_print(f'Summary for {file}', filename)
    write_to_file_and_print(f'Number of Reads with Correct Anchors    : {counter.readC}', filename)
    write_to_file_and_print(f'Number of Reads with MM Anchors         : {counter.editC}', filename)
    write_to_file_and_print(f'Number of Reads with Correct Start Read : {counter.startC}', filename)
    write_to_file_and_print(f'Number of Reads with Correct End Read   : {counter.endC}', filename)
    write_to_file_and_print(f'Number of Reads with 1 Indels Reads     : {counter.indelC}', filename)
    write_to_file_and_print(f'Number of Reads with Wrong Reads        : {counter.missingC}', filename)
    write_to_file_and_print(f'*****************************************************', filename)
    filename.close()

    return aa_dict, indel_dict, reverse

def searchPaired(para, aa_dict, indel_dict, reverse = 0, check = False, number_of_fastq = 2, curr_gz = 1, QC = False):
    gzipped = is_gzipped(para.read1) and is_gzipped(para.read2)
    open_fn = gzip.open if gzipped else open

    reverse = check_reverse(para, para.read1, gzipped) if check else reverse
    logger.info(f'{para.read1} and {para.read2} is reverse: {True if reverse else False}')

    for_anchor = Anchor(para.seq_start, para.seq_end)
    rev_anchor = Anchor(para.seq_end_rev, para.seq_start_rev)
    anchor1 = rev_anchor if reverse else for_anchor
    anchor2 = for_anchor if reverse else rev_anchor

    mm=para.tolerate

    counter = Counter(0,0,0,0,0,0)

    count = 1

    with open_fn(para.read1) as file1, open_fn(para.read2) as file2:
        for line1, line2 in zip_longest(readfq(file1, gzipped), readfq(file2, gzipped)):
            count += 1
            id1, seq1, qual1 = line1
            id2, seq2, qual2 = line2
            assert id1 == id2
            found = 0
            recheck = 0
            for_offset = 0
            rev_offset = 0

            ## Process Read 1
            seq1 = seq1.upper()
            read_start1, read_end1 = seq1.find(anchor1.start.seq), seq1.find(anchor1.end.seq)
            if read_start1 >= 0 and read_end1 >= 0:
                if read_start1 > read_end1:
                    recheck = 1
                else:
                    counter.readC += 1
                    read_end1 += anchor1.end.len

            elif read_start1 >= 0 and read_end1 < 0:
                if MISMATCH < 1: 
                    recheck = 1
                else:
                    temp_seq = seq1[read_start1 + anchor1.start.len:]
                    read_end1, rev_offset = _bitap_search(anchor1.end.bitap, anchor1.end.seq, temp_seq)
                    if read_end1 < 0: 
                        recheck = 1
                    else:
                        counter.readC += 1
                        counter.editC += 1
                        read_end1 = read_end1 + anchor1.end.len + len(temp_seq)

            elif read_start1 < 0 and read_end1 >=0:
                if MISMATCH < 1: 
                    recheck = 1
                else:
                    temp_seq = seq1[:read_end1]
                    read_end1 += anchor1.end.len
                    read_start1, for_offset = _bitap_search(anchor1.start.bitap, anchor1.start.seq, temp_seq)
                    if read_start1 < 0: 
                        recheck = 1
                    else:
                        counter.readC += 1
                        counter.editC += 1
                    
            else:
                recheck = 1

            ## Process Read 2 if required
            if not recheck:
                read_start = read_start1
                read_end = read_end1
                
                #Update Read End -> Assuming Read End may lie in the actual variant 
                '''         
                temp_seq = seq1[read_end:]
                tmp_read_end, rev_offset = _bitap_search(anchor1.end.bitap, anchor1.end.seq, temp_seq, mm=mm)
                if tmp_read_end > 0:
                    read_end = tmp_read_end + anchor1.end.len + len(temp_seq)
                '''

                reverseSeq = 1 if reverse else 0
                sequence = seq1[read_start:read_end]    
                quality = qual1[read_start: read_end]
                paired = 0
                offset = (for_offset, rev_offset)
            else:
                seq2 = seq2.upper()
                read_start2, read_end2 = seq2.find(anchor2.start.seq), seq2.find(anchor2.end.seq)
                if read_start2 >= 0 and read_end2 >= 0:
                    if read_start2 > read_end2:
                        counter.missingC += 1
                        continue
                    counter.readC += 1
                    read_end2 += anchor2.end.len

                elif read_start2 >= 0 and read_end2 < 0:
                    if MISMATCH < 1: 
                        counter.missingC += 1
                        continue

                    temp_seq = seq2[read_start2 + anchor2.start.len:]
                    read_end2, rev_offset = _bitap_search(anchor2.end.bitap, anchor2.end.seq, temp_seq)
                    if read_end2 < 0: 
                        counter.missingC += 1
                        continue
                    
                    counter.readC += 1
                    counter.editC += 1
                    read_end2 = read_end2 + anchor2.end.len + len(temp_seq)

                elif read_start2 < 0 and read_end2 >=0:
                    if MISMATCH < 1: 
                        counter.missingC += 1
                        continue

                    temp_seq = seq2[:read_end2]
                    read_end2 += anchor2.end.len
                    read_start2, for_offset = _bitap_search(anchor2.start.bitap, anchor2.start.seq, temp_seq)
                    if read_start2 < 0: 
                        counter.missingC += 1
                        continue
                    
                    counter.readC += 1
                    counter.editC += 1
                        
                else:
                    counter.missingC += 1
                    continue

                read_start = read_start2
                read_end = read_end2

                '''
                temp_seq = seq2[read_end:]
                tmp_read_end, rev_offset = _bitap_search(anchor2.end.bitap, anchor2.end.seq, temp_seq, mm=mm)
                if tmp_read_end > 0:
                    read_end = tmp_read_end + anchor2.end.len + len(temp_seq)
                '''

                quality = qual2[read_start:read_end]
                reverseSeq = 0 if reverse else 1
                sequence = seq2[read_start:read_end]    
                paired = 1
                offset = (for_offset, rev_offset)
    
            mm = read_end - read_start - para.wt_distance
            if mm < 0: 
                counter.missingC += 1
                continue
            
            if mm > para.indel_range*3: 
                counter.missingC += 1
                continue 
            
            indel_dict[mm] += 1

            if mm > 30: continue
            if mm % 3 != 0: 
                counter.indelC += 1 
                continue

            if reverseSeq: 
                sequence = reverseC(sequence)
                quality = quality[::-1]

            sequence = sequence[para.seq_start.len: -para.seq_end.len] 
            quality = quality[para.seq_start.len: -para.seq_end.len] 
            
            mm_aa = mm/3
            aa_sequence = [table[sequence[i:i+3]] for i in range(0, len(sequence), 3)]
            aa_sequence = ''.join(aa_sequence)
            if aa_sequence not in aa_dict[mm_aa].keys():
                if paired:
                    aa_dict[mm_aa][aa_sequence] = [Variant(aa_sequence) for x in range(number_of_fastq)]
                else:
                    aa_dict[mm_aa][aa_sequence] = [Variant(aa_sequence) for x in range(number_of_fastq)]

            aa_dict[mm_aa][aa_sequence][curr_gz].tot_c += 1
            nt_seq_decoded = encode_DNA(sequence)

            if nt_seq_decoded in aa_dict[mm_aa][aa_sequence][curr_gz].seq_s:
                tmp_qual = aa_dict[mm_aa][aa_sequence][curr_gz].seq_Q[nt_seq_decoded]
                tmp_count = aa_dict[mm_aa][aa_sequence][curr_gz].seq_C[nt_seq_decoded]
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_Q[nt_seq_decoded] = update_qual(tmp_qual, quality, tmp_count)
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_C[nt_seq_decoded] += 1
            else:
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_s.add(nt_seq_decoded)
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_C[nt_seq_decoded] = 1
                aa_dict[mm_aa][aa_sequence][curr_gz].seq_Q[nt_seq_decoded] = quality

    logger.info(f"{count} reads were processed from {para.read1} and {para.read2}")
    filename = open(f'{para.prefix}.summary.log', 'a')
    write_to_file_and_print(f'*****************************************************', filename)
    write_to_file_and_print(f'Summary for {para.read1} and {para.read2}', filename)
    write_to_file_and_print(f'Number of Reads with Correct Anchors    : {counter.readC}', filename)
    write_to_file_and_print(f'Number of Reads with MM Anchors         : {counter.editC}', filename)
    write_to_file_and_print(f'Number of Reads with 1 Indels Reads     : {counter.indelC}', filename)
    write_to_file_and_print(f'Number of Reads with Wrong Reads        : {counter.missingC}', filename)
    write_to_file_and_print(f'*****************************************************', filename)
    filename.close()

    return aa_dict, indel_dict, reverse

def create_file(prefix, placeholder):
    filename = open(f'{prefix}.{placeholder}.log', 'w')
    filename.close()

def write_output(para, aa_dict, number_of_fastq, collapsed = False):
    edit = 'post' if collapsed else 'pre'
    number_of_fastq = 1 if collapsed else number_of_fastq
    all_summary = open(f"{para.prefix}.{edit}collapse.all.log", 'w')
    ind_summary = open(f"{para.prefix}.{edit}collapse.stop.log", 'w')

    total = {}
    for k, v in aa_dict.items():
        curr_total = [0 for i in range(number_of_fastq)]
        if collapsed:
            v = {kx: vx for kx, vx in sorted(v.items(), key=lambda item: item[1].tot_c, reverse=True)}
        else:
            v = {kx: vx for kx, vx in sorted(v.items(), key=lambda item: sum([i.tot_c for i in item[1]]), reverse=True)}
        for k1, v1 in v.items():
            
            if not collapsed: 
                temp_total = sum([i.tot_c for i in v1])

            if collapsed: 
                all_summary.write(f"{k1}:{k}:{v1.tot_c}\n")
            else:
                all_summary.write(f"{k1}:{k}:{','.join([str(i.tot_c) for i in v1])},{temp_total}\n")
        
            if not '_' in k1:
                if collapsed:
                    ind_summary.write(f"{k1}:{k}:{v1.tot_c}\n")
                else:
                    ind_summary.write(f"{k1}:{k}:{','.join([str(i.tot_c) for i in v1])},{temp_total}\n")
            
            if not collapsed:
                for i in range(len(v1)):
                    curr_total[i] += v1[i].tot_c
            else:
                curr_total[0] = v1.tot_c

        all_summary.write('-----\n')
        ind_summary.write('-----\n')
        total[k] = curr_total

    all_summary.close() 
    ind_summary.close()

    if collapsed: return

    with open(f'{para.prefix}.indel-count.log', 'w') as f: 
        f.write(f"Mismatch,{','.join([str(i) for i in list(range(number_of_fastq))])}")
        for k, v in total.items():
            f.write(f"{k},{','.join([str(i) for i in v])}\n")

def _pre_sequence_logo(para, aa_dict, thres = 50, collapsed = True):
    aa = list(set(list(table.values())))
    aa.sort()
    for k, v in aa_dict.items():
        if k <= 0: continue
        if len(v) < 50: continue
        if collapsed:
            v = {kx: vx for kx, vx in v.items() if vx.tot_c > thres}
        else:
            v = {kx: sum([i.tot_c for i in vx]) for kx, vx in v.items() if sum([i.tot_c for i in vx]) > thres}


        logodata = LogoData.from_counts(v)
        
        '''
        temp = {b:{a:0 for a in aa} for b in range(para.indel_range + 2)}
        for c, (seq, v1) in enumerate(v.items()):
            temp_total = v1.tot_c if collapsed else sum([i.tot_c for i in v1])
            length = len(seq)
            for i in range(0, length):
                temp[i][seq[i]] += temp_total
            
        temp = {k:v for k,v in temp.items() if sum(v.values()) > 0}
        temp = np.array([[temp[pos][aa] for key in aa] for pos in range(length)]).T
        '''
        __sequence_plot(para.prefix, temp, length)


def __plot_heatmap(para, pos_dict, aa, length):
    quantity_2D = []
    for k, v in pos_dict.items():
        quantity_2D.append(list(v.values())) 

    quantity_2D = np.array(quantity_2D).T
    log_norm = LogNorm(vmin=quantity_2D.min().min(), vmax=quantity_2D.max().max())
    #cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(quantity_2D.min().min())), 1+math.ceil(math.log10(quantity_2D.max().max())))]
    import matplotlib.colors as mcolors
    plt.switch_backend('agg')
    plt.figure(figsize = (12, 6))
    plt.imshow(quantity_2D, 
               norm=mcolors.LogNorm(), 
               interpolation='nearest', )
    #heatmap = plt.pcolor(quantity_2D)
    aa[-1] = aa[-1].replace('_','*')
    plt.yticks([i for i in range(len(aa))], aa)
    ins = [para.label1] + [f'Ins{i+1}' for i in range(length)] + [para.label2]
    plt.xticks([i + 0.5 for i in range(len(ins))], ins, rotation=90, ha='right')
    plt.colorbar()
    plt.ylabel('Amino Acids')
    plt.xlabel('Peptide Insert Position (nt)')
    plt.title(f'Distribution of in-frame peptide between {para.label1} & {para.label2}')
    plt.savefig(f'{para.prefix}.HeatMap-{length}.png')

def _pre_heatmap(para, aa_dict, collapsed = True):
    #edit This
    aa = list(set(list(table.values())))
    aa.sort()
    for k, v in aa_dict.items():
        if k <= 0: continue
        if len(v) < 50: continue
        temp = {b:{a:0 for a in aa} for b in range(para.indel_range + 2)}
        for seq, v1 in v.items():
            temp_total = v1.tot_c if collapsed else sum([i.tot_c for i in v1])
            length = len(seq)
            for i in range(0, length):
                temp[i][seq[i]] += temp_total
            
        temp = {k:v for k,v in temp.items() if sum(v.values()) > 0}
        diff = para.anchor_dist // 3
        __plot_heatmap(para, temp, aa, length-diff)

def _pre_collapse(aa_list):
    ## Collapse the Dict into a single element array
    tmpVariant = Variant(aa_list[0].aa_seq)
    for v1 in aa_list:
        tmpVariant += v1
    return tmpVariant

def _plot_histo_dist(para, indel_dict):
    indel = [v for k,v in indel_dict.items()]
    total = sum(indel)
    RANGE = list(range(para.indel_range * 3))
    indel_normalized = [round(v/total,5) for v in indel]

    plt.switch_backend('agg')
    plt.figure(figsize = (12, 6))
    plt.bar(RANGE, indel_normalized, color ='maroon')
    plt.xticks(np.arange(min(RANGE), max(RANGE)+1, 1.0))
    plt.xlabel('Insertion Length (nt)')
    plt.ylabel('Probability Density (Normalized Count)')
    plt.title(f'Distribution of Insertion Length Between {para.label1} and {para.label2}')
    plt.savefig(f'{para.prefix}.Normalized_Count.png')

def searchVariant(para):
    create_file(para.prefix, 'summary')
    para._calculate_number_fastq() 

    aa_dict    = {i:{} for i in range(0, INDEL_RANGE)}
    indel_dict = {i:0  for i in range(0, INDEL_RANGE*3)}
    checked    = False
    reverse    = 0
    count      = 0

    if para.read1 and para.read2:
        logger.info("Reading Paired Reads")
        aa_dict, indel_dict, reverse = searchPaired(para, aa_dict, indel_dict, 
                                                    reverse = reverse, check = not checked,
                                                    number_of_fastq = para.fastqC, curr_gz = count, QC = not para.merge)
        checked = True
        count += 1

    if para.readM:
        logger.info("Reading Merged Reads")
        aa_dict, indel_dict, reverse = searchUnpaired(para, para.readM, aa_dict, indel_dict, 
                                                      reverse = reverse, check = not checked, 
                                                      number_of_fastq = para.fastqC, curr_gz = count, QC = not para.merge)
        checked = True
        count += 1
        
    if para.readU:
        for u in para.readU:
            logger.info(f"Reading {u}")
            aa_dict, indel_dict, _ = searchUnpaired(para, u, aa_dict, indel_dict, 
                                                    reverse = reverse, check = True, 
                                                    number_of_fastq = para.fastqC, curr_gz = count, QC = not para.merge)
            count += 1

    logger.info(f"Number of Sets Analysed: {count}")
    aa_dict = {k:v for k,v in aa_dict.items() if v}
    write_output(para, aa_dict, count, collapsed = False)
    if para.plot:
        _plot_histo_dist(para, indel_dict)
        del indel_dict

    if para.collapse: 
        create_file(para.prefix, 'anomaly')
        for k, v in aa_dict.items():
            logger.info(f"Collapsing {k}th insertion!")
            # k: mm_aa
            # v: aa_sequence->dict of List
            hcDict = {}
            lcList = []

            if len(v) <= 10:
                for aa_sequence, value in v.items():
                    value = _pre_collapse(value)
                    hcDict[aa_sequence] = value

                aa_dict[k] = hcDict 
                continue

            for aa_sequence, value in v.items():
                value = _pre_collapse(value)
                total = value.tot_c
                if total > para.threshold:
                    hcDict[aa_sequence] = value

                if total <= para.threshold:
                    lcList.append(value)
            
            logger.debug(f"Number of Low Count Variants: {len(lcList)}")
            updateDct = multiCollapse(para, hcDict, lcList)

            for ku, vu in updateDct.items():
                hcDict[ku].tot_c += vu

            aa_dict[k] = hcDict

        write_output(para, aa_dict, 1, collapsed = True)

    if para.plot:
        _pre_heatmap(para, aa_dict, collapsed = para.collapse)
        #_pre_sequence_logo(para, aa_dict, collapsed = para.collapse)

    del aa_dict
    return

def main():
    args = get_args()
    logging.basicConfig(level=args.loglevel, format=FORMAT)
    if args.usage or not args.config:
        print_usage()
        return
    para = get_toml(args, logger=logger)
   
    logger.info(f"Analysis Started!")
    logger.info(f"Seat back and get your C(K)opi ready!")

    if para.merge: 
        logger.info(f"fastp Analysis!")
        success, gzipped = execute_fastp(para)
        if not success:
            para.merge = False
        else:
            tmp = '.gz' if gzipped else ''
            para.readM = f'{para.prefix}.merged.fastq{tmp}'
    
    start_time = time.time()
    searchVariant(para)

    end_time = time.time()
    total_time = round(end_time - start_time,3)
    logger.info(f"Analysis Ended! in {total_time} seconds!")


if __name__== "__main__":
    main()


'''
def get_parameters():
    ID = 'SCRAMBLER'
    ref_seq = 'GGACATTGAAAAGGTCATGATTACAGACGAAGAGGAAATCAGGACAACCAATCCCGTGGCTACGGAGCAGTATGGTTCTGTAtctACCAACCTCCAGAGAGGCAACAGACAAGCAGCTACCGCAGATGTCAACACACAAGGCGTTCTTCCAGGCATGGTCTGG'
    seq_start = "CAGAGAGGC" #'AGAGCAGCAGCACA' # 14 nt after the first nt.
    seq_end = "CAAGCAGCT" #'GACCCTGCGACCGG' # last 15 nt
    interest = "TEST" #'Heart'#
    read1 = "AmineAAV1-upaired-1.fastq.gz"
    read2 = "AmineAAV1-upaired-2.fastq.gz"
    merged = "AmineAAV1-merged.fastq.gz"

    return Parameters(ref_seq, seq_start,seq_end, read1, read2, merge = False, readMerged = merged)


def get_toml(args):
    pass

def fastp_process(para):
    pass

DESCRIPTION  = 

def get_args():
    parser = argparse.ArgumentParser(
        prog = 'Capture AAV Potential Inserts (CAPI)',
        description = DESCRIPTION
    )
    parser.add_argument('-t', '--toml', dest = 'toml', #required=True,
                        help="Parameters of Analysis in TOML format")
    parser.add_argument('-d', '--debug',
                        help='Print lots of debugging statements',
                        action="store_const",dest="loglevel",const=logging.DEBUG,
                        default=logging.INFO)

    args = parser.parse_args()
    return args

def main():
    args = get_args()
    para = get_parameters() #get_toml(args)
    start_time = time.time()
    logging.basicConfig(level=args.loglevel, format=FORMAT)
    logger.info(f"Analysis Started!")
    logger.info(f"Have you drank your C(K)opi?")

    if para.merge: fastp_process(para)
    searchVariant(para)

    end_time = time.time()
    total_time = round(end_time - start_time,3)
    logger.info(f"Analysis Ended! in {total_time} seconds!")
'''