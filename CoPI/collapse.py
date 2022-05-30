#!/usr/bin/env python3

from bitarray import bitarray
import numpy as np
from CoPI.utils import *
import multiprocessing as mp
import sys

KLEN = 6

import ctypes
alive = mp.Value(ctypes.c_bool, lock=False)
alive.value = True

def create4MerList():
    import itertools
    perm3mer = [''.join(p) for p in itertools.product(DECODING_LST, repeat=KLEN)]
    perm3mer.sort()
    return perm3mer

aa4Lib = create4MerList()
aa4Lib_len = len(aa4Lib)

##-------------------------------------------------------------------------------------------------
## Collapse Modules

def _LFQuery(q, bitArr, length, seeds, mm = 1):
    for i in range(seeds):
        if i>mm: break
        ia = i
        ia *= KLEN
        if (ia+KLEN)>length: 
            ia = length - KLEN

        tmp = q[ia:ia+KLEN]
        index = aa4Lib.index(tmp)
        if bitArr[i][index]: 
            return 0

    return 1

translate = lambda seq: ''.join([table[seq[i:i+3]] for i in range(0, len(seq), 3)])
def _write_anomaly(prefix, lock, aa_seq, tot_c, reason):
    aa_seq = translate(aa_seq)
    with lock:
        with open(f"{prefix}.anomaly.log", 'a') as f:
            f.write(f"{aa_seq}:{tot_c},{reason}\n")

def _check_hits(hits, lc, thres = 5):
    tmp = {}
    for k, v in hits.items():
        if v.seq_C[encode_DNA(k)] > lc*thres:
            tmp[k] = v
    return tmp

def __chkDup(tmp_lst):
    tmp_seq = set()
    for e in tmp_lst:
        if e.p in tmp_seq:
            return True
        else:
            tmp_seq.add(e.p)         
    return False

hamming_distance = lambda x, y: sum([(x1 != y1) for x1, y1 in zip(x, y)])
hamming_position = lambda x, y: [i for i in range(len(x)) if x[i] != y[i]]
def _compare_hits(hits, lc_seq, lc_qual):
    hit = ''
    tmp = []
    cache = namedtuple('cache', 's p q c ql')
    for k, v in hits.items():
        pos = hamming_position(k, lc_seq)
        if len(pos) > 1: continue
        hc_qual = v.seq_Q[encode_DNA(k)]
        tmp_cache = cache(k, pos[0], abs(ord(hc_qual[pos[0]]) - ord(lc_qual[pos[0]])), v.seq_C[encode_DNA(k)], hc_qual)
        tmp.append(tmp_cache)

    if __chkDup(tmp):
        # Case 1: mm at same pos
        if all(e.p == tmp[0].p for e in tmp):
            #print ("Case 1")
            aa = translate(tmp[0].s)
            same = all(aa == translate(e.s) for e in tmp)
            if same: 
                hit = aa

        # Case 2: mm at same pos AND different pos
        else:
            #print ("Case 2")
            length = len(lc_seq)
            tmp_grp = [[]] * length

            for i in range(length):
                tmp_grp[tmp[i].p].append(tmp[i])

            tmp2 = []
            for i in tmp_grp:
                if not i: continue  
                if len(i)>1:
                    aa = translate(i[0].s)
                    same = all(aa == translate(e.s) for e in i)
                    
                    if same: 
                        qual = [qual_set(t.ql, t.c) for t in tmp_grp]
                        qual, qual_sum = update_qual2(qual)
                        qual_dif = qual[i.p] - lc_qual[i.p]
                        new = cache(i.s, i.p, qual_dif, qual_sum, 0)
                        tmp2.append()
                    else:
                        #Undecisive
                        return None
                else:
                    tmp2.append(i)

            min_d = float('inf')
            for e in tmp2:
                if min_d > e.q:
                    min_d = e.q
                    hit = e.s
            hit = translate(hit)

    # Case 3: mm at different pos
    # Lowest Difference (Most likely PCR error, higher confidence in sequencing)
    else:
        #print ("Case 3")
        min_d = float('inf')
        for e in tmp:
            if min_d > e.q:
                min_d = e.q
                hit = e.s
        hit = translate(hit)

    return hit

def _add_hit(tmpDict, hit, count, t = True):
    hit = translate(hit) if t else hit
    if hit in tmpDict.keys():
        tmpDict[hit] += count
    else:
        tmpDict[hit] = count
    return tmpDict

aa_hc = namedtuple('aa_hc', 'seq count')
def processCollapse(para, hcDict, lcList, queue, lock, logger, mm = 1):

    def check_alive(queue):
        queued = False
        while not queued:
            try:
                queued = True
            except:
                pass
            par_alive = mp.parent_process().is_alive()
            if not (par_alive and alive.value):
                alive.value = False
                queue.close()
                sys.exit(1)

    length = lcList[0].aa_len*3
    updateDict = {}
    hcLst = []
    for k, v in hcDict.items():
        # k: aa_seq
        # v: Variant
        hcLst += [decode_DNA(i, length) for i in list(v.seq_s)]

    hcBWA = BWA(''.join(hcLst), length)
    if length > KLEN:
        hcBWA.losslessFiltration(hcLst)

    seedsC = (length // KLEN) + 1
    aa_len = length//3
    for checking, qV in enumerate(lcList):
        check_alive(queue)
        if checking%50 == 0: logger.debug(f"Checking iteration {checking}")
        for q, v in qV.seq_C.items():
            q = decode_DNA(q, length)
            if length > KLEN:
                if _LFQuery(q, hcBWA.bitarr, length, seedsC, mm): 
                    _write_anomaly(para.prefix, lock, q, v, "Not Found (Lossless)")
                    continue
            
            hits = hcBWA.find_match(q, mm)
            hits = [a for a in hits if a % length == 0]
            if not hits: 
                # more than 2 mm
                _write_anomaly(para.prefix, lock, q, v, "Not Found (Aligner)")
                continue

            hits = [hcBWA.extract_sequence(length, a) for a in hits]
            tmp_dict = {h:hcDict[translate(h)] for h in hits}
            tmp_dict = _check_hits(tmp_dict, v)
            if not tmp_dict: 
                _write_anomaly(para.prefix, lock, q, v, 'Lower than threshold')
                continue

            if len(tmp_dict) == 1:
                hit = list(tmp_dict.keys())
                updateDict = _add_hit(updateDict, hit[0], v)
            else:
                hit = _compare_hits(tmp_dict, q, qV.seq_Q[encode_DNA(q)])
                if hit:
                    updateDict = _add_hit(updateDict, hit, v, False)
                else: 
                    _write_anomaly(para.prefix, lock, q, v, "Unsure")

    logger.debug("Completed")
    queue.put([aa_hc(k, v) for k,v in updateDict.items()])          

def queueCollapse(queue, hcList):
    tmp = {}
    while True:
        item = queue.get()
        if item == 'Done':
            hcList += [aa_hc(k, v) for k,v in tmp.items()]
            return
        else:
            for a in item:
                _add_hit(tmp, a.seq, a.count, t = False)

##-------------------------------------------------------------------------------------------------
## Multiprocessing Utils

def partition(lst, n):
    division = len(lst) / n
    return [lst[round(division * i):round(division * (i + 1))] for i in range(n)]

def ranges(N, nb):
    temp_list = [(r.start, r.stop) for r in partition(range(N), nb)]
    for i in temp_list: 
        yield (i[0],i[1])

def multiCollapse(para, hcDict, lcList, logger):
    queue = mp.Queue()
    manager = mp.Manager()
    lock = mp.Lock()
    hcList = manager.list()
    current_processes = []

    for no, i in enumerate(ranges(len(lcList), para.threads)):
        lcListPartitioned = lcList[i[0]:i[1]]
        #para, length, hcDict, lcList, queue, lock
        logger.debug(f"CPU {no} will process {len(lcListPartitioned)} variants")
        pp = mp.Process(target = processCollapse, args = (para, hcDict, lcListPartitioned, queue, lock, logger))
        current_processes.append(pp)

    for pp in current_processes:
        pp.start()

    del lcListPartitioned
    del lcList
    sink = mp.Process(target=queueCollapse, args=(queue, hcList))
    sink.start()

    for pp in current_processes:
        pp.join()

    queue.put('Done')
    sink.join()

    logger.debug("Multi Completed")
    return {a.seq:a.count for a in list(hcList)}

##-------------------------------------------------------------------------------------------------
## DA for BWA & WT & LF (Lossless Filtration)

class Node:

    def __init__(self, value = []):
        self.left  = None #0 
        self.right = None #1
        self.value = bitarray(''.join([str(i) for i in value]))

class Wavelet_Tree():

    def __init__(self, sequence, alphabet):
        self.root = None
        self.length = len(sequence)
        self.alphabet = alphabet + ['$']
        self.CharBit = {k:None for k in self.alphabet}
        
        x = {k:sequence.count(k) for k in self.alphabet}
        self.nt_count = {k:v for k,v in x.items() if v} 
        self.ConstructNode(sequence, keys = self.alphabet)
        #self.printTree()

    def ConstructNode(self, string, node = None, keys = None, charbit = '', root = True):
        if len(keys) == 1:
            self.CharBit[keys[0]] = charbit
            return 

        charbitL = charbitR = charbit
        charbitL += '0'
        charbitR += '1'
        Bitmap = []
        x = {}

        for char in keys:
            if char not in list(self.nt_count.keys()): continue
            x[char] = self.nt_count[char]

        keys.sort(reverse = True)
        #keys = alphabet
        mid = len(keys)//2
        SplitCharL = keys[:mid]
        SplitCharR = keys[mid:]
        Sleft = Sright = ''

        for c in list(string):
            if c in SplitCharL:
                Sleft += c
                Bitmap.append(0)
            else:
                Sright += c
                Bitmap.append(1)

        newNode = Node(Bitmap)
        if not self.root:
            self.root = newNode

        if charbit: 
            current_charbit = charbit[-1]
            if current_charbit == '0':
                node.left = newNode
            elif current_charbit == '1':
                node.right = newNode

        self.ConstructNode(Sleft, newNode, SplitCharL, charbitL, root = False)
        self.ConstructNode(Sright, newNode, SplitCharR, charbitR, root = False)
        if not root: return
        
        self.CharBit = {k:v for k,v in self.CharBit.items() if v}
        self.CharBit_inv = {v:k for k, v in self.CharBit.items()}

    def printTree(self):
        if self.root is not None:
            self._printTree(self.root)

    def _printTree(self, node, left = True):
        if node is not None:
            self._printTree(node.left, left)
            print (left, node.value.__sizeof__(), str(node.value) + ' ')
            self._printTree(node.right, False)

    def _BinaryRank(self, c, p, node):
        return node.value[:p].count(c)

    def Rank(self, c, p):
        ''' 
        c = Character
        p = position/ offset
        '''
        if p < 0:
            return 0
        charbit = list(self.CharBit[c])
        max_counter = len(charbit)
        i = 0
        node = None
        while True: 
            if not node:
                node = self.root
            else:
                if current_charbit == 0:
                    node = node.left
                elif current_charbit == 1:
                    node = node.right
            
            current_charbit = int(charbit[i])
            p = self._BinaryRank(current_charbit, p, node)
            i += 1
            if i == len(charbit):
                break

        return p

    def _BinaryAccess(self, node, o):
        return node.value[o]

    def Access(self, offset):
        charbit = ''
        node = None
        while True:
            if not node:
                node = self.root
            o = int(self._BinaryAccess(node, offset))
            charbit += str(o)

            if charbit in list(self.CharBit_inv.keys()):
                return self.CharBit_inv[charbit]

            r = int(self._BinaryRank(o, offset, node))
            offset = r 
            if o == 0:
                node = node.left
            elif o == 1:
                node = node.right
    
    def Select(self):
        pass

class _Suffix:
    def __init__(self, sequence, position):
        self.seq = sequence
        self.pos = position

    @staticmethod
    def textKey(a):
        return a.seq

class _create_BWA:
    def __init__(self, reference, aa_len):
        reference = reference.upper()
        rotation_list, rotation_list_reverse, suffix_array, bwt, rbwt = [list() for i in range(5)]
        self.ALPHABET = list(set(reference))

        C = {k:1 for k in self.ALPHABET}
        C['$'] = 0
        reference = "%s$" % reference
        reverse_reference = reference[::-1] 

        for i in range(len(reference)):
            new_rotation = "%s%s" % (reference[i:],reference[0:i])
            struct = _Suffix(new_rotation,i)
            rotation_list.append(struct)

            new_rotation_reverse = "%s%s" % (reverse_reference[i:],reverse_reference[0:i])
            struct_rev = _Suffix(new_rotation_reverse,i)
            rotation_list_reverse.append(struct_rev)

            if reference[i]!='$':
                for char in self.ALPHABET:
                    if reference[i] < char:
                        C[char] = C[char] + 1   

        rotation_list.sort(key=_Suffix.textKey)
        rotation_list_reverse.sort(key=_Suffix.textKey)

        for i in rotation_list:
            suffix_array.append(i.pos)
            bwt.append(i.seq[-1:])

        for i in rotation_list_reverse:
            rbwt.append(i.seq[-1:])

        ## Use Wavelet Tree to Rank instead of OCC DS
        self.non_wave = bwt
        self.bwt = Wavelet_Tree(bwt, self.ALPHABET)
        self.rbwt = Wavelet_Tree(rbwt, self.ALPHABET)
        self.n = len(reference)
        self.sampled = 32
        leftover = self.n - 1 - ((self.n - 1) % self.sampled)
        suffix_array = [i if i%self.sampled == 0 or i in range(leftover, self.n) else float('inf') for i in suffix_array]
        
        self.SA_bit = bitarray([1 if a <= self.n else 0 for a in suffix_array])
        self.SA = [a for a in suffix_array if a <= self.n]
        self.reverse_state = True
        self.C = C
        self.ALPHABET = list(self.C.keys()) 

class BWA():
    '''
    A class to process Queries
    '''
    def __init__(self, reference, aa_len):
        bwa = _create_BWA(reference, aa_len)
        
        self.C = bwa.C
        self.ALPHABET = bwa.ALPHABET
        self.bwt = bwa.bwt
        self.rbwt = bwa.rbwt
        self.n = bwa.n

        self.sampled = bwa.sampled
        self.SA_bit = bwa.SA_bit
        self.SA = bwa.SA
        
        self.D = list() 
        self.reverse_state = True

        self.aa_len = aa_len
        self.bitarr = None

    def losslessFiltration(self, hcLst):
        #Assuming Similar Length ie. 7
        aa_len = len(hcLst[0])
        seeds = (aa_len // KLEN) + 1
        bitArr = [bitarray(aa4Lib_len) for i in range(seeds)]
        for b in bitArr:
            b.setall(0)

        for r in hcLst:
            for i in range(seeds):
                ia = i
                ia *= KLEN
                if (ia+KLEN)>aa_len: 
                    ia = aa_len - KLEN

                tmp = r[ia:ia+KLEN]
                index = aa4Lib.index(tmp)
                bitArr[i][index] = 1

        self.bitarr = bitArr
        
    #exact matching - no indels or mismatches allowed
    def exact_match(self, query):
        query = query.upper()
        top = 0
        bot = self.n - 1
        i = len(query) - 1

        while (i>=0 and bot > top):
            c = query[i]
            top = self.C[c] + self.bwt.Rank(c, top)
            bot = self.C[c] + self.bwt.Rank(c, bot)
            i -= 1

        if i < 0: 
            return self.fill_SA(top, bot) 
        else:
            return []

    def find_SA_pos(self, pos):
        '''
        pos - pos with reference to original sequence
        '''
        if pos in self.SA:
            return self.SA.index(pos)

        ## closest_value with reference to SA 
        closest_value = min(self.SA, key=lambda x: abs(x-pos) if x else float('inf'))
        if pos > closest_value:
            closest_value = min([self.n - 1, closest_value + self.sampled])
            if closest_value != (self.n - 1):
                while (closest_value % self.sampled != 0):
                    closest_value -= 1
            
        ## pos - index of element in bwt bit
        while not pos in self.SA:
            temp_pos = self.SA.index(closest_value)
            c = self.bwt.Access(temp_pos)
            temp = self.C[c] + self.bwt.Rank(c, temp_pos)
            if not self.SA[temp]:
                self.SA[temp] = closest_value - 1
            closest_value -= 1

        return self.SA.index(pos)

    def extract_sequence(self, q_len, pos):
        seq = ''
        pos += q_len
        pos = self.find_SA_pos(pos)

        while (len(seq) != q_len):
            c = self.bwt.Access(pos)
            pos = self.C[c] + self.bwt.Rank(c, pos)
            seq += c 

        return seq[::-1]

    def _fill_SA(self, pos):
        c = self.bwt.Access(pos)
        temp = self.C[c] + self.bwt.Rank(c, pos)

        if not self.SA[temp]:
            self._fill_SA(temp)
        
        if self.SA[temp] >= 0:
            self.SA[pos] = self.SA[temp] + 1
            if self.SA[pos] >= self.n:
                self.SA[pos] -= self.n


    def fill_SA(self, start, stop = None):
        if not self.SA[start]:
            self._fill_SA(start)
        
        if not stop: 
            return self.SA[start]
                
        if stop and not self.SA[stop]:
            self._fill_SA(stop)

        if self.SA[start] and self.SA[stop]:
            return self.SA[start:stop]

    def _reverse_SA(self):
        temp = []
        SA = self.SA[::-1]
        for bit in self.SA_bit:
            if int(bit) == 0:
                temp.append(None)
            elif int(bit) == 1:
                curr = SA.pop()
                temp.append(curr)
        self.SA = temp

    # get the position(s) of the query in the reference
    def find_match(self, query, mismatch):
        if self.reverse_state:
            self._reverse_SA()
            self.reverse_state = False
        
        query = query.upper()
        if mismatch: return self.inexact_match(query, mismatch)
        else: return self.exact_match(query) 

    ## inexact matching, z is the max threshold for allowed edits
    def inexact_match(self,query,z):
        self.calculate_d(query)
        SA_indices = self.inexRecur(query, len(query)-1, z, 0, self.n-1)
        
        ## Return the values in the SA
        return [self.fill_SA(x) for x in SA_indices]

    ## recursion function that effectively "walks" through the suffix tree using the SA, BWT, Occ and C datastructures
    def inexRecur(self, query, i, z, top, bot):
        tempset = set()
            
        ## 2 stop conditions, one when too many differences have been encountered, another when the entire query has been matched, terminating in success
        ## Reached the limit of differences at this stage, terminate this path traversal
        if z < self.get_D(i):
            ## Return empty set   
            return set()
        
        ## Empty query string, entire query has been matched, return SA indexes k:l
        if i < 0:
            for m in range(top,bot): tempset.add(m)
            return tempset
            
        result = set()
        ## For each character in the alphabet
        for char in self.ALPHABET:
            ## Find the SA interval for the char
            newK = self.C[char] + self.bwt.Rank(char, top)
            newL = self.C[char] + self.bwt.Rank(char, bot)
            ## If the substring was found
            if newK <= newL:
                ## If the char was correctly aligned, then continue without decrementing z (differences)
                if char == query[i]: 
                    result = result.union(self.inexRecur(query,i-1,z,newK,newL))
                ## Continue but decrement z, to indicate that this was a difference/unalignment
                else: 
                    result = result.union(self.inexRecur(query,i-1,z-1,newK,newL))
        return result

    #calculates the D array for a query, used to prune the tree walk and increase speed for inexact searching
    def calculate_d(self, query):
        top = 0
        bot = self.n - 1
        z = 0

        ## Empty the D array
        self.D = list()
        for i in range(len(query)):
            c = query[i]
            top = self.C[c] + self.rbwt.Rank(c,top) 
            bot = self.C[c] + self.rbwt.Rank(c,bot + 1) - 1
            
            ## If this character has NOT been found
            if top > bot:
                top = 0
                bot = self.n - 1
                z = z + 1
            self.D.append(z)

    ## Gets values from the D array
    ## NOTE D(-1) = 0
    def get_D(self,index):
        if index < 0:
            return 0
        else:
            return self.D[index]