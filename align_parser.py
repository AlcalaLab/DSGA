#!/usr/bin/env python3

from Bio import SeqIO
from operator import itemgetter
from itertools import groupby


def extract_gap_pos(soma_seq):
    gap_pos = [i for i, j in enumerate(soma_seq) if j == '-']

    gap_ranges = []

    for k, g in groupby(enumerate(gap_pos), lambda x:x[0]-x[1]):
        group = list(map(itemgetter(1), g))

        if len(group) > 1:
            gap_ranges.append((group[0],group[-1]))

        else:
            gap_ranges.append((group[0],group[0]))

    return gap_ranges[1:-1]


def check_shared_pointer(germ_pntr, soma_pntr):
    shared_ntds = ''
    for i in range(len(germ_pntr)):
        if germ_pntr[i] == soma_pntr[i]:
            shared_ntds += germ_pntr[i]
        else:
            break
    return len(shared_ntds)


def check_intron_bndry(germ_seq: str, gap_start: int, gap_end: int, max_shift: int = 5):
    intron_pos = None

    intron_len = 1 + gap_end - gap_start

    # Check shifts from -max_shift to +max_shift
    for pos in range(-max_shift, max_shift + 1):
        intron_st = rough_start + pos
        intron_end = intron_st + intron_len

        if intron_st < 0 or intron_end + 2 > len(germ_seq):
            continue

        maybe_dnr = germ_seq[intron_st:intron_st + 2]
        maybe_acc = germ_seq[intron_end - 2:intron_end]
        bndry_sites = [('gt', 'ag'), ('ct','ac')]

        for bndry in bndry_sites:
            if maybe_dnr == bndry[0] and maybe_acc == bndry[1]:
                intron_pos = intron_st
                # print(intron_st, intron_end)
                break

    return intron_pos


def check_pointer_intron(soma_seq, germ_seq, gap_start, gap_end, max_shit: int = 5):
    germ_pntr = germ_seq[gap_start: gap_start + 8]
    soma_pntr = soma_seq[gap_end + 1: gap_end + 9]

    pntr_len = check_shared_pointer(germ_pntr, soma_pntr)

    if pntr_len > 4:
        return 'IES'

    elif check_intron_bndry(
            germ_seq,
            gap_start,
            gap_end,
            max_shift):
        return 'INTRON'

    else:
        return 'IES'


def something(fasta_file):
    soma, germ = [i for i in SeqIO.parse(fasta_file,'fasta')]
    gap_ranges = extract_gap_pos(f'{soma.seq}')
    germ_seq = f'{germ.seq}'
    soma_seq = f'{soma.seq}'


for grange in gap_ranges:
    x = check_pointer_intron(soma_seq, germ_seq, grange[0], grange[1])
    print(grange, x)
