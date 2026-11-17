#!/usr/bin/env python3

import glob, os

from itertools import groupby
from operator import itemgetter
from pathlib import Path

from Bio import SeqIO


def extract_gap_pos(
        soma_seq: str) -> list:
    gap_pos = [i for i, j in enumerate(soma_seq) if j == '-']

    gap_ranges = []

    for k, g in groupby(enumerate(gap_pos), lambda x: x[0]-x[1]):
        group = list(map(itemgetter(1), g))

        if len(group) > 1:
            gap_ranges.append((group[0],group[-1]))

        else:
            gap_ranges.append((group[0],group[0]))

    return gap_ranges[1:-1]


def check_shared_pointer(
        germ_pntr: str,
        soma_pntr: str) -> int:

    shared_ntds = ''

    for i in range(len(germ_pntr)):
        if germ_pntr[i] == soma_pntr[i]:
            shared_ntds += germ_pntr[i]

        else:
            break

    return len(shared_ntds)


def check_intron_bndry(
        germ_seq: str,
        gap_start: int,
        gap_end: int,
        max_shift: int = 5) -> int:

    intron_pos = None

    intron_len = 1 + gap_end - gap_start

    # Check shifts from -max_shift to +max_shift
    for pos in range(-max_shift, max_shift + 1):
        intron_st = gap_start + pos
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


def compare_pointer_intron(
        soma_seq: str,
        germ_seq: str,
        gap_start: int,
        gap_end: int,
        max_shift: int = 5,
        min_gap: int = 25) -> str:

    germ_pntr = germ_seq[gap_start: gap_start + 8]
    soma_pntr = soma_seq[gap_end + 1: gap_end + 9]

    pntr_len = check_shared_pointer(germ_pntr, soma_pntr)

    ies_len_eval = False

    if gap_end - gap_start > min_gap:
        ies_len_eval = True

    if pntr_len > 3 and ies_len_eval:
        return 'IES'

    else:
        intron_pos = check_intron_bndry(
                        germ_seq,
                        gap_start,
                        gap_end,
                        max_shift)

        if intron_pos:
            return 'INTRON'

        else:
            if pntr_len >= 2 and ies_len_eval:
                return 'IES'

            else:
                return 'UNCLEAR'


def eval_align_architecture(
        fasta_file: str,
        max_shift: int = 5,
        min_gap: int = 25):
    soma, germ = [i for i in SeqIO.parse(fasta_file,'fasta')]
    gap_ranges = extract_gap_pos(f'{soma.seq}')
    germ_seq = f'{germ.seq}'
    soma_seq = f'{soma.seq}'

    clear_ranges = []

    for grange in gap_ranges:
        x = compare_pointer_intron(soma_seq, germ_seq, grange[0], grange[1])
        # print(grange, x)
        if x == 'INTRON':
            clear_ranges.append(x)
        elif x == 'UNCLEAR':
            clear_ranges.append(x)
        else:
            clear_ranges.append(grange)

    return clear_ranges


def eval_local_alns(
        local_align_dir: str) -> None:
    """
    Note to self!

    this is just a temporary function!

    The outputs need to be assessed more deeply:
      - update the alignment tables
        - update the locus type
        - update the MDS coordinates for MDSs interrupted by introns
      - check the non-scrambled portions of the scrambled loci
    """

    intronic_dir = f'{local_align_dir.rstrip("/").rpartition("/")[0]}/Single_MDS_With_Intron/'
    intronic_nsc_dir = f'{local_align_dir.rstrip("/").rpartition("/")[0]}/Non_Scrambled_With_Intron/'
    only_nsc_dir = f'{local_align_dir.rstrip("/").rpartition("/")[0]}/Non_Scrambled_Loci/'
    # unclear_dir = f'{local_align_dir.rstrip("/").rpartition("/")[0]}/Unclear_Loci/'

    Path(intronic_dir).mkdir(exist_ok = True, parents = True)
    Path(intronic_nsc_dir).mkdir(exist_ok = True, parents = True)
    Path(only_nsc_dir).mkdir(exist_ok = True, parents = True)
    # Path(unclear_dir).mkdir(exist_ok = True, parents = True)

    for aln in glob.glob(f'{local_align_dir}/*fasta'):
        ies_intron_eval = eval_align_architecture(aln)
        print(ies_intron_eval)
        print('INTRON' not in ies_intron_eval and len(set(ies_intron_eval)) > 0)

        if len(set(ies_intron_eval)) == 1 and 'INTRON' in ies_intron_eval:
            os.system(f'cp {aln} {intronic_nsc_dir}')

        elif len(set(ies_intron_eval)) != 1 and'INTRON' in ies_intron_eval:
            os.system(f'cp {aln} {intronic_dir}')

        elif 'INTRON' not in ies_intron_eval and len(set(ies_intron_eval)) > 0:
            os.system(f'cp {aln} {only_nsc_dir}')

        else:
            pass
            # os.system(f'cp {aln} {unclear_dir}')
