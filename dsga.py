#!/usr/bin/env python3

"""
Code for running the genomic comparisons. Includes the management of the
comparative approach as well as the initial parsing of the genomic architectures.

Dependencies:
- BLAST+
- BioPython

Last updated: 18-09-26
"""

import os, subprocess, sys, warnings

from collections import defaultdict
from pathlib import Path

# just to suppress BioPython's Runtimewarning -- do not want to panic people
warnings.filterwarnings("ignore", category = RuntimeWarning)

from Bio.Seq import Seq
from Bio import SeqIO


def check_exising_files(filename: str):
    return Path(filename).is_file() and Path(filename).stat().st_size > 0


def prepare_fasta(outdir: str, out_name: str, fasta: str, min_len: int = 10000, germ: bool = True):
    fasta_prep_dir = f'{outdir}Assembly_Backup/'
    # print(fasta_prep_dir)

    Path(fasta_prep_dir).mkdir(exist_ok = True, parents = True)

    if germ:
        out_fasta = f'{fasta_prep_dir}{out_name}.Germ.{int(min_len/1000)}Kbp.fasta'
    else:
        out_fasta = f'{fasta_prep_dir}{out_name}.Soma.{min_len}bp.fasta'

    x = []
    nseqs = 1
    for i in SeqIO.parse(fasta, 'fasta'):
        if len(i.seq) < min_len:
            continue
        if germ:
            n = f'{out_name}_XX_Germ_{nseqs}_Len_{len(i.seq)}'
        else:
            n = f'{out_name}_XX_Soma_{nseqs}_Len_{len(i.seq)}'
        i.id = n
        i.description = ''
        i.name = ''
        x.append(i)
        nseqs += 1

    if len(x) == 0:
        if germ:
            print(f'\nERROR: No germline loci found that are longer than {min_len}bp')
        else:
            print(f'\nERROR: No somatic sequences found that are longer than {min_len}bp')
        sys.exit()

    SeqIO.write(x, out_fasta, 'fasta')

    return out_fasta


def prep_blast(outdir: str, out_name: str, germ_fasta: str):
    out_db_dir = f'{outdir}/BLASTN_db/'
    germ_db = f'{out_db_dir}/{out_name}.germdb'

    Path(out_db_dir).mkdir(exist_ok = True, parents = True)

    blastdb_cmd = ['makeblastdb', '-dbtype', 'nucl', '-in', f'{germ_fasta}', '-out', f'{germ_db}']

    if check_exising_files(f'{germ_db}.ndb'):
        return germ_db

    else:
        print('Preparing BLASTN database')

        blastdb_result = subprocess.run(blastdb_cmd, stdout = subprocess.DEVNULL, check = True)

    return germ_db


def blast_germ_soma(
        outdir: str,
        out_name: str,
        germ_fasta: str,
        soma_fasta: str,
        min_germ: int = 10000,
        min_soma: int = 400,
        threads: int = 4) -> str:

    out_tsv_dir = f'{outdir}/SpreadSheets/'

    out_tsv = f'{out_tsv_dir}{out_name}.Germ{int(min_germ/1000)}kbp_Soma{min_soma}bp.BLASTN.tsv'

    Path(out_tsv_dir).mkdir(exist_ok = True, parents = True)

    germ_db = prep_blast(outdir, out_name, germ_fasta)

    blast_cmd = [
        'blastn', '-ungapped', '-outfmt', '6', '-num_threads', f'{threads}', '-db',
        f'{germ_db}', '-query', f'{soma_fasta}', '-out', f'{out_tsv}']

    if check_exising_files(out_tsv):
        return out_tsv

    else:
        print('Running BLASTN -- This may take a while!')
        blastn_result = subprocess.run(blast_cmd, stdout = subprocess.DEVNULL, check = True)

    return out_tsv


def check_overlap(x: tuple, y: tuple, max_overlap: int = 50, thresh: float = 0.5) -> bool:
    overlap_val = len(range(max(x[2], y[2]), min(x[3],y[3])+1))

    min_frag_size = min(max_overlap, min((x[3] - x[2]), (y[3] - y[2])))

    return overlap_val >= int(min_frag_size * thresh)


def drop_overlap_coords(coords_sorted_lst: list, thresh: float = 0.5) -> list:
    coords_keep = [coords_sorted_lst[0]]

    diff = lambda x: abs(x[2]-x[3])

    for x in coords_sorted_lst[1:]:
        if check_overlap(coords_keep[-1], x):
            if diff(coords_keep[-1]) == diff(x):
                coords_keep.append(x)

            else:
                coords_keep[-1] = max(coords_keep[-1], x, key = diff)

        else:
            coords_keep.append(x)

    return coords_keep


def index_ranges(filt_coords, soma: bool = True):
    if soma:
        indexed_ranges = sorted(
                [(j[2], j[3], i) for i, j in enumerate(filt_coords)],
                key = lambda x: x[-1])
    else:
        indexed_ranges = sorted(
                [(j[4], j[5], i) for i, j in enumerate(filt_coords)],
                key = lambda x: x[-1])


    return indexed_ranges


def merge_germ_ranges(germ_ranges):
    sorted_ranges = sorted(germ_ranges, key = lambda x: x[0])
    merged = []
    for r in sorted_ranges:
        if not merged or merged[-1][1] < r[0]:
            merged.append(r)
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], r[1]))
    return merged


def update_filt_coords(filt_coords, filt_germ_ranges):
    updated_coords = []

    for i in filt_germ_ranges:
        soma_st = filt_coords[i[0]][2]
        soma_end = filt_coords[i[1]][3]
        germ_st = filt_coords[i[0]][-2]
        germ_end = filt_coords[i[1]][-1]
        new_coord = (filt_coords[i[0]][0], soma_end - soma_st, soma_st, soma_end, germ_st, germ_end)
        updated_coords.append(new_coord)

    return updated_coords


def check_germ_overlap(filt_coords, min_ies: int = 5):
    germ_overlap_ranges = []

    germ_ranges = index_ranges(filt_coords, False)

    for n in range(len(germ_ranges)-1):
        if germ_ranges[n][0] < germ_ranges[n][1]:
            if germ_ranges[n+1][0] in range(germ_ranges[n][0], germ_ranges[n][1]+1):
                germ_overlap_ranges.append((n, n+1))

            elif abs(germ_ranges[n+1][0] - germ_ranges[n][1]) < min_ies:
                germ_overlap_ranges.append((n, n+1))

            else:
                germ_overlap_ranges.append((n,n))
        else:
            if germ_ranges[n+1][0] in range(germ_ranges[n][1], germ_ranges[n][0]+1):
                germ_overlap_ranges.append((n, n+1))

            elif abs(germ_ranges[n+1][0] - germ_ranges[n][1]) < min_ies:
                germ_overlap_ranges.append((n, n+1))

            else:
                germ_overlap_ranges.append((n,n))

    germ_overlap_ranges.append((n+1, n+1))

    filt_germ_ranges = merge_germ_ranges(germ_overlap_ranges)

    return update_filt_coords(filt_coords, filt_germ_ranges)


def refine_valid_pointers(filt_coords, min_pointer = 2, max_pointer = 25):

    overlapping_pairs = []
    consecutive_mds = []

    soma_ranges = index_ranges(filt_coords, True)

    for i in range(len(soma_ranges)-1):
        mds_st_1, mds_end_1, mds_idx_1 = soma_ranges[i]
        mds_st_2, mds_end_2, mds_idx_2 = soma_ranges[i+1]

        if mds_st_2 > mds_end_1:
            continue

        pointer_size = (mds_end_1 - mds_st_2)+1

        if min_pointer <= pointer_size <= max_pointer:
            overlapping_pairs.append({
                "mds_1": (mds_st_1, mds_end_1),
                "mds_2": (mds_st_2, mds_end_2),
                "original_indices": (mds_idx_1, mds_idx_2),
                "pointer_size": pointer_size})

            consecutive_mds.append((filt_coords[mds_idx_1], filt_coords[mds_idx_2]))

    return consecutive_mds


def eval_sign(coords: list) -> str:
    try:
        return (coords[-1] - coords[-2]) / abs(coords[-1] - coords[-2])
    except ZeroDivisionError:
        return None


def eval_locus_type(coords_lst: list) -> str:
    coords_sorted_lst = sorted(coords_lst, key = lambda x: x[2])

    num_loci = len(set([i[0] for i in coords_sorted_lst]))

    mds_germ_signs = [eval_sign(coords_lst[0])]

    for n in range(1, len(coords_sorted_lst)):
        mds_germ_signs.append(eval_sign(coords_sorted_lst[n]))
        mds_germ_signs.append(eval_sign((coords_sorted_lst[n-1][-1], coords_sorted_lst[n][-2])))

    mds_germ_signs = [i for i in mds_germ_signs if i]

    if len(mds_germ_signs) == 1:
        return 'Single-MDS'

    elif len(set(mds_germ_signs)) == 1 and num_loci == 1:
        return 'Non-Scrambled-Single-Locus'

    elif len(set(mds_germ_signs)) == 1 and num_loci != 1:
        return 'Non-Scrambled-Multi-Locus'

    elif len(set(mds_germ_signs)) != 1 and num_loci == 1:
        return 'Scrambled-Single-Locus'

    elif len(set(mds_germ_signs)) != 1 and num_loci != 1:
        return 'Scrambled-Multi-Locus'

    else:
        return 'This-Should-Not-Work-Contact-Xyrus'


def check_soma_coverage(
        soma_name: str,
        soma_hits: list,
        min_aln_prop: float = 0.6) -> bool:

    min_soma = int(int(soma_name.rpartition("Len_")[-1].partition("_")[0]) * min_aln_prop)

    return sum([i[1] for i in soma_hits]) >= min_soma


def sort_and_check_num_loci(coords_lst):
    num_loci = len(set([i[0] for i in coords_lst]))
    if num_loci > 1:
        sorted_coords = sorted(coords_lst, key = lambda x: (x[0], x[2]))

    else:
        sorted_coords = sorted(coords_lst, key = lambda x:  x[2])

    return sorted_coords, num_loci


def filter_loci(soma_name, sorted_coords, min_aln_prop, min_pointer = 2, max_pointer = 25):
    loci_dict = defaultdict(list)
    init_filt_dict = {}

    skip_soma_germ = []

    best_cover = [0,'']

    soma_len = int(soma_name.rpartition("Soma_")[-1].partition("_")[0])

    for i in sorted_coords:
        loci_dict[i[0]].append(i)

    for k, v in loci_dict.items():
        if not check_soma_coverage(soma_name, v, min_aln_prop):
            skip_soma_germ.append(k)

        elif len(v) == 1:
            init_filt_dict[k] = v

        else:
            filt_coords = drop_overlap_coords(v)
            pointer_eval = refine_valid_pointers(filt_coords, min_pointer = 2, max_pointer = 25)

            # Check note below (so redundant... but to make it hard to miss)
            """
            NOTE FOR MYSELF: there needs to a missing-MDS option!!!!

            add a missing MDS option down the road...
            """

            if pointer_eval:
                init_filt_dict[k] = filt_coords

            else:
                init_filt_dict[k] = filt_coords

    for k, v in init_filt_dict.items():
        if sum([i[1] for i in v]) > best_cover[0]:
            best_cover = [sum([i[1] for i in v]), k]

    if best_cover[-1]:
        return init_filt_dict[best_cover[-1]]

    return None


def filter_hits(
        out_tsv: str,
        min_aln_prop: float = 0.6,
        min_pointer = 2,
        max_pointer = 25,
        multi_filt: bool = True):

    # soma_germ_dict = defaultdict(list)
    skip_soma_germ = []
    soma_germ_summary = []
    soma_germ_dict = defaultdict(list)

    # eval_soma = 'Tintinnopsis_tocantinensis_LKH824_LKH869_XX_Soma_13_Len_1607'

    for line in open(out_tsv).readlines():
        # if eval_soma not in line:
        #     continue
        s, g = line.split('\t')[:2]
        aln = line.split('\t')[3]
        ss, se, gs, ge = line.split('\t')[6:10]

        soma_germ_dict[s].append((g, int(aln), int(ss), int(se), int(gs), int(ge)))


    for k, v in soma_germ_dict.items():
        if not check_soma_coverage(k, v, min_aln_prop):
            skip_soma_germ.append(k)
            continue

        sorted_coords, num_loci = sort_and_check_num_loci(v)

        filt_coords = []
        final_coords = []

        if not multi_filt:
            filt_coords = drop_overlap_coords(sorted_coords)

        else:
            filt_coords = filter_loci(k, sorted_coords, min_aln_prop, min_pointer, max_pointer)

        if not filt_coords:
            continue

        if len(filt_coords) > 1:
            final_coords = check_germ_overlap(filt_coords)
        else:
            final_coords = filt_coords

        if not final_coords:
            skip_soma_germ.append(k)
            continue

        if not check_soma_coverage(k, final_coords):
            skip_soma_germ.append(k)
            continue

        germ_arch_type = eval_locus_type(final_coords)

        mds_num = 1

        for i in final_coords:
            updated_line = '\t'.join(f'{n}' for n in i)

            soma_germ_summary.append(f'{k}\t{updated_line}\t{germ_arch_type}\tMDS-{mds_num}')

            mds_num += 1

    if soma_germ_summary:

        return soma_germ_summary, soma_germ_dict, list(set(skip_soma_germ))

    return None


def refine_nonscrambled():
    pass


def save_summary_tsv(sg_summary: list, outdir: str, out_name: str, multi_filt: bool = True):
    dsga_tsv = f'{outdir}/{out_name}.DSGA'

    if multi_filt:
        dsga_tsv += '_MultiFilt.Summary'

    header = 'Soma\tGerm\tAlignment_Length\tSoma_Start\tSoma_End\tGerm_Start\tGerm_End\tGSA_Type\tMDS_Number\n'

    with open(f'{dsga_tsv}.tsv', 'w+') as w:
        w.write(header)
        w.write('\n'.join(sg_summary))


def eval_germ_soma_arch(
        out_name: str,
        germ_fasta: str,
        soma_fasta: str,
        min_germ: int = 10000,
        min_soma: int = 400,
        min_aln_prop: float = 0.6,
        min_pointer: int = 2,
        max_pointer: int = 25,
        multi_filt: bool = True,
        threads: int = 4):

    out_dir = f'{out_name}_DSGA/'

    germ_filt_fasta = prepare_fasta(out_dir, out_name, germ_fasta, min_germ)
    soma_filt_fasta = prepare_fasta(out_dir, out_name, soma_fasta, min_soma, False)

    out_tsv = blast_germ_soma(
                out_dir,
                out_name,
                germ_filt_fasta,
                soma_filt_fasta,
                min_germ,
                min_soma,
                threads)

    sg_summary, soma_germ_dict, skip_soma_germ = filter_hits(
                                                    out_tsv,
                                                    min_aln_prop,
                                                    min_pointer,
                                                    max_pointer,
                                                    multi_filt)

    if sg_summary:
        save_summary_tsv(sg_summary, out_dir, out_name, multi_filt)


if __name__ == '__main__':
    try:
        soma_fasta = sys.argv[1]
        germ_fasta = sys.argv[2]
        output_name = sys.argv[3]

    except:
        print('\nUsage:\n\n    python3 gsa_genomic_comparison.py [SOMA-FASTA] [GERM-FASTA] [OUTPUT-NAME]\n')
        sys.exit()

    eval_germ_soma_arch(
            output_name,
            germ_fasta,
            soma_fasta,
            min_germ = 10000,
            min_soma = 400,
            min_aln_prop = 0.6,
            min_pointer = 2,
            max_pointer = 25,
            multi_filt = True,
            threads = 4)
