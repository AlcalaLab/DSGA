#!/usr/bin/env python3

"""
Code for running the genomic comparisons. Includes the management of the
comparative approach as well as the initial parsing of the genomic architectures.

Dependencies:
- BLAST+
- BioPython

Last updated: 18-09-26
"""

import glob, os, subprocess, sys, warnings

from collections import defaultdict
from pathlib import Path

# just to suppress BioPython's Runtimewarning -- do not want to panic people
warnings.filterwarnings("ignore", category = RuntimeWarning)

from Bio.Seq import Seq
from Bio import SeqIO




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


def filter_blastn_hits(
        out_tsv: str,
        min_aln_prop: float = 0.6,
        min_pointer = 2,
        max_pointer = 25,
        multi_filt: bool = True,
        id_filt: float = 95.0):

    # soma_germ_dict = defaultdict(list)
    skip_soma_germ = []
    soma_germ_summary = []
    soma_germ_dict = defaultdict(list)

    # eval_soma = 'Tintinnopsis_tocantinensis_LKH824_LKH869_XX_Soma_13_Len_1607'

    for line in open(out_tsv).readlines():
        # if eval_soma not in line:
        #     continue
        if id_filt and float(line.split('\t')[2]) < id_filt:
                continue
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

        else:
            soma_germ_dict[k] = filt_coords

        germ_arch_type = eval_locus_type(final_coords)

        mds_num = 1

        for i in final_coords:
            updated_line = '\t'.join(f'{n}' for n in i)

            soma_germ_summary.append(f'{k}\t{updated_line}\t{germ_arch_type}\tMDS-{mds_num}')

            mds_num += 1

    if soma_germ_summary:

        return soma_germ_summary, soma_germ_dict, list(set(skip_soma_germ))

    return None


def save_summary_tsv(sg_summary: list, outdir: str, out_name: str, multi_filt: bool = True):
    dsga_tsv = f'{outdir}/{out_name}.DSGA'

    if multi_filt:
        dsga_tsv += '_MultiFilt.Summary'

    header = 'Soma\tGerm\tAlignment_Length\tSoma_Start\tSoma_End\tGerm_Start\tGerm_End\tGSA_Type\tMDS_Number\n'

    with open(f'{dsga_tsv}.tsv', 'w+') as w:
        w.write(header)
        w.write('\n'.join(sg_summary))


def nonscram_loci_for_eval(
        sg_summary: list,
        soma_germ_dict: dict
        ) -> dict:

    nsc_eval_dict = {}
    for i in sg_summary:
        if 'Non-Scrambled' in i:

            nsc_locus = sorted(soma_germ_dict[i.split('\t')[0]], key = lambda x: x[2])
            nsc_sst = nsc_locus[0][2]
            nsc_send = nsc_locus[-1][3]

            nsc_gst = nsc_locus[0][-2]
            nsc_gend = nsc_locus[-1][-1]

            frame = 'f'

            if nsc_gst > nsc_gend:
                frame = 'r'

            nsc_eval_dict[i.split('\t')[0]] = [nsc_locus[0][0], nsc_sst, nsc_send, nsc_gst, nsc_gend, frame]

    return nsc_eval_dict


def align_germ_soma(
        out_unaln: str,
        out_aln: str,
        threads: int = 4):

    for f in glob.glob(f'{out_unaln}*fasta'):
        out_f = f.replace("UnAlign","Align").replace(".fasta",".EINSI.fasta")

        mafft_cmd = f'einsi --quiet --thread {threads} {f} > {out_f}'

        mafft_result = subprocess.run(mafft_cmd, shell = True, stdout = subprocess.DEVNULL, check = True)


def prep_nonscram_loci(
        out_dir: str,
        soma_prep_fasta: str,
        germ_prep_fasta: str,
        sg_summary: list,
        soma_germ_dict: dict):

    nsc_eval_dict = nonscram_loci_for_eval(
                        sg_summary,
                        soma_germ_dict)

    out_unaln = f'{out_dir}Soma_Germ_Refinement/UnAligned_Soma_Germ/'
    out_aln = out_unaln.replace("UnAlign","Align")

    Path(out_unaln).mkdir(exist_ok = True, parents = True)
    Path(out_aln).mkdir(exist_ok = True, parents = True)

    germ_loci = list(set([v[0] for v in nsc_eval_dict.values()]))

    soma_seqs = {i.id:f'{i.seq}' for i in SeqIO.parse(soma_prep_fasta, 'fasta') if i.id in nsc_eval_dict.keys()}
    germ_seqs = {i.id:f'{i.seq}' for i in SeqIO.parse(germ_prep_fasta, 'fasta') if i.id in germ_loci}

    for k, v in nsc_eval_dict.items():
        out_fasta = f'{k}_XX_{v[0].rpartition("_XX_")[-1]}.UnAligned.fasta'
        s_st = v[1]-1
        s_end = v[2]
        seqs = f'>{k}_XX_{s_st}_{s_end}\n{soma_seqs[k][s_st:s_end]}\n'

        if v[-1] == 'f':
            g_st = max(v[-3]-10, 0)
            g_end = min(v[-2]+100, int(v[0].rpartition("_")[-1]))
            seqs += f'>{v[0]}_XX_{g_st}_{g_end}\n{germ_seqs[v[0]][g_st:g_end+1]}'

        else:
            g_st = max(v[-2]-100, 0)
            g_end = min(v[-3]+100, int(v[0].rpartition("_")[-1]))
            g_seq = Seq(germ_seqs[v[0]][g_st:g_end+1]).reverse_complement()
            seqs += f'>{v[0]}_XX_{g_end}_{g_st}_RC\n{g_seq}'

        with open(f'{out_unaln}{out_fasta}','w+') as w:
            w.write(seqs)

    align_germ_soma(out_unaln, out_aln)

    return out_aln
