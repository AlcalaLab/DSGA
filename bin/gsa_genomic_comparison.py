#!/usr/bin/env python3

"""
Code for running the genomic comparisons. Includes the management of the
comparative approach as well as the initial parsing of the genomic architectures.

Dependencies:
- BLAST+
- BioPython
- MAFFT

Last updated: 31-07-26
"""

import os, sys

from Bio.Seq import Seq
from Bio import SeqIO
from pathlib import Path

from collections import defaultdict

def prepare_fasta(outdir: str, taxon_name: str,fasta: str, min_len: int = 10000, germ: bool = True):
    fasta_prep_dir = f'{outdir}/Assembly_Backup/'
    # print(fasta_prep_dir)

    Path(fasta_prep_dir).mkdir(exist_ok = True, parents = True)

    if germ:
        out_fasta = f'{fasta_prep_dir}/{taxon_name}.Germ.{int(min_len/1000)}Kbp.fasta'
    else:
        out_fasta = f'{fasta_prep_dir}/{taxon_name}.Soma.{min_len}bp.fasta'

    x = []
    nseqs = 1
    for i in SeqIO.parse(germ_fasta, 'fasta'):
        if len(i.seq) < min_len:
            continue
        if germ:
            n = f'{taxon_name}_XX_Germ_{nseqs}_Len_{len(i.seq)}'
        else:
            n = f'{taxon_name}_XX_Soma_{nseqs}_Len_{len(i.seq)}'
        i.id = n
        i.description = ''
        i.name = ''
        x.append(i)
        nseqs += 1

    SeqIO.write(x, out_fasta, 'fasta')

    return out_fasta


def prep_blast(outdir: str, taxon_name: str, germ_fasta: str):
    out_db_dir = f'{outdir}/BLASTN_db/'
    germ_db = f'{out_db_dir}/{taxon_name}.germdb'

    Path(out_db_dir).mkdir(exist_ok = True, parents = True)

    blastdb_cmd = 'makeblastdb -dbtype nucl ' \
                    f'-in {germ_fasta} ' \
                    f'-out {germ_db}'

    if Path(f'{germ_db}.ndb').is_file() and Path(f'{germ_db}.ndb').stat().st_size > 0:
        return germ_db

    else:
        print('Preparing BLASTN database')
        os.system(blastdb_cmd)

    return germ_db


def blast_germ_soma(
        outdir: str,
        taxon_name: str,
        germ_fasta: str,
        soma_fasta: str,
        min_germ: int = 10000,
        min_soma: int = 400,
        threads: int = 4) -> str:

    out_tsv_dir = f'{outdir}/SpreadSheets/'

    out_tsv = f'{out_tsv_dir}{taxon_name}.Germ{int(min_germ/1000)}kbp_Soma{min_soma}bp.BLASTN.tsv'

    Path(out_tsv_dir).mkdir(exist_ok = True, parents = True)

    germ_db = prep_blast(outdir, taxon_name, germ_fasta)

    blast_cmd = 'blastn -ungapped ' \
                '-outfmt 6 ' \
                f'-num_threads {threads} ' \
                f'-db {germ_db} ' \
                f'-query {soma_fasta} ' \
                f'-out {out_tsv}'

    if Path(out_tsv).is_file() and Path(out_tsv).stat().st_size > 0:
        return out_tsv
    else:
        print('Running BLASTN -- This may take a while!')
        os.system(blast_cmd)

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


def extract_valid_overlaps(sorted_cords, min_pointer = 2, max_pointer = 25):
    # sort by starting position and record the index for easier extraction later
    indexed_ranges = sorted(
            [(r[2], r[3], i) for i, r in enumerate(sorted_coords)],
            key = lambda x: (x[0], x[1]))

    overlapping_pairs = []

    for i in range(len(indexed_ranges)):
        start1, end1, idx1 = indexed_ranges[i]

        for j in range(i + 1, len(indexed_ranges)):
            start2, end2, idx2 = indexed_ranges[j]

            if start2 >= end1:
                break

            # Calculate the explicit overlap bounds
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            overlap_size = overlap_end - overlap_start

            # 3. Filter based on the size constraints
            if min_pointer <= overlap_size <= max_pointer:
                overlapping_pairs.append({
                    "range_1": (start1, end1),
                    "range_2": (start2, end2),
                    "original_indices": (idx1, idx2),
                    "overlap_size": overlap_size
                })


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


def filter_hits(
        out_tsv: str,
        min_aln_prop: float = 0.6):

    # soma_germ_dict = defaultdict(list)
    skip_soma_germ = []
    soma_germ_summary = []
    soma_germ_dict = defaultdict(list)
    # soma_germ_summary = defaultdict(list)

    for line in open(out_tsv).readlines():
        s, g = line.split('\t')[:2]
        aln = line.split('\t')[3]
        ss, se, gs, ge = line.split('\t')[6:10]
        soma_germ_dict[s].append((g, int(aln), int(ss), int(se), int(gs), int(ge)))

        # eval = 'Schmidingerella_arcuata_XX_Soma_1107_Len_1020_Kcov_11.62'

        for k, v in soma_germ_dict.items():
            if not check_soma_coverage(k, v, min_aln_prop):
                skip_soma_germ.append(k)
                continue

            sorted_coords = sorted(v, key = lambda x: x[2])

            print(k)
            filt_coords = drop_overlap_coords(sorted_coords)

            if not check_soma_coverage(k, filt_coords):
                skip_soma_germ.append(k)
                continue

            germ_arch_type = eval_locus_type(filt_coords)
            print(germ_arch_type)
            mds_num = 1

            for i in sorted_coords:
                updated_line = '\t'.join(f'{n}' for n in i)
                soma_germ_summary.append(f'{k}\t{updated_line}\t{germ_arch_type}\tMDS-{mds_num}')
                mds_num += 1

        if soma_germ_summary:
            print(len(soma_germ_dict))
            return soma_germ_summary, soma_germ_dict, list(set(skip_soma_germ))



def refine_nonscrambled():
    pass


def save_summary_tsv(sg_summary: list, outdir: str, taxon_name: str):
    dsga_tsv = f'{outdir}/{taxon_name}.DSGA.Summary.tsv'
    header = 'Soma\tGerm\tAlignment_Length\tSoma_Start\tSoma_End\tGerm_Start\tGerm_End\tGSA_Type\tMDS_Number\n'
    with open(dsga_tsv, 'w+') as w:
        w.write(header)
        w.write('\n'.join(sg_summary))



print('\nworking on the pointer overlap instead of relying on inferences solely based pure overlaps! Do this BY locus pairing?')

print('\nThis code is NOT functional alone, yet!\n')

sys.exit()
"""
for each locus ... find out if it covers a substantial portion

if it doesn't cover at least 20% of the CDS then move along

then, starting with the biggest, identify all the pointer overlaps!

if it doesn't have any, check the next locus...

for all the missing MDSs, search the updated sets of coords for pointers?

if you still can't then just die because this is fucking hard -- likely scrambled in
a complicated up way...
"""


taxon_name = ''

if not taxon_name:
    sys.exit()

out_dir = f'{taxon_name}_DSGA/'

threads = 2

germ_fasta = ''
soma_fasta = ''


germ_filt_fasta = prepare_fasta(out_dir, taxon_name, germ_fasta)
soma_filt_fasta = prepare_fasta(out_dir, taxon_name, soma_fasta, 400, False)

out_tsv = blast_germ_soma(
            out_dir,
            taxon_name,
            germ_filt_fasta,
            soma_filt_fasta,
            10000,
            400,
            threads)

sg_summary, soma_germ_dict, skip_soma_germ = filter_hits(out_tsv)
save_summary_tsv(sg_summary, out_dir, taxon_name)
