#!/usr/bin/env python3

import glob, os, subprocess, sys, warnings

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


def prep_blast(
        outdir: str,
        out_name: str,
        germ_fasta: str,
        quiet: bool = False) -> str:

    out_db_dir = f'{outdir}/BLASTN_db/'
    germ_db = f'{out_db_dir}/{out_name}.germdb'

    Path(out_db_dir).mkdir(exist_ok = True, parents = True)

    blastdb_cmd = ['makeblastdb', '-dbtype', 'nucl', '-in', f'{germ_fasta}', '-out', f'{germ_db}']

    if check_exising_files(f'{germ_db}.ndb'):
        return germ_db

    else:
        if not quiet:
            print('#-------- Preparing BLASTN Database ---------#')

        blastdb_result = subprocess.run(blastdb_cmd, stdout = subprocess.DEVNULL, check = True)

    return germ_db


def blast_germ_soma(
        outdir: str,
        out_name: str,
        germ_fasta: str,
        soma_fasta: str,
        min_germ: int = 10000,
        min_soma: int = 400,
        threads: int = 4,
        quiet: bool = False) -> str:

    out_tsv_dir = f'{outdir}/SpreadSheets/'

    out_tsv = f'{out_tsv_dir}{out_name}.Germ{int(min_germ/1000)}kbp_Soma{min_soma}bp.BLASTN.tsv'

    Path(out_tsv_dir).mkdir(exist_ok = True, parents = True)

    germ_db = prep_blast(
                outdir,
                out_name,
                germ_fasta,
                quiet)

    blast_cmd = [
        'blastn', '-ungapped', '-outfmt', '6', '-num_threads', f'{threads}', '-db',
        f'{germ_db}', '-query', f'{soma_fasta}', '-out', f'{out_tsv}']

    if check_exising_files(out_tsv):
        return out_tsv

    else:
        if not quiet:
            print('#-------- Running BLASTN Comparison ---------#')

        blastn_result = subprocess.run(blast_cmd, stdout = subprocess.DEVNULL, check = True)

    return out_tsv
