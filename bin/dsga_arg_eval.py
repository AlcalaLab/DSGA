#!/usr/bin/env python3

import argparse, sys


def capture_args():
    # create the top-level parser
    gen_descript = 'Delineating Somatic and Germline Architecture\n\nUsage:\n     dsga.py <options>\n\n'

    parser = argparse.ArgumentParser(description = gen_descript,
            usage = argparse.SUPPRESS, add_help = False,
            formatter_class = argparse.RawDescriptionHelpFormatter)

    # parser.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)
    parser.add_argument('--version', action = 'store_true', help = argparse.SUPPRESS)

    sga = parser.add_argument_group('General Options', description = (
        '''--out-dir (-o)        output directory name\n'''
        '''--soma (-s)           input SOMATIC FASTA file\n'''
        '''--germ (-g)           input SOMATIC FASTA file\n'''
        '''--min-soma (-ms)      minimum length of somatic sequences (default = 400 base pairs)\n'''
        '''--min-germ (-mg)      minimum length of germline sequences (default = 10,000 base pairs)\n'''
        '''--threads (-t)        number of CPU threads to use (default = 4)\n'''
        '''--quiet (-q)          no console output\n'''
        '''--gzip (-gz)          tar and gzip output\n'''
        '''--help (-h)           show this help message and exit\n'''))

    adv_sga = parser.add_argument_group('Advanced Options', description = (
        '''--min-aln (-a)        minimum proportion of the somatic sequences length (default = 0.60)\n'''
        '''--min-id (-id)        minimum percent identity for soma-germ alignments (default = 95.0)\n'''
        '''--min-pointer (-p)    minimum length (base pairs) for pointer sequences (default = 2)\n'''
        '''--max-pointer (-P)    maximum length (base pairs) for pointer sequences (default = 25)\n'''
        '''--telomere            single repeat of the telomere sequence\n\n'''))

    sga.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)

    sga.add_argument('--out-dir','-o', action = 'store', metavar = '[output-directory]',
        type = str, help = argparse.SUPPRESS)

    sga.add_argument('--soma', '-s', action = 'store', metavar = '[somatic-fasta]',
        type = str, help = argparse.SUPPRESS)

    sga.add_argument('--germ', '-g', action = 'store', metavar = '[germline-fasta]',
        type = str, help = argparse.SUPPRESS)

    sga.add_argument('--min-soma', '-ms', action = 'store', default = 400,
        metavar = '[min-soma-length]', type = int, help = argparse.SUPPRESS)

    sga.add_argument('--min-germ', '-mg', action = 'store', default = 10000,
        metavar = '[min-germ-length]', type = int, help = argparse.SUPPRESS)

    sga.add_argument('--threads','-t', action = 'store', default = 4,
        metavar = '[Threads]', type = int, help = argparse.SUPPRESS)

    sga.add_argument('--clean', action = 'store_true', help = argparse.SUPPRESS)

    sga.add_argument('--quiet', '-q', action = 'store_true', help = argparse.SUPPRESS)

    sga.add_argument('--gzip', '-gz', action = 'store_true', help = argparse.SUPPRESS)

    adv_sga.add_argument('--min-aln', '-a', action = 'store', default = 0.6,
        metavar = '[min-proportion-soma-aligned]', type = float, help = argparse.SUPPRESS)

    adv_sga.add_argument('--min-id', '-id', action = 'store', default = 95.0,
        metavar = '[min-percent-identity]', type = float, help = argparse.SUPPRESS)

    adv_sga.add_argument('--min-pointer', '-p', action = 'store', default = 2,
        metavar = '[min-pointer-length]', type = int, help = argparse.SUPPRESS)

    adv_sga.add_argument('--max-pointer', '-P', action = 'store', default = 25,
        metavar = '[max-pointer-length]', type = int, help = argparse.SUPPRESS)

    if len(sys.argv) == 1:
#     print(ascii_logo_vsn())
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    args.min_pointer = abs(args.min_pointer)
    args.max_pointer = abs(args.max_pointer)

    args_eval = double_check_args(args)

    if not args_eval:
        sys.exit(1)

    return args_eval[1]


def double_check_args(args):
    req_args = []
    dsga_warnings = []

    if not args.soma:
        req_args.append('    --soma <somatic-sequences-fasta>')

    if not args.germ:
        req_args.append('    --germ <germ-sequences-fasta>')

    if not args.out_dir:
        req_args.append('    --out-dir <output-directory-name>')

    if not 0 < args.min_aln <= 1.0:
        dsga_warnings.append('WARNING: minimum proportion of the soma being aligned is not between 0 and 1')
        dsga_warnings.append('         Defaulting to 0.6')
        args.min_aln = 0.6

    if not 0 < args.min_id <= 100:
        dsga_warnings.append('WARNING: percent identity for parsing germline-soma alignments is not between 0 and 100')
        dsga_warnings.append('         Defaulting to 95.0%')
        args.id = 95.0

    if args.min_pointer >= args.max_pointer:
        dsga_warnings.append('WARNING: minimum pointer length must be smaller than the maximum pointer length')
        dsga_warnings.append('         Defaulting to 2 base pairs')
        args.min_pointer = 2

    if args.max_pointer < args.min_pointer:
        dsga_warnings.append('WARNING: maximum pointer length must be larger than the minimum pointer length')
        dsga_warnings.append('         Defaulting to 25 base pairs')
        args.max_pointer = 25

    if dsga_warnings:
        print('\n#----------------- WARNINGS -----------------#')
        print('\n'.join(dsga_warnings) + '\n')

    if req_args:
        print('\n#----------------- ERRORS ------------------#')
        print('ERROR: Missing the following required arguments:')
        print("\n".join(req_args) + '\n')

        return False

    else:
        return True, args
