##!/usr/bin/env python3


from bin import dsga_arg_eval as ae
from bin import dsga_prep_soma_germ as psg
from bin import dsga_soma_germ_comp as sgc
from bin import dsga_refine_align as rfa

def eval_germ_soma_arch(
        out_name: str,
        germ_fasta: str,
        soma_fasta: str,
        min_germ: int = 10000,
        min_soma: int = 400,
        min_aln_prop: float = 0.6,
        id_filt: float = 95.0,
        min_pointer: int = 2,
        max_pointer: int = 25,
        multi_filt: bool = True,
        threads: int = 4,
        quiet: bool = True):

    if not out_name.rstrip("/").endswith('DSGA'):
        out_dir = f'{out_name}_DSGA/'

    else:
        out_dir = out_name
        out_name.rpartition("DSGA")[0].rstrip(".").rstrip("_")

    if not args.quiet:
        print('\n#--------------------------------------------#')
        print('#-------- Preparing Query Sequences ---------#')
        print('#------- Preparing Germline Sequences -------#')

    germ_filt_fasta = psg.prepare_fasta(
                            out_dir,
                            out_name,
                            germ_fasta,
                            min_germ)

    if not args.quiet:
        print('#------- Preparing Somatic Sequences --------#')

    soma_filt_fasta = psg.prepare_fasta(
                            out_dir,
                            out_name,
                            soma_fasta,
                            min_soma,
                            False)

    if not args.quiet:
        print('\n#--------------------------------------------#')
        print('#-------- Aligning Soma and Germline --------#')

    out_tsv = psg.blast_germ_soma(
                    out_dir,
                    out_name,
                    germ_filt_fasta,
                    soma_filt_fasta,
                    min_germ,
                    min_soma,
                    threads,
                    quiet)


    if not args.quiet:
        print('\n#------- Filtering BLASTN Alignments --------#')

    sg_summary, soma_germ_dict, skip_soma_germ = sgc.filter_blastn_hits(
                                                        out_tsv,
                                                        min_aln_prop,
                                                        min_pointer,
                                                        max_pointer,
                                                        multi_filt,
                                                        id_filt)

    if not args.quiet:
        print('\n#--------------------------------------------#')
        print('#----- Refining Soma-Germline Alignments ----#')
        print('#--------------------------------------------#')

    if sg_summary:
        if not args.quiet:
            print('#------- Generating Local Alignments --------#')

        nsc_align_dir = sgc.prep_nonscram_loci(
                            out_dir,
                            soma_filt_fasta,
                            germ_filt_fasta,
                            sg_summary,
                            soma_germ_dict)

        if not args.quiet:
            print('#-------- Refining Local Alignments ---------#')

        rfa.eval_local_alns(
                nsc_align_dir)

    sgc.save_summary_tsv(sg_summary, out_dir, out_name, multi_filt)


if __name__ == '__main__':

    args = ae.capture_args()
    # print(args)
    """
    NOTE TO SELF!!!

    Currently ONLY WORKS FOR TRANSCRIPTOME TO GERMLINE comparisons
      - Skip the intron inferences when working with SOMATIC GENOMES
      - Need to handle/remove the telomeric repeats from SOMATIC GENOME SEQUENCES

    Will need to update this later
    """

    eval_germ_soma_arch(
            args.out_dir,
            args.germ,
            args.soma,
            args.min_germ,
            args.min_soma,
            args.min_aln,
            args.min_id,
            args.min_pointer,
            args.max_pointer,
            True,
            args.threads,
            args.quiet)

    # if args.command == 'add-taxa':
    #     if args.quiet:
    #         print('\n#-------------------------------------------#')
    #         print('#------------ Adding New Taxon -------------#')
    #         print('#-------------------------------------------#')
    #
    #     add_new_taxa(args)
    #
    #
    # if args.command == 'msa':
    #     if args.quiet:
    #         print('\n#-------------------------------------------#')
    #         print('#------------- Generating MSAs -------------#')
    #         print('#-------------------------------------------#')
    #
    #     msa_processing(args)


    #
    # if args.command == 'contam':
    #     if args.quiet:
    #         print('\n#-------------------------------------------#')
    #         print('#----------- Assessing By-Catch ------------#')
    #         print('#-------------------------------------------#')
    #
    #     by_catch_assess(args)
    #
    # if args.command == 'genetic-code':
    #     if args.quiet:
    #         print('\n#-------------------------------------------#')
    #         print('#--------- Inferring Genetic Code ----------#')
    #         print('#-------------------------------------------#')
    #
    #     genetic_code_eval(args)

    # print(args)

    # if args.quiet:
    #     print('\nThanks for using PhyG!')

#
# if __name__ == '__main__':
#     try:
#         soma_fasta = sys.argv[1]
#         germ_fasta = sys.argv[2]
#         output_name = sys.argv[3]
#
#     except:
#         print('\nUsage:\n\n    python3 gsa_genomic_comparison.py [SOMA-FASTA] [GERM-FASTA] [OUTPUT-NAME]\n')
#         sys.exit()
#
