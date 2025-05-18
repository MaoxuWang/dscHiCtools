import sys
from argparse import ArgumentParser
import pkg_resources
from argparse import RawTextHelpFormatter
from dscHiCtools import (
    mapBarcode, tag2bam, filterBarcode, splitBarcode, sub2cool, splitContacts, scAB
)
import os 


def parse_args():
    version = pkg_resources.require("dscHiCtools")[0].version

    barcodes_dir = os.path.join(os.path.dirname(__file__), "barcodes")
    CpG_dir = os.path.join(os.path.dirname(__file__), "cpg_freq")

    print(f'{barcodes_dir}', file = sys.stderr)

    parser = ArgumentParser(
        prog="dscHiCtools",
        description=(
            "Tools for demultiplexed single-cell HiC data processing\n"
            "Author: Maoxu Wang, Xie Lab, Peking University\n"
            "Version: {}".format(version)
        ),
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + str(version)
    )
    subparsers = parser.add_subparsers(title="Subcommands")


    # mapBarcode
    parser_mapBarcode = subparsers.add_parser(
        "mapBarcode", description=(
            "Add cell barcodes from I2 to read name of R1 & R2 and mapp them to the given cellbarcode whitelist.\n"
            "Output: barcoded fastq file of R1 and R2. Cell barcodes will be converted into RNA cell barcodes in 10X multiome."
        ),
        formatter_class=RawTextHelpFormatter
    )

    parser_mapBarcode.add_argument(
                        "--read1",
                        required=True,
                        help="read1 input fastq format file containing barcode reads (.gz)")
    parser_mapBarcode.add_argument(
                        "--index",
                        required=True,
                        help="fastq format file containing barcode reads (.gz)")
    parser_mapBarcode.add_argument("--read2",
                        required=True,
                        help="read2 input fastq format file containing barcode reads (.gz)")
    parser_mapBarcode.add_argument(
                        "-o", 
                        "--outdir",
                        type=str,
                        required=True,
                        help="The output directory for stat and barcoded fastq file. stat: cellbarcode_count.txt; fastq: barcoded.fastq.gz ")
    parser_mapBarcode.add_argument(
                        "--sampleName",
                        type=str,
                        default="out",
                        help="sample name")
    parser_mapBarcode.add_argument(
                        "--threads",
                        type=int,
                        default=10,
                        help="threads")
    parser_mapBarcode.add_argument(
                        "--start_idx",
                        type=int,
                        default=0,
                        help="start index of cellbarcode in index read")
    parser_mapBarcode.add_argument(
                        "--end_idx",
                        type=int,
                        default=16,
                        help="end index of cellbarcode in index read")
    parser_mapBarcode.add_argument(
                        "--min_mismatch",
                        type=int,
                        default=1,
                        help="The minimum mismatch between fastq barcode and ref barcode")
    parser_mapBarcode.add_argument(
                        "--n_lines",
                        type=int,
                        default=-1,
                        help="Total lines of input fastq file to be processed. Default: all (-1)")
    parser_mapBarcode.add_argument(
                        "--ref_barcode",
                        type=str,
                        default=os.path.join(barcodes_dir, "737K-cratac-v1.txt.gz"),
                        help="The 10X reference barcode file")
    parser_mapBarcode.add_argument(
                        "--multiome_RNA_barcode",
                        type=str,
                        default=os.path.join(barcodes_dir, "737K-arc-RNA-v1.txt.gz"),
                        help="The 10X multiome reference RNA barcode file")
    parser_mapBarcode.add_argument(
                        "--multiome_DNA_barcode",
                        type=str,
                        default=os.path.join(barcodes_dir, "737K-arc-DNA-v1.txt.gz"),
                        help="The 10X multiome reference DNA barcode file")
    parser_mapBarcode.add_argument(
                    "--chemistry",
                    type=str,
                    default="atac",
                    choices=["atac", "multiome"],
                    help="library type")   
    parser_mapBarcode.add_argument(
                        "--suffix",
                        type=str,
                        help="identified added as suffix of cell barcode. default: None")
    parser_mapBarcode.set_defaults(func=mapBarcode.main)

    # addTags2bam
    parser_addTag = subparsers.add_parser(
        "addTag", description=(
            "Add cell barcode tags to bam file for 10X dscHi-C data\n"
            "Output: bam file with CB tags"
        ),
        formatter_class=RawTextHelpFormatter
    )

    parser_addTag.add_argument("--input_bam",
                        required=True,
                        help="input bam (.bam)")
    parser_addTag.add_argument("--output_bam",
                        required=True,
                        help="output_bam bam (.bam)")
    parser_addTag.add_argument("--threads",
                        type=int,
                        default=10,
                        help="threads")
    parser_addTag.add_argument(
                        "--samtools",
                        type=str,
                        default="samtools",
                        help="samtools executable path")
    parser_addTag.set_defaults(func=tag2bam.main)


    # filterBarcode
    parser_filterBarcode = subparsers.add_parser(
        "filterBarcode", description=(
            "Filter cell barcodes, only save barcodes listed to bam to file.\n"
            "Output: bam file with listed barcodes"
        ),
        formatter_class=RawTextHelpFormatter
    )

    parser_filterBarcode.add_argument("--input_bam",
                        required=True,
                        help="input bam (.bam)")
    parser_filterBarcode.add_argument("--output_bam",
                        required=True,
                        help="output_bam bam (.bam)")
    parser_filterBarcode.add_argument(
                        "--barcode_file",
                        type=str,
                        required=True,
                        help="barcodes to save (.tsv or .tsv.gz)")
    parser_filterBarcode.add_argument(
                        "--threads",
                        type=int,
                        default=10,
                        help="threads")
    parser_filterBarcode.add_argument(
                        "--samtools",
                        type=str,
                        default="samtools",
                        help="samtools executable path")
    parser_filterBarcode.set_defaults(func=filterBarcode.main)

    # splitBarcode
    parser_splitBarcode = subparsers.add_parser(
        "splitBarcode", description=(
            "split bam into different files by cell barcodes.\n"
            "Output: single cell bam"
        ),
        formatter_class=RawTextHelpFormatter
    )

    parser_splitBarcode.add_argument("--input_bam",
                        required=True,
                        help="input bam (.bam)")
    parser_splitBarcode.add_argument("--outdir",
                        required=True,
                        help="output directory")
    parser_splitBarcode.add_argument("--n_reads",
                        type=int,
                        help="number of reads in sam.gz to process")
    parser_splitBarcode.add_argument(
                        "--barcode_file",
                        type=str,
                        required=True,
                        help="barcodes to save (.tsv or .tsv.gz)")
    parser_splitBarcode.add_argument(
                        "--prefix",
                        type=str,
                        help="Identifier appended to output name")
    parser_splitBarcode.add_argument(
                        "--threads",
                        type=int,
                        default=10,
                        help="threads")
    parser_splitBarcode.add_argument(
                        "--samtools",
                        type=str,
                        default="samtools",
                        help="samtools executable path")
    parser_splitBarcode.add_argument(
                        "--header",
                        type=str,
                        help="header for sam.gz")
    parser_splitBarcode.set_defaults(func=splitBarcode.main)


    # sub2cool
    parser_sub2cool = subparsers.add_parser(
        "sub2cool", description=(
            "extract sub-cluster cells and convert to .cool file\n"
            "Output: .cool file"
        ),
        formatter_class=RawTextHelpFormatter
    )

    parser_sub2cool.add_argument("--input_pairs",
                        required=True,
                        help="input contacts.pairs.txt(.gz)")
    parser_sub2cool.add_argument("--outdir",
                        required=True,
                        help="output directory")
    parser_sub2cool.add_argument(
                        "--metadata",
                        type=str,
                        required=True,
                        help="file that specifies barcodes and sub-group")
    parser_sub2cool.add_argument(
                        "--threads",
                        type=int,
                        default=10,
                        help="threads")
    parser_sub2cool.add_argument(
                        "--resolution",
                        nargs='+', 
                        default=["1000", "2000", "5000", "10000", "25000", "50000", "100000", "250000", "500000", "1000000", "2500000"],
                        help="cooler resolution sperated by ' '")
    parser_sub2cool.add_argument(
                        "--cooler",
                        type=str,
                        default="/share/home/wangmx/anaconda3/envs/cooler/bin/cooler",
                        help="cooler executable path")
    parser_sub2cool.add_argument(
                        "--hic2cool",
                        type=str,
                        default="/share/home/wangmx/anaconda3/envs/cooler/bin/hic2cool",
                        help="hic2cool executable path")
    parser_sub2cool.add_argument(
                        "--juicer",
                        type=str,
                        default="java -Xmx256G -jar /share/home/wangmx/software/juicer/scripts/common/juicer_tools_1.22.01.jar",
                        help="juicer executable path")

    parser_sub2cool.add_argument(
                        "--species",
                        type=str,
                        required=True,
                        help="mm10 or hg38")
    parser_sub2cool.add_argument(
                        "-n", "--normalize",
                        action='store_true',
                        help="if set, normalize the matrices")
    parser_sub2cool.add_argument(
                        "--suffix",
                        type=str,
                        default=".contacts.pairs.txt.gz",
                        help="suffix of contact pairs file name")
    parser_sub2cool.set_defaults(func=sub2cool.main)

    # scAB value calculation
    parser_scAB = subparsers.add_parser(
        "scAB", description=(
            "calculate scA/B matrix for dscHiC data\n"
            "Output: matrix file(.txt)"
        ),
        formatter_class=RawTextHelpFormatter
    )
    parser_scAB.add_argument("--input_contact_pairs",
                        required=True,
                        help="contact pairs.txt(.gz), with cell barcode information")
    parser_scAB.add_argument("--output",
                        required=True,
                        help="output matrix file (.txt)")
    parser_scAB.add_argument("--ref_CpG",
                        default=os.path.join(CpG_dir, "hg38.cpg.1m.txt"),
                        help="reference CpG frequency file generated from genome.fa (default: hg38.cpg.1m.txt)")
    parser_scAB.add_argument("--res",
                        type=int,
                        default=1000000,
                        help="resolution of bin")
    parser_scAB.add_argument("--cores",
                        type=int,
                        default=10,
                        help="multi-cores used for running")
    parser_scAB.add_argument("--barcode_file",
                        help="barcodes to save (.tsv or .tsv.gz)")
    parser_scAB.add_argument("--weighted",
                        type=bool,
                        default=True,
                        help="whether use weighted mean rather than mean")
    parser_scAB.add_argument("--input_format",
                        type=str,
                        choices=["scVI-3d", "raw", "higashi_imputed"],
                        default="raw",
                        help="Input format specified (default: raw e.g. as Tan 2018 described)")
    parser_scAB.add_argument("--cellbarcode_dict",
                        type=str,
                        help="convert 0,1,2... from higashi format to cellbarcode")
    parser_scAB.add_argument("--cis_only",
                        action='store_true',
                        help="if set, only use cis contacts")
    parser_scAB.set_defaults(func=scAB.main)

    # splitContacts
    parser_splitContacts = subparsers.add_parser(
        "splitContacts", description=(
            "split contacts pairs file into different files by cell barcodes.\n"
            "Output: single cell contacts pairs file"
        ),
        formatter_class=RawTextHelpFormatter
    )

    parser_splitContacts.add_argument("--input_contacts",
                        required=True,
                        help="input contacts pairs txt(.gz) or directory for h5py file")
    parser_splitContacts.add_argument("--outdir",
                        required=True,
                        help="output directory")
    parser_splitContacts.add_argument(
                        "--threads",
                        type=int,
                        default=10,
                        help="threads")
    parser_splitContacts.add_argument(
                        "--n_lines",
                        type=int,
                        default=-1,
                        help="Total lines to process. default: all lines")
    parser_splitContacts.add_argument(
                        "--format",
                        type=str,
                        choices=["pairs", "h5py"],
                        default="pairs",
                        help="Pairs of h5py file to split. default: pair file format")
    parser_splitContacts.add_argument(
                        "--id_convert_file",
                        type=str,
                        help="Cell id to cell barcode (higashi required)")
    parser_splitContacts.add_argument(
                        "--header",
                        action='store_true',
                        help="Whether to output pairs header. default: False")
    parser_splitContacts.add_argument(
                        "--neighbour",
                        type=int,
                        default=5,
                        help="whether to use neighbor information to impute (higashi parameters)")
    parser_splitContacts.add_argument("--barcode_file",
                        help="barcodes to save (.tsv or .tsv.gz) default: save all")
    parser_splitContacts.add_argument("--resolution",
                        default=1000000,
                        type=int,
                        help="resolution for each bin default: 1Mb")
    parser_splitContacts.set_defaults(func=splitContacts.main)

    return parser


def main():
    parser = parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    else:
        options = parser.parse_args()
        options.func(options)

if __name__ == '__main__':
    main()
