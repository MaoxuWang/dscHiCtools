from dscHiCtools import cli
import sys
from argparse import ArgumentParser
import pkg_resources
from argparse import RawTextHelpFormatter

version = pkg_resources.require("dscHiCtools")[0].version
parser = ArgumentParser(
    prog="dscHiC-Tools",
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
                    default="/share/home/wangmx/software/pipelines/HT-HiC/src_1.0/data/737K-cratac-v1.txt.gz",
                    help="The 10X reference barcode file")
parser_mapBarcode.add_argument(
                    "--multiome_RNA_barcode",
                    type=str,
                    default="/share/home/wangmx/software/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz",
                    help="The 10X multiome reference RNA barcode file")
parser_mapBarcode.add_argument(
                    "--multiome_DNA_barcode",
                    type=str,
                    default="/share/home/wangmx/software/cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz",
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
parser_mapBarcode.set_defaults(func=cli.run_mapBarcode)


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
parser_addTag.set_defaults(func=cli.run_addTag)


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
parser_filterBarcode.set_defaults(func=cli.run_filterBarcode)

# filterBarcode
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
parser_splitBarcode.add_argument(
                    "--barcode_file",
                    type=str,
                    required=True,
                    help="barcodes to save (.tsv or .tsv.gz)")
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
parser_splitBarcode.set_defaults(func=cli.run_splitBarcode)


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
                    help="cooler resolution")
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
                    "-n",
                    action='store_true',
                    help="Don't normalize the matrices")
parser_sub2cool.set_defaults(func=cli.run_sub2cool)


def main():
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    else:
        options = parser.parse_args()
        sys.setrecursionlimit(200000)
        options.func(options)

if __name__ == "__main__":
    main()
