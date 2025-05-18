import pysam
from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from dscHiCtools.chunkFiles import chunk_bam
import time 
from itertools import islice
from dscHiCtools.chunkFiles import (
    get_n_lines, 
    chunk_reads
)
import gzip
import sys 
    
    
def splitBarcodeSam(
    input_sam,
    outdir,
    barcodes,
    indexes,
    index,
    header_f,
    prefix
):
    """
        sam.gz as input 
        sam is not sorted
    """
    n = 1
    t = time.time()
    outputBam_dict = {}
    for barcode in barcodes:
        output_sam = os.path.join(outdir, prefix + barcode + "_" + str(index) + ".sam.gz")
        outputBam_dict[barcode] = gzip.open(output_sam, "wb")
        with gzip.open(header_f, 'rt') as header:
            for line in header.readlines():
                outputBam_dict[barcode].write(line.encode())

    with gzip.open(input_sam, "rt") as textfile1:
        totallines = islice(
            textfile1, indexes[0] * 1, indexes[1] * 1
        )
        while True:
            try:
                map_info = next(totallines)

                if n % 1000000 == 0:
                    utils.eprint(
                        "[splitBarcode::] Processed 1,000,000 reads in {}. Total "
                        "reads: {:,} in child {}".format(
                            utils.secondsToText(time.time() - t), n, os.getpid()
                        )
                    )
                    sys.stdout.flush()
                    t = time.time()  
                n += 1

                cell_barcode = map_info.split("\t")[0].split(":")[0]
                if cell_barcode in barcodes:
                    outputBam_dict[cell_barcode].write(map_info.encode())
                else:
                    continue
            except StopIteration as e:
                break

    for barcode in barcodes:
        outputBam_dict[barcode].close()


def splitBarcodeBam(
    input_bam,
    outdir,
    barcodes,
    index,
    interval,
    prefix
):
    """
        map cell barcodes through intervals of input bam file (for multi-thread running)
        Args:
            input_bam: input bam file path
            outdir: output directory 
            barcode_tree: BKtree for reference barcode file
            index: specify different prefix for temparary sc bam
            interval: interval of bam files
            min_mismatch: min mismatch bases allowed for mapping cell barcode to reference
        Returns:
            tagged.bam
            mapped.mtx: mapped cell barcode - readsN
            markDuplicates.mtx: mapped cell barcode - readName - start - end 
            unmapped.mtx: unmapped cell barcode - readsN - sequence
    """
    ## count the processed reads 
    n = 1
    t = time.time()
    outputBam_dict = {}

    inputBam = pysam.AlignmentFile(input_bam, "rb")
    header = inputBam.header
    for barcode in barcodes:
        output_bam = os.path.join(outdir, prefix + barcode + "_" + str(index) + ".bam")
        outputBam_dict[barcode] = pysam.AlignmentFile(output_bam, "wb", header=header)
   
    for i in interval:
        for read in inputBam.fetch(i[0], i[1], i[2]):
            if n % 1000000 == 0:
                utils.eprint(
                    "[splitBarcode::] Processed 1,000,000 reads in {}. Total "
                    "reads: {:,} in child {}".format(
                        utils.secondsToText(time.time() - t), n, os.getpid()
                    )
                )
                sys.stdout.flush()
                t = time.time()          
            try:
                cell_barcode = read.get_tag("CB")
                n += 1
            except Exception as e:
                cell_barcode = read.qname.split(":")[0]
            if cell_barcode in barcodes:
                outputBam_dict[cell_barcode].write(read)
            else:
                continue
    for barcode in barcodes:
        outputBam_dict[barcode].close()
    inputBam.close()


@utils.log_info
def main(args):
    """
        split cell barcode to single bam 
        Args:
            input_bam: bwa-men generated Hi-C mapped mode bam
            outdir: output directory that stores single cell bam
            threads: >= 1
            barcode_f: barcode file 
            samtools: path to samtools
        Returns:
            cellbarcode.bam
    """
    
    barcodes = utils.readBarcode(args.barcode_file)
    threads = 10 if args.threads > 10 else args.threads
    prefix = utils.getPrefix(args.prefix)
    utils.eprint(f"[splitBarcode::] Prefix {prefix} will be added...")

    if args.input_bam[-3::] == "bam":

        inputBam = pysam.AlignmentFile(args.input_bam, "rb")
        idx_file = args.input_bam + ".bai"
        if not os.path.exists(idx_file):
            utils.eprint("[splitBarcode::] Indexing input bam file...")
            Args_m = [f'{args.samtools} index -@ {threads} {args.input_bam}']
            subprocess.check_call(Args_m, shell=True)

        intervals = chunk_bam(inputBam, threads) 
        inputBam.close()
        
        ## parallel
        if threads <= 1:
            splitBarcodeBam(
                input_bam=args.input_bam,
                outdir=args.outdir,
                barcodes=barcodes,
                index = 1,
                interval=intervals[1]
            )
            utils.eprint("[splitBarcode::] Cell barcodes are filtered.")
        else:
            i = 0
            utils.eprint(f"[splitBarcode::] Filtering cell barcodes is running with {threads} cores.")
            p = multiprocess.Pool(processes=threads)
            
            for interval in intervals.values():
                i += 1
                p.apply_async(
                    splitBarcodeBam,
                    args=(
                        args.input_bam,
                        args.outdir,
                        barcodes,
                        i,
                        interval
                    ),
                    error_callback=utils.print_error,
                )

            p.close()
            p.join()
            utils.eprint("[splitBarcode::] Filtering done.")
            utils.eprint("[splitBarcode::] Merging results...")
        
            ## merge bam
            for barcode in barcodes:
                output_bam = os.path.join(args.outdir, prefix + barcode + ".bam")
                tempFiles = os.path.join(args.outdir, prefix + barcode + "_*.bam")
                Args_m = [f'{args.samtools} merge -fpc -@ {threads} {output_bam} {tempFiles}']
                subprocess.check_call(Args_m, shell=True)

                Args_m = [f'rm {tempFiles}']
                subprocess.check_call(Args_m, shell=True)

    elif args.input_bam[-6::] == "sam.gz":
        if args.n_reads is not None:
            n_lines = args.n_reads
        else:
            n_lines = get_n_lines(args.input_bam, "sam")
        utils.eprint(f"[splitBarcode::] Total input sam has {n_lines} reads")

        ## parallel
        if threads <= 1:
            splitBarcodeSam(
                input_bam=args.input_bam,
                outdir=args.outdir,
                barcodes=barcodes,
                indexes=[0, n_lines],
                index = "",
                header_f=args.header,
                prefix=prefix
            )
            utils.eprint("[splitBarcode::] Cell barcodes are filtered.")
        else:
            chunk_indexes = chunk_reads(n_lines, threads)
            utils.eprint(f"[splitBarcode::] Filtering cell barcodes is running with {threads} cores.")
            p = multiprocess.Pool(processes=threads)
            utils.eprint(f"[splitBarcode::] Index are {chunk_indexes}")

            for i, indexes in enumerate(chunk_indexes):
                p.apply_async(
                    splitBarcodeSam,
                    args=(
                        args.input_bam,
                        args.outdir,
                        barcodes,
                        indexes,
                        str(i),
                        args.header,
                        prefix
                    ),
                    error_callback=utils.print_error,
                )

            p.close()
            p.join()
            utils.eprint("[splitBarcode::] Filtering done.")
            utils.eprint("[splitBarcode::] Merging results...")
        
            ## merge bam
            for barcode in barcodes:
                output_bam = os.path.join(args.outdir, prefix + barcode + ".sam.gz")
                tempFiles = os.path.join(args.outdir, prefix + barcode + "_*.sam.gz")
                Args_m = [f'{args.samtools} merge -fpc -@ {threads} {output_bam} {tempFiles}']
                subprocess.check_call(Args_m, shell=True)

                Args_m = [f'rm {tempFiles}']
                subprocess.check_call(Args_m, shell=True)

    else:
        raise("Input Format Error! (.bam or .sam.gz)")
   
