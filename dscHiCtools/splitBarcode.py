import pysam
from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from dscHiCtools.chunkFiles import chunk_bam
import time 


def splitBarcodeBam(
    input_bam,
    outdir,
    barcodes,
    index,
    interval
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
        output_bam = os.path.join(outdir, barcode + "_" + str(index) + ".bam")
        outputBam_dict[barcode] = pysam.AlignmentFile(output_bam, "wb", header=header)
   
    for i in interval:
        for read in inputBam.fetch(i[0], i[1], i[2]):
            if n % 1000000 == 0:
                utils.eprint(
                    "[splitBarcode::] Processed 1,000,000 reads in {}. Total "
                    "reads: {:,} in child {}".format(
                        secondsToText(time.time() - t), n, os.getpid()
                    )
                )
                sys.stdout.flush()
                t = time.time()          
            try:
                cell_barcode = read.get_tag("CB")
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
    inputBam = pysam.AlignmentFile(args.input_bam, "rb")
    idx_file = args.input_bam + ".bai"
    if not os.path.exists(idx_file):
        utils.eprint("[splitBarcode::] Indexing input bam file...")
        Args_m = [f'{args.samtools} index -@ {args.threads} {args.input_bam}']
        subprocess.check_call(Args_m, shell=True)

    intervals = chunk_bam(inputBam, args.threads)
    barcodes = utils.readBarcode(args.barcode_file)
    inputBam.close()
    ## parallel
    if args.threads <= 1:
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
        utils.eprint(f"[splitBarcode::] Filtering cell barcodes is running with {args.threads} cores.")
        p = multiprocess.Pool(processes=args.threads)
        
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
            output_bam = os.path.join(args.outdir, barcode + ".bam")
            tempFiles = os.path.join(args.outdir, barcode + "_*.bam")
            Args_m = [f'{args.samtools} merge -fpc -@ {args.threads} {output_bam} {tempFiles}']
            subprocess.check_call(Args_m, shell=True)

            Args_m = [f'rm {tempFiles}']
            subprocess.check_call(Args_m, shell=True)
