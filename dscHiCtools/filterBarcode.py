import pysam
from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from dscHiCtools.chunkFiles import chunk_bam
import time 


def filterBarcodeBam(
    input_bam,
    output_bam,
    barcodes,
    interval
):
    """
        map cell barcodes through intervals of input bam file (for multi-thread running)
        Args:
            input_bam: input bam file path
            output_bam: output bam file path
            barcode_tree: BKtree for reference barcode file
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

    inputBam = pysam.AlignmentFile(input_bam, "rb")
    header = inputBam.header
    outputBam = pysam.AlignmentFile(output_bam, "wb", header=header)
   
    for i in interval:
        for read in inputBam.fetch(i[0], i[1], i[2]):
            if n % 1000000 == 0:
                utils.eprint(
                    "[filterBarcode::] Processed 1,000,000 reads in {}. Total "
                    "reads: {:,} in child {}".format(
                        utils.secondsToText(time.time() - t), n, os.getpid()
                    )
                )
                sys.stdout.flush()
                t = time.time() 
            try:
                cell_barcode = read.get_tag("CB")
            except Exception as e:
                cell_barcode = read.qname.split(":")[0]

            if cell_barcode in barcodes:
                outputBam.write(read)
            else:
                continue
    outputBam.close()
    inputBam.close()


@utils.log_info
def main(args):
    """
        count cell barcode, add cell barcode to read name
        Args:
            input_bam: bwa-men generated Hi-C mapped mode bam
            output_bam: cell barcode mapped bam with CB tag (.cellbarcoded.bam)
            threads: >= 1
            ref_barcode: reference cell barcode file (.txt) 
            min_mismatch: minimal mismatch allowed for cell barcode mapping
            outMarkDuplicates: whether to print cellbarcode-fragment infomation to mark duplicates
        Returns:
            tagged.bam
            mapped.mtx: mapped cell barcode - readsN
            markDuplicates.mtx: mapped cell barcode - readName - start - end 
            unmapped.mtx: unmapped cell barcode - readsN - sequence
            if outMarkDuplicates true: bam.markduplicates.txt.gz
    """
    inputBam = pysam.AlignmentFile(args.input_bam, "rb")
    idx_file = args.input_bam + ".bai"
    if not os.path.exists(idx_file):
        utils.eprint("[filterBarcode::] Indexing input bam file...")
        Args_m = [f'{args.samtools} index -@ {args.threads} {args.input_bam}']
        subprocess.check_call(Args_m, shell=True)
    intervals = chunk_bam(inputBam, args.threads)

    barcodes = utils.readBarcode(args.barcode_file)
    inputBam.close()
    ## parallel
    if args.threads <= 1:
        filterBarcodeBam(
            input_bam=args.input_bam,
            output_bam=args.output_bam,
            barcodes=barcodes,
            interval=intervals[1]
        )
        utils.eprint("[filterBarcode::] Cell barcodes are filtered.")
    else:
        i = 0
        utils.eprint(f"[filterBarcode::] Filtering cell barcodes is running with {args.threads} cores.")
        p = multiprocess.Pool(processes=args.threads)
        
        bamTempFiles = []
        for interval in intervals.values():
            i += 1
            output_bam_temp = args.output_bam + str(i) + ".bam"
            bamTempFiles.append(output_bam_temp)
            p.apply_async(
                filterBarcodeBam,
                args=(
                    args.input_bam,
                    output_bam_temp,
                    barcodes,
                    interval
                ),
                error_callback=utils.print_error,
            )

        p.close()
        p.join()
        utils.eprint("[filterBarcode::] Filtering done.")
        utils.eprint("[filterBarcode::] Merging results...")
    
        ## merge bam
        tempFiles = " ".join(bamTempFiles)
        Args_m = [f'{args.samtools} merge -fpc -@ {args.threads} {args.output_bam} {tempFiles}']
        subprocess.check_call(Args_m, shell=True)

        if os.path.exists(args.output_bam):
            [os.remove(i) for i in bamTempFiles]
        else:
            raise Exception("[filterBarcode::] samtools merge failed, temp files not deleted")