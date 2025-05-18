import pybktree
import Levenshtein
from dscHiCtools import utils 
import os 
import sys 
import time 
from collections import Counter
import gzip 
import multiprocess
import subprocess
from itertools import islice
from dscHiCtools.write import (
    write_unmapped, 
    write_mapped
)
from dscHiCtools.chunkFiles import (
    get_n_lines, 
    chunk_reads
)

## mapping:     
def makeBKtree(barcode_f, RNA_to_DNA = None):
    """
        make the BKtree for reference barcode file;
        This is for fast matching
        Returns:
            a BKtree
    """
    if RNA_to_DNA is None:
        barcode_tree = pybktree.BKTree(Levenshtein.hamming, barcode_f)
    else:
        DNA_barcode_f = []
        for barcode in barcode_f:
            DNA_barcode_f.append(RNA_to_DNA.get(barcode))
        barcode_tree = pybktree.BKTree(Levenshtein.hamming, DNA_barcode_f)

    utils.eprint("[mapBarcode::] Generated barcode tree from reference barcode file")
    return(barcode_tree)


def mapCellBarcode(
    read1,
    index,
    read2,
    start_idx,
    end_idx,
    indexes,
    out_r1,
    out_r2,
    reverse,
    barcode_tree,
    min_mismatch,
    suffix,
    chemistry,
    DNA_to_RNA
):
    """
        Add cell barcodes to reads through indexes of input file (for multi-thread running)
        Args:
            input: input R1.fastq R2.fastq R3.fastq file containing barcode reads
            index: Pair of first and last index for islice
            out_r1: output ; read name added by cell barcode
            out_r2
            barcode_tree: BKtree for reference barcode file
            min_mismatch: min mismatch bases allowed for mapping cell barcode to reference
            suffix: additional identifier
            chemistry: library type
            DNA_to_RNA: convert mapped DNA barcodes to RNA barcodes
        Returns:
            mapped.mtx: mapped cell barcode - readsN
            unmapped.mtx: unmapped cell barcode - readsN - sequence
            barcoded reads
    """
    ## count the processed reads 
    n = 1
    t = time.time()
    mapped = Counter()
    unmapped = Counter()

    ## must output to a .gz file
    out_r1_file = gzip.open(out_r1, "wb")
    out_r2_file = gzip.open(out_r2, "wb")

    with gzip.open(read1, "rt") as textfile1, gzip.open(index, "rt") as textfile2, gzip.open(
        read2, "rt"
    ) as textfile3:
        totallines = islice(
            zip(textfile1, textfile2, textfile3), indexes[0] * 4, indexes[1] * 4
        )
        while True:
            try:
                readNames = next(totallines)
                seqs = next(totallines)
                _char = next(totallines)
                qualities = next(totallines)

                # Progress info
                if n % 10000000 == 0:
                    utils.eprint(
                        "[mapBarcode::] Processed 10,000,000 reads in {}. Total "
                        "reads: {:,} in child {}".format(
                            utils.secondsToText(time.time() - t), n, os.getpid()
                        )
                    )
                    sys.stdout.flush()
                    t = time.time()
                n += 1
                read = seqs[1].strip()
                if reverse:
                    cell_barcode_string = utils.reverseRead(read[start_idx:end_idx])
                else:
                    cell_barcode_string = read[start_idx:end_idx]
                if len(cell_barcode_string) != 16:
                    utils.eprint("[mapBarcode::] wrong cell barcode length found!")
                    continue

                ### unmapped
                cell_barcode = utils.alignReference(cell_barcode_string, barcode_tree, min_mismatch)
                if( not cell_barcode):
                    unmapped[cell_barcode_string] = unmapped.get(cell_barcode_string, 0) + 1 
                    continue
                ### mapped
                else:
                    if chemistry == "multiome":
                        cell_barcode = DNA_to_RNA.get(cell_barcode)
                    mapped[cell_barcode] = mapped.get(cell_barcode, 0) + 1 

                    if suffix is not None:
                        cell_barcode = cell_barcode + "-" + suffix
                    utils.writeFile(out_r1_file, 
                        cell_barcode,
                        readNames,
                        seqs,
                        _char,
                        qualities,
                        0)
                    utils.writeFile(out_r2_file, 
                        cell_barcode,
                        readNames,
                        seqs,
                        _char,
                        qualities,
                        2)
                   
            except StopIteration as e:
                break
    
    out_r1_file.close()
    out_r2_file.close()
    utils.eprint(
        "[mapBarcode::] Mapping done for process {}. Processed {:,} reads".format(os.getpid(), n - 1)
    )
    sys.stdout.flush()
    return mapped, unmapped


@utils.log_info
def main(args):
    """
        add cell barcode to read name of R1, R3
        Args:
            read1: read1 file from 10X data (.fastq.gz)
            index: index file from 10X data (.fastq.gz)
            read2: read2 file from 10X data (.fastq.gz)
            outdir: directory for output
            threads: >= 1
            sampleName: output prefix that specifys sample
            ref_barcode: reference cell barcode file (.txt) 
            min_mismatch: minimal mismatch allowed for cell barcode mapping
            chemistry: forward, reverse, auto-detect (default: auto-detect)
            suffix: additional identifier 
            n_lines: file line count
            multiome_RNA_barcode: mapping multiome DNA barcodes to RNA barcodes 
            multiome_DNA_barcode: mapping multiome DNA barcodes to RNA barcodes 
        Returns:
            barcoded.fastq.gz (using RNA barcode)
            QC/ mapped.tsv.gz: mapped cell barcode - readsN
            QC/ unmapped.tsv.gz: unmapped cell barcode - readsN - sequence
    """
    if args.n_lines is None or args.n_lines <= 0 :
        args.n_lines = get_n_lines(args.read1)
    n_reads = int(args.n_lines / 4)
    utils.eprint("[mapBarcode::] Start adding cell barcodes...")
    utils.eprint(f"[mapBarcode::] Total input fastq files have {args.n_lines} lines.\nProcessing {n_reads} reads")

    RNA_to_DNA = None 
    DNA_to_RNA = None 
    if args.chemistry == "multiome":
        utils.eprint(f"[mapBarcode::] Chemistry is set to be {args.chemistry}")
        utils.eprint(f"[mapBarcode::] Creating DNA <-> RNA barcodes dictionary...")

        RNA_to_DNA, DNA_to_RNA = utils.multiomeBarcodesDict(args.multiome_RNA_barcode, args.multiome_DNA_barcode)
        ## create barcde tree from RNA barcodes (all from arc)
        RNA_barcodes = utils.readBarcode(args.multiome_RNA_barcode)
        barcode_tree = makeBKtree(RNA_barcodes, RNA_to_DNA)
    else:
        DNA_barcodes = utils.readBarcode(args.ref_barcode)
        barcode_tree = makeBKtree(DNA_barcodes)
    if not os.path.isfile(args.read1) :
        raise Exception(f"{args.read1} not found !")

    if not os.path.isfile(args.index) :
        raise Exception(f"{args.index} not found !")

    if not os.path.isfile(args.read2) :
        raise Exception(f"{args.read2} not found !")

    ## matrix out to QC subdirectory
    outfolder = os.path.join(os.path.dirname(args.outdir), "QC")

    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    ## detect strand orientation
    reverse = utils.autoDetectChemistry(
            args.index,
            barcode_tree,
            args.min_mismatch
        )

    if args.threads <= 1:
        out_r1 = os.path.join(args.outdir, args.sampleName + ".R1.barcoded.fastq.gz")
        out_r2 = os.path.join(args.outdir, args.sampleName + ".R2.barcoded.fastq.gz")

        mappedMtx, unmappedMtx = mapCellBarcode(
            read1=args.read1,
            index=args.index,
            read2=args.read2,
            start_idx=args.start_idx,
            end_idx=args.end_idx,
            indexes=[0, n_reads],
            out_r1=out_r1,
            out_r2=out_r2,
            reverse=reverse,
            barcode_tree=barcode_tree,
            min_mismatch=args.min_mismatch,
            suffix=args.suffix,
            chemistry=args.chemistry,
            DNA_to_RNA=DNA_to_RNA
        )
        utils.eprint("[mapBarcode::] Cell barcodes are added.")
    else:
        utils.eprint(f"[mapBarcode::] Adding cell barcodes is running with {args.threads} cores.")
        p = multiprocess.Pool(processes=args.threads)
        chunk_indexes = chunk_reads(n_reads, args.threads)

        i = 0
        r1_single_files = []
        r2_single_files = []
        parallel_results = []

        ## multi-processing
        for indexes in chunk_indexes:
            i += 1
            out_r1 = os.path.join(args.outdir, args.sampleName + "_" + str(i) + ".R1.barcoded.fastq.gz")
            out_r2 = os.path.join(args.outdir, args.sampleName + "_" + str(i) + ".R2.barcoded.fastq.gz")
            r1_single_files.append(out_r1)
            r2_single_files.append(out_r2)

            p.apply_async(
                mapCellBarcode,
                args=(
                    args.read1,
                    args.index,
                    args.read2,
                    args.start_idx,
                    args.end_idx,
                    indexes,
                    out_r1,
                    out_r2,
                    reverse,
                    barcode_tree,
                    args.min_mismatch,
                    args.suffix,
                    args.chemistry,
                    DNA_to_RNA
                ),
                error_callback=utils.print_error,
                callback=parallel_results.append,
            )

        p.close()
        p.join()
        utils.eprint("[mapBarcode::] Merging results...")
        mappedMtx, unmappedMtx = utils.merge_results(parallel_results)

        ## merge fastq.gz subprocess
        out_final_r1 = os.path.join(args.outdir, args.sampleName + ".R1.barcoded.fastq.gz")
        out_final_r2 = os.path.join(args.outdir, args.sampleName + ".R2.barcoded.fastq.gz")

        Args_m = [f'cat {" ".join(r1_single_files)} > {out_final_r1}']
        subprocess.check_call(Args_m, shell=True)

        Args_m = [f'cat {" ".join(r2_single_files)} > {out_final_r2}']
        subprocess.check_call(Args_m, shell=True)

        ## remove temporary files
        for my_file in r1_single_files:
            if os.path.exists(my_file):
                os.remove(my_file)
        for my_file in r2_single_files:
            if os.path.exists(my_file):
                os.remove(my_file)
        utils.eprint("[mapBarcode::] Cell barcodes are added and mapped!")

    ## save to stat
    unmapped_file = args.sampleName + ".unmapped.tsv.gz"
    mapped_file = args.sampleName + ".mapped.tsv.gz"

    with open(os.path.join(outfolder, args.sampleName + ".total_readsN.txt"), 'w') as read_f:
        read_f.write(str(n_reads))

    ### mapped barcode
    write_mapped(mappedMtx, outfolder, mapped_file)

    ### unmapped barcode
    write_unmapped(
        merged_no_match=unmappedMtx,
        outfolder=outfolder,
        filename=unmapped_file,
    )