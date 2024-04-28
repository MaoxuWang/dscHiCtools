import os
from dscHiCtools.utils import *
from dscHiCtools.chunkFiles import *
from dscHiCtools.write import *
import subprocess
import sys
import time
import multiprocess
import functools
import pysam
import re

def mapBarcode(
    read1,
    index,
    read2,
    start_idx,
    end_idx,
    outdir,
    threads,
    sampleName,
    ref_barcode,
    min_mismatch,
    chemistry,
    suffix,
    n_lines,
    multiome_RNA_barcode,
    multiome_DNA_barcode
):
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
    if n_lines is None or n_lines <= 0 :
        n_lines = get_n_lines(read1)
    n_reads = int(n_lines / 4)
    print("Start adding cell barcodes...")
    print(f"Total input fastq files have {n_lines} lines.\nProcessing {n_reads} reads")

    RNA_to_DNA = None 
    DNA_to_RNA = None 
    if chemistry == "multiome":
        print(f"Chemistry is set to be {chemistry}")
        print(f"Creating DNA <-> RNA barcodes dictionary...")

        RNA_to_DNA, DNA_to_RNA = multiomeBarcodesDict(multiome_RNA_barcode, multiome_DNA_barcode)
        ## create barcde tree
        RNA_barcodes = readBarcode(ref_barcode)
        barcode_tree = makeBKtree(RNA_barcodes, RNA_to_DNA)
    else:
        DNA_barcodes = readBarcode(ref_barcode)
        barcode_tree = makeBKtree(DNA_barcodes)
    if not os.path.isfile(read1) :
        raise Exception(f"{read1} not found !")

    if not os.path.isfile(index) :
        raise Exception(f"{index} not found !")

    if not os.path.isfile(read2) :
        raise Exception(f"{read2} not found !")

    ## matrix out to QC subdirectory
    outfolder = os.path.join(os.path.dirname(outdir), "QC")

    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    ## detect strand orientation
    reverse = autoDetectChemistry(
            index,
            barcode_tree,
            min_mismatch
        )

    if threads <= 1:
        out_r1 = outdir + "/" + sampleName + ".R1.barcoded.fastq.gz"
        out_r2 = outdir + "/" + sampleName + ".R2.barcoded.fastq.gz"

        mappedMtx, unmappedMtx = mapCellBarcode(
            read1=read1,
            index=index,
            read2=read2,
            start_idx=start_idx,
            end_idx=end_idx,
            indexes=[0, n_reads],
            out_r1=out_r1,
            out_r2=out_r2,
            reverse=reverse,
            barcode_tree=barcode_tree,
            min_mismatch=min_mismatch,
            suffix=suffix,
            chemistry=chemistry,
            DNA_to_RNA=DNA_to_RNA
        )
        print("Cell barcodes are added.")
    else:
        print(f"Adding cell barcodes is running with {threads} cores.")
        p = multiprocess.Pool(processes=threads)
        chunk_indexes = chunk_reads(n_reads, threads)

        i = 0
        r1_single_files = []
        r2_single_files = []
        parallel_results = []

        ## multi-processing
        for indexes in chunk_indexes:
            i += 1
            out_r1 = outdir + "/" + sampleName + "_" + str(i) + ".R1.barcoded.fastq.gz"
            out_r2 = outdir + "/" + sampleName + "_" + str(i) + ".R2.barcoded.fastq.gz"
            r1_single_files.append(out_r1)
            r2_single_files.append(out_r2)

            p.apply_async(
                mapCellBarcode,
                args=(
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
                ),
                error_callback=print_error,
                callback=parallel_results.append,
            )

        p.close()
        p.join()
        print("Merging results...")
        mappedMtx, unmappedMtx = merge_results(parallel_results)

        ## merge fastq.gz subprocess
        out_final_r1 = outdir + "/" + sampleName + ".R1.barcoded.fastq.gz"
        out_final_r2 = outdir + "/" + sampleName + ".R2.barcoded.fastq.gz"

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
        print("Cell barcodes are added and mapped!")

    ## save to stat
    unmapped_file = sampleName + ".unmapped.tsv.gz"
    mapped_file = sampleName + ".mapped.tsv.gz"

    ### unmapped barcode
    write_unmapped(
        merged_no_match=unmappedMtx,
        outfolder=outfolder,
        filename=unmapped_file,
    )

    ### mapped barcode
    write_mapped(mappedMtx, outfolder, mapped_file)



def addTags2bam(
    input_bam,
    output_bam,
    threads,
    samtools
):
    """
        add cell barcode tags to bam file
        input:
            bam file, the header of read name is barcode 
        output:
            bam file, with "CB:" tag added
    """
    inputBam = pysam.AlignmentFile(input_bam, "rb")
    intervals = chunk_bam(inputBam, threads)

    # assert inputBam.check_index(),  "ErrorType: %s file does not have index file." % input_bam
    if inputBam.check_index() is not True:
        Args_m = [f'{samtools} index -@ {threads} {inputBam}']
        subprocess.check_call(Args_m, shell=True)

    inputBam.close()
    p = multiprocess.Pool(threads)

    print("Adding CB tags to bam file...")
    print(f"Multi-cores: {threads}")
    i = 0
    bamTempFiles = []
    for interval in intervals.values():
        i += 1
        output_bam_temp = output_bam + str(i) + ".bam"
        bamTempFiles.append(output_bam_temp)
        p.apply_async(
            addCBtag,
            args=(
                interval,
                input_bam,
                output_bam_temp
            ),
            error_callback=print_error,
        )

    p.close()
    p.join()
    print("Filtering done.")
    print("Merging results...")
    
    tempFiles = " ".join(bamTempFiles)
    Args_m = [f'{samtools} merge -fpc -@ {threads} {output_bam} {tempFiles}']
    subprocess.check_call(Args_m, shell=True)

    if os.path.exists(output_bam):
        [os.remove(i) for i in bamTempFiles]
    else:
        raise Exception("samtools merge failed, temp files not deleted")


def filterBarcode(
    input_bam,
    output_bam,
    threads,
    barcode_f,
    samtools
):
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
    inputBam = pysam.AlignmentFile(input_bam, "rb")
    intervals = chunk_bam(inputBam, threads)

    if (not inputBam.check_index()):
        Args_m = [f'{samtools} index -@ {threads} {input_bam}']
        subprocess.check_call(Args_m, shell=True)

    # assert inputBam.check_index(),  "ErrorType: %s file does not have index file." % input_bam
    if inputBam.check_index() is not True:
        Args_m = [f'{samtools} index -@ {threads} {inputBam}']
        subprocess.check_call(Args_m, shell=True)

    barcodes = readBarcode(barcode_f)
    inputBam.close()
    ## parallel
    if threads <= 1:
        filterBarcodeBam(
            input_bam=input_bam,
            output_bam=output_bam,
            barcodes=barcodes,
            interval=intervals[1]
        )
        print("Cell barcodes are filtered.")
    else:
        i = 0
        print(f"Filtering cell barcodes is running with {threads} cores.")
        p = multiprocess.Pool(processes=threads)
        
        bamTempFiles = []
        for interval in intervals.values():
            i += 1
            output_bam_temp = output_bam + str(i) + ".bam"
            bamTempFiles.append(output_bam_temp)
            p.apply_async(
                filterBarcodeBam,
                args=(
                    input_bam,
                    output_bam_temp,
                    barcodes,
                    interval
                ),
                error_callback=print_error,
            )

        p.close()
        p.join()
        print("Filtering done.")
        print("Merging results...")
    
        ## merge bam
        tempFiles = " ".join(bamTempFiles)
        Args_m = [f'{samtools} merge -fpc -@ {threads} {output_bam} {tempFiles}']
        subprocess.check_call(Args_m, shell=True)

        if os.path.exists(output_bam):
            [os.remove(i) for i in bamTempFiles]
        else:
            raise Exception("samtools merge failed, temp files not deleted")


def splitBarcode(
    input_bam,
    outdir,
    threads,
    barcode_f,
    samtools
):
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
    inputBam = pysam.AlignmentFile(input_bam, "rb")
    intervals = chunk_bam(inputBam, threads)

    if (not inputBam.check_index()):
        Args_m = [f'{samtools} index -@ {threads} {input_bam}']
        subprocess.check_call(Args_m, shell=True)

    if inputBam.check_index() is not True:
        Args_m = [f'{samtools} index -@ {threads} {inputBam}']
        subprocess.check_call(Args_m, shell=True)

    barcodes = readBarcode(barcode_f)
    inputBam.close()
    ## parallel
    if threads <= 1:
        splitBarcodeBam(
            input_bam=input_bam,
            outdir=outdir,
            barcodes=barcodes,
            index = 1,
            interval=intervals[1]
        )
        print("Cell barcodes are filtered.")
    else:
        i = 0
        print(f"Filtering cell barcodes is running with {threads} cores.")
        p = multiprocess.Pool(processes=threads)
        
        for interval in intervals.values():
            i += 1
            p.apply_async(
                splitBarcodeBam,
                args=(
                    input_bam,
                    outdir,
                    barcodes,
                    i,
                    interval
                ),
                error_callback=print_error,
            )

        p.close()
        p.join()
        print("Filtering done.")
        print("Merging results...")
    
        ## merge bam
        for barcode in barcodes:
            output_bam = os.path.join(outdir, barcode + ".bam")
            tempFiles = os.path.join(outdir, barcode + "_*.bam")
            Args_m = [f'{samtools} merge -fpc -@ {threads} {output_bam} {tempFiles}']
            subprocess.check_call(Args_m, shell=True)

            Args_m = [f'rm {tempFiles}']
            subprocess.check_call(Args_m, shell=True)


def sub2cool(
    input_pairs,
    outdir,
    metadata,
    threads,
    resolution,
    cooler,
    juicer,
    hic2cool,
    species,
    normalize
):
    """
        extract contacts pairs of sub-cluster-cells and convert to .cool
        Args:
            input_pairs: raw dscHiC contact pairs file
            outdir: output directory 
            metadata: specify cellbarcodes and corresponding sub-cluster
            threads: >= 1
            resolution: list different resolution size
            cooler: path
            juicer: path
            hic2cool: path
            species: mm10 or hg38
            normalize: whether to normalize matrics in juicer pre
        Returns:
            .sub_cluster.cool

    """
    cluster_dict = pd.read_csv(metadata,
                        sep="\t")   
    groups = set(list(cluster_dict['group']))
    
    df_iter = chunk_dataframe(input_pairs, threads)
    parallel_results = []
    p = multiprocess.Pool(processes=threads)
    for chunk in df_iter:
        p.apply_async(
            sub_cluster2cool,
                args=(
                chunk,
                groups,
                cluster_dict,
            ),
            callback=parallel_results.append,
            error_callback=print_error,
        ) 
    p.close()
    p.join()
    print("[::dscHiCtools] Merging results...")

    for _group in groups:
        print(f"[::dscHiCtools] converting _group...")

        outfile_contacts_pairs = os.path.join(outdir, _group + ".contacts.pairs.txt.gz")
        wirte_all_group(
            outfile_contacts_pairs,
            _group,
            parallel_results)

        outfile_cooler = os.path.join(outdir, _group + "." + ".cool")

        ## convert
        resolution_symbol = ",".join(resolution)
        ## juicer .hic
        print("[::juicer] convert to .hic...")
        if not normalize:
            print("Normalizing matrics...")
            Args_m = [f'{juicer} pre --threads {threads} -r {resolution_symbol} {outfile_contacts_pairs} {outdir + "/" + _group}.hic {species}']
        else:
            Args_m = [f'{juicer} pre -n --threads {threads} -r {resolution_symbol} {outfile_contacts_pairs} {outdir + "/" + _group}.hic {species}']
        subprocess.check_call(Args_m, shell=True)

        ## convert to .mcool
        print("[::hic2cool] convert to .mcool...")

        Args_m = [f'{hic2cool} convert -r 0 {outdir + "/" + _group}.hic {outdir + "/" + _group}.mcool']
        subprocess.check_call(Args_m, shell=True)
    
        print("[::cooler] balancing...")
        for _resolution in resolution:
            Args_m = [f'{cooler} balance -p {threads} {outdir + "/" + _group}.mcool::/resolutions/{_resolution}']
            subprocess.check_call(Args_m, shell=True)


        # Args_m = [f'{cooler} cload pairs -c1 1 -p1 2 -c2 3 -p2 4 {chromSize}:{resolution} {outfile_contacts_pairs} {outfile_cooler}']
        # subprocess.check_call(Args_m, shell=True)