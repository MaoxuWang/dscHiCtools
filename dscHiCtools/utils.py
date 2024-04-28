import gzip
import pybktree
from collections import Counter
from collections import defaultdict
import Levenshtein
from itertools import islice
import functools
import sys
from dscHiCtools.secondsToText import secondsToText
import random
import string
import pysam
import time
import os 
from dscHiCtools.write import write_4DN_contacts
import pandas as pd 

def print_error(value):
    print("[Error::] ", value)


def log_info(func):
    """
        Decorator that prints function arguments and runtime
    """

    @functools.wraps(func)
    def wrapper(args):
        print(
            "Function {} called with the following arguments:\n".format(func.__name__),
            file=sys.stderr
        )
        for arg in vars(args):
            print(str(arg) + "\t" + str(getattr(args, arg)), file=sys.stderr)
        start_time = time.time()
        func(args)
        elapsed = secondsToText(time.time() - start_time)
        print(f'\nFunction completed in  {elapsed}\n',
             file=sys.stderr)

    return wrapper


## get reverse reads 5'-3'
def reverseRead(read):
    """
        get reverse reads ; direction: 5'-3'
        Returns: 
            new read
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    new_read = ''
    for base in read:
        new_read = complement[base] + new_read
    return new_read


def readBarcode(ref_barcode):
    """
        read cellbarcodes file
        one cellbarcode per row
        auto detect gzipped or not
        Args:
            ref_barcode
        Returns:
            (set) barcodes
    """
    barcode_f = []
    if ref_barcode[-2:]  == "gz":
        with gzip.open(ref_barcode, "rt") as ref:
            while True:
                barcode = ref.readline().rstrip()
                if(len(barcode) == 0):
                    break
                barcode_f.append(barcode)
    else:
        with open(ref_barcode, "rt") as ref:
            while True:
                barcode = ref.readline().rstrip()
                if(len(barcode) == 0):
                    break
                barcode_f.append(barcode)
    return set(barcode_f)


def multiomeBarcodesDict(multiome_RNA_barcode, multiome_DNA_barcode):
    """
        create dictionary that maps RNA barcodes to DNA barcodes or DNA barcodes to RNA barcodes
    """
    RNA_to_DNA = defaultdict()
    DNA_to_RNA = defaultdict()

    with gzip.open(multiome_RNA_barcode, "rt") as RNA_barcodes, gzip.open(multiome_DNA_barcode, "rt") as DNA_barcodes:
        while True:
            RNA = RNA_barcodes.readline().strip()
            if len(RNA) == 0:
                break 
            DNA = DNA_barcodes.readline().strip()
            RNA_to_DNA[RNA] = DNA
            DNA_to_RNA[DNA] = RNA 
    return RNA_to_DNA, DNA_to_RNA


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

    print("Generated barcode tree from reference barcode file")
    return(barcode_tree)


## convert base quality
def phred33toQ(qual):
    """
        convert to quality
    """
    return ord(qual) - 33


def findOffset(seq_barcode, ref_barcode, qual):
    """
        find the offset of mismatch 
        Returns: 
            the sequencing quality of this base
    """
    offset = "Na"
    for i in range(len(ref_barcode)):
        if seq_barcode[i] == ref_barcode[i]:
            continue
        offset = i
    if offset == "Na":
        return 0
    else:
        return (phred33toQ(qual[i]))


# merge
def merge_results(parallel_results):
    """
        Merge chunked results from parallel processing.
        Args: 
            parallel_results (list): List of dict with mapping results.
        Returns: 
            merged_results (Counter): Total reads per cell as a Counter
    """
    mapped_results = Counter()
    unmapped_results = Counter()

    for chunk in parallel_results:
        mapped = chunk[0]
        unmapped = chunk[1]

        mapped_results.update(mapped)
        unmapped_results.update(unmapped) 
    
    return (mapped_results, unmapped_results)


def writeFile(
    outFile,
    cellBarcode,
    readNames,
    seqs,
    chars,
    qualities,
    i
):
    readName = "@" + cellBarcode + ":" + readNames[i][1:]
 
    outFile.write(readName.encode())
    outFile.write(seqs[i].encode())
    outFile.write(chars[i].encode())
    outFile.write(qualities[i].encode())


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
                    print(
                        "Processed 10,000,000 reads in {}. Total "
                        "reads: {:,} in child {}".format(
                            secondsToText(time.time() - t), n, os.getpid()
                        )
                    )
                    sys.stdout.flush()
                    t = time.time()
                n += 1
                read = seqs[1].strip()
                if reverse:
                    cell_barcode_string = reverseRead(read[start_idx:end_idx])
                else:
                    cell_barcode_string = read[start_idx:end_idx]
                if len(cell_barcode_string) != 16:
                    print("wrong cell barcode length found!")
                    continue

                ### unmapped
                cell_barcode = alignReference(cell_barcode_string, barcode_tree, min_mismatch)
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
                    writeFile(out_r1_file, 
                        cell_barcode,
                        readNames,
                        seqs,
                        _char,
                        qualities,
                        0)
                    writeFile(out_r2_file, 
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
    print(
        "Mapping done for process {}. Processed {:,} reads".format(os.getpid(), n - 1)
    )
    sys.stdout.flush()
    return mapped, unmapped



def addCBtag(intervals, input_bam, output_bam):
    """
        add CB tag to bam file
        Args:
            intervals: slice of bam file that is used for parallel computing
            input_bam: input bam file
            output_bam: CB tag added bam file
        Returns:
            CB tag added bam file
    """
    samfile = pysam.AlignmentFile(input_bam, "rb")
    header = samfile.header

    outputBam = pysam.AlignmentFile(output_bam, "wb", header=header)
    for i in intervals:
        for read in samfile.fetch(i[0], i[1], i[2]):
            cell_barcode = read.qname.split(":")[0]
            tags = read.tags
            tags.append(('CB', cell_barcode))
            read.set_tags(tags)
            outputBam.write(read)
    outputBam.close()
    samfile.close()


def alignReference(
    cell_barcode_string,
    barcode_tree,
    min_mismatch
    ):
    """
        if multiple hits occured under min mismatch, discard it
        Args:
            cell_barcode_string: cell barcoed from R3
            barcode_tree: BKtree generated by reference file
            min_mismatch: min mismatch allowed 
        Returns:
            valid cell barcode
    """
    hits = list(barcode_tree.find(cell_barcode_string, min_mismatch))
    if(len(hits) == 1):
        valid_barcode = list(hits[0])[1]
        return(valid_barcode)
    elif(len(hits) > 1):
        ## discard multiple hits
        return 0
    else:
        ## discard barcode with no hits
        return 0


def autoDetectChemistry(
    index,
    barcode_tree,
    min_mismatch,
    n = 10000):
    """
        auto detect whether to reverse the i2 index read
        Args:
            input_bam: input bam file path
            barcode_tree: BKtree for reference barcode file
            min_mismatch: min mismatch bases allowed for mapping cell barcode to reference
        Returns: 
            reverse: whether to reverse barcode (associated with chemistry)
    """
    print("Auto detecting orientation...")
    i = 0
    fd_result =0
    rs_result = 0

    with gzip.open(index, 'rt') as textfile1:
        totallines = islice(
            textfile1, 0 * 4, n * 4
        )
        while True:
            try:    
                readNames = next(totallines)
                seq = next(totallines).strip()
                _char = next(totallines)
                qualities = next(totallines)

            except StopIteration as e:
                break

            cell_barcode_string = seq
            cell_barcode_string_reverse = reverseRead(seq)
            if alignReference(cell_barcode_string, barcode_tree, min_mismatch) != 0:
                fd_result += 1
            if alignReference(cell_barcode_string_reverse, barcode_tree, min_mismatch) != 0:
                rs_result += 1
            i += 1

    print(f'Forward: {fd_result} hits of {i} reads ({round(fd_result / i, 4) * 100}% rate)')
    print(f'Reverse: {rs_result} hits of {i} reads ({round(rs_result / i, 4) * 100}% rate)')
    if fd_result > rs_result:
        print("orientation sets to be forward")
        return False
    else:
        print("orientation sets to be reverse")
        return True 


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
                print(
                    "Processed 1,000,000 reads in {}. Total "
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
                outputBam.write(read)
            else:
                continue
    outputBam.close()
    inputBam.close()


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
                print(
                    "Processed 1,000,000 reads in {}. Total "
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


def cluster2dict(metadata):
    """
        create dict that map barcode to cluster
    """
    n = 0
    cluster_dict = {}
    if metadata[-2:] == 'gz':
        with gzip.open(metadata, 'rt') as file:
            while True:
                line = file.readline().rstrip()
                if len(line) == 0:
                    break
                n += 1
                barcode, cluster = line.split("\t")
                cluster_dict[barcode] = cluster
    else:
        with open(metadata, 'rt') as file:
            while True:
                line = file.readline().rstrip()
                if len(line) == 0:
                    break
                n += 1
                barcode, cluster = line.split("\t")
                cluster_dict[barcode] = cluster
    print(f"Read {n} cell barcodes.", file=sys.stderr)
    return cluster_dict


def wirte_all_group(
    outfile_contacts_pairs, 
    _group, 
    parallel_results):
    """
        iterate parallel results to merge contacts paris
    """
    df_group = pd.DataFrame()
    for chunk in parallel_results:
        df_group = pd.concat([df_group, chunk[_group]])
    ## write
    df_group['readID'] = "."
    df_group['strand1'] = "."
    df_group['strand2'] = "."

    df_group = df_group.sort_values(by=['chr1', 'chr2', 'pos1', 'pos2'])
    write_4DN_contacts(df_group.loc[:,  ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]],
        outfile_contacts_pairs
    )


def sub_cluster2cool(
    df,
    groups,
    cluster_dict):
    """
        extract contacts of sub-cluster-cellbarcode (for multi-process running)
        Args:
            df: contact pairs file
            cluster_dict: df that maps cell barcode to cluster
        Returns:
           contacts pairs dataframe for sub-cluster
    """
    sub_file_dict = {}

    df.columns = ["cellbarcode", "chr1", "pos1", "chr2", "pos2"]
    df_clustered = df.set_index("cellbarcode").join(cluster_dict.set_index("cellbarcode"), on="cellbarcode", how="inner")

    for _group in groups:
        sub_file_dict[_group] = df_clustered[df_clustered['group'] == _group]

    return sub_file_dict