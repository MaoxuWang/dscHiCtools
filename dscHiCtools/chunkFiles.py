import gzip
import sys
from math import floor
from math import ceil
import pysam
import pandas as pd 

## calculate Total lines of input file

def get_n_lines(file_path, _format = "fq"):
    """
        Determines how many lines have to be processed
        depending on options and number of available lines.
        Checks that the number of lines is a multiple of 4.
        Args:
            file_path (string): Path to a fastq.gz file

        Returns:
            n_lines (int): Number of lines in the file
    """
    print("Counting number of reads...")
    n_lines = 0
    with gzip.open(file_path, "rt", encoding="utf-8", errors="ignore") as f:
        while True:
            line = f.readline().rstrip()
            if(len(line) == 0):
                break
            n_lines += 1 
    if _format == "fq":
        if n_lines % 4 != 0:
            sys.exit(
                "{}'s number of lines is not a multiple of 4. The file "
                "might be corrupted.\n Exiting".format(file_path)
            )
    return n_lines

## assign index of input file for multi-thread running

def get_indexes(start_index, chunk_size, nth):
    """
        Creates indexes from a reference index, a chunk size an nth number

        Args:
            start_index (int): first position
            chunk_size (int): Chunk size
            nth (int): The nth number
        
        Returns:
            list: First and last position of indexes
    """
    start_index = nth * chunk_size
    stop_index = chunk_size + nth * chunk_size
    return [start_index, stop_index]

## assign reads of input file for multi-thread running

def chunk_reads(n_reads, n):
    """
        Creates a list of indexes for the islice iterator from the map_reads function.
        For .fastq file
        Args:
            n_reads (int): Number of reads to split
            n (int): How many buckets for the split.
        Returns:
            indexes (list(list)): Each entry contains the first and the last index for a read.
    """
    indexes = list()
    if n_reads % n == 0:
        chunk_size = int(n_reads / n)
        rest = 0
    else:
        chunk_size = floor(n_reads / n)
        rest = n_reads - (n * chunk_size)
    for i in range(0, n):
        indexes.append(get_indexes(i, chunk_size, i))
    indexes[-1][1] += rest
    return indexes


def add_start_coords(intervals, chrom_lengths, bamfile):
    """
        given the intervals that will be handled by each core,
        break into genomic regions (removing chromosome-spanning intervals)
    """
    intervals = [[1, 0]] + intervals
    ranges = [intervals[x - 1] + intervals[x] for x in range(1, len(intervals))]
    d = {}
    x = 0
    for i in ranges:
        x += 1
        if i[0] == i[2]:
            d[x] = [(bamfile.get_reference_name(i[0] - 1), i[1], i[3])]
        else:
            d[x] = [
                (bamfile.get_reference_name(i[0] - 1), i[1], chrom_lengths[i[0] - 1])
            ]
            nchrom = i[2] - i[0]
            for y in range(nchrom - 1):
                d[x].append(
                    (bamfile.get_reference_name(i[0] + y), 0, chrom_lengths[i[0] + y])
                )
            d[x].append((bamfile.get_reference_name(i[2] - 1), 0, i[3]))
    return d


def find_chromosome_break(position, chromosomes, current_chrom):
    """
        find position that form the interval 
    """
    assert position <= sum(chromosomes), "position past end of genome"
    if position <= chromosomes[current_chrom]:
        return [current_chrom + 1, position]
    else:
        position = position - chromosomes[current_chrom]
        return find_chromosome_break(position, chromosomes, current_chrom + 1)

def chunk_bam(bamfile, threads):
    """
        chunk file into n chunks for multi-processing
    """
    chrom_lengths = bamfile.lengths
    chunksize = sum(chrom_lengths) / int(threads)
    intervals = []
    for x in range(1, threads + 1):
        position = chunksize * x
        intervals.append(find_chromosome_break(position, chrom_lengths, 0))
    return add_start_coords(intervals, chrom_lengths, bamfile)


def chunk_dataframe(
    mtx,
    threads
):
    n_line = get_n_lines(mtx, "pairs")
    chunkSize = ceil(n_line / threads)

    if(mtx[-2::] == "gz"):
        df_iter = pd.read_csv(mtx,
                        sep="\t",
                        compression = "gzip",
                        iterator=True,
                        on_bad_lines='warn',
                        chunksize=chunkSize) 
    else:
        df_iter = pd.read_csv(mtx,
                        sep="\t", 
                        header=None, 
                        iterator=True,
                        on_bad_lines='warn',
                        chunksize=chunkSize)
    return df_iter