import gzip
from collections import Counter
from collections import defaultdict
from itertools import islice
import functools
import sys
import random
import string
import pysam
import time
import os 
from dscHiCtools.write import write_4DN_contacts
import pandas as pd 

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def print_error(value):
    eprint("[Error::] ", value)

def log_info(func):
    """
        Decorator that prints function arguments and runtime
    """

    @functools.wraps(func)
    def wrapper(args):
        eprint(
            "[dscHiCtools::] Function {} called with the following arguments:\n".format(func.__name__)
        )
        for arg in vars(args):
            eprint(str(arg) + "\t" + str(getattr(args, arg)))
        start_time = time.time()
        func(args)
        elapsed = secondsToText(time.time() - start_time)
        eprint(f'\nFunction completed in  {elapsed}\n')

    return wrapper


def getPrefix(prefix_s):
    if prefix_s is not None:
        if prefix_s[-1::] != "_":
            prefix = prefix_s + "_"
        else:
            prefix = prefix_s
    else:
        prefix = ""
    return prefix


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


## time record
def secondsToText(secs, lang="EN"):
	"""
	Converts datetime to human readable hours, minutes, secondes format.

	Args:
		secs (float): Secondes
		lang (string): Language
	
	Returns:
		string: Human readable datetime format.
	"""
	days = secs//86400
	hours = (secs - days*86400)//3600
	minutes = (secs - days*86400 - hours*3600)//60
	seconds = secs - days*86400 - hours*3600 - minutes*60

	if lang == "ES":
		days_text = "día{}".format("s" if days!=1 else "")
		hours_text = "hora{}".format("s" if hours!=1 else "")
		minutes_text = "minuto{}".format("s" if minutes!=1 else "")
		seconds_text = "segundo{}".format("s" if seconds!=1 else "")
	elif lang == "DE":
		days_text = "Tag{}".format("e" if days!=1 else "")
		hours_text = "Stunde{}".format("n" if hours!=1 else "")
		minutes_text = "Minute{}".format("n" if minutes!=1 else "")
		seconds_text = "Sekunde{}".format("n" if seconds!=1 else "")
	elif lang == "RU":
		days_text = pluralizeRussian(days, "день", "дня", "дней")
		hours_text = pluralizeRussian(hours, "час", "часа", "часов")
		minutes_text = pluralizeRussian(minutes, "минута", "минуты", "минут")
		seconds_text = pluralizeRussian(seconds, "секунда", "секунды", "секунд")
	else:
		#Default to English
		days_text = "day{}".format("s" if days!=1 else "")
		hours_text = "hour{}".format("s" if hours!=1 else "")
		minutes_text = "minute{}".format("s" if minutes!=1 else "")
		seconds_text = "second{}".format("s" if seconds!=1 else "")

	result = ", ".join(filter(lambda x: bool(x),[
	"{0} {1}".format(days, days_text) if days else "",
	"{0} {1}".format(hours, hours_text) if hours else "",
	"{0} {1}".format(minutes, minutes_text) if minutes else "",
	"{0:.4} {1}".format(seconds, seconds_text) if seconds else ""
	]))
	return result


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
    eprint("[dscHiCtools::] Auto detecting orientation...")
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

    eprint(f'[dscHiCtools::] Forward: {fd_result} hits of {i} reads ({round(fd_result / i, 4) * 100}% rate)')
    eprint(f'[dscHiCtools::] Reverse: {rs_result} hits of {i} reads ({round(rs_result / i, 4) * 100}% rate)')
    if fd_result > rs_result:
        eprint("[dscHiCtools::] orientation sets to be forward")
        return False
    else:
        eprint("[dscHiCtools::] orientation sets to be reverse")
        return True 


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
