from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from itertools import islice
from dscHiCtools.chunkFiles import (
    get_n_lines, 
    chunk_reads
)
import time 
import sys 
import gzip
import h5py
import numpy as np
import pandas as pd 
from math import floor

def splitContacts(
    input_contacts,
    outdir,
    indexes,
    process_info
):
    """
        split contacts file based on cell barcode 
        Args:
            input_contacts
            outdir
            indexes: chunk for multiprocess
            process_info: number of process
        Returns:
            write to cellbarcode.temp.contacts.pairs.txt.gz
            ## pairs format v1.0
            #columns: readID chr1 position1 chr2 position2 strand1 strand2
    """
    ## count the processed pairs 
    n = 1
    t = time.time()
    file_dict = {}
    if input_contacts[-2::] == "gz":
        textfile1 = gzip.open(input_contacts, "rt")
    else:
        textfile1 = open(input_contacts, "rt")
    totallines = islice(
        textfile1, indexes[0] * 1, indexes[1] * 1
    )
    
    line = next(totallines)
    cell_barcode = line.rstrip().split("\t")[0]
    output_file = os.path.join(outdir, cell_barcode + "_" + process_info)
    file_dict[cell_barcode] = output_file
    out_contacts_temp = open(output_file, 'w')
    _, chrA, binA, chrB, binB = line.rstrip().split("\t")
    new_line = ['.', chrA, str(binA), chrB, str(binB), '.', '.']
    out_contacts_temp.write("\t".join(new_line) + "\n")
    while True:
        try:
            line = next(totallines)
            if cell_barcode == line.rstrip().split("\t")[0]:
                _, chrA, binA, chrB, binB = line.rstrip().split("\t")
                new_line = ['.', chrA, str(binA), chrB, str(binB), '.', '.']
                out_contacts_temp.write("\t".join(new_line) + "\n")
            else:
                out_contacts_temp.close()
                cell_barcode = line.rstrip().split("\t")[0]
                output_file = os.path.join(outdir, cell_barcode + "_" + process_info)
                file_dict[cell_barcode] = output_file
                out_contacts_temp = open(output_file, 'wt')
                _, chrA, binA, chrB, binB = line.rstrip().split("\t")
                new_line = ['.', chrA, str(binA), chrB, str(binB), '.', '.']
                out_contacts_temp.write("\t".join(new_line) + "\n")

            # Progress info
            if n % 10000000 == 0:
                utils.eprint(
                    "[splitContacts::] Processed 10,000,000 pairs in {}. Total "
                    "pairs: {:,} in child {}".format(
                        utils.secondsToText(time.time() - t), n, os.getpid()
                    )
                )
                sys.stdout.flush()
                t = time.time()
            n += 1
        except StopIteration as e:
            break

    utils.eprint(
        "[splitContacts::] Mapping done for process {}. Processed {:,} pairs".format(os.getpid(), n - 1)
    )
    sys.stdout.flush()
    textfile1.close()
    return file_dict


def id2barcode(id_convert_file):
    id_convert = {}
    if id_convert_file[-2:]  == "gz":
        ref = gzip.open(id_convert_file, "rt")
    else:
        ref = open(id_convert_file, "rt")

    while True:
        line = ref.readline().rstrip()
        if(len(line) == 0):
            break
        _barcode, _id = line.split("\t")
        id_convert[_id] = _barcode
    
    return id_convert


def divide_chunks(l, n): 
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n] 



def splitHDF5(
    input_h5py_file,
    outdir,
    indexes,
    convert_id_dict,
    resolution
):
    """
        split HDF5 file into single contacts pairs file for scGAD computation
        Args:
            input_h5py_file: .hdf5 file (by chrom)
            outdir
            indexes: chunk for multiprocess
            convert_id_dict: dictionary map id to cellbarcode
        Returns:
            write to cell.contacts.pairs.txt
            #columns: chrA, binA, chrB, binB, count
    """
    indexes = list(map(str, indexes))
    file_record = {}
    _chr = os.path.basename(input_h5py_file).split("_")[0]
    with h5py.File(input_h5py_file) as impute_f:
        base_df = pd.DataFrame(impute_f['coordinates'])
        base_df.columns = ["binA", "binB"]
        base_df['binA'] = base_df['binA'] * resolution
        base_df['binB'] = base_df['binB'] * resolution
        base_df['chrA'] = _chr 
        base_df['chrB'] = _chr 
        for _id in indexes:
            outFile = os.path.join(outdir, f'{convert_id_dict[_id]}.{_chr}.contacts.pairs.txt')
            file_record[convert_id_dict[_id]] = outFile

            base_df['count'] = np.array(impute_f[f"cell_{_id}"])
            base_df.loc[:, ['chrA', 'binA', 'chrB', 'binB', 'count']].to_csv(
                outFile, 
                header=True,
                sep='\t', 
                index = False
            )
    return file_record


def checkSorted(file, top=100):
    n = 1
    if file[-2::] == "gz":
        with gzip.open(file, 'rt') as file:  
            line = file.readline().rstrip()
            cb, _, _, _, _ = line.split("\t")
            while n <= top:
                line = file.readline().rstrip()
                if line[0] != cb:
                    return False 
                cb, _, _, _, _ = line.split("\t")
                n += 1
    else:
        with open(file, 'rt') as file:  
            line = file.readline().rstrip()
            cb, _, _, _, _ = line.split("\t")
            while n <= top:
                line = file.readline().rstrip()
                if line[0] != cb:
                    return False
                cb, _, _, _, _ = line.split("\t")
                n += 1
    return True 


def mergeDict(parallel_results):
    final_dict = {}
    for _dict in parallel_results:
        for cell_barcode in _dict.keys():
            if final_dict.get(cell_barcode, None) is None:
                final_dict[cell_barcode] = [_dict[cell_barcode]]
            else:
                final_dict[cell_barcode].append(_dict[cell_barcode])

    return final_dict


@utils.log_info
def main(args):
    """
        split contacts file based on cell barcode in multiprocess
        Args:
            input_contacts: contacts pairs file 
            outdir: output directory that stores single cell bam
            threads: >= 1
            n_lines: number of lines needs to be processed 
        Returns:
            cellbarcode.bam
    """
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    if args.format == "pairs":
        utils.eprint(f"[splitContacts::] Pairs format mode.")
        if args.n_lines is None or args.n_lines <= 0 :
            args.n_lines = get_n_lines(args.input_contacts, "contacts")
        
        chunk_indexes = chunk_reads(args.n_lines, args.threads)

        sorted_ = checkSorted(args.input_contacts)

        if not sorted_:
            utils.eprint(f"[splitContacts::] {args.input_contacts} is not sorted. Sorting by cell barcode...")
            if args.input_contacts[-2::] == "gz":
                Args_m = [f'zcat {args.input_contacts} | sort -k1,1 | pigz > {args.input_contacts + "temp"}']
                subprocess.check_call(Args_m, shell=True)
                Args_m = [f'mv {args.input_contacts + "temp"} {args.input_contacts}']
                subprocess.check_call(Args_m, shell=True)
            else:
                Args_m = [f'sort -k1,1 {args.input_contacts} > {args.input_contacts + "temp"}']
                subprocess.check_call(Args_m, shell=True)
                Args_m = [f'mv {args.input_contacts + "temp"} {args.input_contacts}']
                subprocess.check_call(Args_m, shell=True)
        else:
            utils.eprint(f"[splitContacts::] {args.input_contacts} is sorted.")

        utils.eprint(f"[splitContacts::] Running with {args.threads} cores.")
        p = multiprocess.Pool(processes=args.threads)

        parallel_results = []
        for i, indexes in enumerate(chunk_indexes):
            p.apply_async(
                splitContacts,
                args=(
                    args.input_contacts,
                    args.outdir,
                    indexes,
                    str(i)
                ),
                error_callback=utils.print_error,
                callback=parallel_results.append,
            )

        p.close()
        p.join()
       
    elif args.format == "h5py":
        utils.eprint(f"[splitContacts::] h5py format mode.")
        convert_id_dict = id2barcode(args.id_convert_file)
        cell_num = len(convert_id_dict.keys())
        input_h5py_files = [os.path.join(args.input_contacts, file) for file in os.listdir(args.input_contacts) if f"{args.neighbour}_impute.hdf5" in file ]
        parallel_results = []
        
        p = multiprocess.Pool(processes=args.threads)
        for input_h5py_file in input_h5py_files:
            for indexes in divide_chunks(list(range(cell_num)), floor(cell_num / args.threads)):
                p.apply_async(
                    splitHDF5,
                    args=(
                        input_h5py_file,
                        args.outdir,
                        indexes,
                        convert_id_dict,
                        args.resolution
                    ),
                    error_callback=utils.print_error,
                    callback=parallel_results.append,
                )
        p.close()
        p.join()
    utils.eprint("[splitContacts::] Merging results...")
    files_dict = mergeDict(parallel_results)    

    ## output single files
    if args.header:
        with open(os.path.join(args.outdir, "header"), 'wt') as header:
            header.write('## pairs format v1.0\n')
            header.write('#columns: readID chr1 position1 chr2 position2 strand1 strand2\n')

    if args.barcode_file is not None:
        barcodes_saved = utils.readBarcode(args.barcode_file)
    else:
        barcodes_saved = list(files_dict.keys())

    for barcode in barcodes_saved:
        output_contacts_file = os.path.join(args.outdir, barcode + ".contacts.pairs.txt.gz")
        if args.header:
            tempFiles = [os.path.join(args.outdir, "header")]
        else:
            tempFiles = []
        for _file in files_dict[barcode]:
            tempFiles.append(_file)
        
        Args_m = [f'cat {" ".join(tempFiles)} | pigz > {output_contacts_file}']
        subprocess.check_call(Args_m, shell=True)

        if args.header:
            for _file in tempFiles[1:]:
                os.remove(_file)
        else:
            for _file in tempFiles:
                os.remove(_file)
    utils.eprint("[splitContacts::] Merging done.")
