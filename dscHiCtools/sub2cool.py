from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from itertools import islice
from dscHiCtools.chunkFiles import (
    get_n_lines, 
    chunk_reads
)
import gzip
import pandas as pd 
from dscHiCtools.write import write_4DN_contacts


def sub_cluster2cool(
    input_contact_pairs,
    groups,
    index,
    chunk_id,
    outdir,
    cluster_dict):
    """
        extract contacts of sub-cluster-cellbarcode (for multi-process running)
        Args:
            input_contact_pairs: contact pairs file
            index: for multi-processing
            chunk_id: process id
            outdir: output path for result
            cluster_dict: df that maps cell barcode to cluster
        Returns:
           contacts pairs dataframe for sub-cluster
    """
    file_record = {}
    if input_contact_pairs[-2::] == "gz":
        textfile1 = gzip.open(input_contact_pairs, "rt")
    else:
        textfile1 = open(input_contact_pairs, "rt")
    
    totallines = islice(
        textfile1, index[0] * 1, index[1] * 1
    )
    out_file_dict = {}
    for _group in groups:
        output_file = os.path.join(outdir, _group + "_" + str(chunk_id) + ".contact_pairs")
        file_record[_group] = output_file
        out_file_dict[_group] = open(output_file, "w")
   
    while True:
        try:
            line = next(totallines)
            cellbarcode, chrA, binA, chrB, binB = line.rstrip().split("\t")
            if cellbarcode in cluster_dict:
                row = [".", chrA, binA, chrB, binB, ".", "."]
                out_file_dict[cluster_dict[cellbarcode]].write(
                    "\t".join(list(map(str, row))) + "\n"
                )
            else:
                continue
        except StopIteration as e:
            break
    
    for _group in out_file_dict:
        out_file_dict[_group].close()
    textfile1.close()
    return file_record


def barcode2dict(metadata,
    indir,
    suffix):
    #suffix: ".contacts.pairs.txt.gz" or ".contact.rmbkl.tsv.gz"
    n = 0
    cluster_dict = {}
    cluster_ = set()
    if metadata[-2:] == 'gz':
        file = gzip.open(metadata, 'rt')
    else:
        file = open(metadata, 'rt')
    
    ## remove header
    line = file.readline().rstrip()

    while True:
        line = file.readline().rstrip()
        if len(line) == 0:
            break
        n += 1
        barcode, cluster = line.split("\t")
        file_candidate = os.path.join(indir, barcode + suffix)

        if not os.path.exists(file_candidate):
            utils.eprint(f"[sub2cool::] {file_candidate} not found.")
            continue
        cluster_.add(cluster)
        if cluster not in cluster_dict:
            cluster_dict[cluster] = [ file_candidate ]
        else:
            cluster_dict[cluster].append(file_candidate)

    file.close()

    utils.eprint(f"[sub2cool::] Read {n} cell barcodes and {len(cluster_)} groups.")
    return cluster_dict, cluster_


def cluster2dict(metadata):
    """
        create dict that map barcode to cluster
    """
    n = 0
    cluster_dict = {}
    cluster_ = set()
    if metadata[-2:] == 'gz':
        file = gzip.open(metadata, 'rt')
    else:
        file = open(metadata, 'rt')
    
    ## remove header
    line = file.readline().rstrip()

    while True:
        line = file.readline().rstrip()
        if len(line) == 0:
            break
        n += 1
        barcode, cluster = line.split("\t")
        cluster_.add(cluster)
        cluster_dict[barcode] = cluster
    file.close()

    utils.eprint(f"[sub2cool::] Read {n} cell barcodes and {len(cluster_)} groups.")
    return cluster_dict, cluster_


def mergeResults(parallel_results):
    final_dict = {}
    for item in parallel_results:
        for _group in item.keys():
            if _group not in final_dict:
                final_dict[_group] = [item[_group]]
            else:
                final_dict[_group].append(item[_group])
    
    return final_dict


def writeFile(
    files_write,
    output_contacts_file,
    header,
    _format,
    remove
):     
    # _format: merged ; sc
    tempFiles = []
    for _file in files_write:
        tempFiles.append(_file)
        if _format == "merged":
            Args_m = [f'cat {" ".join(tempFiles)} | sort -k2,2d -k4,4d -k3,3n -k5,5n | pigz > {output_contacts_file + ".temp"}']
        else:
            # Args_m = [f'zcat {" ".join(tempFiles)} | grep -v "#" | sort -k2,2d -k4,4d -k3,3n -k5,5n | pigz > {output_contacts_file + ".temp"}']
            Args_m = [f'zcat {" ".join(tempFiles)} | sort -k2,2d -k4,4d -k3,3n -k5,5n | pigz > {output_contacts_file + ".temp"}']

        utils.eprint(f"[sub2cool::] executing {Args_m}...")
        subprocess.check_call(Args_m, shell=True)
    if remove:
        for _file in tempFiles[1:]:
            os.remove(_file)
    
    Args_m = [f'cat {header} {output_contacts_file + ".temp"} > {output_contacts_file}']
    subprocess.check_call(Args_m, shell=True)

    Args_m = [f'rm {output_contacts_file + ".temp"}']
    subprocess.check_call(Args_m, shell=True)


def convertFormat(
    _group,
    final_dict,
    outdir,
    header_f,
    resolution,
    normalize,
    juicer,
    threads,
    species,
    hic2cool,
    cooler,
    _format="merged"
):
    """
        format: sc or merged
        sc: a path stores for single cell file
        merged: only one merged contacts pairs file   
    """
    utils.eprint(f"[sub2cool::] converting {_group}...")
    outfile_contacts_pairs = os.path.join(outdir, _group + ".contacts.pairs.txt.gz")
    
    writeFile(
        final_dict[_group],
        outfile_contacts_pairs,
        header_f,
        _format,
        _format == "merged")
    outfile_cooler = os.path.join(outdir, _group + "." + ".cool")
    ## convert
    resolution_symbol = ",".join(resolution)
    ## juicer .hic
    utils.eprint("[sub2cool::] convert to .hic...")
    if normalize:
        utils.eprint("[sub2cool::] Normalizing matrics...")
        Args_m = [f'{juicer} pre --threads {threads} -r {resolution_symbol} {outfile_contacts_pairs} {outdir + "/" + _group}.hic {species}']
    else:
        Args_m = [f'{juicer} pre -n --threads {threads} -r {resolution_symbol} {outfile_contacts_pairs} {outdir + "/" + _group}.hic {species}']
    subprocess.check_call(Args_m, shell=True)

    ## convert to .mcool
    utils.eprint("[sub2cool::] convert to .mcool...")

    Args_m = [f'{hic2cool} convert -r 0 {outdir + "/" + _group}.hic {outdir + "/" + _group}.mcool']
    subprocess.check_call(Args_m, shell=True)

    utils.eprint("[sub2cool::] balancing...")
    for _resolution in resolution:
        Args_m = [f'{cooler} balance -p {threads} {outdir + "/" + _group}.mcool::/resolutions/{_resolution}']
        subprocess.check_call(Args_m, shell=True)
    utils.eprint(f"[sub2cool::] {_group} finished ...")



@utils.log_info
def main(args):
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
    header_f = os.path.join(args.outdir, "header.gz")
    with gzip.open(header_f, 'wb') as header:
        header.write('## pairs format v1.0\n'.encode())
        header.write('#columns: readID chr1 position1 chr2 position2 strand1 strand2\n'.encode())

    ## if input_pairs is directory that stores single cell contacts pairs
    if args.input_pairs[-2:] != "gz":
        utils.eprint("[sub2cool::] merge single cell contact pairs file directly ...")

        ## merge path directly
        barcode_file_dict, groups = barcode2dict(
            args.metadata,
            args.input_pairs,
            args.suffix)

        for _group in groups:
            p = multiprocess.Pool(processes=args.threads)
            p.apply_async(
                convertFormat,
                    args=(
                        _group,
                        barcode_file_dict,
                        args.outdir,
                        header_f,
                        args.resolution,
                        args.normalize,
                        args.juicer,
                        args.threads // len(groups),
                        args.species,
                        args.hic2cool,
                        args.cooler,
                        "sc",
                ),
                error_callback=utils.print_error,
            )
        p.close()
        p.join()

        utils.eprint("[sub2cool::] Finished.")

        return 

    cluster_dict, groups = cluster2dict(args.metadata)
    
    args.n_lines = get_n_lines(args.input_pairs, "contacts")
        
    chunk_indexes = chunk_reads(args.n_lines, args.threads)

    parallel_results = []
    p = multiprocess.Pool(processes=args.threads)
    for chunk_id, indexes in enumerate(chunk_indexes):
        p.apply_async(
            sub_cluster2cool,
                args=(
                    args.input_pairs,
                    groups,
                    indexes,
                    chunk_id,
                    args.outdir,
                    cluster_dict,
            ),
            callback=parallel_results.append,
            error_callback=utils.print_error,
        ) 
    p.close()
    p.join()
    utils.eprint("[sub2cool::] Merging results...")
    final_dict = mergeResults(parallel_results)
    
    del parallel_results
    
    p = multiprocess.Pool(processes=len(final_dict.keys()))
    for _group in final_dict.keys():
        p.apply_async(
            convertFormat,
                args=(
                    _group,
                    final_dict,
                    args.outdir,
                    header_f,
                    args.resolution,
                    args.normalize,
                    args.juicer,
                    args.threads // len(final_dict.keys()),
                    args.species,
                    args.hic2cool,
                    args.cooler
            ),
            error_callback=utils.print_error,
        )
    
    p.close()
    p.join()