import pandas as pd
from math import floor
from collections import defaultdict
import numpy as np
import time
import sys
import multiprocess
import gzip
from math import ceil
import cooler
from dscHiCtools.chunkFiles import get_n_lines
from dscHiCtools import utils


"""
    calculate scA/B matrix for dscHiC data
    Param:
        Input contact pairs file
        reference CpG frequency file of certain resolution
        weighed: bool, mean of other bins by weighted contacts
        resolution: resolution for binning 
        input_format: whether from Tan described method, scVI-3d or higashi imputed result
    Returns:
        output format:
        cellbarcode       bin_index A/B
        AAACGAAGTCCAGACC	chr4_0	0.005192
        AAACTGCTCACATTCT	chr4_0	0.021393
        AACAAAGCACAAACAA	chr4_0	0.005112
        AACCAACCAGGTCTGC	chr4_0	0.0287
        AACGGGATCAAACCAC	chr4_0	0.008948

        bin_index: chr_start of bin
"""



def cpg2ref(ref):
    cpg_dict = {}
    with open(ref, "r") as f:
        for line in f:
            chrom, bin_index, cpg_freq = line.rstrip().split('\t')
            cpg_dict[(chrom, bin_index)] = float(cpg_freq)
    return cpg_dict


def contacts2bin(df, res, weighted):
    """
        make bins for contact pairs
    """
    if(weighted):
        df["bin1"] = [str(floor(i)*res) for i in df["bin1"] / res]
        df["bin2"] = [str(floor(i)*res) for i in df["bin2"] / res]       
    else:
        df["bin1"] = [str(floor(i)*res) for i in df["bin1"] / res]
        df["bin2"] = [str(floor(i)*res) for i in df["bin2"] / res]   
        df.drop_duplicates(subset=["cell_barcode", "chr1","bin1", "chr2", "bin2"],
                           keep='first',
                           inplace=True)
    return(df)


def smoothContacts(df, 
    weighted, 
    output,
    cpg_dict,
    cis_only):
    """
        used for raw AB contact pairs as input
        e.g. chr1 bin1 chr2 bin2
        weighted mean is calculated by adding value one by one 
    """
    if(weighted):
        color_data = defaultdict(lambda: [])
        for _, row in df.iterrows():
            if cis_only:
                if( row["chr1"] != row["chr2"] ):
                    continue
            if( (row["chr1"], row["bin1"]) == (row["chr2"], row["bin2"]) ):
                continue
            else:
                if (cpg_dict.get((row["chr2"], row["bin2"])) is not None):

                    color_data[(row["cell_barcode"], row["chr1"], row["bin1"])].append(cpg_dict[(row["chr2"], row["bin2"])])

                if (cpg_dict.get((row["chr1"], row["bin1"])) is not None):

                    color_data[(row["cell_barcode"], row["chr2"], row["bin2"])].append(cpg_dict[(row["chr1"], row["bin1"])])
        with gzip.open(output, "wb") as f:
            for cell_barcode,chrom, bin_index in sorted(color_data.keys()):
                cpg_freq_averaged = np.mean(color_data[(cell_barcode,chrom, bin_index)])
                f.write('{}\t{}\t{}\n'.format(cell_barcode, chrom + "_" +str(bin_index), str(cpg_freq_averaged)).encode())

    else:
        color_data = defaultdict(lambda: [])
        for _, row in df.iterrows():
            if cis_only:
                if( row["chr1"] != row["chr2"] ):
                    continue
            if( (row["chr1"], row["bin1"]) == (row["chr2"], row["bin2"]) ):
                continue
            else:
                color_data[(row["cell_barcode"], row["chr1"], row["bin1"])].append(cpg_dict[(row["chr2"], row["bin2"])])
                color_data[(row["cell_barcode"], row["chr2"], row["bin2"])].append(cpg_dict[(row["chr1"], row["bin1"])])

        with gzip.open(output, "wb") as f:
            for cell_barcode,chrom, bin_index in sorted(color_data.keys()):
                cpg_freq_averaged = np.mean(color_data[(cell_barcode,chrom, bin_index)])
                f.write('{}\t{}\t{}\n'.format(cell_barcode, chrom + "_" +str(bin_index), str(cpg_freq_averaged)).encode())
    return(True)


def smoothContactsMC(df, 
    weighted, 
    cpg_dict,
    res,
    barcode,
    cis_only):
    """
        Multi-core running for raw AB contact pairs as input
        e.g. chr1 bin1 chr2 bin2
        weighted mean is calculated by adding value one by one 
    """
    df.columns = ["cell_barcode", "chr1","bin1", "chr2", "bin2"]
    df = contacts2bin(df, res, weighted)

    if(weighted):
        color_data = defaultdict(lambda: [])
        for _, row in df.iterrows():
            ## filter barcode
            if( barcode is not None and row["cell_barcode"] not in barcode ):
                continue
            if cis_only:
                if( row["chr1"] != row["chr2"] ):
                    continue
            if( (row["chr1"], row["bin1"]) == (row["chr2"], row["bin2"]) ):
                continue
            else:
                if (cpg_dict.get((row["chr2"], row["bin2"])) is not None):
                    color_data[(row["cell_barcode"], row["chr1"], row["bin1"])].append(cpg_dict[(row["chr2"], row["bin2"])])
                if (cpg_dict.get((row["chr1"], row["bin1"])) is not None):

                    color_data[(row["cell_barcode"], row["chr2"], row["bin2"])].append(cpg_dict[(row["chr1"], row["bin1"])])
    else:
        color_data = defaultdict(lambda: [])
        for _, row in df.iterrows():
            if( barcode is not None and row["cell_barcode"] not in barcode ):
                continue
            if cis_only:
                if( row["chr1"] != row["chr2"] ):
                    continue
            if( (row["chr1"], row["bin1"]) == (row["chr2"], row["bin2"]) ):
                continue
            else:
                color_data[(row["cell_barcode"], row["chr1"], row["bin1"])].append(cpg_dict[(row["chr2"], row["bin2"])])
                color_data[(row["cell_barcode"], row["chr2"], row["bin2"])].append(cpg_dict[(row["chr1"], row["bin1"])])
    return(color_data)


def smoothContacts_scVI(
    df, 
    output,
    cpg_dict,
    cis_only):
    """
        df: scVI-3d input format : cellbarcode chrA binA chrB binB count
        output: AB_value.txt 
        cpg_dict: dictionary that stores CpG frequency for each bin
        weighted value is calculated by color * weight (normalized counts vector)
    """
    color_data = defaultdict(lambda: [])
    weights_data = defaultdict(lambda: [])

    for _, row in df.iterrows():
        if cis_only:
            if( row["chr1"] != row["chr2"] ):
                continue
        if( (row["chr1"], row["bin1"]) == (row["chr2"], row["bin2"]) ):
            continue
        else:
            if (cpg_dict.get((row["chr2"], str(row["bin2"]))) is not None):
                color_data[(row["cell_barcode"], row["chr1"], str(row["bin1"]))].append(cpg_dict[(row["chr2"], str(row["bin2"]))])
                weights_data[(row["cell_barcode"], row["chr1"], str(row["bin1"]))].append(row["count"])
            if (cpg_dict.get((row["chr1"], str(row["bin1"]))) is not None):
                color_data[(row["cell_barcode"], row["chr2"], str(row["bin2"]))].append(cpg_dict[(row["chr1"], str(row["bin1"]))])
                weights_data[(row["cell_barcode"], row["chr2"], str(row["bin2"]))].append(row["count"])

    with gzip.open(output, "wb") as f:
        for cell_barcode,chrom, bin_index in sorted(color_data.keys()):
            cpg_freq_averaged = np.average(a=color_data[(cell_barcode,chrom, bin_index)],
                                          weights=weights_data[(cell_barcode,chrom, bin_index)]
                                )
            f.write('{}\t{}\t{}\n'.format(cell_barcode, chrom + "_" +str(bin_index), str(cpg_freq_averaged)).encode())

    return(True)


def smoothContacts_scVIMC(
    df, 
    cpg_dict,
    barcode):
    """
        multi-core
        df: scVI-3d input format : cellbarcode chrA binA chrB binB count
        cpg_dict: dictionary that stores CpG frequency for each bin
        weighted value is calculated by color * weight (normalized counts vector)
    """
    color_data = defaultdict(lambda: [[], []])

    for _, row in df.iterrows():
        ## filter barcode
        if cis_only:
            if( row["chr1"] != row["chr2"] ):
                continue
        if( (row["chr1"], row["bin1"]) == (row["chr2"], row["bin2"]) ):
            continue
        else:
            if (cpg_dict.get((row["chr2"], str(row["bin2"]))) is not None):
                color_data[(row["cell_barcode"], row["chr1"], str(row["bin1"]))][0].append(cpg_dict[(row["chr2"], str(row["bin2"]))])
                color_data[(row["cell_barcode"], row["chr1"], str(row["bin1"]))][1].append(row["count"])
            if (cpg_dict.get((row["chr1"], str(row["bin1"]))) is not None):
                color_data[(row["cell_barcode"], row["chr2"], str(row["bin2"]))][0].append(cpg_dict[(row["chr1"], str(row["bin1"]))])
                color_data[(row["cell_barcode"], row["chr2"], str(row["bin2"]))][1].append(row["count"])
    return(color_data)


def mergeResults(
    parallel_results,
    input_format,
    output):

    if (input_format != "raw"):
        final_results = defaultdict(lambda: [[], []])
        for chunks in parallel_results:
            for key in chunks.keys():
                final_results[key][0].extend(chunks[key][0])
                final_results[key][1].extend(chunks[key][1])
        ## outPut file 
       
        with gzip.open(output, "wb") as f:
            for cell_barcode,chrom, bin_index in sorted(final_results.keys()):
                if (np.sum(final_results[(cell_barcode,chrom, bin_index)][1]) > 0):
                    cpg_freq_averaged = np.average(a=final_results[(cell_barcode,chrom, bin_index)][0],
                                                weights=final_results[(cell_barcode,chrom, bin_index)][1]
                                        )
                    f.write('{}\t{}\t{}\n'.format(cell_barcode, chrom + "_" +str(bin_index), str(cpg_freq_averaged)).encode())
                else:
                    utils.eprint(f'[scAB::] {cell_barcode}:{chrom}-{bin_index} weights NA')
    else:
        final_results = defaultdict(lambda: [])
        for chunks in parallel_results:
            for key in chunks.keys():
                final_results[key].extend(chunks[key])
    
        with gzip.open(output, "wb") as f:
                for cell_barcode,chrom, bin_index in sorted(final_results.keys()):
                    cpg_freq_averaged = np.mean(final_results[(cell_barcode,chrom, bin_index)])
                    f.write('{}\t{}\t{}\n'.format(cell_barcode, chrom + "_" +str(bin_index), str(cpg_freq_averaged)).encode())


@utils.log_info
def main(args):
    start_time = time.time()

    cpg_dict = cpg2ref(args.ref_CpG)

    ## single core running
    if (args.cores == 1):
        utils.eprint("[scAB::] scAB value calculating is in single core running mode.")
        if ( args.input_format == "scVI-3d" ):
            utils.eprint("[scAB::] scVI-3d input accessed...")
            if(args.input_contact_pairs[-2::] == "gz"):
                df = pd.read_csv(args.input_contact_pairs, sep="\t", header=None, compression = "gzip") 
            else:
                df = pd.read_csv(args.input_contact_pairs, sep="\t", header=None) 
            df.columns = ["cell_barcode", "chr1","bin1", "chr2", "bin2", "count"]
            utils.eprint(df.head())
            _ = smoothContacts_scVI(df, args.output, cpg_dict, args.cis_only)
            spent_time = utils.secondsToText(time.time() - start_time)
            utils.eprint(f'[scAB::] scA/B matrix created using {spent_time}.')            
        else:
            if(args.input_contact_pairs[-2::] == "gz"):
                df = pd.read_csv(args.input_contact_pairs, sep="\t", header=None, compression = "gzip") 
            else:
                df = pd.read_csv(args.input_contact_pairs, sep="\t", header=None) 
            df.columns = ["cell_barcode", "chr1","bin1", "chr2", "bin2"]
            if ( args.barcode_file ):
                utils.eprint(f"[scAB::] filtering cell barcode by {args.barcode_file}.")
                barcode = pd.read_csv(args.barcode_file, sep="\t", header=None)  
                barcode.columns = ["cell_barcode"]
                df = barcode.merge(df, on="cell_barcode")
            else:
                utils.eprint("[scAB::] keeping all cell barcodes.")
            df = contacts2bin(df, args.res, args.weighted)
            _ = smoothContacts(df, args.weighted, args.output, cpg_dict, args.cis_only)
            spent_time = utils.secondsToText(time.time() - start_time)
            utils.eprint(f'[scAB::] scA/B matrix created using {spent_time}.')
    ## multi-core running
    else:
        if ( args.input_format == "higashi_imputed" ):
            if (args.cellbarcode_dict is None):
                raise Exception("[scAB::] Please specify converting dictionary for higashi...") 
            
            higashi_barcode_dict = {}
            with open(args.cellbarcode_dict, 'r') as f:
                for line in f.readlines():
                    line_elements = line.strip().split()
                    higashi_barcode_dict[line_elements[1]] = line_elements[0]

            totals_idx = np.array_split(range(len(higashi_barcode_dict)), args.cores, axis=0)
            p = multiprocess.Pool(processes=args.cores)
            parallel_results = []
            for idx_array in totals_idx:
                chunk = pd.DataFrame()
                for i in idx_array:
                    c1 = cooler.Cooler(f'{args.input_contact_pairs}::cells/cell_{i}')
                    temp_1 = c1.matrix(balance=False, as_pixels=True, join=True)[:,:]
                    temp_1 = temp_1.drop(["end1","end2"], axis=1)
                    temp_1.columns = ["chr1", "bin1", "chr2", "bin2", "count"]
                    temp_1['cell_barcode'] = higashi_barcode_dict[str(i)]
                    chunk = pd.concat([chunk, temp_1])
                
                p.apply_async(
                    smoothContacts_scVIMC,
                    args=(
                        chunk,
                        cpg_dict,
                        list(higashi_barcode_dict.values()),
                        args.cis_only
                    ),
                    callback=parallel_results.append,
                    error_callback=utils.print_error,
                    ) 
            p.close()
            p.join()
            utils.eprint("[scAB::] Multi-process finished")
            utils.eprint("[scAB::] Merging results...")

        else:
            if ( args.barcode_file ):
                utils.eprint(f"[scAB::] filtering cell barcode by {args.barcode_file}.")
                barcode = pd.read_csv(args.barcode_file, sep="\t", header=None) 
                barcode.columns = ["cell_barcode"]
                barcode = list(barcode['cell_barcode'])
            else:
                barcode = None
                utils.eprint("[scAB::] Keeping all barcodes")
            n_line = get_n_lines(args.input_contact_pairs, "pairs")
            chunkSize = ceil(n_line / args.cores)
            utils.eprint(f"[scAB::] scAB value calculating is running with {args.cores} cores.")
            p = multiprocess.Pool(processes=args.cores)


            ## read file first
            if(args.input_contact_pairs[-2::] == "gz"):
                df_iter = pd.read_csv(args.input_contact_pairs,
                                sep="\t", header=None,
                                compression = "gzip",
                                iterator=True,
                                chunksize=chunkSize) 
            else:
                df_iter = pd.read_csv(args.input_contact_pairs,
                                sep="\t", 
                                header=None, 
                                iterator=True,
                                chunksize=chunkSize)
            df_iter.columns = ["cell_barcode", "chr1","bin1", "chr2", "bin2", "count"]
            parallel_results = []
            if ( args.input_format == "scVI-3d" ):
                for chunk in df_iter:
                    p.apply_async(
                        smoothContacts_scVIMC,
                            args=(
                                chunk,
                                cpg_dict,
                                barcode,
                                args.cis_only
                            ),
                            callback=parallel_results.append,
                            error_callback=utils.print_error,
                    ) 
                p.close()
                p.join()
                utils.eprint("[scAB::] Mapping finished")
                utils.eprint("[scAB::] Merging results")
            else:
                for chunk in df_iter:
                    p.apply_async(
                        smoothContactsMC,
                            args=(
                                chunk,
                                args.weighted,
                                cpg_dict,
                                args.res,
                                barcode,
                                args.cis_only
                            ),
                            callback=parallel_results.append,
                            error_callback=utils.print_error,
                    ) 
                p.close()
                p.join()
                utils.eprint("[scAB::] Mapping finished")
                utils.eprint("[scAB::] Merging results")

        mergeResults(parallel_results, args.input_format, args.output)
        spent_time = utils.secondsToText(time.time() - start_time)
        utils.eprint(f'[scAB::] scA/B matrix created using {spent_time}.')

if __name__ == "__main__":
    main()
