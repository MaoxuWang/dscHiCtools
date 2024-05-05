from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from dscHiCtools.chunkFiles import chunk_dataframe
import gzip
import pandas as pd 


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
    utils.eprint(f"[sub2cool::] Read {n} cell barcodes.")
    return cluster_dict


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
    cluster_dict = pd.read_csv(args.metadata,
                        sep="\t")   
    groups = set(list(cluster_dict['group']))
    
    df_iter = chunk_dataframe(args.input_pairs, args.threads)
    parallel_results = []
    p = multiprocess.Pool(processes=args.threads)
    for chunk in df_iter:
        p.apply_async(
            sub_cluster2cool,
                args=(
                chunk,
                groups,
                cluster_dict,
            ),
            callback=parallel_results.append,
            error_callback=utils.print_error,
        ) 
    p.close()
    p.join()
    utils.eprint("[sub2cool::] Merging results...")

    for _group in groups:
        utils.eprint(f"[sub2cool::] converting _group...")

        outfile_contacts_pairs = os.path.join(args.outdir, _group + ".contacts.pairs.txt.gz")
        wirte_all_group(
            outfile_contacts_pairs,
            _group,
            parallel_results)

        outfile_cooler = os.path.join(args.outdir, _group + "." + ".cool")

        ## convert
        resolution_symbol = ",".join(args.resolution)
        ## juicer .hic
        utils.eprint("[sub2cool::] convert to .hic...")
        if not args.normalize:
            utils.eprint("[sub2cool::] Normalizing matrics...")
            Args_m = [f'{args.juicer} pre --threads {args.threads} -r {resolution_symbol} {outfile_contacts_pairs} {args.outdir + "/" + _group}.hic {args.species}']
        else:
            Args_m = [f'{args.juicer} pre -n --threads {args.threads} -r {resolution_symbol} {outfile_contacts_pairs} {args.outdir + "/" + _group}.hic {args.species}']
        subprocess.check_call(Args_m, shell=True)

        ## convert to .mcool
        utils.eprint("[sub2cool::] convert to .mcool...")

        Args_m = [f'{args.hic2cool} convert -r 0 {args.outdir + "/" + _group}.hic {args.outdir + "/" + _group}.mcool']
        subprocess.check_call(Args_m, shell=True)
    
        utils.eprint("[sub2cool::] balancing...")
        for _resolution in args.resolution:
            Args_m = [f'{args.cooler} balance -p {args.threads} {args.outdir + "/" + _group}.mcool::/resolutions/{_resolution}']
            subprocess.check_call(Args_m, shell=True)