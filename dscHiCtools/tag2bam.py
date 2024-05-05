import pysam
from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from dscHiCtools.chunkFiles import chunk_bam


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


@utils.log_info
def main(args):
    """
        add cell barcode tags to bam file
        input:
            bam file, the header of read name is barcode 
        output:
            bam file, with "CB:" tag added
    """
    inputBam = pysam.AlignmentFile(args.input_bam, "rb")
    idx_file = args.input_bam + ".bai"
    if not os.path.exists(idx_file):
        utils.eprint("[addTag::] Indexing input bam file...")
        Args_m = [f'{args.samtools} index -@ {args.threads} {args.input_bam}']
        subprocess.check_call(Args_m, shell=True)
    intervals = chunk_bam(inputBam, args.threads)

    inputBam.close()
    p = multiprocess.Pool(args.threads)

    utils.eprint("[addTag::] Adding CB tags to bam file...")
    utils.eprint(f"[addTag::] Multi-cores: {args.threads}")
    i = 0
    bamTempFiles = []
    for interval in intervals.values():
        i += 1
        output_bam_temp = args.output_bam + str(i) + ".bam"
        bamTempFiles.append(output_bam_temp)
        p.apply_async(
            addCBtag,
            args=(
                interval,
                args.input_bam,
                output_bam_temp
            ),
            error_callback=utils.print_error,
        )

    p.close()
    p.join()
    utils.eprint("[addTag::] Filtering done.")
    utils.eprint("[addTag::] Merging results...")
    
    tempFiles = " ".join(bamTempFiles)
    Args_m = [f'{args.samtools} merge -fpc -@ {args.threads} {args.output_bam} {tempFiles}']
    subprocess.check_call(Args_m, shell=True)

    if os.path.exists(args.output_bam):
        [os.remove(i) for i in bamTempFiles]
    else:
        raise Exception("[addTag::] samtools merge failed, temp files not deleted")