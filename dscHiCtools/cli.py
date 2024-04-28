from dscHiCtools.processing import (
    mapBarcode,
    filterBarcode,
    splitBarcode,
    sub2cool,
    addTags2bam
)
from dscHiCtools import utils


@utils.log_info
def run_mapBarcode(options):
    """
        Wraps the mapBarcode function for use on the command line
    """
    mapBarcode(
        read1=options.read1,
        index=options.index,
        read2=options.read2,
        start_idx=options.start_idx,
        end_idx=options.end_idx,
        outdir=options.outdir,
        threads=options.threads,
        sampleName=options.sampleName,
        ref_barcode=options.ref_barcode,
        min_mismatch=options.min_mismatch,
        chemistry=options.chemistry,
        suffix=options.suffix,
        n_lines=options.n_lines,
        multiome_RNA_barcode=options.multiome_RNA_barcode,
        multiome_DNA_barcode=options.multiome_DNA_barcode
    )


@utils.log_info
def run_addTag(options):
    """
        Wraps the addTags2bam function for use on the command line
    """
    addTags2bam(
        input_bam=options.input_bam,
        output_bam=options.output_bam,
        threads=options.threads,
        samtools=options.samtools
    )


@utils.log_info
def run_filterBarcode(options):
    """
        Wraps the filterBarcode function for use on the command line
    """
    filterBarcode(
        input_bam=options.input_bam,
        output_bam=options.output_bam,
        threads=options.threads,
        barcode_f=options.barcode_file,
        samtools=options.samtools
    )


@utils.log_info
def run_splitBarcode(options):
    """
        Wraps the splitBarcode function for use on the command line
    """
    splitBarcode(
        input_bam=options.input_bam,
        outdir=options.outdir,
        threads=options.threads,
        barcode_f=options.barcode_file,
        samtools=options.samtools
    )


@utils.log_info
def run_sub2cool(options):
    """
        Wraps the splitBarcode function for use on the command line
    """
    sub2cool(
        input_pairs=options.input_pairs,
        outdir=options.outdir,
        metadata=options.metadata,
        threads=options.threads,
        resolution=options.resolution,
        cooler=options.cooler,
        hic2cool=options.hic2cool,
        juicer=options.juicer,
        species=options.species,
        normalize=options.n
    )