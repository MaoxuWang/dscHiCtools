## Tools for demultiplexed single-cell HiC data
### Author: Maoxu Wang, Xie Lab, Peking University

### 1. Install
The following packages are required:
- python >=3.6
- python-levenshtein>=0.12.0
- pysam>=0.16.0.1 
- multiprocess>=0.70.6.1
- pandas>=0.23.4
- pybktree==1.1
```shell
    python setup.py install
```

### 2. Add barcodes to read name
R2 file (or I2 file) records the cell barcode information, add it to the read name of R1 and R3
```shell
    dscHiCtools addBarcode \
        --read1 $Read1 \
        --read2 $Read2 \
        --read3 $Read3 \
        -o $outdir \
        --sampleName $sampleid \
        --threads n
```

#### outputs:
- .R1.barcoded.fastq.gz
- .R2.barcoded.fastq.gz
- .R3.barcoded.fastq.gz

### 3. Map barcodes from bam file
map cell barcodes to reference (true) and add "CB" tag
```shell
    dscHiCtools mapBarcode \
        --input_bam ${sampleid}.sorted.bam \
        --output_bam ${sampleid}.cellbarcoded.bam \
        --outMarkDuplicates \ ## whether to print fragment information to mark duplicates
        --min_mismatch 1 \ ## min mismatch allowed (hamming distance)
        --ref_barcode FILE \ ## barcode whitelist
        --threads n
```

#### outputs:
- output_bam
- QC/.mapped.tsv.gz : (tag/tcount)
- QC/.unmapped.tsv.gz
- [QC/.MarkDuplicates.txt.gz]

### 4. Filter barcodes from bam file 
only save barcodes listed to bam to file
```shell
    dscHiCtools filterBarcode \
        --input_bam output/test.cellbarcoded.bam  \
        --output_bam output/test.filtered.bam \
        --barcode_file barcodes.tsv.gz \
        --threads n
```
#### outputs
    .filtered.bam


### 5. Split into single cells
split bam into different files by cell barcodes

Set `ulimit -n 65535` if you encounter error of "Too many open files"
```shell
    dscHiCtools splitBarcode \
        --input_bam INPUT_BAM \
        --outdir "." \
        --barcode_file cellbarcode.txt \ ## cellbarcode
        --threads n 
```
#### outputs
    .sc.bam


### 6. Merge sub-cluster single cells to pseudo-bulk
merge and then convet to .hic and .mcool file 
```shell

    dscHiCtools sub2cool \
        --input_pairs mESC.contacts.pairs.txt.gz \
        --outdir "." \
        --metadata cellbarcode_group.metadata.txt \ ## No-header: cellbarcode, sub-cluster
        --threads n \
        --species mm10 \
        --resolution 10000 20000 50000 \ ## separated by " " 
        -n ## if set, don't normalize matrices 
```
default resolution choice: ["1000", "2000", "5000", "10000", "25000", "50000", "100000", "250000", "500000", "1000000", "2500000"]
#### outputs
- sub_cluster.contacts.pairs.txt.gz (4DN format)
- sub_cluster.hic
- sub_cluster.mcool
