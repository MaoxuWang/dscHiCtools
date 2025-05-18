import os
import gzip
import shutil
import pandas as pd
from scipy import io


def write_to_files(sparse_matrix, top_cells, ordered_tags_map, data_type, outfolder):
    """Write the umi and read sparse matrices to file in gzipped mtx format.

    Args:
        sparse_matrix (dok_matrix): Results in a sparse matrix.
        top_cells (set): Set of cells that are selected for output.
        ordered_tags_map (dict): Tags in order with indexes as values.
        data_type (string): A string definning if the data is umi or read based.
        outfolder (string): Path to the output folder.
    """
    prefix = os.path.join(outfolder,data_type + '_count')
    os.makedirs(prefix, exist_ok=True)
    io.mmwrite(os.path.join(prefix,'matrix.mtx'),sparse_matrix)
    with gzip.open(os.path.join(prefix,'barcodes.tsv.gz'), 'wb') as barcode_file:
        for barcode in top_cells:
            barcode_file.write('{}\n'.format(barcode).encode())
    with gzip.open(os.path.join(prefix,'features.tsv.gz'), 'wb') as feature_file:
        for feature in ordered_tags_map:
            feature_file.write('{}\n'.format(feature).encode())
    with open(os.path.join(prefix,'matrix.mtx'),'rb') as mtx_in:
        with gzip.open(os.path.join(prefix,'matrix.mtx') + '.gz','wb') as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove(os.path.join(prefix,'matrix.mtx'))


def write_dense(sparse_matrix, index, columns, outfolder, filename):
    """
    Writes a dense matrix in a csv format
    
    Args:
       sparse_matrix (dok_matrix): Results in a sparse matrix.
       index (list): List of TAGS
       columns (set): List of cells
       outfolder (str): Output folder
       filename (str): Filename
    """
    prefix = os.path.join(outfolder)
    os.makedirs(prefix, exist_ok=True)
    pandas_dense = pd.DataFrame(sparse_matrix.todense(), columns=columns, index=index)
    pandas_dense.to_csv(os.path.join(outfolder,filename), sep='\t')


def write_mapped(
    merged_match, 
    outfolder, 
    filename):
    """
    Writes a list of top unmapped sequences

    Args:
        merged_no_match (Counter): Counter of unmapped sequences
        outfolder (string): Path of the output folder
        filename (string): Name of the output file
        top_unknowns (int): Number of unmapped sequences to output
    Returns:
        file outpu: sample.unmapped.tsv.gz
    """
    _num = 0
    for element in merged_match:
        _num += merged_match[element]

    with gzip.open(os.path.join(outfolder, filename),'wb') as known_file:
        for element in merged_match:
            known_file.write(f'{element}\t{str(merged_match[element])}\n'.encode())
    return _num


def write_unmapped(
    merged_no_match, 
    outfolder, 
    filename,
    top_unknowns=500):
    """
    Writes a list of top unmapped sequences

    Args:
        merged_no_match (Counter): Counter of unmapped sequences
        outfolder (string): Path of the output folder
        filename (string): Name of the output file
        top_unknowns (int): Number of unmapped sequences to output
    Returns:
        file outpu: sample.unmapped.tsv.gz
    """
    _num = 0
    for element in merged_no_match:
        _num += merged_no_match[element]

    top_unmapped = merged_no_match.most_common(top_unknowns)

    with gzip.open(os.path.join(outfolder, filename),'wb') as unknown_file:
        for element in top_unmapped:
            unknown_file.write(f'{element[0]}\t{str(element[1])}\n'.encode())
    return _num


def write_4DN_contacts(
    df, 
    outfile_contacts_pairs):
    """
    Writes a list of top unmapped sequences

    Args:
        df (DataFrame): dataframe 
        outfile_contacts_pairs (string): Path of the output file
    Returns:
        file outpu: sample.contacts.pairs.txt.gz
    """
    with gzip.open(outfile_contacts_pairs,'wb') as contacts_file:
        contacts_file.write('## pairs format v1.0\n'.encode())
        contacts_file.write('#columns: readID chr1 position1 chr2 position2 strand1 strand2\n'.encode())
        for _, row in df.iterrows():
            line = "\t".join(list(map(str, row)))
            contacts_file.write(f'{line}\n'.encode())
