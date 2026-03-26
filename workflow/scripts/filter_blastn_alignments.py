#!/usr/bin/env python3

import click
import pandas as pd
#from Bio import SeqIO 

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n python filter_blastn_alignment.py -i <BLASTn results> -o <Output file name>'
)
@click.option('-i',
    '--input',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='BLASTn results'
)
@click.option('-o',
    '--output',
    default="kept_alignment.txt",
    type=click.File("w"),
    show_default=True,
    help=('Output file name')
)

def process_data(input, output):
    """
    keep contigs with >=95% of sequences alinged to the genome with >=98% identity
    """
 
    # define the column names for blastn results
    # replaced qseqid with contig, and sseqid with chr for plotting
    columns = [
        'contig', 'chr', 'pident', 'length',
        'mismatch', 'gapopen', 'qstart', 'qend', 
        'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen'
    ]
    
    # read in the blastn results
    df = pd.read_csv(input, sep='\t', header=None, names=columns).sort_values(by=['contig','qstart'])

    # calculate hit_length/qlen
    df['hit_coverage'] = df['length']/df['qlen']

    # keep contigs with hit_coverage >= 0.95 and pident >= 98
    filtered_df = df[(df['hit_coverage'] >= 0.95) & (df['pident'] >= 98)].copy()

    # sstart and send might be swapped, find min as start and max as end
    filtered_df['start'] = filtered_df[['sstart', 'send']].min(axis=1)
    filtered_df['end'] = filtered_df[['sstart', 'send']].max(axis=1)
    filtered_df_selected =  filtered_df[['chr', 'start', 'end', 'contig']].sort_values(by='start')

    # save the filtered df
    filtered_df_selected.to_csv(output, sep='\t', index=False, header=True)

if __name__ == '__main__':
    process_data()
