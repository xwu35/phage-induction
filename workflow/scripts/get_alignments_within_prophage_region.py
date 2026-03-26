#!/usr/bin/env python3

import click
import pandas as pd

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n python filter_blastn_alignment.py -b <Filtered blastn hits> -p <Prophage region coordinates> -o <Output file name>'
)
@click.option('-b',
    '--blastn',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='Filtered blastn hits'
)
@click.option('-p',
    '--prophage',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='Prophage region coordinates'
)
@click.option('-o',
    '--output',
    default="blastn_alignments_within_prophage_regions.txt",
    type=click.File("w"),
    show_default=True,
    help=('Output file name')
)

def process_data(blastn, prophage, output):
    """
    Keep alignments within prophage region
    """

    # read in filtered blastn hits
    blastn = pd.read_csv(blastn, sep='\t')
    
    # read in prophage region coordinates, bed file used
    columns = ['chr', 'start_p', 'end_p', 'region']
    prophage = pd.read_csv(prophage, sep='\t', header=None, names=columns)
   
    # merge dfs
    merged = pd.merge(blastn, prophage, on = 'chr', how = 'left')
    # if the alignment is within prophage region, then the start and end has to be within porphage coordinates
    # since bed file is used, +1 to convert to 1-based coordinates
    filtered_merged = merged[(merged['start'] >= merged['start_p']+1) & (merged['end'] <= merged['end_p'])].copy()

    # save the filtered df
    filtered_merged.to_csv(output, sep='\t', index=False, header=True)

if __name__ == '__main__':
    process_data()
