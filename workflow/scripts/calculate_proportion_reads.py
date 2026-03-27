#!/usr/bin/env python3

import os
import click
import pandas as pd

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n python calculate_proportion_reads.py -i <Read counts by prophage region> -t <Trimmed read counts> -o <Output file name>'
)
@click.option('-i',
    '--input',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='Read counts by prophage region'
)
@click.option('-t',
    '--trimmed',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='Number of trimmed reads'
)
@click.option('-o',
    '--output',
    default="region_counts_proportion.txt",
    type=click.File("w"),
    show_default=True,
    help=('Output file name')
)

def process_data(input, trimmed, output):
    """
    Calculate the proportion of reads mapped to each prophage region by dividing by the total number of trimmed reads
    """
    
    # read in the number of trimmed reads
    trimmed = pd.read_csv(trimmed, sep='\t')[['file', 'num_seqs']]
    trimmed['sample'] = trimmed['file'].apply(lambda x: os.path.basename(x).split('_R1.trimmomatic.fastq.gz')[0])
    trimmed['total_reads'] = trimmed['num_seqs']*2
    total_reads_dict = dict(zip(trimmed['sample'], trimmed['total_reads']))

    # read in read counts by region
    columns = ["chr", "start", "end", "region", "counts"]
    region = pd.read_csv(input, sep='\t', header=None, names=columns)
    samplenames = os.path.basename(input).split('_counts.txt')[0]

    # calculate proportion of reads
    region['pct_reads'] = round(region['counts']/total_reads_dict[samplenames], 3)

    # save the filtered df
    region.to_csv(output, sep='\t', index=False, header=True)

if __name__ == '__main__':
    process_data()
