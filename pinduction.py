#!/usr/bin/env python3

import sys
import os
import subprocess
import click
from datetime import datetime

version = "1.1.0"
@click.version_option(version, "--version", "-v")

def validate_test_run(ctx, param, value):
    """
    A callback to make --reads_dir, --sample_info, --genome_dir and --prophage_dir required if --test is not specified.
    """
    if not ctx.params.get('test') and value is None:
        raise click.BadParameter(f"Option '{param.name}' is required unless --test is used.")
    return value

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n python pinduction.py --sample_info <sample information table> '
    '--reads_dir <reads directory> --genome_dir <genome sequences directory> '
    '--prophage_dir <genome prophage prediction results directory> '
    '--output_dir <output directory>'
)
@click.option(
    '--sample_info',
    callback=validate_test_run,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help='Sample information table (tab separated).' 
    ' The table must contain six columns (sample, R1, R2, genome, seq_name, prophage_file)'
)
@click.option(
    '--reads_dir',
    callback=validate_test_run,
    type=click.Path(dir_okay=True, exists=True, resolve_path=True),
    help='Reads directory'
)
@click.option(
    '--genome_dir',
    callback=validate_test_run,
    type=click.Path(dir_okay=True, exists=True, resolve_path=True),
    help='Genome sequences directory'
)
@click.option(
    '--prophage_dir',
    callback=validate_test_run,
    type=click.Path(dir_okay=True, exists=True, resolve_path=True),
    help='Genome prophage prediction results directory'
)
@click.option(
    '--output_dir',
    default="OUTPUT",
    type=click.Path(dir_okay=True, resolve_path=True),
    show_default=True,
    help=('Output directory')
)
@click.option(
    '--adapter',
    default='',
    type=str,
    show_default=True,
    help=('Adapter sequences. By default, the NexteraPE-PE.fa file within `phage-induction/db/adapters` is used')
)
@click.option(
    '--mapper',
    default="bowtie2",
    type=str,
    show_default=True,
    help=('Mapping software; available options are: bowtie2, minimap2')
)
@click.option(
    '--step',
    default='fastqc',
    type=str,
    show_default=True,
    help=('Steps to run; available options are: fastqc, trimming, mapping, assemble, identification, alignment, annotation')
)
@click.option(
    '--test',
    is_flag=True,
    default=False,
    show_default=True,
    help='Test run'
)
@click.option(
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce'
)
@click.option(
    '--conda_envs',
    default='',
    show_default=True,
    help='Directory to store conda environments.'
    ' By default, the "conda_env" directory within `pinduction` is used'
)
@click.option(
    '--profile',
    default='slurm',
    show_default=True,
    help='Snakemake profile for cluster execution'
)

def run_pinduction(sample_info, reads_dir, genome_dir, prophage_dir, output_dir, adapter, mapper, step, test, dryrun, conda_envs, profile):
            
    # get snakefile, default adapter and default conda envs path
    script_dir=os.path.dirname(os.path.abspath(__file__))
    snakefile=os.path.join(script_dir, "workflow", "Snakefile")
    default_adapters=os.path.join(script_dir, "db", "adapters", "NexteraPE-PE.fa")
    default_envs=os.path.join(script_dir, "conda_envs")

    if test:
        sample_info=os.path.join(script_dir, "test_data", "sample_genome_prophage_info.tsv")
        reads_dir=os.path.join(script_dir, "test_data", "reads")
        genome_dir=os.path.join(script_dir, "test_data", "genome_sequences")
        prophage_dir=os.path.join(script_dir, "test_data", "prophage_prediction")
        output_dir="test_output"

    # write run log if it is not a dry run
    if not dryrun:
        os.makedirs(output_dir, exist_ok=True)
        logfile = os.path.join(output_dir, f"{os.path.basename(output_dir)}_run.log")
        with open(logfile, "w") as log:
            log.write("================pInduction run log==============\n")
            log.write(f"Start time: {datetime.now()}\n")
            log.write(f"pInduction version: {version}\n")
            log.write(f"Sample table: {sample_info}\n")
            log.write(f"Raw reads direcotry: {reads_dir}\n")
            log.write(f"Genome sequence direcotry: {genome_dir}\n")
            log.write(f"Genome prophage prediction direcotry: {prophage_dir}\n")
            log.write(f"Results directory: {output_dir}\n")
            log.write(f"Mapping software: {mapper}")
          
    cmd = (
        'snakemake --snakefile {snakefile} '
        '--use-conda --conda-frontend mamba '
        '--conda-prefix {envs} '
        '--profile {profile} --rerun-incomplete ' 
        '--printshellcmds --nolock --show-failed-logs '
        '{dryrun} '
        '--config sample_info={sample_info} reads_dir={reads} '
        'genome_dir={genome} prophage_dir={prophage} '
        'results_dir={results} adapter={adapter} '
        'mapper={mapper} step={step}'
        ).format(
            snakefile=snakefile,
            envs=default_envs if conda_envs=='' else conda_envs,
            profile=profile,
            dryrun='--dryrun' if dryrun else '',
            sample_info=sample_info,
            reads=reads_dir,
            genome=genome_dir,
            prophage=prophage_dir,
            results=output_dir,
            adapter=default_adapters if adapter=='' else adapter,
            mapper=mapper,
            step=step
            )

    # run snakemake with command-line config
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError:
        print("Snakemake failed. see log for details.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
  run_pinduction()
