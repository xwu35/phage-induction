#!/usr/bin/env python

import os
import pandas as pd

# function to get variables and reads path
def parse_samples_and_reads(metadata, reads_dir=""):
    R1_map = {}
    R2_map = {}
    df = pd.read_table(metadata).set_index("sample")
    samples = df.index.tolist()
    for sample in samples:
      sample_R1 = df.loc[sample, "R1"]
      sample_R2 = df.loc[sample, "R2"]
      R1_map[sample] = os.path.join(reads_dir, sample_R1)
      R2_map[sample] = os.path.join(reads_dir, sample_R2)
    return samples, R1_map, R2_map 

def parse_sample_genome_and_sequences(metadata, genome_dir="", prophage_dir=""):
    genome_seq_map = {}
    sample_genome_dict = {}
    df = pd.read_table(metadata).set_index("sample")
    # link sample to genome seq and prophage file, consider only one genome per sample for now: this is used for mapping
    for sample in df.index:
      sample_genome_dict[sample] = [df.loc[sample, "genome"], os.path.join(genome_dir, df.loc[sample, "seq_name"]), os.path.join(prophage_dir, df.loc[sample, "prophage_file"])]
    # link genome to its sequence: this is used for blastn alignment
    # select genome and seq_name column and remove duplicate rows
    new_df = df[['genome', 'seq_name']].drop_duplicates(subset=['genome', 'seq_name']).set_index("genome")
    genomes = new_df.index.tolist()
    for genome in genomes:
      seq = new_df.loc[genome, "seq_name"]
      genome_seq_map[genome] = os.path.join(genome_dir, seq)
    return sample_genome_dict, genomes, genome_seq_map


def parse_metadata(metadata, reads_dir="", genome_dir="", prophage_dir=""):
    R1_map = {}
    R2_map = {}
    sample_genome_dict = {}
    genome_seq_map = {}
    df = pd.read_table(metadata).set_index("sample")
    samples = df.index.tolist()
    for sample in samples:
      # link sample to its sequence: this is used for trimming
      sample_R1 = df.loc[sample, "R1"]
      sample_R2 = df.loc[sample, "R2"]
      R1_map[sample] = os.path.join(reads_dir, sample_R1)
      R2_map[sample] = os.path.join(reads_dir, sample_R2)
      # link sample to genome seq and prophage file, consider only one genome per sample for now: this is used for mapping
      sample_genome_dict[sample] = [df.loc[sample, "genome"], os.path.join(genome_dir, df.loc[sample, "seq_name"]), os.path.join(prophage_dir, df.loc[sample, "prophage_file"])]
    # link genome to its sequence: this is used for blastn alignment
    # select genome and seq_name column and remove duplicate rows
    new_df = df[['genome', 'seq_name']].drop_duplicates(subset=['genome', 'seq_name']).set_index("genome")
    genomes = new_df.index.tolist()
    for genome in genomes:
      seq = new_df.loc[genome, "seq_name"]
      genome_seq_map[genome] = os.path.join(genome_dir, seq)
    return samples, R1_map, R2_map, sample_genome_dict, genomes, genome_seq_map