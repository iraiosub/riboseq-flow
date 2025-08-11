#!/usr/bin/env python3

# Author: Oscar Wilkins

import pysam
import argparse
import re
import pandas as pd
import gzip
import csv
import os


def find_all_orfs(tx_fasta_d, info_dict, max_tx=1_000_000):
    """
    Finds all ORFs with canonical start codon and stops codons. Values are 1 based!
    """
    counter = 0
    df_list = []
    for tx_id, tx_seq in tx_fasta_d.items():
        cds_start = int(info_dict[tx_id][3])  # 1 based
        starts = [m.start() for m in re.finditer('ATG', tx_seq)]  # 0 based
        starts1 = [a + 1 for a in starts]

        stops = [m.start() + 1 for m in re.finditer('TGA', tx_seq)]  # first nt of stop codon
        stops += [m.start() + 1 for m in re.finditer('TAA', tx_seq)]
        stops += [m.start() + 1 for m in re.finditer('TAG', tx_seq)]

        start_frames_list = [(a - cds_start) % 3 for a in starts1]
        start_frames = {}
        for start, frame in zip(starts1, start_frames_list):
            start_frames[start] = frame

        stop_frames_list = [(a - cds_start) % 3 for a in stops]
        stop_frames = {}
        for stop, frame in zip(stops, stop_frames_list):
            stop_frames[stop] = frame

        first_stop = {}
        for start in starts1:
            z = [100_000_000]
            try:
                first_stop[start] = min(z + [a for a in stops if a > start and start_frames[start] == stop_frames[a]])
            except KeyError:
                print(start)
                print([a for a in stops if a not in stop_frames.keys()])
                print(start_frames)
                print(stop_frames)
                assert 1 == 0

        all_linked_starts = {}
        furthest_starts = {}
        for start, stop in first_stop.items():
            try:
                all_linked_starts[stop].append(start)
            except KeyError:
                all_linked_starts[stop] = [start]

        frames = []
        annotated = []
        for stop, starts in all_linked_starts.items():
            if stop != 100_000_000:
                furthest_starts[stop] = min(starts)
                frames.append((min(starts) - cds_start) % 3)
                annotated.append(min(starts) == cds_start)

        this_df = pd.DataFrame.from_dict({"transcript_id": tx_id, "orf_start": furthest_starts.values(),
                                          "orf_stop": furthest_starts.keys(),
                                          "orf_frame": frames, "annotated": annotated})

        df_list.append(this_df)

        counter += 1

        if counter >= max_tx:
            break

        if counter % 1000 == 0:
            print(str(counter) + " transcripts searched for ORFs")

    orf_df = pd.concat(df_list)

    return orf_df



def load_tsv_into_dict(filename, key_column):
    """
    Loads TSV data into a dictionary for O(1) lookup based on a specified key,
    skipping the header line.

    Args:
        filename: The path to the TSV file.
        key_column: The name of the column to use as the dictionary key.

    Returns:
        A dictionary where the keys are values from the key_column and the
        values are the entire rows of data (as lists).
    """
    data_dict = {}

    open_func = gzip.open if filename.endswith(".gz") else open
    with open_func(filename, 'rt') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)  # Extract the header line
        for row in reader:
            key = row[header.index(key_column)]
            data_dict[key] = row
    return data_dict


def read_fasta(filename):
  """
  This function reads a FASTA file and stores the sequences in a dictionary.

  Args:
      filename: The path to the FASTA file.

  Returns:
      A dictionary where the key is the sequence name and the value is the sequence.
  """
  sequences = {}
  current_seq_name = None
  current_seq = ""

  open_func = gzip.open if filename.endswith(".gz") else open

  with open_func(filename, "rt") as f:
    for line in f:
      line = line.rstrip()  # Remove trailing whitespace
      if line.startswith(">"):
        # New sequence
        if current_seq_name:
          sequences[current_seq_name] = current_seq
        current_seq_name = line[1:]  # Get sequence name without '>'
        current_seq = ""
      else:
        # Add line to current sequence
        current_seq += line

  if current_seq_name:
    sequences[current_seq_name] = current_seq
  return sequences





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", "-b")
    parser.add_argument("--fasta")
    parser.add_argument("--info")
    parser.add_argument("--output")
    parser.add_argument('--flank', type=int, default=300)
    args = parser.parse_args()


    tx_fasta_d = read_fasta(args.fasta)

    info_dict = load_tsv_into_dict(args.info, 'transcript_id')

    orf_df = find_all_orfs(tx_fasta_d, info_dict)

    # Save ORFs using the basename of the --info file
    info_basename = os.path.splitext(os.path.basename(args.info))[0]
    orf_df.to_csv(f"{info_basename}.orf_predictions.csv.gz", compression='gzip', index=False)


    x = 0
    bam_tx_list = []
    bam_pos_list = []
    bam_footprint_list = []

    chunk_size = 100
    bam_path = args.bam

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        transcripts = bam.references

        all_tx = []
        all_pos = []
        all_footprints = []

        for i in range(0, len(transcripts), chunk_size):
            chunk = transcripts[i:i+chunk_size]

            # For this chunk, temporary lists
            chunk_tx_list = []
            chunk_pos_list = []
            chunk_footprint_list = []

            for tx in chunk:
                if tx not in info_dict:
                    continue

                for record in bam.fetch(tx):
                    if record.is_reverse:
                        continue

                    footprint = str(record.query_sequence)
                    footprint_start = record.reference_start
                    footprint_length = len(footprint)
                    tx_info = info_dict[tx]
                    annotated_start = int(tx_info[3])
                    footprint_frame = (footprint_start - annotated_start) % 3

                    if footprint[0] == tx_fasta_d[tx][footprint_start]:
                        mismatch = 'MM'
                    else:
                        mismatch = 'm'

                    estimated_A_site = footprint_start + footprint_length - 12
                    # estimated_A_site = footprint_start + footprint_length - 12 - (footprint_length - 28)/2

                    chunk_tx_list.append(tx)
                    chunk_pos_list.append(estimated_A_site)
                    chunk_footprint_list.append(f"{footprint_length}_{footprint_frame}_{mismatch}")

            # Append chunk data to all data
            all_tx.extend(chunk_tx_list)
            all_pos.extend(chunk_pos_list)
            all_footprints.extend(chunk_footprint_list)

        # After all chunks processed, build full DataFrame once
        footprint_df = pd.DataFrame({
            'transcript_id': all_tx,
            'A_site_estimate': all_pos,
            'footprint_type': all_footprints
        })

        # Proceed with merging and filtering as before:
        joint_df = pd.merge(footprint_df, orf_df, on='transcript_id', how='inner')

        joint_df = joint_df[
            (joint_df['orf_start'] - args.flank <= joint_df['A_site_estimate']) &
            (joint_df['A_site_estimate'] <= joint_df['orf_stop'] + args.flank)
        ]

        joint_df.to_csv(args.output + '.riboloco.csv.gz', compression='gzip', index=False)

        joint_df = joint_df[
            (joint_df['orf_start'] <= joint_df['A_site_estimate']) &
            (joint_df['A_site_estimate'] <= joint_df['orf_stop'])
        ]

        joint_df['n'] = joint_df.groupby(
            ['transcript_id', 'orf_start', 'orf_stop', 'footprint_type']
        )['footprint_type'].transform('size')

        joint_df = joint_df.reset_index().drop_duplicates(
            subset=['transcript_id', 'footprint_type', 'annotated', 'orf_start', 'orf_stop', 'n']
        )

        joint_df = joint_df.drop(['A_site_estimate', 'index'], axis=1)

        # Output the final summary CSV
        joint_df.to_csv(args.output + '.riboloco_summary.csv.gz', compression='gzip', index=False)




    # annotated_only_df = joint_df[joint_df['annotated'] == 1]
    # annotated_only_df['n2'] = annotated_only_df.reset_index().groupby('footprint_type')['n'].transform('sum')
    # print(annotated_only_df)




main()





