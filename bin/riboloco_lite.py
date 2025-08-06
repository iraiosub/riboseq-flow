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
    counter = 0
    df_list = []
    for tx_id, tx_seq in tx_fasta_d.items():
        cds_start = int(info_dict[tx_id][3])  # 1 based
        starts = [m.start() for m in re.finditer('ATG', tx_seq)]  # 0 based
        starts1 = [a + 1 for a in starts]

        stops = [m.start() + 1 for m in re.finditer('TGA', tx_seq)]
        stops += [m.start() + 1 for m in re.finditer('TAA', tx_seq)]
        stops += [m.start() + 1 for m in re.finditer('TAG', tx_seq)]

        start_frames_list = [(a - cds_start) % 3 for a in starts1]
        start_frames = {start: frame for start, frame in zip(starts1, start_frames_list)}
        stop_frames_list = [(a - cds_start) % 3 for a in stops]
        stop_frames = {stop: frame for stop, frame in zip(stops, stop_frames_list)}

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
            all_linked_starts.setdefault(stop, []).append(start)

        frames = []
        annotated = []
        for stop, starts in all_linked_starts.items():
            if stop != 100_000_000:
                min_start = min(starts)
                furthest_starts[stop] = min_start
                frames.append((min_start - cds_start) % 3)
                annotated.append(min_start == cds_start)

        this_df = pd.DataFrame.from_dict({
            "transcript_id": tx_id,
            "orf_start": furthest_starts.values(),
            "orf_stop": furthest_starts.keys(),
            "orf_frame": frames,
            "annotated": annotated
        })

        df_list.append(this_df)
        counter += 1
        if counter >= max_tx:
            break
        if counter % 1000 == 0:
            print(str(counter) + " transcripts searched for ORFs")

    return pd.concat(df_list)


def load_tsv_into_dict(filename, key_column):
    data_dict = {}
    open_func = gzip.open if filename.endswith(".gz") else open
    with open_func(filename, 'rt') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)
        for row in reader:
            key = row[header.index(key_column)]
            data_dict[key] = row
    return data_dict


def read_fasta(filename):
    sequences = {}
    current_seq_name = None
    current_seq = ""
    open_func = gzip.open if filename.endswith(".gz") else open
    with open_func(filename, "rt") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_seq_name:
                    sequences[current_seq_name] = current_seq
                current_seq_name = line[1:]
                current_seq = ""
            else:
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

    info_basename = os.path.splitext(os.path.basename(args.info))[0]
    orf_df.to_csv(f"{info_basename}.orf_predictions.csv.gz", compression='gzip', index=False)

    # BAM processing
    chunk_size = 100
    bam_path = args.bam
    all_tx, all_pos, all_footprints, all_read_starts, all_read_names = [], [], [], [], []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        transcripts = bam.references
        for i in range(0, len(transcripts), chunk_size):
            chunk = transcripts[i:i + chunk_size]
            for tx in chunk:
                if tx not in info_dict:
                    continue
                tx_info = info_dict[tx]
                annotated_start = int(tx_info[3])
                tx_seq = tx_fasta_d[tx]
                for record in bam.fetch(tx):
                    if record.is_reverse:
                        continue
                    footprint = str(record.query_sequence)
                    footprint_start = record.reference_start
                    footprint_length = len(footprint)
                    footprint_frame = (footprint_start - annotated_start) % 3
                    mismatch = 'MM' if footprint[0] == tx_seq[footprint_start] else 'm'
                    # estimated_A_site = footprint_start + footprint_length - 12 - (footprint_length - 28) / 2
                    estimated_A_site = footprint_start + footprint_length - 12
                    
                    all_tx.append(tx)
                    all_pos.append(estimated_A_site)
                    all_footprints.append(f"{footprint_length}_{footprint_frame}_{mismatch}")
                    all_read_starts.append(footprint_start)
                    all_read_names.append(record.query_name)

    footprint_df = pd.DataFrame({
        'transcript_id': all_tx,
        'A_site_estimate': all_pos,
        'footprint_type': all_footprints,
        'read_start': all_read_starts,
        'read_name': all_read_names
    })

    joint_df = pd.merge(footprint_df, orf_df, on='transcript_id', how='inner')
    joint_df = joint_df[
        (joint_df['orf_start'] - args.flank <= joint_df['A_site_estimate']) &
        (joint_df['A_site_estimate'] <= joint_df['orf_stop'] + args.flank)
    ]

    # Reorder columns: existing + new read info at the end
    output_cols = ['transcript_id', 'A_site_estimate', 'footprint_type',
                   'orf_start', 'orf_stop', 'orf_frame', 'annotated',
                   'read_start', 'read_name']
    joint_df = joint_df[output_cols]

    joint_df.to_csv(args.output + '.riboloco.csv.gz', compression='gzip', index=False)

    # Summary
    joint_df_filtered = joint_df[
        (joint_df['orf_start'] <= joint_df['A_site_estimate']) &
        (joint_df['A_site_estimate'] <= joint_df['orf_stop'])
    ].copy()

    joint_df_filtered['n'] = joint_df_filtered.groupby(
        ['transcript_id', 'orf_start', 'orf_stop', 'footprint_type']
    )['footprint_type'].transform('size')

    joint_df_filtered = joint_df_filtered.reset_index().drop_duplicates(
        subset=['transcript_id', 'footprint_type', 'annotated', 'orf_start', 'orf_stop', 'n']
    )

    joint_df_filtered = joint_df_filtered.drop(['A_site_estimate', 'index'], axis=1)

    joint_df_filtered.to_csv(args.output + '.riboloco_summary.csv.gz', compression='gzip', index=False)


if __name__ == "__main__":
    main()
