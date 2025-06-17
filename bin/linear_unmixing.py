#!/usr/bin/env python3

# Author: Oscar Wilkins

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import argparse


def convert_frame_of_dict(in_frame_dict, frame_shift):
    d = {}
    for key, value in in_frame_dict.items():
        current_frame = int(key.split("_")[1])
        new_frame = (current_frame + frame_shift) % 3

        if 'm' in key or 'MM' in key:
            new_key = key.split("_")[0] + "_" + str(new_frame) + "_" + key.split("_")[2]
        else:
            new_key = key.split("_")[0] + "_" + str(new_frame)

        d[new_key] = value
    return d


def build_dist_vector(all_keys, frame_dict):
    dist_vector = []
    for key in all_keys:
        if key in frame_dict.keys():
            dist_vector.append(frame_dict[key])
        else:
            dist_vector.append(0) # no coverage in this frame
    return np.array(dist_vector)


def build_dist_vector(all_keys, frame_dict):
    dist_vector = []
    for key in all_keys:
        if key in frame_dict.keys():
            dist_vector.append(frame_dict[key])
        else:
            dist_vector.append(0) # no coverage in this frame
    return np.array(dist_vector)


def fit_weights(D, y):
    def objective(w, D, y):
        return np.linalg.norm(D @ w - y)**2

    x0 = np.array([1/3, 1/3, 1/3])
    constraints = [{'type': 'eq', 'fun': lambda w: np.sum(w) - 1}]
    bounds = [(0, 1)] * 3

    result = minimize(objective, x0, args=(D, y), bounds=bounds, constraints=constraints)
    return result.x


def kl_divergence(p, q, eps=1e-12):
    # Ensure no zero entries
    p = np.clip(p, eps, 1)
    q = np.clip(q, eps, 1)
    return np.sum(p * np.log(p / q))

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--expected_dist', '-e',
        type=str,
        required=True,
        help='Path to the expected distribution file'
    )

    parser.add_argument(
        '--all_footprints', '-a',
        type=str,
        required=True,
        help='Path to the file with all footprints across all ORFs'
    )

    parser.add_argument(
        '--output', '-o',
        type=str,
        required=True,
        help='Path to the output file'
    )

    parser.add_argument(
        '--bootstrap_number', '-b',
        type=int,
        required=False,
        default=5,
        help='Number of bootstraps to do'
    )

    parser.add_argument(
        '--min_count', '-m',
        type=int,
        required=False,
        default=1,
        help='Minimum footprints in ORF'
    )

    args = parser.parse_args()

    # Read in data and process

    actual_df = pd.read_csv(args.all_footprints)

    if actual_df.empty:
        print("No ORF footprint data found in input. Exiting.")
        exit(0)


    orf_ids = list(set(actual_df['orf_id']))

    expected_df = pd.read_csv(args.expected_dist)

    frame0_dict = dict(zip(expected_df['footprint_type'], expected_df['frac']))
    frame1_dict = convert_frame_of_dict(frame0_dict, 1)
    frame2_dict = convert_frame_of_dict(frame0_dict, 2)

    all_keys = set(frame1_dict) | set(frame2_dict) | set(frame0_dict)
    all_keys = list(all_keys) 

    frame0_dist = build_dist_vector(all_keys, frame0_dict)
    frame1_dist = build_dist_vector(all_keys, frame1_dict)
    frame2_dist = build_dist_vector(all_keys, frame2_dict)

    d1 = frame0_dist.astype(float)
    d2 = frame1_dist.astype(float)
    d3 = frame2_dist.astype(float)

    d1 /= d1.sum()
    d2 /= d2.sum()
    d3 /= d3.sum()

    D = np.column_stack([d1, d2, d3])

    # Process each orf
    all_results = {}

    grouped = actual_df.groupby('orf_id')

    for i, this_orf in enumerate(orf_ids):

        if i % 1000 == 0:
            print(i)

        df = grouped.get_group(this_orf)

        this_orf_dict = dict(zip(df['footprint_type'], df['footprint_type_n']))

        y = build_dist_vector(all_keys, this_orf_dict)

        y = y.astype(float)

        total_y = int(y.sum())

        if total_y < args.min_count:
            continue

        y /= total_y

        n_bootstrap = args.bootstrap_number
        weights_samples = []
        kl_divergences = []

        probabilities = y

        for _ in range(n_bootstrap):
            # Sample counts with replacement from bins according to y probabilities
            resampled_indices = np.random.choice(len(y), size=total_y, p=probabilities)

            y_resampled = np.random.multinomial(total_y, probabilities).astype(float)
            y_resampled /= y_resampled.sum()


            w = fit_weights(D, y_resampled)
            weights_samples.append(w)

            y_fit = D @ w
            y_fit /= y_fit.sum()  # normalize just in case

            kl = kl_divergence(y_resampled, y_fit)
            kl_divergences.append(kl)

        weights_samples = np.array(weights_samples)
        weights_mean = weights_samples.mean(axis=0)
        weights_std = weights_samples.std(axis=0)

        result_dict = {
            'mean_weight_frame_0': weights_mean[0],
            'mean_weight_frame_1': weights_mean[1],
            'mean_weight_frame_2': weights_mean[2],
            'std_weight_frame_0': weights_std[0],
            'std_weight_frame_1': weights_std[1],
            'std_weight_frame_2': weights_std[2],
            'kl_mean': np.mean(kl_divergences),
            'kl_std': np.std(kl_divergences),
        }

        all_results[this_orf] = result_dict


    df_results = pd.DataFrame.from_dict(all_results, orient='index')

    df_results.index.name = 'orf_id'       # name the index
    df_results = df_results.reset_index()  # turn index into a column

    df_results.to_csv(args.output)

if __name__ == "__main__":
    main()
