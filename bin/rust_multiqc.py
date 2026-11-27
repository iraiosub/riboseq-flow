#!/usr/bin/env python3

"""Build compact RUST-inspired QC views from full profiles and read-length support."""

import argparse
import csv
import gzip
from html import escape
import json
from pathlib import Path


MIN_PSITES = 1000
PARENT = {
    "parent_id": "rust_qc", "parent_name": "RUST-inspired QC",
    "section_href": "https://doi.org/10.1038/ncomms12915",
    "parent_description": "A local RUST-inspired reimplementation of codon-context sequence-bias QC, "
                          "using read-count-weighted codon frequencies. Original method: "
                          "<a href='https://doi.org/10.1038/ncomms12915'>O'Connor, Andreev and Baranov (2016)</a>. "
                          "Displayed lengths require at least 1,000 CDS P-sites after subsampling and "
                          "dominant-frame selection. Full PDF and TSV outputs retain all analysed lengths.",
}


def read_table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_report(qc_summary, profiles, lengths):
    counts = {(row["sample"], int(row["length"])): int(row["psites_used"])
              for path in lengths for row in read_table(path)}
    data = {}
    for path in profiles:
        for row in read_table(path):
            key = (row["sample"], int(row["length"]))
            if counts.get(key, 0) >= MIN_PSITES:
                data.setdefault(key, {})[int(row["offset_codon"])] = float(row["kld"])

    def write_section(identifier, title, kind, values, description, xcats=None, ycats=None, **config):
        content = dict(PARENT, id=identifier, section_name=title, plot_type=kind,
                       data=values, description=description, xcats=xcats, ycats=ycats,
                       pconfig=dict(id=identifier + "_plot", title="RUST-inspired QC: " + title, **config))
        Path(identifier + "_mqc.json").write_text(json.dumps(content, indent=2) + "\n")

    if not data:
        write_section("rust_support", "Read-length support", "html",
                      "<p>No read length meets the 1,000-P-site display minimum.</p>", "")
        return

    qc = {row["name"]: row for row in read_table(qc_summary)}
    frequent_profiles = {}
    for identifier, title, column, description in (
        ("rust_frequent", "Most frequent read length", "most_abundant_length",
         "The most frequent useful read length per sample, as selected in Ribo-seq QC."),
        ("rust_periodic", "Most periodic read length", "most_periodic_length",
         "The most periodic expected read length from Ribo-seq QC, with at least 10% of the "
         "count at the most frequent useful length."),
    ):
        selected = {}
        unavailable = []
        for sample in sorted(qc):
            value = qc[sample][column]
            length = int(float(value)) if value not in ("", "NA") else None
            key = (sample, length)
            if key in data:
                selected[sample + " | " + str(length) + " nt"] = data[key]
            else:
                unavailable.append(escape(sample) + (": " + str(length) + " nt" if length is not None else ": no selection"))
        if identifier == "rust_frequent":
            frequent_profiles = selected
        elif selected and selected == frequent_profiles:
            description += " These selections coincide with the most frequent lengths for the reported samples."
        if unavailable:
            description += " Unavailable after RUST display filtering: " + "; ".join(unavailable) + "."
        if selected:
            write_section(identifier, title, "linegraph", selected, description,
                          xlab="Codon position relative to P-site", ylab="Sequence bias (KLD)",
                          xDecimals=False, height=320, tt_label="{point.x} codons: {point.y:.6f} KLD")
        else:
            write_section(identifier, title, "html", "<p>" + description + "</p>", "")

    keys = sorted(data)
    offsets = sorted({position for points in data.values() for position in points})
    matrix = [[data[key].get(position) for position in offsets] for key in keys]
    sample_labels = {sample: "S" + str(index + 1) for index, sample in enumerate(sorted({key[0] for key in keys}))}
    labels = [sample_labels[sample] + " | " + str(length) + " nt | " + sample for sample, length in keys]
    legend = "; ".join(label + " = " + escape(sample) for sample, label in sample_labels.items())
    write_section("rust_all_lengths", "All supported read lengths", "heatmap", matrix,
                  "Each row is one sample and read length; colour shows KLD at each codon position. "
                  "Values use the same scale across rows. " + legend + ". Hover for full labels and values.",
                  xcats=list(map(str, offsets)), ycats=labels,
                  xTitle="Codon position relative to P-site", yTitle={"text": "Sample | read length (nt)"},
                  height=max(400, 120 + 22 * len(keys)), square=False,
                  xcats_samples=False, ycats_samples=True, decimalPlaces=6)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc-summary", required=True)
    parser.add_argument("--profiles", nargs="+", required=True)
    parser.add_argument("--lengths", nargs="+", required=True)
    args = parser.parse_args()
    build_report(args.qc_summary, args.profiles, args.lengths)
