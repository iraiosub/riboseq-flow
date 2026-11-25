#!/usr/bin/env python3

"""Summarise the retained riboWaltz P-sites for MultiQC."""

import argparse
import csv
import gzip
from collections import Counter
from html import escape
import json
from pathlib import Path


REGIONS = {"5utr": "5' UTR", "cds": "CDS", "3utr": "3' UTR"}
WINDOWS = {"start": (-25, 40), "stop": (-40, 25)}
PARENT = {
    "parent_id": "psite_qc", "parent_name": "riboWaltz: P-site QC",
    "section_href": "https://github.com/LabTranslationalArchitectomics/riboWaltz",
    "parent_description": "P-site positioning and diagnostics from "
                          "<a href='https://github.com/LabTranslationalArchitectomics/riboWaltz'>riboWaltz</a> "
                          "(<a href='https://doi.org/10.1371/journal.pcbi.1006169'>Lauria et al., 2018</a>). "
                          "Plots pool the read lengths retained after length and periodicity filtering for each sample.",
}


def summarise(paths):
    samples = {}
    for path in sorted(paths):
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = row["sample"]
                if sample not in samples:
                    samples[sample] = {
                        "lengths": set(), "transcripts": set(), "regions": Counter(),
                        "frames": {region: Counter() for region in REGIONS},
                        "start": Counter(), "stop": Counter(),
                    }
                counts = samples[sample]
                counts["lengths"].add(int(row["length"]))
                if row.get("transcript"):
                    counts["transcripts"].add(row["transcript"])
                region = row["psite_region"]
                if region not in REGIONS:
                    continue
                start = int(row["psite_from_start"])
                stop = int(row["psite_from_stop"])
                counts["regions"][region] += 1
                counts["frames"][region][start % 3] += 1
                for end, position in (("start", start), ("stop", stop)):
                    lower, upper = WINDOWS[end]
                    if lower <= position <= upper:
                        counts[end][position] += 1
    return dict(sorted(samples.items()))


def write_plot(identifier, title, kind, data, description, parent, categories=None, **config):
    if not data:
        return
    content = dict(parent, id=identifier, section_name=title, description=description,
                   plot_type=kind, data=data,
                   pconfig=dict(id=identifier + "_plot", title="P-sites: " + title, **config))
    if categories is not None:
        content["categories"] = categories
    Path(identifier + "_mqc.json").write_text(json.dumps(content, indent=2) + "\n")


def main(paths, transcripts="representative"):
    samples = summarise(paths)
    retained = "".join("<li><strong>{}</strong>: {} nt{}</li>".format(
        escape(sample), ", ".join(map(str, sorted(counts["lengths"]))),
        "; {} transcripts with retained P-sites".format(len(counts["transcripts"])) if counts["transcripts"] else "")
        for sample, counts in samples.items())
    scope = "All BAM-observed transcripts with an annotated CDS" if transcripts == "all_cds" else "Representative/custom transcript reference"
    parent = dict(PARENT, extra="<p><strong>Transcript selection:</strong> " + scope + ".</p>"
                               "<p><strong>Read lengths retained for riboWaltz</strong></p><ul>" + retained + "</ul>")
    if transcripts == "all_cds":
        parent["parent_description"] += " Counts represent transcript assignments; genomic tracks count each unambiguous read once."

    regions = {sample: {label: counts["regions"][region] for region, label in REGIONS.items()}
               for sample, counts in samples.items() if sum(counts["regions"].values())}
    write_plot("psite_regions", "Transcript regions", "bargraph", regions,
               "Proportion of assigned P-sites in each transcript region. Switch to counts "
               "to view the underlying totals.", parent,
               categories={label: {} for label in REGIONS.values()},
               stacking="percent", cpswitch=True, cpswitch_c_active=False, ylab="P-sites (%)",
               height=max(240, 120 + 28 * len(regions)))

    frames = {sample + " | " + label: {"Frame " + str(frame): counts["frames"][region][frame]
                                      for frame in range(3)}
              for sample, counts in samples.items() for region, label in REGIONS.items()
              if sum(counts["frames"][region].values())}
    write_plot("psite_frames", "Reading frames", "bargraph", frames,
               "P-site frames relative to the annotated CDS start, shown separately for "
               "each sample and region. Regions without assigned P-sites are omitted.", parent,
               categories={"Frame " + str(frame): {} for frame in range(3)},
               stacking="percent", cpswitch=True, cpswitch_c_active=False, ylab="P-sites (%)",
               sort_samples=False, height=max(280, 120 + 26 * len(frames)))

    for end, (lower, upper) in WINDOWS.items():
        profiles = {}
        for sample, counts in samples.items():
            total = sum(counts[end].values())
            if total:
                profiles[sample] = {position: counts[end][position] / total
                                    for position in range(lower, upper + 1)}
        reference = "CDS start" if end == "start" else "CDS end"
        write_plot("psite_" + end, "Around " + end + " codon", "linegraph", profiles,
                   "Each sample is normalised within this window. Zero marks the annotated "
                   + reference + ". Samples with no P-sites here are omitted.", parent,
                   xlab="Distance from " + reference + " (nt)", ylab="P-site frequency",
                   xmin=lower, xmax=upper, ymin=0, xDecimals=False, height=340)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transcripts", choices=("representative", "all_cds"), default="representative")
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args()
    main(args.paths, args.transcripts)
