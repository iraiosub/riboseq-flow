#!/usr/bin/env python3

"""Collate nf-core-style versions.yml files for publication and native MultiQC reporting."""

from pathlib import Path
import sys

import yaml


def collect_versions(paths, nextflow_version, pipeline_version):
    versions = {}
    for path in paths:
        with open(path) as handle:
            for process, tools in yaml.load(handle, Loader=yaml.BaseLoader).items():
                process = process.split(":")[-1]
                for tool, version in tools.items():
                    versions.setdefault(process, {}).setdefault(tool, set()).add(version)
    versions["Workflow"] = {"Nextflow": {nextflow_version}, "riboseq-flow": {pipeline_version}}
    return {process: {tool: sorted(values)[0] if len(values) == 1 else sorted(values)
                      for tool, values in sorted(tools.items())}
            for process, tools in sorted(versions.items())}


if __name__ == "__main__":
    versions = collect_versions(sys.argv[3:], sys.argv[1], sys.argv[2])
    Path("software_versions.yml").write_text(yaml.safe_dump(versions, sort_keys=True))
