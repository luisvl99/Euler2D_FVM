"""Read the files written by euler_cli and by the GUI export."""

import json
import re
from pathlib import Path

import pandas as pd

_KEY_VALUE = re.compile(r"(\w+)=(\S+)")


def _number(text):
    for convert in (int, float):
        try:
            return convert(text)
        except ValueError:
            pass
    return text


def read_csv(path):
    """Return (DataFrame, metadata) for one exported CSV file.

    The metadata dict holds every key=value pair of the '#' header, for
    example case, Nx, Lx, CFL, Scheme, Flux, Reconstruction, iter, sim_time.
    """
    meta = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                break
            meta.update((key, _number(value)) for key, value in _KEY_VALUE.findall(line))
    return pd.read_csv(path, comment="#"), meta


def read_run(folder):
    """Return the run.json summary of a euler_cli output folder."""
    with open(Path(folder) / "run.json", encoding="utf-8") as f:
        return json.load(f)
