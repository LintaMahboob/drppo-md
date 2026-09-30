"""Shared helpers for the DrPPO MD pipeline."""
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = REPO / "config" / "proteins.yaml"


def load_config(path=None):
    """Read the YAML configuration file."""
    path = Path(path) if path else DEFAULT_CONFIG
    with open(path) as fh:
        return yaml.safe_load(fh)


def get_protein(cfg, name):
    """Return the settings of one protein or exit with a helpful message."""
    try:
        return cfg["proteins"][name]
    except KeyError:
        known = ", ".join(cfg["proteins"])
        sys.exit(f"Unknown protein '{name}'. Defined in config: {known}")
