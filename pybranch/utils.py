# pybranch/utils.py
"""A collection of utility functions for the PyBranch project.

This module provides helper functions that are used across different parts of
the analysis package, such as loading configuration files.
"""

import yaml
from typing import Dict, Any

def load_config(config_path: str = 'config.yaml') -> Dict[str, Any]:
    """Loads and parses the YAML configuration file.

    This function opens the specified YAML file, safely loads its contents,
    and returns them as a Python dictionary.

    Args:
        config_path (str): The path to the YAML configuration file.
            Defaults to 'config.yaml'.

    Returns:
        Dict[str, Any]: A dictionary containing the configuration settings.
        Returns an empty dictionary if the file is empty.

    Raises:
        FileNotFoundError: If the configuration file cannot be found at the
            specified path.
        yaml.YAMLError: If the configuration file contains invalid YAML syntax
            and cannot be parsed.
    """
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        if config is None:
            return {}
        return config
    except FileNotFoundError:
        print(f"Error: Configuration file not found at '{config_path}'")
        raise
    except yaml.YAMLError as e:
        print(f"Error: Could not parse YAML configuration file: {e}")
        raise