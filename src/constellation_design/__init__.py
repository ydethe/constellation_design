# -*- coding: utf-8 -*-
"""

.. include:: ../../README.md

# Testing

## Run the tests

To run tests, just run:

    pdm run pytest

## Test reports

[See test report](../tests/report.html)

[See test results](../tests/results/fig_comparison.html)

[See coverage](../coverage/index.html)

.. include:: ../../CHANGELOG.md

"""
import os
import logging

from rich.logging import RichHandler
from setuptools_scm import get_version  # type: ignore


# création de l'objet logger qui va nous servir à écrire dans les logs
logger = logging.getLogger("constellation_design_logger")
logger.setLevel(os.environ.get("LOGLEVEL", "INFO").upper())

stream_handler = RichHandler()
logger.addHandler(stream_handler)
