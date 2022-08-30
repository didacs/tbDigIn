"""Main entry point for all PYTOMEBIO tools."""

import logging
import sys
from typing import Callable
from typing import Dict
from typing import List
from typing import Optional
from typing import Type

import defopt

from pytomebio.tools.delim_filter import delim_filter
from pytomebio.tools.delim_filter import ColumnFilter

TOOLS: List[Callable] = sorted(
    [delim_filter],
    key=lambda f: f.__name__,
)


PARSERS: Dict[Type, Callable] = {ColumnFilter: ColumnFilter.build}


def main(argv: Optional[List[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    logger = logging.getLogger(__name__)
    if len(argv) != 0 and all(arg not in argv for arg in ["-h", "--help"]):
        logger.info("Running command: tomebio-tools " + " ".join(argv))
    try:
        defopt.run(funcs=TOOLS, argv=argv, parsers=PARSERS)
        logger.info("Completed successfully.")
    except Exception as e:
        logger.info("Failed on command: " + " ".join(argv))
        raise e
