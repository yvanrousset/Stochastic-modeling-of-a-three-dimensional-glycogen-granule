from datetime import datetime
from . import core
from . import utils
from . import algorithm

import logging

# set up logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.WARNING)
# log to file
fh = logging.FileHandler(
    "{}_{:%Y-%m-%d}.log".format(__name__, datetime.now()))
fh.setLevel(logging.DEBUG)
# log to terminal
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add handlers
#logger.addHandler(fh)
logger.addHandler(ch)
