# -*- coding: utf_8 -*-

from collections import OrderedDict

import re
import sqlite3

import cc_proteins.utils as utils

LOG = utils.get_logger()

REGEX_ENSEMBL_MOUSE_ID = re.compile("ENSMUS([EGTP])[0-9]{11}", re.IGNORECASE)
REGEX_ENSEMBL_HUMAN_ID = re.compile("ENS([EGTP])[0-9]{11}", re.IGNORECASE)
REGEX_MGI_ID = re.compile("MGI:[0-9]{1,}", re.IGNORECASE)
REGEX_REGION = re.compile("(CHR|)*\s*([0-9]{1,2}|X|Y|MT)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)


class Region:
    """Encapsulates a genomic region.

    Attributes:
        chromosome (str): chromosome name
        start_position (int): start position
        end_position (int): end position
    """
    def __init__(self):
        """Initialization."""
        self.chromosome = ''
        self.start_position = None
        self.end_position = None

    def __str__(self):
        """Return string representing this region.

        Returns:
            str: in format of chromosome:start_position-end_position
        """
        return '{}:{}-{}'.format(self.chromosome,
                                 self.start_position,
                                 self.end_position)

    def __repr__(self):
        """Internal representation.

        Returns:
            dict: keys being the attributes
        """
        return {'chromosome': self.chromosome,
                'start_position': self.start_position,
                'end_position': self.end_position}


def dictify_row(cursor, row):
    """Turns the given row into a dictionary where the keys are the column names

    Args:
        cursor: database cursor
        row: database row

    Returns: ``OrderedDict``
    """
    d = OrderedDict()
    for i, col in enumerate(cursor.description):
        d[col[0]] = row[i]
    return d


def connect_to_database(version=None):
    """Connect to the ensimpl database

    Args:
        version (int): the Ensembl version number

    Returns:
        a connection to the database
    """
    try:
        database = config.get_ensimpl_db(version)
        return sqlite3.connect(database)
    except Exception as e:
        LOG.error('Error connecting to database: {}'.format(str(e)))
        raise e


def nvl(value, default):
    """Returns value if value has a value, else default.

    Args:
        value: the evalue to evaluate
        default: the default value

    Returns:
        value or default
    """
    return value if value else default


def nvli(value, default):
    """Returns value as an int if value can be converted, else default.

    Args:
        value: the value to evaluate and convert
        default: the default value

    Returns:
        value or default
    """
    ret = default
    if value:
        try:
            ret = int(value)
        except ValueError:
            pass
    return ret


def get_multiplier(factor):
    """Get multiplying factor.

    The factor value should be 'mb', 'm', or 'k' and the correct multiplier
    will be returned.

    Args:
        factor (str): 'mb', 'm', or 'k'

    Returns:
        int: the multiplying value
    """
    if factor:
        factor = factor.lower()

        if factor == 'mb':
            return 10000000
        elif factor == 'm':
            return 1000000
        elif factor == 'k':
            return 1000

    return 1


def str_to_region(location):
    """Parse a string into a genomic location.

    Args:
        location (str): genomic location (range)

    Returns:
        ``ensimpl.search.search_database.Region``: region object

    Raises:
        ValueError: if location is invalid
    """
    if not location:
        raise ValueError('No location specified')

    valid_location = location.strip()

    if len(valid_location) <= 0:
        raise ValueError('Empty location')

    match = REGEX_REGION.match(valid_location)

    if not match:
        raise ValueError('Invalid location string')

    loc = Region()
    loc.chromosome = match.group(2)
    loc.start_position = match.group(4)
    loc.end_position = match.group(7)
    multiplier_one = match.group(5)
    multiplier_two = match.group(8)

    loc.start_position = int(loc.start_position)
    loc.end_position = int(loc.end_position)

    if multiplier_one:
        loc.start_position *= get_multiplier(multiplier_one)

    if multiplier_two:
        loc.end_position *= get_multiplier(multiplier_two)

    return loc


def validate_ensembl_id(ensembl_id):
    """Validate an id to make sure it conforms to the convention.

    Args:
        ensembl_id (str): the ensembl identifer

    Returns:
        the valid ensembl id

    Raises:
        ValueError: if the Ensembl ID is invalid
    """
    if not ensembl_id:
        raise ValueError('No Ensembl ID')

    valid_id = ensembl_id.strip()

    if len(valid_id) <= 0:
        raise ValueError('Empty Ensembl ID')

    if REGEX_ENSEMBL_HUMAN_ID.match(ensembl_id):
        return valid_id
    elif REGEX_ENSEMBL_MOUSE_ID.match(ensembl_id):
        return valid_id

    raise ValueError('Invalid Ensembl ID')


