# -*- coding: utf-8 -*-
"""Useful generic utilities for the package.
"""
from collections import OrderedDict
from functools import cmp_to_key
from operator import itemgetter as ig

import bz2
import gzip
import logging
import os
import random
import string

from urllib.request import urlopen

logging.basicConfig(format='[cc_proteins] [%(asctime)s] %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p')


def format_time(start, end):
    """
    Format length of time between start and end.

    Args:
        start (long): the start time
        end (long): the end time

    Returns:
        str: a formatted string of hours, minutes, and seconds
    """
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)


class ReverseProxied(object):
    """
    Wrap the application in this middleware and configure the
    front-end server (Apaccto add these headers, to let you quietly bind
    this to a URL other than / and to an HTTP scheme that is 
    different than what is used locally.

    http://flask.pocoo.org/snippets/35/

    In apache:

    <Proxy http://server:port>
        Header add X_FORWADED_PATH "/myprefix"
        RequestHeader set X_FORWARDED_PATH "/myprefix"
    </Proxy>

    ProxyPass /myprefix http://server:port disablereuse=On
    ProxyPassReverse /myprefix http://server:port

    In nginx:

    location /myprefix {
        proxy_pass http://192.168.0.1:5001;
        proxy_set_header Host $host;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Scheme $scheme;
        proxy_set_header X-Script-Name /myprefix;
        }

    """
    def __init__(self, app):
        """
        Constructor.

        Args:
            app: the application
        """
        self.app = app

    def __call__(self, environ, start_response):
        """
        Middleware call.

        Args:
            environ (dict): application environment
            start_response (Response.start_response): response method

        Returns:
            application
        """
        script_name = environ.get('HTTP_X_FORWARDED_PATH', '')

        if script_name:
            environ['SCRIPT_NAME'] = script_name
            path_info = environ['PATH_INFO']
            if path_info.startswith(script_name):
                environ['PATH_INFO'] = path_info[len(script_name):]

        scheme = environ.get('HTTP_X_SCHEME', '')

        if scheme:
            environ['wsgi.url_scheme'] = scheme

        return self.app(environ, start_response)


def get_logger():
    """Get the logger.

    Returns:
        logger: the logging object
    """
    return logging.getLogger(__name__)


def configure_logging(level=0):
    """Configure the logger with the specified level.

        0 = logging.WARNING is informational
        1 = logging.INFO is user debug
        2> = logging.DEBUG is developer debug

    Args:
        level (int): the logging level; defaults to 0
    """
    if level == 0:
        get_logger().setLevel(logging.WARN)
    elif level == 1:
        get_logger().setLevel(logging.INFO)
    elif level > 1:
        get_logger().setLevel(logging.DEBUG)


def dictify_row(cursor, row):
    """Turns the given row into a dictionary where the keys are the column names.

    Args:
        cursor (sqlite3.Cursor): the database cursor
        row (dict): the current row

    Returns:
        OrderedDict: dictionary with keys being column names
    """
    d = OrderedDict()
    for i, col in enumerate(cursor.description):
        d[col[0]] = row[i]
    return d


def dictify_cursor(cursor):
    """
    Converts all cursor rows into dictionaries where keys are the column names.

    Args:
        cursor (sqlite3.Cursor): the database cursor

    Returns:
        list: list of ``OrderedDict`` objects with keys being column names
    """
    return [dictify_row(cursor, row) for row in cursor]


def format_time(start, end):
    """Format length of time between start and end.

    Args:
        start (float): the start time
        end (float): the end time

    Returns:
        str: a formatted string of hours, minutes, and seconds
    """
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    return "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)


def str2bool(val):
    """Convert a string into a boolean.  Valid strings that return ``True``
    are: ``true``, ``1``, ``t``, ``y``, ``yes``

    Canse-sensitivity does NOT matter.  ``yes`` is the same as ``YeS``.

    Args:
        val (str): A string representing the boolean value of ``True``.

    Returns:
        bool: ``True`` if `val` represents a boolean ``True``.
    """
    return str(val).lower() in ['true', '1', 't', 'y', 'yes']


def open_resource(resource, mode='rb'):
    """Open different types of files and return the handle.

    Args:
        resource (str): a string representing a file or url
        mode (str): standard file open modes

    Returns:
        a handle to the opened ``resource``
    """
    if not resource:
        return None

    if not isinstance(resource, str):
        return resource

    if resource.endswith(('.gz', '.Z', '.z')):
        return gzip.open(resource, mode)
    elif resource.endswith(('.bz', '.bz2', '.bzip2')):
        return bz2.BZ2File(resource, mode)
    elif resource.startswith(('http://', 'https://', 'ftp://')):
        return urlopen(resource)
    else:
        return open(resource, mode)


def get_file_name(url, directory=None):
    """Get the last part of the file name.

    Args:
        url (str): url of file
        directory: (optional) directory name

    Returns:
        str: full path of file
    """
    download_file_name = url.split('/')[-1]
    local_directory = directory if directory else os.getcwd()
    return os.path.abspath(os.path.join(local_directory, download_file_name))


def is_url(url):
    """Check if this is a URL or not by just checking the protocol

    Args:
        url (str): a string in the form of a url

    Returns:
        bool: True if it has a valid protocol
    """
    if url and url.startswith(('http://', 'https://', 'ftp://')):
        return True

    return False


def delete_file(file_name):
    """ Delete specified ``file_name``.

    Args:
        file_name (str): name of file to delete
    """
    try:
        os.remove(file_name)
    except:
        pass


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy.

    Args:
        x (dict): a dictionary
        y (dict): a dictionary

    Returns:
        dict: merged dictionary
    """
    z = x.copy()
    z.update(y)
    return z


def create_random_string(size=6, chars=string.ascii_uppercase + string.digits):
    """Generate a random string.

    Args:
        size (int): defaults to 6 characters
        chars (list): ``list`` of characters to use

    Returns:
        str: a random generated string
    """
    return ''.join(random.choice(chars) for _ in range(size))


def is_sorted(l):
    """Check to see if l is ascending order.

    Args:
        l (list): ``list`` of values

    Returns:
        bool: True if sorted, False otherwise
    """
    return all(l[i] <= l[i+1] for i in iter(range(len(l)-1)))


def bin_search(a, from_index, to_index, key):
    """
    Searches a range of the specified list (a) for the specified value (key)
    using the binary search algorithm. The range must be sorted prior to making
    this call. If it is not sorted, the results are undefined.  If the range
    contains multiple elements with the specified value, there is no guarantee
    which one will be found.

    Take from Java Arrays.binarySearch source code and adapted.

    Args:
        a (list): ``list`` of values
        from_index (int): index of first element
        to_index (int): index of last element
        key: value to find

    Returns:
        int: index of the search key, if it is contained in the ``list`` within
        the specified range; otherwise, (-(insertion point) - 1). The insertion
        point is defined as the point at which the key would be inserted into
        the array: the index of the first element in the range greater than the
        key, or to_index if all elements in the range are less than the
        specified key. Note that this guarantees that the return value will be
        >= 0 if and only if the key is found.
    """
    low = from_index
    high = to_index - 1

    if key > a[high]:
        return -(to_index + 1)

    while low <= high:
        mid = (low + high)
        mid_val = a[mid]

        if mid_val < key:
            low = mid + 1
        elif mid_val > key:
            high = mid - 1
        else:
            # key found
            return mid

    return -(low + 1)


def binary_search(a, key):
    """
    Searches a range of the specified list (a) for the specified value (key)
    using the binary search algorithm. The range must be sorted prior to making
    this call. If it is not sorted, the results are undefined.  If the range
    contains multiple elements with the specified value, there is no guarantee
    which one will be found.

    Take from Java Arrays.binarySearch source code and adapted.

    Args:
        a (list): ``list`` of values
        key: value to find

    Returns:
        int: index of the search key, if it is contained in the ``list``;
        otherwise, (-(insertion point) - 1). The insertion point is defined as
        the point at which the key would be inserted into the array: the index
        of the first element in the range greater than the key, or to_index if
        all elements in the range are less than the specified key. Note that
        this guarantees that the return value will be >= 0 if and only if the
        key is found.
    """
    return bin_search(a, 0, len(a), key)


def multikeysort(items, columns):
    """Sort a ``list`` of ``dicts`` by multiple keys in ascending or descending
    order. To sort in descending order, prepend a '-' (minus sign) on the
    column name.

    Pulled from: https://stackoverflow.com/questions/1143671/python-sorting-list-of-dictionaries-by-multiple-keys

    Examples:
        >>> my_list = [
            {'name': 'apple', 'count': 10, 'price': 1.00},
            {'name': 'banana', 'count': 5, 'price': 1.00},
            {'name': 'orange', 'count': 20, 'price': 2.00},
        ]

        >>> multikeysort(my_list, ['-name', 'count'])
        [{'count': 20, 'name': 'orange', 'price': 2.0},
         {'count': 5, 'name': 'banana', 'price': 1.0},
         {'count': 10, 'name': 'apple', 'price': 1.0}]

        >>> multikeysort(my_list, ['-price', 'count'])
        [{'count': 20, 'name': 'orange', 'price': 2.0},
         {'count': 5, 'name': 'banana', 'price': 1.0},
         {'count': 10, 'name': 'apple', 'price': 1.0}]

    Args:
        items (list): The ``list`` of ``dict`` objects.
        columns (list): A ``list`` of columns names to sort `items`.

    Returns:
        list: The sorted ``list``.

    """
    comparers = [
        ((ig(col[1:].strip()), -1) if col.startswith('-') else (
        ig(col.strip()), 1))
        for col in columns
    ]

    def comparer(left, right):
        comparer_iter = (
            cmp(fn(left), fn(right)) * mult
            for fn, mult in comparers
        )
        return next((result for result in comparer_iter if result), 0)

    return sorted(items, key=cmp_to_key(comparer))


def download_file(url, file_name=None, directory=None):
    """Download ``url`` and save to ``file_name`` if specified.

     Args:
         url (str): url to download
         file_name (str): the file name to save to, defaults to ``url`` filename
         directory (str): the directory tp save to

     Returns:
         str: full path to downloaded file
     """

    CHUNK_SIZE = 1024

    download_file_name = url.split('/')[-1]

    local_file_name = file_name if file_name else download_file_name

    local_directory = directory if directory else os.getcwd()

    local_path = os.path.abspath(os.path.join(local_directory, local_file_name))
    '''

    with open(local_path, 'wb') as f:
        r = requests.get(url, stream=True)
        total_length = r.headers.get('content-length')
        if total_length is None:  # no content length header
            f.write(r.content)
        else:
            total_length = int(total_length)
            with click.progressbar(length=total_length,
                                   label=download_file_name) as bar:
                for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                    f.write(chunk)
                    f.flush()
                    bar.update(len(chunk))
    '''
    print(local_path)
    return local_path

