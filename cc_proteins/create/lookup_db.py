# -*- coding: utf-8 -*-
"""This module is specific to ensimpl db operations.

Todo:
    * better documentation
"""


import sqlite3
import time

from collections import OrderedDict

import cc_proteins.utils as utils

LOG = utils.get_logger()


def initialize(db):
    """Initialize the lookup database.

    Args:
        db (str): the name of the database file
    """
    LOG.info('Initializing database: {}'.format(db))

    start = time.time()
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    LOG.info('Generating tables...')
    for sql in SQL_CREATE_TABLES:
        LOG.debug(sql)
        cursor.execute(sql)

    cursor.close()
    conn.commit()
    conn.close()

    LOG.info('Database initialized in: {}'.format(
        utils.format_time(start, time.time())))


def insert_ids(db, ids, chunk_size=50000):
    """Insert ids into the database.

    Args:
        db (str): the name of the database file
        ids (list): all the ids
        chunk_size (int): how many ids to insert at a time
    """
    LOG.info('Inserting IDs into database: {}'.format(db))

    start = time.time()
    conn = sqlite3.connect(db)

    sql_ids_insert = 'INSERT INTO id_lookup VALUES (?, ?, ?, ?, ?)'

    # chunk
    prev_x = 0
    for x in range(chunk_size, len(ids), chunk_size):
        cursor = conn.cursor()
        LOG.debug('Inserting ids {:,} -> {:,}'.format(prev_x, x - 1))
        cursor.executemany(sql_ids_insert, ids[prev_x:x])
        cursor.close()
        conn.commit()
        prev_x = x

    cursor = conn.cursor()
    LOG.debug('Inserting ids {:,} -> {:,}'.format(prev_x, len(ids)))
    cursor.executemany(sql_ids_insert, ids[prev_x:])
    cursor.close()
    conn.commit()

    LOG.info('IDs inserted in: {}'.format(
        utils.format_time(start, time.time())))


def insert_gtp(db, gtp, chunk_size=50000):
    """Insert the gene, transcript, protein information into the database.

    Args:
        db (str): the name of the database file
        gtp (list): a list of the gene, transcript, protein information
        chunk_size (int): how many ids to insert at a time
    """
    LOG.info('Inserting transcripts, proteins into database: {}'.format(db))
    start = time.time()
    conn = sqlite3.connect(db)

    sql_gtpe_insert = ('INSERT INTO ensembl_gtp '
                       'VALUES (?, ?, ?, ?, ?, ?, ?, ?, '
                       '?, ?, ?, ?, ?, ?, ?, ?, ?)')

    LOG.info('Inserting gtp data into database: {}'.format(db))

    # chunk
    prev_x = 0
    for x in range(chunk_size, len(gtp), chunk_size):
        cursor = conn.cursor()
        LOG.debug('Inserting rows {:,} -> {:,}'.format(prev_x, x - 1))
        cursor.executemany(sql_gtpe_insert, gtp[prev_x:x])
        cursor.close()
        conn.commit()
        prev_x = x

    cursor = conn.cursor()
    LOG.debug('Inserting rows {:,} -> {:,}'.format(prev_x, len(gtp)))
    cursor.executemany(sql_gtpe_insert, gtp[prev_x:])
    cursor.close()
    conn.commit()

    LOG.info('GTP data inserted in: {}'.format(
        utils.format_time(start, time.time())))


def finalize(db):
    """Finalize the database.  Move everything to where it needs to be and
    create the necessary indices.

     Args:
        db (str): the name of the database file
     """
    start = time.time()
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    LOG.info("Finalizing database....")

    for sql in SQL_INDICES:
        LOG.debug(sql)
        cursor.execute(sql)

    conn.commit()
    conn.close()

    LOG.info("Finalizing complete: {0}".format(
        utils.format_time(start, time.time())))


def dictify_row(cursor, row):
    """Turns the given row into a dictionary where the keys are the column
    names.

    Args:
        cursor: database cursor
        row: database row

    Returns: ``OrderedDict``
    """
    d = OrderedDict()
    for i, col in enumerate(cursor.description):
        d[col[0]] = row[i]
    return d


def lookup_gene_id_mgi(db, mgi_id):
    """Search for th Ensembl gene id by the MGI id.

    Args:
        db (str): the database
        mgi_id (str): the MGI identifier in form of MGI:NNNNNNN

    Returns:
        the Ensembl ID if found, None if not
    """
    sql = '''
        select distinct gene_id ensembl_id 
          from ensembl_gtp 
         where mgi_id = :mgi_id
           and ensembl_id like 'ENSMUSG%'
    '''

    sql_params = {'mgi_id': mgi_id}

    result = None

    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    for row in cursor.execute(sql, sql_params):
        result = row['ensembl_id']

    cursor.close()

    return result


def lookup_gene_id_protein(db, protein_id):
    """Attempt to find the Ensembl gene id by the Ensembl protein id

    Args:
        db (str): the database
        protein_id (str): the Ensembl protein identifier

    Returns:
        the Ensembl ID if found, None if not
    """
    sql = '''
        select distinct ensembl_id ensembl_id 
          from id_lookup 
         where ensembl_id like 'ENSMUSG%'
           and strain_ensembl_id in (
               select protein_id 
                 from ensembl_gtp 
                where gene_id in (
                    select gene_id 
                      from ensembl_gtp 
                     where protein_id = :protein_id
           ))
    '''

    sql_params = {'protein_id': protein_id}

    result = None

    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    for row in cursor.execute(sql, sql_params):
        result = row['ensembl_id']

    cursor.close()

    return result


def lookup_gene_id(db, protein_id, mgi_id):
    """Wrapper around both `lookup_gene_id_mgi` and `lookup_gene_id_protein`.

    Args:
        db (str): the database
        protein_id (str): the Ensembl protein identifier
        mgi_id (str): the MGI identifier in form of MGI:NNNNNNN

    Returns:
        the Ensembl ID if found, None if not
    """
    ensembl_id1 = lookup_gene_id_protein(db, protein_id)
    ensembl_id2 = None

    if mgi_id:
        ensembl_id2 = lookup_gene_id_mgi(db, mgi_id)

    if ensembl_id2:
        return ensembl_id2
    else:
        return ensembl_id1


def lookup_gene_name(db, ensembl_id):
    """Attempt to find the symbol or name of a gene.

    Args:
        db (str): the database
        ensembl_id (str): the Ensembl gene identifier

    Returns:
        the name of the gene if found, None if not
    """
    sql = '''
        SELECT distinct gene_name gene_name
          FROM ensembl_gtp
         WHERE gene_id = :gene_id
    '''

    sql_params = {'gene_id': ensembl_id}

    result = None

    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    for row in cursor.execute(sql, sql_params):
        result = row['gene_name']

    cursor.close()

    return result


def lookup_gtp(db, strain_id, strain_ensembl_id):
    """Get the mapping of gene, transcipts and proteins.

    Args:
        db (str): the database
        strain_id (str): the Ensembl strain identifier
        strain_ensembl_id (str): the Ensembl protein identifier

    Returns:
        list: a ``list`` of Ensembl genes, transcripts and proteins
    """
    sql = '''
        SELECT *
          FROM ensembl_gtp
         WHERE strain_id = :strain_id
           AND protein_id = :protein_id
    '''

    sql_params = {'strain_id': strain_id,
                  'protein_id': strain_ensembl_id}

    results = []

    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    for row in cursor.execute(sql, sql_params):
        results.append(dictify_row(cursor, row))

    cursor.close()

    return results


def lookup_id(db, strain_id, strain_ensembl_id):
    """Get the identifier information for a given protein id.

    Args:
        db (str): the database
        strain_id (str): the Ensembl strain identifier
        strain_ensembl_id (str): the Ensembl protein identifier

    Returns:
        list: a ``list`` of Ensembl identifiers
    """
    sql = '''
        SELECT *
          FROM id_lookup
         WHERE strain_id = :strain_id
           AND strain_ensembl_id = :protein_id
    '''

    sql_params = {'strain_id': strain_id,
                  'protein_id': strain_ensembl_id}

    results = []

    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    for row in cursor.execute(sql, sql_params):
        results.append(dictify_row(cursor, row))

    cursor.close()

    return results


SQL_CREATE_TABLES = ['''
    CREATE TABLE IF NOT EXISTS id_lookup (
       id_lookup_key INTEGER,
       ensembl_id TEXT NOT NULL,
       strain_ensembl_id TEXT NOT NULL,
       strain_id TEXT NOT NULL,
       source_name TEXT NOT NULL,
       PRIMARY KEY (id_lookup_key)
    );
''', '''
    CREATE TABLE IF NOT EXISTS ensembl_gtp (
        gtp_key INTEGER PRIMARY KEY,
        gene_id TEXT NOT NULL,
        gene_version TEXT,
        gene_name TEXT,
        mgi_id TEXT,
        gene_chrom TEXT,
        gene_start INTEGER,
        gene_end INTEGER,
        gene_strand INTEGER,
        transcript_id TEXT,
        transcript_version TEXT,
        transcript_name TEXT,
        transcript_start INTEGER,
        transcript_end INTEGER,
        protein_id TEXT,
        protein_version TEXT,
        strain_id TEXT
    );
''']

SQL_INDICES = [
    'CREATE INDEX IF NOT EXISTS idx_gtp_strain_id ON ensembl_gtp(strain_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtp_gene_id ON ensembl_gtp(gene_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtp_mgi_id ON ensembl_gtp(mgi_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtp_transcript_id ON ensembl_gtp(transcript_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtp_protein_id ON ensembl_gtp(protein_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_id_lookup_strain_id ON id_lookup (strain_id ASC);',
    'CREATE INDEX IF NOT EXISTS idx_id_lookup_ensembl ON id_lookup (ensembl_id ASC);',
    'CREATE INDEX IF NOT EXISTS idx_id_lookup_strain_ensembl_id ON id_lookup (strain_ensembl_id ASC);'
]

