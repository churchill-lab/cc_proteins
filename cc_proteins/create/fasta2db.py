# -*- coding: utf-8 -*-

import codecs
import gzip
import json
import os
import sqlite3
import time

from Bio import SeqIO

import cc_proteins.utils as utils
import cc_proteins.create.lookup_db as lookup_db
import cc_proteins.create.ensembl_db as ensembl_db


LOG = utils.get_logger()

DEFAULT_CONFIG = 'ftp://ftp.jax.org/churchill-lab/ccproteins/config.json'

FIELDS = ['strain_id',
          'protein_id',
          'protein_chrom',
          'protein_start',
          'protein_end',
          'protein_strand',
          'transcript_id',
          'transcript_name',
          'transcript_start',
          'transcript_end',
          'gene_id',
          'gene_name',
          'gene_start',
          'gene_end',
          'mgi_id',
          'ref_protein_id',
          'ref_transcript_id',
          'ref_gene_id',
          'seq']


def initialize(db):
    """Initialize the Fasta database

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

    LOG.info('Inserting values...')
    for sql in SQL_INSERT_DATA:
        LOG.debug(sql)
        cursor.execute(sql)

    cursor.close()
    conn.commit()
    conn.close()

    LOG.info('Database initialized in: {}'.format(
        utils.format_time(start, time.time())))


def rec_to_list(rec):
    """Convert a dict to list.

    Args:
        rec (dict): what to convert

    Returns:
        list: list of variables

    """
    line = []
    for f in FIELDS:
        line.append(rec.get(f, ''))
    return line


def insert_fasta(database, strain_id, fasta_file, db_lookup):
    """Insert a Fasta record into the database.

    Args:
        database (str): Fasta database
        strain_id (str): strain identifier
        fasta_file (str): the Fasta file to parse
        db_lookup (str): lookup database
    """
    LOG.debug('Parsing strain: {}, Fasta: {}'.format(strain_id, fasta_file))

    start = time.time()

    good_ids = 0
    bad_ids = 0
    fixed_ids = 0
    good_gtp = 0
    bad_gtp = 0
    records = []
    fields = ['?'] * len(FIELDS)
    fields = ','.join(fields)

    try:
        with gzip.open(fasta_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                elems = record.description.split()
                protein_id = elems[0].split('.')[0]
                protein_location = elems[2]

                loc = protein_location.split(':')
                protein_chrom = loc[2]
                protein_start = int(loc[3])
                protein_end = int(loc[4])
                protein_strand = int(loc[5])

                gene_id = elems[3].split(':')[1]
                transcript_id = elems[4].split(':')[1]

                gene_id = gene_id.split('.')[0]
                transcript_id = transcript_id.split('.')[0]

                mgi_id = None
                mgi_idx = elems[-1].find('MGI:')
                if mgi_idx > -1:
                    # 'Symbol;Acc:MGI:2145458]']
                    mgi_id = elems[-1][mgi_idx:-1]

                rec = {'protein_id': protein_id,
                       'protein_chrom': protein_chrom,
                       'protein_start': protein_start,
                       'protein_end': protein_end,
                       'protein_strand': protein_strand,
                       'transcript_id': transcript_id,
                       'transcript_start': None,
                       'transcript_end': None,
                       'transcript_name': None,
                       'gene_id': gene_id,
                       'gene_start': None,
                       'gene_end': None,
                       'gene_name': None,
                       'ref_gene_id': None,
                       'ref_transcript_id': None,
                       'ref_protein_id': None,
                       'strain_id': strain_id,
                       'mgi_id': mgi_id,
                       'seq': str(record.seq)}

                if strain_id == 'mus_musculus':
                    # reference strain
                    rec['ref_gene_id'] = gene_id
                    rec['ref_transcript_id'] = transcript_id
                    rec['ref_protein_id'] = protein_id
                else:

                    #print('------------------------------------------------------')
                    #print(rec)

                    ids = lookup_db.lookup_id(db_lookup, strain_id, protein_id)

                    #print('lookup_db.lookup_id(' ,db_lookup, strain_id, protein_id, ',')
                    #print('lookup_db.lookup_gtp(', db_lookup, strain_id, protein_id, ',')

                    if ids:
                        good_ids += 1
                        for i in ids:
                            if i['ensembl_id'][:7] == 'ENSMUSG':
                                rec['ref_gene_id'] = i['ensembl_id']
                            elif i['ensembl_id'][:7] == 'ENSMUSG':
                                rec['ref_transcript_id'] = i['ensembl_id']
                            elif i['ensembl_id'][:7] == 'ENSMUSP':
                                rec['ref_protein_id'] = i['ensembl_id']

                        if rec['ref_gene_id'] is None or len(rec['ref_gene_id']) == 0:
                            rec['ref_gene_id'] = lookup_db.lookup_gene_id(db_lookup,
                                                                          protein_id,
                                                                          mgi_id)

                            if rec['ref_gene_id']:
                                #print('FIXED: ', rec['ref_gene_id'])
                                #print(record)
                                fixed_ids += 1
                            else:
                                #print('NOT FIXED: ', str(rec))
                                #print(record)
                                LOG.debug("Cannot fix the following record: {}".format(protein_id))



                    else:
                        #LOG.debug('no_ids')
                        bad_ids += 1
                        #print('lookup_db.lookup_gene_id(', db_lookup, protein_id, mgi_id, ')')
                        rec['ref_gene_id'] = lookup_db.lookup_gene_id(db_lookup, protein_id, mgi_id)

                        if rec['ref_gene_id']:
                            #print('FIXED: ', rec['ref_gene_id'])
                            #print(record)
                            fixed_ids += 1
                        else:
                            #print('NOT FIXED: ', str(rec))
                            #print(record)

                            LOG.debug("Cannot fix the following record: {}".format(protein_id))

                gtp = lookup_db.lookup_gtp(db_lookup, strain_id, protein_id)
                if gtp:
                    good_gtp += 1
                    if len(gtp) > 1:
                        LOG.error('Too many records for {}'.format(protein_id))
                        LOG.error((len(gtp)))
                        break

                    rec['transcript_id'] = gtp[0]['transcript_id']
                    rec['transcript_start'] = gtp[0]['transcript_start']
                    rec['transcript_end'] = gtp[0]['transcript_end']
                    rec['transcript_name'] = gtp[0]['transcript_name']
                    rec['gene_id'] = gtp[0]['gene_id']
                    rec['gene_start'] = gtp[0]['gene_start']
                    rec['gene_end'] = gtp[0]['gene_end']
                    rec['gene_name'] = gtp[0]['gene_name']

                    if not rec['gene_name']:
                        rec['gene_name'] = lookup_db.lookup_gene_name(db_lookup, rec['ref_gene_id'])

                else:
                    bad_gtp += 1

                #if rec['ref_gene_id'] is None or len(rec['ref_gene_id']) == 0:
                #    print('----------------------------------')
                #    print('record=', record)
                #    print('UGH:', rec)

                records.append(rec_to_list(rec))

        conn = sqlite3.connect(database)
        sql_fasta_insert = ('INSERT INTO fasta_record '
                            'VALUES (null, {})'.format(fields))

        cursor = conn.cursor()
        LOG.debug('Inserting {:,} records'.format(len(records)))
        cursor.executemany(sql_fasta_insert, records)
        cursor.close()
        conn.commit()
    except Exception as e:
        LOG.error(e)

    LOG.debug('Good Identifiers: {}'.format(good_ids))
    LOG.debug('Bad Identifiers: {}'.format(bad_ids))
    LOG.debug('Fixed Identifiers: {}'.format(fixed_ids))
    LOG.debug('Good GTP: {}'.format(good_gtp))
    LOG.debug('Bad GTP: {}'.format(bad_gtp))

    LOG.info('Transcripts, proteins inserted in: {}'.format(
        utils.format_time(start, time.time())))


def finalize(Database):
    """Finalize the database.  Move everything to where it needs to be and
    create the necessary indices.

     Args:
        Database (str): the name of the database file
     """
    start = time.time()
    conn = sqlite3.connect(Database)
    cursor = conn.cursor()

    LOG.info("Finalizing database....")

    for sql in SQL_INDICES:
        LOG.debug(sql)
        cursor.execute(sql)

    conn.commit()
    conn.close()

    LOG.info("Finalizing complete: {0}".format(
        utils.format_time(start, time.time())))


def generate_database(database, resource, version):
    """Create fasta database.

    Args:
        database (str): database to create
        resource (str): configuration to parse for Ensembl information
        version (int): Ensembl version to use
    """
    database = os.path.abspath(database)

    LOG.debug('Database: {}'.format(database))
    LOG.debug('Resource: {}'.format(resource))
    LOG.debug('Ensembl Version: {}'.format(version))

    config = {}

    try:
        reader = codecs.getreader("utf-8")
        with utils.open_resource(resource) as fd:
            config = json.load(reader(fd))
    except IOError as io_error:
        LOG.error('Unable to access resource: {}'.format(resource))
        LOG.debug(io_error)
    except Exception as exc:
        LOG.error('Unable to parse resource: {}'.format(resource))
        LOG.debug(exc)
        return

    max_version = max([int(v) for (v, val) in config['configurations'].items()])

    if not version:
        version = str(max_version)

    if not config:
        LOG.error('Unable to determine the Ensembl releases and locations '
                  'for download')
        LOG.error('Please make sure that the resource "{}" is accessible '
                  'and in the correct format'.format(resource))
        raise Exception("Unable to create lookup databases")

    if version not in config['configurations']:
        LOG.error('Ensembl version {} is unknown'.format(version))
        LOG.error('Please make sure that the resource "{}" is accessible '
                  'and in the correct format'.format(resource))
        raise Exception("Unable to create lookup databases")

    ensembl_version = config['configurations'][version]

    dir_name = os.path.dirname(database)
    tmp_lookup_db = 'lookup.{}.db3'.format(version)
    tmp_lookup_db = os.path.join(dir_name, tmp_lookup_db)
    #utils.delete_file(tmp_lookup_db)

    LOG.info('Creating lookup database {} ...'.format(tmp_lookup_db))

    try:
        '''
        lookup_db.initialize(tmp_lookup_db)

        gtp = ensembl_db.extract_ensembl_gtp(ensembl_version,
                                             ensembl_version['reference'])
        lookup_db.insert_gtp(tmp_lookup_db, gtp)

        for strain in ensembl_version['strains']:
            compara_ids = ensembl_db.extract_compara_ids(ensembl_version, strain)
            lookup_db.insert_ids(tmp_lookup_db, compara_ids)
            gtp = ensembl_db.extract_ensembl_gtp(ensembl_version, strain)
            lookup_db.insert_gtp(tmp_lookup_db, gtp)

        LOG.info('Finalizing lookup db...')

        lookup_db.finalize(tmp_lookup_db)
        '''

        LOG.info('Creating Fasta db ...')
        utils.delete_file(database)

        initialize(database)

        for strain in ensembl_version['strains']:
            fasta_file = utils.download_file(strain['fasta_file'],
                                             directory=dir_name)

            insert_fasta(database, strain['strain_id'], fasta_file, tmp_lookup_db)


        fasta_file = utils.download_file(ensembl_version['reference']['fasta_file'],
                                         directory=dir_name)

        insert_fasta(database, ensembl_version['reference']['strain_id'], fasta_file, tmp_lookup_db)

        finalize(database)

        #utils.delete_file(tmp_lookup_db)
    except sqlite3.Error as e:
        LOG.error('Database Error: {}'.format(e))


SQL_CREATE_TABLES = ['''
    CREATE TABLE IF NOT EXISTS fasta_record (
        fasta_record_id INTEGER PRIMARY KEY,
        strain_id TEXT,
        protein_id TEXT,
        protein_chrom TEXT,
        protein_start INTEGER,
        protein_end INTEGER,
        protein_strand INTEGER,
        transcript_id TEXT,
        transcript_name TEXT,
        transcript_start TEXT,
        transcript_end TEXT,
        gene_id TEXT,
        gene_name TEXT,
        gene_start TEXT,
        gene_end TEXT,
        mgi_id TEXT,
        ref_protein_id TEXT,
        ref_transcript_id TEXT,
        ref_gene_id TEXT,
        seq TEXT
    );
''', '''
    CREATE TABLE IF NOT EXISTS strains
    (
        strain_num INTEGER PRIMARY KEY,
        strain_id TEXT,
        strain_key TEXT,
        strain_desc TEXT,
        strain_order INTEGER
    );
''']

SQL_INSERT_DATA = [
    'INSERT INTO strains VALUES (null, "mus_musculus_aj", "A", "A/J", 1);',
    'INSERT INTO strains VALUES (null, "mus_musculus", "B", "C57BL/6J", 2);',
    'INSERT INTO strains VALUES (null, "mus_musculus_129s1svimj", "C", "129S1/SvImJ", 3);',
    'INSERT INTO strains VALUES (null, "mus_musculus_nodshiltj", "D", "NOD/ShiLtJ", 4);',
    'INSERT INTO strains VALUES (null, "mus_musculus_nzohlltj", "E", "NZO/H1LtJ", 5);',
    'INSERT INTO strains VALUES (null, "mus_musculus_casteij", "F", "CAST/EiJ", 6);',
    'INSERT INTO strains VALUES (null, "mus_musculus_pwkphj", "G", "PWK/PhJ", 7);',
    'INSERT INTO strains VALUES (null, "mus_musculus_wsbeij", "H", "WSB/EiJ", 8);'
]


SQL_INDICES = [
    'CREATE INDEX IF NOT EXISTS idx_fasta_record_strain_idx ON fasta_record (strain_id ASC);',
    'CREATE INDEX IF NOT EXISTS idx_strains_strain_idx ON strains (strain_id ASC);'
]





