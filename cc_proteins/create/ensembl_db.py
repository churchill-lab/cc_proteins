# -*- coding: utf-8 -*-
"""This module is specific to Ensembl database operations utilizing their
public MySQL server.

Todo:
    * better documentation
"""

import pymysql

import cc_proteins.utils as utils

LOG = utils.get_logger()

COMPARA_DB_PREFIX = 'ensembl_compara_'

SQL_ENSEMBL_GET_IDS_GENE = '''
SELECT smpsi.source_stable_id ensembl_id, ss.stable_id strains_ensembl_id, gd.name strain_name, ss.source_name
  FROM seq_member ss,
       genome_db gd,
       seq_member_projection_stable_id smpsi
 WHERE ss.genome_db_id = gd.genome_db_id
   AND ss.seq_member_id = smpsi.target_seq_member_id
   AND gd.name = %s
   AND smpsi.source_stable_id like 'ENSMUS%%'
 ORDER BY smpsi.source_stable_id
'''

SQL_ENSEMBL_GET_IDS_PROTEIN = '''
SELECT ss.stable_id ensembl_id, st.stable_id strains_id, gt.name strain_name, ss.source_name 
  FROM seq_member ss,
       genome_db gd,
       seq_member_projection smp,
       seq_member st,
       genome_db gt
 WHERE ss.genome_db_id = gd.genome_db_id
   AND ss.seq_member_id = smp.source_seq_member_id 
   AND smp.target_seq_member_id = st.seq_member_id
   AND st.genome_db_id = gt.genome_db_id
   AND gt.name = %s
   AND ss.stable_id like 'ENSMUS%%'
 ORDER BY ss.stable_id   
'''

SQL_ENSEMBL_STRAIN_GTP = '''
SELECT g.stable_id gene_id,
       g.version gene_version,
       (select display_label
          from xref
         where xref.xref_id = g.display_xref_id) gene_name,
       (select x1.dbprimary_acc
          from xref x1
         where g.display_xref_id = x1.xref_id
           and x1.external_db_id = 1400) mgi_id,
       s.name gene_chrom,
       g.seq_region_start gene_start,
       g.seq_region_end gene_end,
       g.seq_region_strand gene_strand,
       t.stable_id transcript_id,
       t.version transcript_version,
       (select display_label
          from xref
         where xref.xref_id = t.display_xref_id) transcript_name,
       t.seq_region_start transcript_start,
       t.seq_region_end transcript_end,
       (select stable_id 
          from translation 
         where translation.transcript_id = t.transcript_id) protein_id,
       (select version
          from translation 
         where translation.transcript_id = t.transcript_id) protein_version
  FROM gene g,
       transcript t,
       seq_region s
 WHERE g.gene_id = t.gene_id
   AND g.seq_region_id = s.seq_region_id
   AND s.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                   '11',  '12', '13', '14', '15', '16', '17', '18', '19', '20',
                   '21', '22', 'X', 'Y', 'MT')                  
 order by g.stable_id, t.stable_id
'''


def connect_to_database(server, port, user_id, password, db):
    """Connect to Ensembl.

    Args:
        server (str): the ip address or server name
        port (int): port number
        user_id (str): Ensembl MySQL userid
        password (str): Ensembl MySQL password
        db (str): the database to use

    Returns:
        a connection to the database

    Raises:
        `pymysql.Error`
    """
    try:
        LOG.debug('Connecting to {} ...'.format(server))
        connection = pymysql.connect(host=server,
                                     port=int(port),
                                     user=user_id,
                                     password=password,
                                     db=db)
                                     #cursorclass=pymysql.cursors.DictCursor)

        LOG.debug('Connected')
        return connection
    except pymysql.Error as e:
        LOG.error('Unable to connect to Ensembl: {}'.format(e))
        raise e


def extract_compara_ids(ensembl_version, strain):
    """Extract the ids from Ensembl Compara.

    Args:
        ensembl_version
        strain_info (:obj:`cc_protein_viewer.create.db.StrainInfo`):
            contains information about the Ensembl strain

    Returns:
        list: a ``list`` of ``dict`` representing ids
    """
    strain_ids = []

    try:
        conn = connect_to_database(ensembl_version['db_server'],
                                   ensembl_version['db_port'],
                                   ensembl_version['db_user'],
                                   ensembl_version['db_password'],
                                   ensembl_version['compara_db'])

        strain_id = strain['strain_id']

        LOG.debug('Strain: {}'.format(strain))
        LOG.debug('Extracting gene ids...')
        with conn.cursor() as cursor:
            num_rows = cursor.execute(SQL_ENSEMBL_GET_IDS_GENE, (strain_id, ))
            LOG.debug('{:,} records returned'.format(num_rows))

            for row in cursor:
                strain_ids.append((None, row[0], row[1], row[2], row[3]))

            LOG.debug('{:,} ids extracted'.format(len(strain_ids)))

        LOG.debug('Extracting protein ids...')
        with conn.cursor() as cursor:
            num_rows = cursor.execute(SQL_ENSEMBL_GET_IDS_PROTEIN, (strain_id, ))

            LOG.debug('{:,} records returned'.format(num_rows))

            for row in cursor:
                strain_ids.append((None, row[0], row[1], row[2], row[3]))

            LOG.debug('{:,} ids extracted'.format(len(strain_ids)))
    except pymysql.Error as e:
        LOG.error('Unable to extract ids from strain: {}').format(e)
        return None

    return strain_ids


def extract_ensembl_gtp(ensembl_version, strain):
    """Extract the gene, transcript and protein information from Ensembl.

    Args:
        ensembl_version
        strain_info (:obj:`cc_protein_viewer.create.db.StrainInfo`):
            contains information about the Ensembl strain

    Returns:
        list: a list of the gene, transcript and protein information
    """
    gtp = []
    genes = {}

    try:
        conn = connect_to_database(ensembl_version['db_server'],
                                   ensembl_version['db_port'],
                                   ensembl_version['db_user'],
                                   ensembl_version['db_password'],
                                   strain['db'])

        LOG.debug('Strain: {}'.format(strain))
        LOG.debug('Extracting transcript, protein information...')
        strain_id = strain['strain_id']
        with conn.cursor() as cursor:
            count = 0
            num_rows = cursor.execute(SQL_ENSEMBL_STRAIN_GTP)
            LOG.debug('{:,} records returned'.format(num_rows))

            for row in cursor:
                d = [None]
                d.extend((row[x] for x in range(0, 15)))
                d.append(strain_id)
                gtp.append(d)
                genes[row[0]] = 1

                if count and count % 10000 == 0:
                    LOG.debug('{:,} rows, {:,} genes extracted'.format(count, len(genes)))
                count += 1

            LOG.debug('{:,} genes extracted'.format(len(genes)))
    except pymysql.Error as e:
        LOG.error('Unable to extract genes from ensembl: {}').format(e)
        return None

    return gtp




