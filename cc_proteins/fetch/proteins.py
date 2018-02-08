# -*- coding: utf_8 -*-

import sqlite3

from collections import OrderedDict

import cc_proteins.utils as utils


LOG = utils.get_logger()

PROTEINS = OrderedDict()
PROTEINS['A'] = 'alanine'
PROTEINS['B'] = 'aspartate or asparagine'
PROTEINS['C'] = 'cysteine'
PROTEINS['D'] = 'aspartate'
PROTEINS['E'] = 'glutamate'
PROTEINS['F'] = 'phenylalanine'
PROTEINS['G'] = 'glycine'
PROTEINS['H'] = 'histidine'
PROTEINS['I'] = 'isoleucine'
PROTEINS['K'] = 'lysine'
PROTEINS['L'] = 'leucine'
PROTEINS['M'] = 'methionine'
PROTEINS['N'] = 'asparagine'
PROTEINS['O'] = 'pyrrolysine'
PROTEINS['P'] = 'proline'
PROTEINS['Q'] = 'glutamine'
PROTEINS['R'] = 'arginine'
PROTEINS['S'] = 'serine'
PROTEINS['T'] = 'threonine'
PROTEINS['U'] = 'selenocysteine'
PROTEINS['V'] = 'valine'
PROTEINS['W'] = 'tryptophan'
PROTEINS['Y'] = 'tyrosine'
PROTEINS['Z'] = 'glutamic acid or glutamine'


STRAINS = {'mus_musculus_aj': {'key':'A', 'color':"#F0F000", 'name':"A/J"},
           'mus_musculus': {'key':'B', 'color':"#808080", 'name':"C57BL/6J"},
           'mus_musculus_129s1svimj': {'key':'C', 'color':"#F08080", 'name':"129S1/SvImJ"},
           'mus_musculus_nodshiltj': {'key':'D', 'color':"#1010F0", 'name':"NOD/ShiLtJ"},
           'mus_musculus_nzohlltj': {'key':'E', 'color':"#00A0F0", 'name':"NZO/H1LtJ"},
           'mus_musculus_casteij': {'key':'F', 'color':"#00A000", 'name':"CAST/EiJ"},
           'mus_musculus_pwkphj': {'key':'G', 'color':"#F00000", 'name':"PWK/PhJ"},
           'mus_musculus_wsbeij': {'key':'H', 'color':"#9000E0", 'name':"WSB/EiJ"}}


class ProteinException(Exception):
    """Protein exception class."""
    pass


def search_proteins(database, sequence):
    """Search for proteins

    Args:
        database (str): the sqlite database
        sequence (str): the sequence to find

    Returns:
        list: of matching gene information

    Raises:
        ProteinException: when sqlite error or other error occurs
    """
    results = OrderedDict()

    try:
        sql_search = ('SELECT fr.* '
                      '  FROM fasta_record fr, '
                      '       strains s '
                      ' WHERE fr.strain_id = s.strain_id '
                      '   AND fr.ref_gene_id IS NOT NULL '
                      '   AND fr.seq LIKE :sequence '
                      ' ORDER BY fr.ref_gene_id, s.strain_order ')

        conn = sqlite3.connect(database)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        for row in cursor.execute(sql_search, {'sequence': sequence}):
            ref_gene_id = row['ref_gene_id']
            print(row)

            gene = results.get(ref_gene_id, {'ensembl_id': ref_gene_id,
                                             'symbol': row['gene_name'],
                                             'proteins': {}})

            proteins = gene['proteins']

            ref_protein_id = row['ref_protein_id']
            if not ref_protein_id:
                ref_protein_id = row['protein_id']

            protein = proteins.get(ref_protein_id, {'A': 0,
                                                    'B': 0,
                                                    'C': 0,
                                                    'D': 0,
                                                    'E': 0,
                                                    'F': 0,
                                                    'G': 0,
                                                    'H': 0})

            strain_key = STRAINS[row['strain_id']]['key']
            protein[strain_key] += 1
            proteins[ref_protein_id] = protein
            gene['proteins'] = proteins
            results[ref_gene_id] = gene

    except sqlite3.Error as e:
        raise ProteinException(e)

    return list(results.values())


def search_genes(database, term):
    """Search for proteins by gene.

    Args:
        database (str): The sqlite database.
        term (str): The search term.

    Returns:
        list: Matching gene information.

    Raises:
        GeneException: when sqlite error or other error occurs
    """
    results = OrderedDict()

    try:
        sql_search = ('SELECT fr.* '
                      '  FROM fasta_record fr, '
                      '       strains s '
                      ' WHERE fr.strain_id = s.strain_id '
                      '   AND (fr.ref_gene_id like :term '
                      '    OR fr.gene_name like :term )'
                      ' ORDER BY fr.gene_name, s.strain_order')

        conn = sqlite3.connect(database)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        search_term = term

        while search_term[0] == '%':
            search_term = search_term[1:]
        while search_term[-1] == '%':
            search_term = search_term[:-1]

        search_term = '%{}%'.format(search_term)

        for row in cursor.execute(sql_search, {'term': search_term}):
            ref_gene_id = row['ref_gene_id']
            print(row['ref_gene_id'], row['ref_protein_id'], row['protein_id'])

            gene = results.get(ref_gene_id, {'ensembl_id': ref_gene_id,
                                             'symbol': row['gene_name'],
                                             'strains': {}})


            strains = gene['strains']
            strain_key = STRAINS[row['strain_id']]['key']

            strain = strains.get(strain_key, {'count': 0,
                                              'ref_protein_id': None})

            if row['ref_protein_id']:
                ref_protein_id = row['ref_protein_id']
                if strain['ref_protein_id']:
                    if ref_protein_id not in strain['ref_protein_id']:
                        strain['ref_protein_id'].append(ref_protein_id)
                else:
                    strain['ref_protein_id'] = [ref_protein_id]

            strain['count'] = strain['count'] + 1

            strains[strain_key] = strain
            gene['strains'] = strains
            results[ref_gene_id] = gene

    except sqlite3.Error as e:
        raise ProteinException(e)

    print('returning=', list(results.values()))

    return list(results.values())


def gene_fasta(database, ensembl_id):
    """Search for proteins by gene.

    Args:
        database (str): The sqlite database.
        term (str): The search term.

    Returns:
        list: Matching gene information.

    Raises:
        GeneException: when sqlite error or other error occurs
    """
    results = OrderedDict()

    try:
        sql_search = ('SELECT fr.* '
                      '  FROM fasta_record fr, '
                      '       strains s '
                      ' WHERE fr.strain_id = s.strain_id '
                      '   AND fr.ref_gene_id = :ensembl_id '
                      ' ORDER BY s.strain_order')

        conn = sqlite3.connect(database)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        for row in cursor.execute(sql_search, {'ensembl_id': ensembl_id}):
            results[row['protein_id']] = row['seq']
    except sqlite3.Error as e:
        raise ProteinException(e)

    return results

import zlib
import struct
import time



def gene_fasta_yield(database, ensembl_id, protein_nums=None):
    """Search for proteins by gene.

    Args:
        database (str): The sqlite database.
        term (str): The search term.

    Returns:
        list: Matching gene information.

    Raises:
        GeneException: when sqlite error or other error occurs
    """
    try:
        sql_search = ('SELECT fr.* '
                      '  FROM fasta_record fr, '
                      '       strains s '
                      ' WHERE fr.strain_id = s.strain_id '
                      '   AND fr.ref_gene_id = :ensembl_id '
                      ' ORDER BY s.strain_order')

        conn = sqlite3.connect(database)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        yield (
                b'\037\213\010\000' +  # Gzip file, deflate, no filename
                struct.pack('<L', int(time.time())) +  # compression start time
                b'\002\377'  # maximum compression, no OS specified
        )

        # bookkeeping: the compression state, running CRC and total length
        compressor = zlib.compressobj(
                9, zlib.DEFLATED, -zlib.MAX_WBITS, zlib.DEF_MEM_LEVEL, 0)
        crc = zlib.crc32(b"")
        length = 0

        if protein_nums is None:
            protein_nums = []

        print(protein_nums)

        idx = 0

        for row in cursor.execute(sql_search, {'ensembl_id': ensembl_id}):
            if idx in protein_nums:
                data = '>{}\n{}\n'.format(row['protein_id'], row['seq']).encode()
                chunk = compressor.compress(data)
                if chunk:
                    yield chunk
                crc = zlib.crc32(data, crc) & 0xffffffff
                length += len(data)
            idx += 1

        # Finishing off, send remainder of the compressed data, and CRC and length
        yield compressor.flush()
        yield struct.pack(b"<2L", crc, length & 0xffffffff)

        cursor.close()
        conn.close()
    except sqlite3.Error as e:
        raise ProteinException(e)





