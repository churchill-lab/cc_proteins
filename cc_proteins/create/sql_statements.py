# -*- coding: utf-8 -*-

SQL_ENSEMBL_65_GENE_SELECT = '''
SELECT g.gene_id,
       s.name seq_id,
       g.seq_region_start,
       g.seq_region_end,
       g.seq_region_strand,
       g.stable_id ensembl_id,
       g.version ensembl_id_version,
       x1.display_label symbol,
       x1.description,
       x2.xref_id,
       x2.dbprimary_acc external_id,
       x2.display_label external_label,
       x2.version external_version,
       ed2.db_name,
       ed2.db_display_name
  FROM gene g,
       xref x1,
       xref x2,
       object_xref o,
       external_db ed1,
       external_db ed2,
       seq_region s
 WHERE g.display_xref_id = x1.xref_id
   AND g.seq_region_id = s.seq_region_id
   AND o.ensembl_id = g.gene_id
   AND o.xref_id = x2.xref_id
   AND x1.external_db_id = ed1.external_db_id
   AND x2.external_db_id = ed2.external_db_id
   AND ed2.db_name IN ('EntrezGene', 'MGI', 'UniGene', 'EntrezGene', 'HGNC')
'''

SQL_ENSEMBL_65_GENE_SELECT_WHERE = '''
   AND g.stable_id = ?
'''

SQL_ENSEMBL_65_GENE_SELECT_ORDER_BY = '''
 ORDER BY g.stable_id 
'''

SQL_ENSEMBL_65_SELECT_GTPE = '''
SELECT g.stable_id gene_id,
       g.version gene_version,
       x1.display_label gene_name,
       s.name gene_chrom,
       g.seq_region_start gene_start,
       g.seq_region_end gene_stop,
       g.seq_region_strand gene_strand,
       t.stable_id transcript_id,
       t.version transcript_version,
       x2.display_label transcript_name,
       s.name transcript_chrom,
       t.seq_region_start transcript_start,
       t.seq_region_end transcript_stop,
       t.seq_region_strand transcript_strand,
       translation.stable_id protein_id,
       e.stable_id ensembl_exon_id,
       e.version ensembl_exon_version,
       s.name exon_chrom,
       e.seq_region_start exon_start,
       e.seq_region_end exon_stop,
       e.seq_region_strand exon_strand,
       et.rank exon_number
  FROM gene g,
       transcript t
  LEFT JOIN translation ON t.transcript_id = translation.transcript_id,
       exon e,
       exon_transcript et,
       seq_region s,
       xref x1,
       xref x2
 WHERE g.gene_id = t.gene_id
   AND g.seq_region_id = s.seq_region_id
   AND t.transcript_id = et.transcript_id
   AND et.exon_id = e.exon_id
   AND g.display_xref_id = x1.xref_id
   AND t.display_xref_id = x2.xref_id
'''

SQL_ENSEMBL_65_SELECT_GTPE_WHERE = '''
   AND g.stable_id = ?
'''

SQL_ENSEMBL_65_SELECT_GTPE_ORDER_BY = '''
 ORDER BY g.stable_id, t.stable_id, et.rank
'''


SQL_ENSEMBL_65_SELECT_SYNONYMS = '''
SELECT * 
  FROM external_synonym 
 WHERE xref_id IN ?
'''
















SQL_CREATE_GTF_TABLE = '''
    CREATE TABLE IF NOT EXISTS gtf (
        gtf_key INTEGER PRIMARY KEY,
        species_id TEXT NOT NULL,
        gene_id TEXT  NOT NULL,
        transcript_id TEXT,
        ensembl_id TEXT NOT NULL,
        ensembl_id_version INTEGER,
        ensembl_name TEXT,
        seqid TEXT NOT NULL,
        start INTEGER,
        end INTEGER,
        strand INTEGER,
        source_key INTEGER NOT NULL,
        type_key INTEGER NOT NULL,
        exon_number INTEGER
    )
'''

SQL_CREATE_GTF_TYPES_TABLE = '''
    CREATE TABLE IF NOT EXISTS gtf_types (
        _key INTEGER PRIMARY KEY,
        gtf_type TEXT NOT NULL
    )
'''

SQL_CREATE_GTF_SOURCES_TABLE = '''
    CREATE TABLE IF NOT EXISTS gtf_sources (
        _key INTEGER PRIMARY KEY,
        gtf_source TEXT NOT NULL
    )
'''

SQL_CREATE_GENES_TABLE = '''
    CREATE TABLE IF NOT EXISTS ensembl_genes (
       _ensembl_genes_key INTEGER,
       ensembl_id TEXT NOT NULL,
       ensembl_version TEXT,
       species_id TEXT NOT NULL,
       symbol TEXT,
       name TEXT,
       synonyms TEXT,
       external_ids TEXT,
       chromosome TEXT NOT NULL,
       start_position INTEGER NOT NULL,
       end_position INTEGER NOT NULL,
       strand INTEGER NOT NULL,
       complete INTEGER NOT NULL,
       PRIMARY KEY (_ensembl_genes_key)
    )
'''

SQL_CREATE_GENES_LOOKUP_TABLE = '''
    CREATE TABLE IF NOT EXISTS ensembl_genes_lookup (
       _ensembl_genes_lookup_key INTEGER,
       ensembl_gene_id TEXT,
       lookup_value TEXT COLLATE NOCASE,
       ranking_id TEXT,
       species_id TEXT,
       PRIMARY KEY (_ensembl_genes_lookup_key)
    );
'''


SQL_CREATE_GENES_LOOKUP_TMP_TABLE = '''
    CREATE TABLE IF NOT EXISTS ensembl_genes_lookup_tmp (
       ensembl_gene_id TEXT,
       lookup_value TEXT,
       ranking_id TEXT,
       species_id TEXT
    )
'''

SQL_CREATE_ENSEMBL_SEARCH_TABLE = '''
    CREATE VIRTUAL TABLE IF NOT EXISTS ensembl_search
        USING fts4(_ensembl_genes_lookup_key,
        lookup_value)
'''

SQL_INSERT_GTF_TABLE = '''
    INSERT
      INTO gtf
    VALUES (null, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
'''

SQL_INSERT_GTF_TYPES_TABLE = '''
    INSERT
      INTO gtf_types
    VALUES (?, ?)
'''

SQL_INSERT_GTF_SOURCES_TABLE = '''
    INSERT
      INTO gtf_sources
    VALUES (?, ?)
'''

SQL_INSERT_GENES_TABLE = '''
    INSERT
      INTO ensembl_genes
    VALUES (null, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
'''

SQL_INSERT_GENES_LOOKUP_TMP_TABLE = '''
    INSERT
      INTO ensembl_genes_lookup_tmp
    VALUES (?, ?, ?, ?)
'''

SQL_INSERT_GENES_LOOKUP_TABLE = '''
INSERT
  INTO ensembl_genes_lookup
SELECT distinct null,
       ensembl_gene_id,
       lookup_value,
       ranking_id,
       species_id TEXT
  FROM ensembl_genes_lookup_tmp
 WHERE lookup_value is not null
 ORDER BY ensembl_gene_id, lookup_value, ranking_id
'''

SQL_INSERT_GENES_LOOKUP_TABLE_TMP_IDS = '''
    INSERT INTO ensembl_genes_lookup_tmp
    SELECT distinct gene_id, ensembl_id, ?, species_id
      FROM gtf
     WHERE type_key = ?
'''

SQL_INSERT_ENSEMBL_SEARCH_TABLE = '''
    INSERT
      INTO ensembl_search
    SELECT _ensembl_genes_lookup_key, lookup_value
      FROM ensembl_genes_lookup
'''

SQL_UPDATE_GENES_TABLE_EXT = '''
    UPDATE ensembl_genes
       SET name = ?,
           synonyms = ?,
           external_ids = ?,
           complete = 1
     WHERE ensembl_id = ?
'''

SQL_SELECT_GENES_TABLE = '''
    SELECT *
      FROM genes
     ORDER BY _ensembl_genes_key
'''

SQL_SELECT_GENES_TABLE_BY_SPECIES = '''
    SELECT *
      FROM ensembl_genes
     WHERE species_id = :species_id
     ORDER BY _ensembl_genes_key
'''

SQL_SELECT_GTF_TABLE_BY_SPECIES = '''
    SELECT *
      FROM gtf
     WHERE species_id = :species_id
     ORDER BY gtf_key
'''

SQL_SELECT_GTF_TYPES_TABLE = '''
    SELECT *
      FROM gtf_types
'''

SQL_SELECT_GTF_SOURCES_TABLE = '''
    SELECT *
      FROM gtf_sources
'''

SQL_INDICES_GTF = [
    'CREATE INDEX IF NOT EXISTS idx_gtf_species_id ON gtf(species_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_gene_id ON gtf(gene_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_transcript_id ON gtf(transcript_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_ensembl_id ON gtf(ensembl_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_ensembl_id_version ON gtf(ensembl_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_seqid ON gtf(seqid ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_start ON gtf(start ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_end ON gtf(end ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_source_key ON gtf(source_key ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_type_key ON gtf(type_key ASC)',
    'CREATE INDEX IF NOT EXISTS idx_gtf_exon_number ON gtf(exon_number ASC)'
]

SQL_INDICES_GTF_TYPES = [
    'CREATE INDEX IF NOT EXISTS idx_gtf_types_gtf_type ON gtf_types(gtf_type ASC)',
]

SQL_INDICES_GTF_SOURCES = [
    'CREATE INDEX IF NOT EXISTS idx_gtf_sources_gtf_source ON gtf_sources(gtf_source ASC)',
]

SQL_INDICES_GENE = [
    'CREATE INDEX IF NOT EXISTS idx_ensembl_gene_id ON ensembl_genes (ensembl_id ASC)',
]

SQL_INDICES_GENE_LOOKUP = [
    'CREATE INDEX IF NOT EXISTS idx_lookup_ensembl_gene_id ON ensembl_genes_lookup (ensembl_gene_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_lookup_value ON ensembl_genes_lookup (lookup_value ASC)',
    'CREATE INDEX IF NOT EXISTS idx_lookup_id ON ensembl_genes_lookup (ranking_id ASC)',
    'CREATE INDEX IF NOT EXISTS idx_lookup_species_id ON ensembl_genes_lookup (species_id ASC)',
]

SQL_DROP_GENES_LOOKUP_TMP_TABLE = '''
    DROP TABLE ensembl_genes_lookup_tmp
'''


