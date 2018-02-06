# -*- coding: utf-8 -*-
import os
import time

from Bio.Align.Applications import ClustalOmegaCommandline
from flask import current_app

from cc_proteins.app import create_celery_app
from cc_proteins.fetch import proteins
import cc_proteins.utils as utils

celery = create_celery_app()


@celery.task(ignore_result=True, bind=True)
def search_proteins(self, sequences):
    """Celery task to search for peptide strings in proteins.

    Args:
        self: Celery reference
        sequences (list): peptide strings to search for

    Returns:
        dict: a dictionary specifying the state and results
    """
    results = {}
    total_sequences = len(sequences)

    for i, sequence in enumerate(sequences):
        print(sequence)
        temp_seq = sequence.replace('*', '%').replace('?', '_')
        print(temp_seq)
        print(current_app.config['FASTA_DB'])
        response = proteins.search_proteins(current_app.config['FASTA_DB'], temp_seq)
        results[sequence] = response
        self.update_state(state='PROGRESS',
                          meta={'current': i+1, 'total': total_sequences}),
        time.sleep(1)

    self.update_state(state='SUCCESS',
                      meta={'current': total_sequences, 'total': total_sequences, 'result':results}),
    return {'current': total_sequences, 'total': total_sequences, 'result':results}



@celery.task(ignore_result=True, bind=True)
def clustalo_gene(self, gene_id):
    """Celery task to search for peptide strings in proteins.

    Args:
        self: Celery reference
        sequences (list): peptide strings to search for

    Returns:
        dict: a dictionary specifying the state and results
    """
    meta = {'gene': gene_id}

    print('in clustalo_gene')

    self.update_state(state='PROGRESS', meta=meta)

    try:
        # get config
        tmp_dir = '/app/cc_proteins/tmp'
        db = current_app.config['FASTA_DB']

        # get the sequences for the gene
        sequences = proteins.gene_fasta(db, gene_id)

        # dump to file
        fname = utils.create_random_string(20)
        input_file = os.path.join(tmp_dir, '{}.fasta'.format(fname))
        output_file = os.path.join(tmp_dir, '{}.clustal'.format(fname))

        print('input_file=', input_file)
        print('output_file=', output_file)


        num_sequences = 0
        with open(input_file, "w") as fd:
            for gene, sequence in sequences.items():
                fd.write(">{}\n".format(gene))
                fd.write("{}\n".format(sequence))
                num_sequences += 1

        meta['num_sequences'] = num_sequences
        self.update_state(state='PROGRESS', meta=meta),

        # perform clustalo
        clustalomega_cline = ClustalOmegaCommandline(infile=input_file,
                                                     outfile=output_file,
                                                     outfmt='clustal',
                                                     verbose=True,
                                                     auto=True)

        clustalomega_cline()

        print('Clustalo done')

        # Initial implementation had state='FAILURE', but ran into several issues.
        #
        # store_errors_even_if_ignored=True (IF using 'FAILURE' state)
        #
        # Essentially the task is completing even on request errors, it's just
        # getting bad responses from the request.
        #
        # Instead of FAILURE, we will set SUCCESS and handle it on the client side.

        meta['output_id'] = fname
    except Exception as e:
        print('Exception in call_api task: ', str(e))
        meta['error'] = 'Error'
        meta['error_message'] = str(e)

    self.update_state(state='SUCCESS', meta=meta)

    return meta

