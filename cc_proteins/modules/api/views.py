# -*- coding: utf_8 -*-
import os

from functools import wraps

from flask import Blueprint
from flask import Response
from flask import current_app
from flask import jsonify
from flask import request
from flask import send_file
from flask import stream_with_context

from cc_proteins.fetch import proteins
from cc_proteins import utils


api = Blueprint('api', __name__, template_folder='templates', url_prefix='/api')


def support_jsonp(func):
    """Wraps JSONified output for JSONP requests."""

    @wraps(func)
    def decorated_function(*args, **kwargs):
        callback = request.args.get('callback', False)
        if callback:
            resp = func(*args, **kwargs)
            resp.set_data('{}({})'.format(str(callback),
                                          resp.get_data(as_text=True)))
            resp.mimetype = 'application/javascript'
            return resp
        else:
            return func(*args, **kwargs)

    return decorated_function


@api.route('/search/genes/<term>', methods=['GET'])
def search_genes(term):

    print('calling search...')

    results = proteins.search_genes(current_app.config['FASTA_DB'], term)

    print('search was called')

    return jsonify(results)


@api.route('/search/proteins', methods=['POST'])
def search_proteins():
    import cc_proteins.modules.api.tasks as protein_tasks

    json_data = request.get_json()
    sequences = json_data['sequences']
    print('calling search...')
    task = protein_tasks.search_proteins.apply_async(args=[sequences])
    print('search was called')
    return jsonify({'task_id': task.task_id})


@api.route('/search/proteins/status/<task_id>')
def search_proteins_status(task_id):
    """
    Query Celery for task status based on its taskid/jobid

    Celery status can be one of:
    PENDING - Job not yet run or unknown status
    PROGRESS - Job is currently running
    SUCCESS - Job completed successfully
    FAILURE - Job failed
    REVOKED - Job get canceled
    """
    import cc_proteins.modules.api.tasks as protein_tasks

    task = protein_tasks.search_proteins.AsyncResult(task_id)
    print('/clustalo/proteins/status/', task, task.state)

    if task.state == 'PENDING':
        # job did not start yet
        response = {
            'state': task.state,
            'current': 0,
            'total': 1,
            'status': 'Pending...'
        }
    elif task.state != 'FAILURE':
        response = {
            'state': task.state,
            'current': task.info.get('current', 0),
            'total': task.info.get('total', 1),
            'status': task.info.get('status', '')
        }
        print('TASK_INFO=', str(task.info))
        if 'result' in task.info:
            response['result'] = task.info['result']
            response['status'] = 'DONE'
    else:
        # something went wrong in the background job
        response = {
            'state': task.state,
            'current': 1,
            'total': 1,
            'status': str(task.info),  # this is the exception raised
        }
    return jsonify(response)


@api.route('/clustalo/gene', methods=['POST'])
def clustalo_gene():
    import cc_proteins.modules.api.tasks as protein_tasks
    print('calling search...')
    ensembl_id = request.values.get('ensemblID')
    task = protein_tasks.clustalo_gene.apply_async(args=[ensembl_id])
    print('search was called')
    return jsonify({'task_id': task.task_id})


@api.route('/clustalo/gene/status/<task_id>')
def clustalo_gene_status(task_id):
    """
    Query Celery for task status based on its taskid/jobid

    Celery status can be one of:
    PENDING - Job not yet run or unknown status
    PROGRESS - Job is currently running
    SUCCESS - Job completed successfully
    FAILURE - Job failed
    REVOKED - Job get canceled
    """
    import cc_proteins.modules.api.tasks as protein_tasks

    task = protein_tasks.clustalo_gene.AsyncResult(task_id)

    print('/clustalo/gene/status/', task, task.state)

    if task.state == 'PENDING':
        # job did not start yet
        response = {
            'state': task.state,
            'status': 'Pending...'
        }
    elif task.state != 'FAILURE':
        response = {
            'state': task.state,
            'status': 'RUNNING'
        }
        if 'output_id' in task.info:
            response['status'] = 'DONE'
    else:
        # something went wrong in the background job
        response = {
            'state': task.state,
            'status': str(task.info),  # this is the exception raised
        }
    return jsonify(response)


@api.route('/clustalo/gene/data')
def clustalo_gene_data():
    """
    Query Celery for task status based on its taskid/jobid

    Celery status can be one of:
    PENDING - Job not yet run or unknown status
    PROGRESS - Job is currently running
    SUCCESS - Job completed successfully
    FAILURE - Job failed
    REVOKED - Job get canceled
    """
    import cc_proteins.modules.api.tasks as protein_tasks

    task_id = request.values.get('task_id', '')

    try:
        task = protein_tasks.clustalo_gene.AsyncResult(task_id)
        print('/clustalo/gene/data/', task, task.state)
        if 'output_id' in task.info:
            tmp_dir = '/app/cc_proteins/tmp'

            # try to keep the temp directory cleaned
            fname = '{}.fasta'.format(task.info['output_id'])
            fasta_file = os.path.join(tmp_dir, fname)
            utils.delete_file(fasta_file)

            fname = '{}.clustal'.format(task.info['output_id'])
            output_file = os.path.join(tmp_dir, fname)
            return send_file(output_file, mimetype="application/octet-stream")
    except Exception as e:
        print(str(e))
        return 'Error: '  + str(e)


@api.route("/gene/fasta", methods=['GET'])
def streamed_fasta_response():
    ensembl_id = request.values.get('ensemblID', '')
    response = Response(stream_with_context(proteins.gene_fasta_yield(current_app.config['FASTA_DB'], ensembl_id)), mimetype='application/gzip')
    response.headers['Content-Disposition'] = 'attachment; filename={}.fasta.gz'.format(ensembl_id)
    return response