import time

import click

import cc_proteins.create.fasta2db as fasta2db
import cc_proteins.utils as utils



@click.command('generate', options_metavar='<options>',
             short_help='create an annotation database')
@click.argument('database', metavar='<database>',
                type=click.Path(file_okay=True, exists=False,
                                resolve_path=True, dir_okay=False))
@click.option('-r', '--resource', default=fasta2db.DEFAULT_CONFIG)
@click.option('--version', default=None)
@click.option('-v', '--verbose', count=True)
def cli(database, resource, version, verbose):
    """
    Creates a new protein sequence database
    """
    utils.configure_logging(verbose)
    LOG = utils.get_logger()
    LOG.info("Creating database...")

    tstart = time.time()
    fasta2db.generate_database(database, resource, version)
    tend = time.time()

    LOG.info("Creation time: {}".format(utils.format_time(tstart, tend)))


