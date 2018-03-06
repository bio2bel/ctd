# -*- coding: utf-8 -*-

"""Run this script with :code:`python3 -m bio2bel_ctd`"""

import logging

import click

from .constants import DEFAULT_CACHE_CONNECTION
from .manager import Manager
from .models import Chemical, Gene


@click.group()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.pass_context
def main(ctx, connection):
    """Convert CTD to BEL"""
    logging.basicConfig(level=10, format="%(asctime)s - %(levelname)s - %(message)s")
    ctx.obj = Manager(connection=connection)


@main.command()
@click.pass_obj
def populate(manager):
    """Populates the database"""
    manager.populate()


@main.command()
@click.option('-y', '--yes', is_flag=True)
@click.pass_obj
def drop(manager, yes):
    """Drops database"""
    if yes or click.confirm('Drop everything?'):
        manager.drop_all()


@main.command()
@click.pass_obj
def summarize(manager):
    """Summarizes the database"""
    for k, v in manager.summarize().items():
        click.echo('{}: {}'.format(k.replace('_', ' ').title(), v))


@main.command()
@click.option('-v', '--debug', is_flag=True)
@click.option('-h', '--host')
@click.option('-p', '--port', type=int)
@click.pass_obj
def web(manager, debug, host, port):
    """Run the web app"""
    from .web import get_app
    app = get_app(connection=manager, url='/')
    app.run(debug=debug, host=host, port=port)


@main.group()
def manage():
    pass


@manage.group()
def chemicals():
    pass


@chemicals.command()
@click.argument('mesh_id')
@click.pass_obj
def get(manager, mesh_id):
    """Get a chemical by its MeSH identifier. Try MESH:C490728 for lapatinib"""
    chemical = manager.get_chemical_by_mesh(mesh_id)

    if chemical is None:
        click.echo('Not found: {}'.format(mesh_id))
    else:
        click.echo('MeSH Identifier: {}'.format(chemical.chemical_id))
        click.echo('Name: {}'.format(chemical.chemical_name))
        if chemical.definition:
            click.echo('Definition: {}'.format(chemical.definition))


@chemicals.command()
@click.option('--limit', type=int, default=5)
@click.option('--offset', type=int)
@click.pass_obj
def ls(manager, limit, offset):
    query = manager.session.query(Chemical)

    if limit > 0:
        query = query.limit(limit)

    if offset is not None:
        query = query.offset(offset)

    click.echo('\t'.join(('MeSH', 'Name', 'Definition', 'Parents')))

    for chemical in query:
        click.echo('\t'.join(map(str, (
            chemical.chemical_id,
            chemical.chemical_name,
            chemical.definition,
            '|'.join(map(str, chemical.parent_ids))
        ))))


@manage.group()
def genes():
    pass


@genes.command()
@click.argument('entrez_id')
@click.pass_obj
def get(manager, entrez_id):
    """Get a gene by its Entrez Gene identifier"""
    gene = manager.get_gene_by_entrez_id(entrez_id)

    if gene is None:
        click.echo('Not found: {}'.format(entrez_id))
    else:
        click.echo('Entrez Gene Identifier: {}'.format(gene.gene_id))
        click.echo('Name: {}'.format(gene.gene_name))
        click.echo('Symbol: {}'.format(gene.gene_symbol))

        for ixn in gene.chemical_interactions.limit(5):
            click.echo(ixn)


@genes.command()
@click.option('--limit', type=int, default=5)
@click.option('--offset', type=int)
@click.pass_obj
def ls(manager, limit, offset):
    query = manager.session.query(Gene)

    if limit > 0:
        query = query.limit(limit)

    if offset is not None:
        query = query.offset(offset)

    click.echo('\t'.join(('EGID', 'Name', 'Symbol')))

    for gene in query:
        click.echo('\t'.join(map(str, (
            gene.gene_id,
            gene.gene_name,
            gene.gene_symbol,
        ))))


if __name__ == '__main__':
    main()
