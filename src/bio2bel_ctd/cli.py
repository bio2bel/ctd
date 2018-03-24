# -*- coding: utf-8 -*-

"""Run this script with :code:`python3 -m bio2bel_ctd`"""

import logging
import sys

import click

from .constants import DEFAULT_CACHE_CONNECTION
from .manager import Manager
from .models import Action, Chemical, Gene


@click.group()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.pass_context
def main(ctx, connection):
    """Convert CTD to BEL"""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    logging.getLogger('bio2bel.utils').setLevel(logging.ERROR)
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
    """Manage chemicals"""


def _echo_chemical(chemical, interaction_limit=None):
    click.echo('ID: {}'.format(chemical.id))
    click.echo('Name: {}'.format(chemical.chemical_name))
    click.echo('MeSH Identifier: {}'.format(chemical.chemical_id))
    if chemical.cas_rn:
        click.echo('CAS Registry Number: {}'.format(chemical.cas_rn))
    if chemical.definition:
        click.echo('Definition: {}'.format(chemical.definition))

    interactions = chemical.gene_interactions
    if interaction_limit is not None and interaction_limit > 0:
        interactions = interactions.limit(interaction_limit)

    if interactions:
        click.echo('Interactions (sample)')
        for interaction in interactions:
            click.echo(' - ({}) {}'.format(interaction.id, interaction))


@chemicals.command()
@click.argument('mesh_id')
@click.option('-i', '--interaction-limit', type=int, default=5)
@click.pass_obj
def get(manager, mesh_id, interaction_limit):
    """Get a chemical by its MeSH identifier. Try MESH:C490728 for lapatinib"""
    chemical = manager.get_chemical_by_mesh(mesh_id)

    if chemical is None:
        click.echo('MeSH Identifier not found: {}'.format(mesh_id))
        sys.exit(0)

    _echo_chemical(chemical, interaction_limit=interaction_limit)


@chemicals.command()
@click.argument('cas_rn')
@click.option('-i', '--interaction-limit', type=int, default=5)
@click.pass_obj
def get_cas(manager, cas_rn, interaction_limit):
    """Get a chemical by its CAS Registry Number. Try 55-18-5 for Diethylnitrosamine"""
    chemical = manager.get_chemical_by_cas(cas_rn)

    if chemical is None:
        click.echo('CAS Registry Number not found: {}'.format(cas_rn))
        sys.exit(0)

    _echo_chemical(chemical, interaction_limit=interaction_limit)


def _echot(*t):
    click.echo('\t'.join(map(str, t)))


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

    _echot('MeSH', 'Name', 'Definition', 'Parents')

    for chemical in query:
        _echot(
            chemical.chemical_id,
            chemical.chemical_name,
            chemical.definition,
            '|'.join(map(str, chemical.parent_ids))
        )


@manage.group()
def genes():
    """Manage genes"""


@genes.command()
@click.argument('entrez_id')
@click.pass_obj
def get(manager, entrez_id):
    """Get a gene by its Entrez Gene identifier"""
    gene = manager.get_gene_by_entrez_id(entrez_id)

    if gene is None:
        click.echo('Not found: {}'.format(entrez_id))
        sys.exit(0)

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

    _echot('EGID', 'Name', 'Symbol')

    for gene in query:
        _echot(
            gene.gene_id,
            gene.gene_name,
            gene.gene_symbol,
        )


@manage.group()
def ixns():
    """Manage chemical-gene interactions"""


@ixns.command()
@click.argument('ixn_id')
@click.pass_obj
def get(manager, ixn_id):
    """Get a interaction by its database identifier"""
    ixn = manager.get_interaction_by_id(ixn_id)

    if ixn is None:
        click.echo('Interaction not found: {}'.format(ixn_id))
        sys.exit(0)

    click.echo('Chemical: {}'.format(ixn.chemical))
    click.echo('Gene: {}'.format(ixn.gene))
    click.echo('Interaction: {}'.format(ixn.interaction))

    for action in ixn.interaction_actions:
        click.echo('Action: {}'.format(action))

    for gene_form in ixn.gene_forms:
        click.echo('Gene Form: {}'.format(gene_form))


@manage.group()
def actions():
    """Manage chemical-gene interaction actions"""


@actions.command()
@click.pass_obj
def ls(manager):
    _echot('Type Name', 'Code', 'Description', 'Parent Code')
    for action in manager.session.query(Action).all():
        _echot(
            action.type_name,
            action.code,
            action.description,
            action.parent_code
        )


@main.command()
@click.pass_obj
def upload_bel(manager):
    """Upload BEL to database"""
    manager.upload_bel_graph()


if __name__ == '__main__':
    main()
