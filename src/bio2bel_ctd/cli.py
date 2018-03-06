# -*- coding: utf-8 -*-

"""Run this script with :code:`python3 -m bio2bel_ctd`"""

import logging

import click

from bio2bel_ctd import Manager
from bio2bel_ctd.constants import DEFAULT_CACHE_CONNECTION


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


if __name__ == '__main__':
    main()
