# -*- coding: utf-8 -*-

from .manager import Manager

__all__ = [
    'enrich_chemicals',
]


def enrich_chemicals(graph, connection=None):
    """Enriches chemicals in the graph

    :param pybel.BELGraph graph: A BEL graph
    :type connection: str or bio2bel_ctd.Manager
    """
    m = Manager.ensure(connection=connection)
    m.enrich_chemicals(graph)
