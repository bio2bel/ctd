# -*- coding: utf-8 -*-

import pyctd.manager.database
from bio2bel.utils import get_connection
from pyctd.manager.database import DbManager
from pyctd.manager.query import QueryManager
from .constants import DATA_DIR, MODULE_NAME
from .models import ChemGeneIxn, Chemical, Disease, Gene, Pathway

__all__ = [
    'Manager'
]


def _get_connection_string(connection):
    return get_connection(module_name=MODULE_NAME, connection=connection)


# Monkey patch PyCTD connection loader
pyctd.manager.database.get_connection_string = _get_connection_string


class Manager(QueryManager, DbManager):
    # Override the directory in which data gets stored
    pyctd_data_dir = DATA_DIR

    def populate(self, *args, **kwargs):
        """Populates the database"""
        self.db_import(*args, **kwargs)

    def _count_model(self, model):
        return self.session.query(model).count()

    def count_genes(self):
        """Counts the genes in the database

        :rtype: int
        """
        return self._count_model(Gene)

    def count_chemicals(self):
        """Counts the chemicals in the database

        :rtype: int
        """
        return self._count_model(Chemical)

    def count_chemical_gene_interactions(self):
        """Counts the chemical-gene interactions in the database

        :rtype: int
        """
        return self._count_model(ChemGeneIxn)

    def count_pathways(self):
        return self._count_model(Pathway)

    def count_diseases(self):
        return self._count_model(Disease)

    def summarize(self):
        """Returns a summary dictionary of the database

        :rtype: dict[str,int]
        """
        return dict(
            chemicals=self.count_chemicals(),
            genes=self.count_genes(),
            chemical_gene_interactions=self.count_chemical_gene_interactions(),
            diseases=self.count_diseases(),
            pathways=self.count_pathways(),
        )

    def enrich_chemicals(self, graph):
        """Finds chemicals that can be mapped and enriched with the CTD

        :param pybel.BELGraph graph:
        """
        raise NotImplementedError

    @staticmethod
    def ensure(connection=None):
        """Checks and allows for a Manager to be passed to the function.

        :param connection: can be either a already build manager or a connection string to build a manager with.
        """
        if connection is None or isinstance(connection, str):
            return Manager(connection=connection)

        if isinstance(connection, Manager):
            return connection

        raise TypeError
