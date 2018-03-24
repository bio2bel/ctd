# -*- coding: utf-8 -*-

import pyctd.manager.database
from bio2bel.abstractmanager import AbstractManager
from bio2bel.utils import get_connection
from pyctd.manager.database import DbManager
from pyctd.manager.models import Base
from pyctd.manager.query import QueryManager
from .constants import DATA_DIR, MODULE_NAME
from .enrichment_utils import add_chemical_gene_interaction
from .models import ChemGeneIxn, Chemical, Disease, Gene, Pathway

__all__ = [
    'Manager'
]


def _get_connection_string(connection):
    return get_connection(module_name=MODULE_NAME, connection=connection)


# Monkey patch PyCTD connection loader
pyctd.manager.database.get_connection_string = _get_connection_string


class _PyCTDManager(QueryManager, DbManager):
    # Override the directory in which data gets stored
    pyctd_data_dir = DATA_DIR


class Manager(AbstractManager, _PyCTDManager):
    module_name = MODULE_NAME

    @property
    def base(self):
        return Base

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

    def get_chemical_by_mesh(self, mesh_id):
        """Gets a chemical by its MeSH identifier

        :param str mesh_id: A MeSH identifier of a chemical
        :rtype: Optional[Chemical]
        """
        return self.session.query(Chemical).filter(Chemical.chemical_id == mesh_id).one_or_none()

    def get_gene_by_entrez_id(self, entrez_id):
        """Gets a gene by its Entrez Gene identifier

        :param str entrez_id: An Entrez Gene identifier of a gene
        :rtype: Optional[Gene]
        """
        return self.session.query(Gene).filter(Gene.gene_id == entrez_id).one_or_none()

    def enrich_graph_chemical(self, graph, mesh_id):
        """Enriches the BEL graph with chemical-gene interactions for the given chemical

        :param pybel.BELGraph graph: A BEL graph
        :param mesh_id: A MeSH identifier of a chemical
        """
        chemical = self.get_chemical_by_mesh(mesh_id)
        if chemical is None:
            return

        for ixn in chemical.gene_interactions:
            add_chemical_gene_interaction(graph, ixn)

    def enrich_graph_gene(self, graph, entrez_id):
        """Enriches the BEL graph with chemical-gene interactions for the given gene

        :param pybel.BELGraph graph: A BEL graph
        :param str entrez_id: An Entrez Gene identifier of a gene
        """
        gene = self.get_gene_by_entrez_id(entrez_id)
        if gene is None:
            return

        for ixn in gene.chemical_interactions:
            add_chemical_gene_interaction(graph, ixn)

    def enrich_chemicals(self, graph):
        """Finds chemicals that can be mapped and enriched with the CTD

        :param pybel.BELGraph graph:
        """
        raise NotImplementedError
