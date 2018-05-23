# -*- coding: utf-8 -*-

"""Bio2BEL CTD Manager."""

import logging

from tqdm import tqdm

import pybel
import pyctd
import pyctd.manager
import pyctd.manager.database
from bio2bel import AbstractManager
from bio2bel.utils import get_connection
from pybel import BELGraph
from pybel.constants import IDENTIFIER, NAME, NAMESPACE
from pyctd.manager.database import DbManager
from pyctd.manager.query import QueryManager
from pyctd.manager.table import get_table_configurations
from .constants import DATA_DIR, MODULE_NAME
from .enrichment_utils import add_chemical_gene_interaction
from .models import Base, ChemGeneIxn, Chemical, Disease, Gene, Pathway

__all__ = [
    'Manager'
]

log = logging.getLogger(__name__)


def _get_connection_string(connection):
    return get_connection(module_name=MODULE_NAME, connection=connection)


# Monkey patch PyCTD connection loader
pyctd.manager.database.get_connection_string = _get_connection_string

_exclude_tables = {
    'exposure_event'
}


class _PyCTDManager(QueryManager, DbManager):
    # Override the directory in which data gets stored
    pyctd_data_dir = DATA_DIR


def _get_urls():
    return [
        pyctd.manager.defaults.url_base + pyctd.manager.table_conf.tables[model]['file_name']
        for model in pyctd.manager.table_conf.tables
    ]


class Manager(AbstractManager, _PyCTDManager):
    """Bio2BEL manager for the CTD."""

    module_name = MODULE_NAME

    # Compensate for some weird structuring of PyCTD code
    tables = get_table_configurations()

    @property
    def _base(self):
        return Base

    def is_populated(self):
        """Check if the database is already populated.

        :rtype: bool
        """
        return 0 < self.count_chemical_gene_interactions()

    def populate(self, urls=None, force_download=False, only_tables=None, exclude_tables=None):
        """Updates the CTD database

        1. downloads all files from CTD
        2. drops all tables in database
        3. creates all tables in database
        4. import all data from CTD files

        :param iter[str] urls: An iterable of URL strings
        :param bool force_download: force method to download
        """
        if not urls:
            urls = _get_urls()

        log.info('Update CTD database from %s', urls)

        self.drop_all()
        self.download_urls(urls=urls, force_download=force_download)
        self.create_all()
        self.import_tables(only_tables=only_tables, exclude_tables=(exclude_tables or _exclude_tables))
        self.session.close()

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

    def list_chemicals(self):
        """List all chemicals.

        :rtype: list[Chemical]
        """
        return self._list_model(Chemical)

    def count_chemical_gene_interactions(self):
        """Counts the chemical-gene interactions in the database

        :rtype: int
        """
        return self._count_model(ChemGeneIxn)

    def count_pathways(self):
        """Counts the pathways in the database

        :rtype: int
        """
        return self._count_model(Pathway)

    def count_diseases(self):
        """Counts the diseases in the database

        :rtype: int
        """
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

    def get_chemical_by_cas(self, cas_rn):
        """Gets a chemical by its CAS Registry Number

        :param str cas_rn: A CAS Registry Number
        :rtype: Optional[Chemical]
        """
        return self.session.query(Chemical).filter(Chemical.cas_rn == cas_rn).one_or_none()

    def get_gene_by_entrez_id(self, entrez_id):
        """Gets a gene by its Entrez Gene identifier

        :param str entrez_id: An Entrez Gene identifier of a gene
        :rtype: Optional[Gene]
        """
        return self.session.query(Gene).filter(Gene.gene_id == entrez_id).one_or_none()

    def get_interaction_by_id(self, ixn_id):
        """Gets an interaction by its database identifier

        :param int ixn_id: An interaction database identifier
        :rtype: Optional[ChemGeneIxn]
        """
        return self.session.query(ChemGeneIxn).filter(ChemGeneIxn.id == ixn_id).one_or_none()

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

    def enrich_graph_genes(self, graph):
        """Enriches the BEL graph with chemical-gene interactions for all Entrez genes

        :param pybel.BELGraph graph: A BEL graph
        """
        for gene_node, data in graph.nodes(data=True):
            namespace = data.get(NAMESPACE)
            if namespace not in {'EG', 'EGID', 'ENTREZ'}:
                continue

            identifier = data.get(IDENTIFIER)
            name = data.get(NAME)

            if identifier is not None:
                self.enrich_graph_gene(graph, identifier)
            elif name is not None:
                self.enrich_graph_gene(graph, name)
            else:
                raise KeyError

    def enrich_chemicals(self, graph):
        """Finds chemicals that can be mapped and enriched with the CTD

        :param pybel.BELGraph graph:
        """
        raise NotImplementedError

    def to_bel_graph(self):
        """Converts all possible aspects of the database to BEL.

        :rtype: pybel.BELGraph

        .. warning:: Not complete!

        To do:

        - add namespaces
        - use cursors
        - multiprocessing
        """
        graph = BELGraph(name='CTD', version='1.0.0')

        for ixn in tqdm(self.session.query(ChemGeneIxn).all(), total=self.count_chemical_gene_interactions()):
            add_chemical_gene_interaction(graph, ixn)

        return graph
