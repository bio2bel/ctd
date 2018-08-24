# -*- coding: utf-8 -*-

"""Bio2BEL CTD Manager."""

import logging
from typing import List, Mapping, Optional

import pyctd
import pyctd.manager
import pyctd.manager.database
from pyctd.manager.database import DbManager
from pyctd.manager.query import QueryManager
from pyctd.manager.table import get_table_configurations
from sqlalchemy.ext.declarative import DeclarativeMeta
from tqdm import tqdm

import bio2bel_mesh
import pybel
from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from bio2bel.utils import get_connection
from pybel import BELGraph
from pybel.constants import IDENTIFIER, NAME, NAMESPACE
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

    'exposure_event',
}


class _PyCTDManager(QueryManager, DbManager):
    # Override the directory in which data gets stored
    pyctd_data_dir = DATA_DIR


def _get_urls():
    return [
        pyctd.manager.defaults.url_base + pyctd.manager.table_conf.tables[model]['file_name']
        for model in pyctd.manager.table_conf.tables
    ]


class Manager(AbstractManager, BELManagerMixin, FlaskMixin, _PyCTDManager):
    """Bio2BEL manager for the CTD."""

    module_name = MODULE_NAME
    flask_admin_models = [Gene, Chemical, Disease, Pathway, ChemGeneIxn]
    # Compensate for some weird structuring of PyCTD code
    tables = get_table_configurations()

    @property
    def _base(self) -> DeclarativeMeta:
        return Base

    def is_populated(self) -> bool:
        """Check if the database is already populated."""
        return 0 < self.count_chemical_gene_interactions()

    def populate(self, urls=None, force_download=False, only_tables=None, exclude_tables=None) -> None:
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

        log.info('downloading CTD database from %s', urls)
        self.download_urls(urls=urls, force_download=force_download)
        log.info('importing tables')
        self.import_tables(only_tables=only_tables, exclude_tables=(exclude_tables or _exclude_tables))

    def count_genes(self) -> int:
        """Count the genes in the database."""
        return self._count_model(Gene)

    def list_chemicals(self) -> List[Chemical]:
        """List all chemicals."""
        return self._list_model(Chemical)

    def count_chemicals(self) -> int:
        """Count the chemicals in the database."""
        return self._count_model(Chemical)

    def list_chemical_gene_interactions(self) -> List[ChemGeneIxn]:
        """List all chemical-gene interactions."""
        return self._list_model(ChemGeneIxn)

    def count_chemical_gene_interactions(self) -> int:
        """Count the chemical-gene interactions in the database."""
        return self._count_model(ChemGeneIxn)

    def count_pathways(self) -> int:
        """Count the pathways in the database."""
        return self._count_model(Pathway)

    def count_diseases(self) -> int:
        """Count the diseases in the database."""
        return self._count_model(Disease)

    def summarize(self) -> Mapping[str, int]:
        """Return a summary dictionary of the database."""
        return dict(
            chemicals=self.count_chemicals(),
            genes=self.count_genes(),
            chemical_gene_interactions=self.count_chemical_gene_interactions(),
            diseases=self.count_diseases(),
            pathways=self.count_pathways(),
        )

    def get_chemical_by_mesh(self, mesh_id: str) -> Optional[Chemical]:
        """Get a chemical by its MeSH identifier, if it exists.

        :param mesh_id: A MeSH identifier of a chemical
        """
        return self.session.query(Chemical).filter(Chemical.chemical_id == mesh_id).one_or_none()

    def get_chemical_by_cas(self, cas_rn: str) -> Optional[Chemical]:
        """Get a chemical by its CAS Registry Number, if it exists.

        :param str cas_rn: A CAS Registry Number
        :rtype: Optional[Chemical]
        """
        return self.session.query(Chemical).filter(Chemical.cas_rn == cas_rn).one_or_none()

    def get_gene_by_entrez_id(self, entrez_id: str) -> Optional[Gene]:
        """Get a gene by its Entrez Gene identifier, if it exists.

        :param entrez_id: An Entrez Gene identifier of a gene
        :rtype: Optional[Gene]
        """
        return self.session.query(Gene).filter(Gene.gene_id == entrez_id).one_or_none()

    def get_interaction_by_id(self, ixn_id: int) -> Optional[ChemGeneIxn]:
        """Get an interaction by its database identifier

        :param ixn_id: An interaction database identifier
        """
        return self.session.query(ChemGeneIxn).filter(ChemGeneIxn.id == ixn_id).one_or_none()

    def enrich_graph_chemical(self, graph: BELGraph, mesh_id: str) -> None:
        """Enrich the BEL graph with chemical-gene interactions for the given chemical.

        :param graph: A BEL graph
        :param mesh_id: A MeSH identifier of a chemical
        """
        chemical = self.get_chemical_by_mesh(mesh_id)
        if chemical is None:
            return

        for ixn in chemical.gene_interactions:
            add_chemical_gene_interaction(graph, ixn)

    def enrich_graph_gene(self, graph: BELGraph, entrez_id: str) -> None:
        """Enrich the BEL graph with chemical-gene interactions for the given gene.

        :param graph: A BEL graph
        :param entrez_id: An Entrez Gene identifier of a gene
        """
        gene = self.get_gene_by_entrez_id(entrez_id)
        if gene is None:
            return

        for ixn in gene.chemical_interactions:
            add_chemical_gene_interaction(graph, ixn)

    def enrich_graph_genes(self, graph: BELGraph) -> None:
        """Enrich the BEL graph with chemical-gene interactions for all Entrez genes.

        :param graph: A BEL graph
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

    def enrich_chemicals(self, graph: BELGraph) -> None:
        """Find chemicals that can be mapped and enriched with the CTD.

        :param pybel.BELGraph graph: A BEL graph
        """
        for chemical_node, data in graph.nodes(data=True):
            namespace = data.get(NAMESPACE)
            if namespace not in {'MESHC', 'MESH'}:
                continue

            identifier = data.get(IDENTIFIER)
            name = data.get(NAME)

            if identifier is not None:
                self.enrich_graph_chemical(graph, identifier)
            elif name is not None:
                self.enrich_graph_chemical(graph, name)
            else:
                raise KeyError

    def to_bel(self) -> BELGraph:
        """Convert all possible aspects of the database to BEL.

        .. warning:: Not complete!

        To do:

        - add namespaces
        - use cursors
        - multiprocessing
        """
        graph = BELGraph(name='CTD', version='1.0.0')

        mesh_manager = bio2bel_mesh.Manager(engine=self.engine, session=self.session)
        mesh_manager.add_namespace_to_graph(graph)

        for chem_gene_ixn in tqdm(self.list_chemical_gene_interactions(),
                                  total=self.count_chemical_gene_interactions()):
            add_chemical_gene_interaction(graph, chem_gene_ixn)

        return graph
