# -*- coding: utf-8 -*-

import logging
import os

from bio2bel.testing import AbstractTemporaryCacheClassMixin
from bio2bel_ctd import Manager

log = logging.getLogger(__name__)

dir_path = os.path.dirname(os.path.realpath(__file__))
resources_dir = os.path.join(dir_path, 'resources')

chemicals_url = os.path.join(resources_dir, 'CTD_chemicals.tsv.gz')
genes_url = os.path.join(resources_dir, 'CTD_genes.tsv.gz')
chemical_gene_interaction_types_url = os.path.join(resources_dir, 'CTD_chem_gene_ixn_types.tsv')
chemical_gene_interactions_url = os.path.join(resources_dir, 'CTD_chem_gene_ixns.tsv.gz')

_urls = [
    chemicals_url,
    genes_url,
    chemicals_url,
    genes_url
]

_only_tables = [
    'action',
    'chemical',
    'gene',
    'chem_gene_ixn',
]


class _TestManager(Manager):
    pyctd_data_dir = resources_dir


class TemporaryCacheClassMixin(AbstractTemporaryCacheClassMixin):
    Manager = _TestManager


class PopulatedDatabaseMixin(TemporaryCacheClassMixin):

    @classmethod
    def populate(cls):
        cls.manager.populate(urls=_urls, only_tables=_only_tables)
