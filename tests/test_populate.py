# -*- coding: utf-8 -*-

"""Tests the database gets populated"""

import logging

from tests.constants import PopulatedDatabaseMixin, TemporaryCacheClassMixin, resources_dir

log = logging.getLogger(__name__)


class TestConfig(TemporaryCacheClassMixin):
    def test_config(self):
        self.assertEqual(resources_dir, self.manager.pyctd_data_dir)


class TestImport(PopulatedDatabaseMixin):
    def test_count_genes(self):
        self.assertEqual(3, self.manager.count_genes())
