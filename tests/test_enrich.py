# -*- coding: utf-8 -*-

"""

Example page: http://ctdbase.org/detail.go?type=relationship&ixnId=3020460

"""

import unittest

from bio2bel_ctd import enrich_chemicals
from pybel import BELGraph
from pybel.constants import (
    CITATION, CITATION_REFERENCE, CITATION_TYPE, CITATION_TYPE_PUBMED, DECREASES, OBJECT,
    RELATION, SUBJECT,
)
from pybel.dsl import abundance, rna
from tests.constants import PopulatedDatabaseMixin

ex_mesh_name = abundance(namespace='MESH', name='Diethylnitrosamine')
ex_mesh_id = abundance(namespace='MESH', identifier='D004052')

ex_chebi_name = abundance(namespace='CHEBI', name='N-nitrosodiethylamine')
ex_chebi_id = abundance(namespace='CHEBI', identifier='34873')

abcc6 = rna(namespace='ENTREZ', name='ABCC6', identifier='368')
abcc6_tuple = abcc6.as_tuple()


class TestStuff(PopulatedDatabaseMixin):
    def help_test_graph(self, graph, c_tuple):
        """Diethylnitrosamine results in decreased expression of ABCC6 mRNA"""
        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        enrich_chemicals(graph, connection=self.connection)

        self.assertIn(abcc6_tuple, graph)
        # check that certain reactions are there?

        self.assertEqual(2, graph.number_of_nodes())

        self.assertIn(abcc6_tuple, graph.edge[c_tuple])

        self.assertEqual(1, graph.number_of_edges())

        key = list(graph.edge[c_tuple][abcc6_tuple])[0]
        data = graph.edge[c_tuple][abcc6_tuple][key]

        self.assertIn(RELATION, data)
        self.assertEqual(DECREASES, data[RELATION])

        self.assertNotIn(SUBJECT, data)
        self.assertNotIn(OBJECT, data)

        self.assertIn(CITATION, data)

        self.assertIn(CITATION_TYPE, data[CITATION])
        self.assertEqual(CITATION_TYPE_PUBMED, data[CITATION][CITATION_TYPE])

        self.assertIn(CITATION_REFERENCE, data[CITATION])
        self.assertEqual('19638242', data[CITATION][CITATION_REFERENCE])

        # TODO future PyBEL might have an evidence code system

    def test_enrich_mesh_name(self):
        graph = BELGraph()
        c_tuple = graph.add_node_from_data(ex_mesh_name)
        self.help_test_graph(graph, c_tuple)

    def test_enrich_mesh_id(self):
        graph = BELGraph()
        c_tuple = graph.add_node_from_data(ex_mesh_id)
        self.help_test_graph(graph, c_tuple)

    def test_enrich_chebi_id(self):
        graph = BELGraph()
        c_tuple = graph.add_node_from_data(ex_chebi_id)
        self.help_test_graph(graph, c_tuple)

    def test_enrich_chebi_name(self):
        graph = BELGraph()
        c_tuple = graph.add_node_from_data(ex_chebi_name)
        self.help_test_graph(graph, c_tuple)


if __name__ == '__main__':
    unittest.main()
