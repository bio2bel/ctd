# -*- coding: utf-8 -*-

import logging
from typing import Set

from pybel import BELGraph
from pybel.constants import DECREASES, INCREASES, REGULATES
from pybel.dsl import (
    CentralDogma, abundance as abundance_dsl, activity, complex_abundance as complex_abundance_dsl, fragment,
    gene as gene_dsl, gmod, pmod, protein as protein_dsl, reaction, rna as rna_dsl, translocation,
)
from .constants import MODULE_NAME
from .models import ChemGeneIxn

log = logging.getLogger('bio2bel_ctd')


def get_dsl_chemical(ixn: ChemGeneIxn) -> abundance_dsl:
    """Return a PyBEL DSL object for the chemical from the interaction."""
    return abundance_dsl(
        namespace='mesh',
        name=str(ixn.chemical.chemical_name),
        identifier=str(ixn.chemical.chemical_id)
    )


def get_dsl_gene(ixn: ChemGeneIxn) -> CentralDogma:
    """Return a PyBEL DSL object for the gene from the interaction."""
    gene_forms = list(ixn.gene_forms)

    # do checking here too?
    gene_form = gene_forms[0].gene_form

    if gene_form == 'gene':
        dsl = gene_dsl
    elif gene_form == 'mRNA':
        dsl = rna_dsl
    elif gene_form == 'protein':
        dsl = protein_dsl
    else:
        raise ValueError('unknown form: {}'.format(gene_form))

    return dsl(
        namespace='ncbigene',
        name=str(ixn.gene.gene_symbol),
        identifier=str(ixn.gene.gene_id)
    )


def _ixn_is_changes_entity(ixn: ChemGeneIxn, interaction_action: str, gene_form: str):
    """
    :param ixn: A chemical-gene interaction
    :param interaction_action:
    :param gene_form:
    :rtype: bool
    """
    actions = list(ixn.interaction_actions)
    forms = list(ixn.gene_forms)

    return (
            len(actions) == 1 and
            actions[0].interaction_action == interaction_action and
            len(forms) == 1 and
            forms[0].gene_form == gene_form
    )


def _ixn_is_changes_gene(ixn: ChemGeneIxn, interaction_action: str) -> bool:
    """
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :param str interaction_action:
    """
    return _ixn_is_changes_entity(ixn, interaction_action, 'gene')


def _ixn_is_changes_mrna(ixn: ChemGeneIxn, interaction_action: str) -> bool:
    """
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :param str interaction_action:
    :rtype: bool
    """
    return _ixn_is_changes_entity(ixn, interaction_action, 'mRNA')


def _ixn_is_changes_protein(ixn: ChemGeneIxn, interaction_action: str) -> bool:
    """
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :param str interaction_action:
    :rtype: bool
    """
    return _ixn_is_changes_entity(ixn, interaction_action, 'protein')


def ixn_is_increases_mrna(ixn: ChemGeneIxn) -> bool:
    """Checks if the interaction results in the increase of the expression of the mRNA of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_mrna(ixn, 'increases^expression')


def ixn_is_decreases_mrna(ixn: ChemGeneIxn):
    """Checks if the interaction results in the decrease of the expression of the mRNA of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_mrna(ixn, 'decreases^expression')


def ixn_is_regulates_mrna(ixn: ChemGeneIxn):
    """Checks if the interaction results in the regulation of the expression of the mRNA of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_mrna(ixn, 'affects^expression')


def ixn_is_increases_protein(ixn: ChemGeneIxn):
    """Checks if the interaction results in the increase of the expression of the protein of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_protein(ixn, 'increases^expression')


def ixn_is_decreases_protein(ixn: ChemGeneIxn):
    """Checks if the interaction results in the decrease of the expression of the protein of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_protein(ixn, 'decreases^expression')


def ixn_is_regulates_protein(ixn: ChemGeneIxn):
    """Checks if the interaction results in the decrease of the expression of the protein of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_protein(ixn, 'affects^expression')


def ixn_is_increases_activity(ixn: ChemGeneIxn):
    """Checks if the interaction results in the decrease of the activity of the protein of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_protein(ixn, 'increases^activity')


def ixn_is_decreases_activity(ixn: ChemGeneIxn):
    """Checks if the interaction results in the decrease of the activity of the protein of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_protein(ixn, 'decreases^activity')


def ixn_is_increases_phosphorylation(ixn: ChemGeneIxn):
    return _ixn_is_changes_protein(ixn, 'increases^phosphorylation')


def ixn_is_decreases_phosphorylation(ixn: ChemGeneIxn):
    return _ixn_is_changes_protein(ixn, 'decreases^phosphorylation')


def ixn_is_increases_hydroxylation(ixn: ChemGeneIxn):
    return _ixn_is_changes_protein(ixn, 'increases^hydroxylation')


def ixn_is_decreases_hydroxylation(ixn: ChemGeneIxn):
    return _ixn_is_changes_protein(ixn, 'decreases^hydroxylation')


def _add_ixn_causal_expression(graph: BELGraph, ixn: ChemGeneIxn, edge_type) -> Set[str]:
    """Adds an interaction to the graph after

    :param str edge_type: Either :data:`pybel.constants.INCREASES` or :data:`pybel.constants.DECREASES`
    :return: The hash of the added edge
    :rtype: str
    """
    gene = get_dsl_gene(ixn)
    chemical = get_dsl_chemical(ixn)

    return {
        graph.add_qualified_edge(
            chemical,
            gene,
            edge_type,
            evidence=ixn.interaction,
            citation=str(reference.pubmed_id),
            annotations={
                'Species': str(ixn.organism_id),
                'bio2bel': MODULE_NAME
            }
        )
        for reference in ixn.pubmed_ids
    }


def add_ixn_increases_expression(graph, ixn):
    """

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_causal_expression(graph, ixn, INCREASES)


def add_ixn_decreases_expression(graph, ixn):
    """

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_causal_expression(graph, ixn, DECREASES)


def add_ixn_regulates_expression(graph, ixn):
    """

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_causal_expression(graph, ixn, REGULATES)


def _add_ixn_causal_actvity(graph: BELGraph, ixn: ChemGeneIxn, edge_type) -> Set[str]:
    """Adds an interaction to the graph after

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :param str edge_type: Either :data:`pybel.constants.INCREASES` or :data:`pybel.constants.DECREASES`
    :return: The hash of the added edge
    :rtype: str
    """
    gene = get_dsl_gene(ixn)
    chemical = get_dsl_chemical(ixn)

    return {
        graph.add_qualified_edge(
            chemical,
            gene,
            edge_type,
            evidence=ixn.interaction,
            citation=str(reference.pubmed_id),
            annotations={
                'Species': str(ixn.organism_id),
                'bio2bel': MODULE_NAME
            },
            object_modifier=activity(),
        )
        for reference in ixn.pubmed_ids
    }


def add_ixn_increases_activity(graph: BELGraph, ixn: ChemGeneIxn):
    """Add an interaction that represents the chemical increasing the activity of a protein."""
    return _add_ixn_causal_actvity(graph, ixn, INCREASES)


def add_ixn_decreases_activity(graph: BELGraph, ixn: ChemGeneIxn):
    """Add an interaction that represents the chemical decreasing the activity of a protein."""
    return _add_ixn_causal_actvity(graph, ixn, DECREASES)


def _add_ixn_changes_pmod(graph: BELGraph, ixn: ChemGeneIxn, edge_type, pmod_name):
    chemical = get_dsl_chemical(ixn)
    protein = get_dsl_gene(ixn).with_variants(pmod(name=pmod_name))

    return {
        graph.add_qualified_edge(
            chemical,
            protein,
            edge_type,
            evidence=ixn.interaction,
            citation=str(reference.pubmed_id),
            annotations={
                'Species': str(ixn.organism_id)
            },
        )
        for reference in ixn.pubmed_ids
    }


def _add_ixn_changes_gmod(graph: BELGraph, ixn: ChemGeneIxn, edge_type: str, pmod_name: str) -> Set[str]:
    chemical = get_dsl_chemical(ixn)
    gene = get_dsl_gene(ixn).with_variants(gmod(name=pmod_name))

    return {
        graph.add_qualified_edge(
            chemical,
            gene,
            edge_type,
            evidence=ixn.interaction,
            citation=str(reference.pubmed_id),
            annotations={
                'Species': str(ixn.organism_id)
            },
        )
        for reference in ixn.pubmed_ids
    }


def add_ixn_increases_methylation(graph, ixn):
    """Adds an interaction that represents the chemical increasing the methylation of a gene

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_changes_gmod(graph, ixn, INCREASES, 'Me')


def add_ixn_decreases_methylation(graph, ixn):
    """Adds an interaction that represents the chemical decreasing the methylation of a gene

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_changes_gmod(graph, ixn, DECREASES, 'Me')


def add_ixn_regulates_methylation(graph, ixn):
    """Adds an interaction that represents the chemical regulation of the methylation of a gene

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_changes_gmod(graph, ixn, REGULATES, 'Me')


def add_ixn_increases_phosphorylation(graph, ixn):
    """Adds an interaction that represents the chemical increasing the phosphorylation of a protein

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_changes_pmod(graph, ixn, INCREASES, 'Ph')


def add_ixn_decreases_phosphorylation(graph, ixn):
    """Adds an interaction that represents the chemical decreasing the phosphorylation of a protein

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_changes_pmod(graph, ixn, DECREASES, 'Ph')


def add_ixn_increases_hydroxylation(graph, ixn):
    """Adds an interaction that represents the chemical increasing the hydroxylation of a protein

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_changes_pmod(graph, ixn, INCREASES, 'Hy')


def add_ixn_decreases_hydroxylation(graph, ixn):
    """Adds an interaction that represents the chemical decreasing the hydroxylation of a protein

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_changes_pmod(graph, ixn, DECREASES, 'Hy')


def add_ixn_increases_oxidation(graph, ixn):
    """Adds an interaction that represents the chemical increasing the oxidation of a protein

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_changes_pmod(graph, ixn, INCREASES, 'Ox')


def add_ixn_decreases_oxidation(graph, ixn):
    """Adds an interaction that represents the chemical decreasing the hydroxylation of a protein

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_changes_pmod(graph, ixn, DECREASES, 'Ox')


def ixn_is_binding(ixn):
    """Checks if the interaction represents the binding of a chemical to a protein

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool

    Example from CTD: `10-(fluoroethoxyphosphinyl)-N-(biotinamidopentyl)decanamide binds to ACTA1 protein
    <http://ctdbase.org/detail.go?type=relationship&ixnId=2714292>`_
    """
    return _ixn_is_changes_protein(ixn, 'affects^binding')


def add_ixn_binding(graph: BELGraph, ixn: ChemGeneIxn) -> Set[str]:
    """Add an interaction that represents the chemical binding to a protein.

    :param graph: A BEL graph
    :param ixn: A chemical-gene interaction
    :return: The hash of the added edge
    """
    chemical = get_dsl_chemical(ixn)
    protein = get_dsl_gene(ixn)
    chemical_protein_complex = complex_abundance_dsl(members=[
        chemical,
        protein
    ])

    return {
        graph.add_increases(
            chemical,
            chemical_protein_complex,
            evidence=ixn.interaction,
            citation=str(reference.pubmed_id),
            annotations={
                'Species': str(ixn.organism_id)
            },
        )
        for reference in ixn.pubmed_ids
    }


def ixn_is_decreasing_expression_and_activity(ixn: ChemGeneIxn):
    """

    :param ixn: A chemical-gene interaction
    :rtype: bool

    Chemical: 103D5R
    Gene: hypoxia inducible factor 1 alpha subunit
    Interaction: 103D5R results in decreased expression of and results in decreased activity of HIF1A protein
    Action: decreases^activity
    Action: decreases^expression
    Gene Form: protein
    """
    raise NotImplementedError  # TODO


def ixn_is_binding_and_increasing(ixn: ChemGeneIxn):
    """

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool

    Chemical: 10-hydroxywarfarin
    Gene: nuclear receptor subfamily 1 group I member 2
    Interaction: 10-hydroxywarfarin binds to and results in increased activity of NR1I2 protein
    Action: affects^binding
    Action: increases^activity
    Gene Form: protein
    """
    raise NotImplementedError  # TODO


def add_ixn_affect_localization(graph: BELGraph, ixn: ChemGeneIxn) -> Set[str]:
    """"""
    chemical = get_dsl_chemical(ixn)
    protein = get_dsl_gene(ixn)

    return {
        graph.add_qualified_edge(
            chemical,
            protein,
            REGULATES,
            evidence=ixn.interaction,
            citation=str(reference.pubmed_id),
            annotations={
                'Species': str(ixn.organism_id),
                'bio2bel': MODULE_NAME
            },
            object_modifier=translocation({}, {})
        )
        for reference in ixn.pubmed_ids
    }


def add_ixn_increases_cleavage(graph: BELGraph, ixn: ChemGeneIxn) -> Set[str]:
    """"""
    chemical = get_dsl_chemical(ixn)
    protein = get_dsl_gene(ixn)
    cleaved_protein = protein.with_variants(fragment())

    cleavage_reaction = reaction(
        reactants=[protein],
        products=[cleaved_protein]
    )
    return {
        graph.add_increases(
            chemical,
            cleavage_reaction,
            evidence=ixn.interaction,
            citation=str(reference.pubmed_id),
            annotations={
                'Species': str(ixn.organism_id),
                'bio2bel': MODULE_NAME
            },
        )
        for reference in ixn.pubmed_ids
    }


def add_ixn_increases_chemical_synthesis(graph: BELGraph, ixn: ChemGeneIxn):
    """"""
    chemical = get_dsl_chemical(ixn)
    protein = get_dsl_gene(ixn)

    return [
        graph.add_qualified_edge(
            protein,
            chemical,
            INCREASES,
            evidence=ixn.interaction,
            citation=pubmed.pubmed_id,
            annotations={
                'Species': str(ixn.organism_id)
            },
        )
        for pubmed in ixn.pubmed_ids
    ]


def add_chemical_gene_interaction(graph, ixn: ChemGeneIxn):
    """Adds a chemical-gene interaction to the BEL graph

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    """
    if ixn_is_increases_mrna(ixn) or ixn_is_increases_protein(ixn):
        return add_ixn_increases_expression(graph, ixn)

    if ixn_is_decreases_mrna(ixn) or ixn_is_decreases_protein(ixn):
        return add_ixn_decreases_expression(graph, ixn)

    if ixn_is_regulates_mrna(ixn) or ixn_is_regulates_protein(ixn):
        return add_ixn_regulates_expression(graph, ixn)

    if ixn_is_increases_activity(ixn):
        return add_ixn_increases_activity(graph, ixn)

    if ixn_is_decreases_activity(ixn):
        return add_ixn_decreases_activity(graph, ixn)

    if ixn_is_increases_phosphorylation(ixn):
        return add_ixn_increases_phosphorylation(graph, ixn)

    if ixn_is_decreases_phosphorylation(ixn):
        return add_ixn_decreases_phosphorylation(graph, ixn)

    if ixn_is_increases_hydroxylation(ixn):
        return add_ixn_increases_hydroxylation(graph, ixn)

    if ixn_is_decreases_hydroxylation(ixn):
        return add_ixn_decreases_hydroxylation(graph, ixn)

    if _ixn_is_changes_protein(ixn, 'increases^oxidation'):
        return add_ixn_increases_oxidation(graph, ixn)

    if _ixn_is_changes_protein(ixn, 'decreases^oxidation'):
        return add_ixn_decreases_oxidation(graph, ixn)

    if _ixn_is_changes_gene(ixn, 'increases^methylation'):
        return add_ixn_increases_methylation(graph, ixn)

    if _ixn_is_changes_gene(ixn, 'decreases^methylation'):
        return add_ixn_decreases_methylation(graph, ixn)

    if _ixn_is_changes_gene(ixn, 'affects^methylation'):
        return add_ixn_regulates_methylation(graph, ixn)

    if ixn_is_binding(ixn):
        return add_ixn_binding(graph, ixn)

    if _ixn_is_changes_protein(ixn, 'affects^localization'):
        return add_ixn_affect_localization(graph, ixn)

    if _ixn_is_changes_protein(ixn, 'increases^cleavage'):
        return add_ixn_increases_cleavage(graph, ixn)

    if _ixn_is_changes_protein(ixn, 'increases^chemical synthesis'):
        return add_ixn_increases_chemical_synthesis(graph, ixn)

    if len(list(ixn.interaction_actions)) > 1:
        return

    log.debug('did not map (%d) %s', ixn.id, ixn.interaction)


def test(size=None, limit=1000):
    from bio2bel_ctd import Manager
    from bio2bel_ctd.models import ChemGeneIxn
    from pybel import BELGraph

    m = Manager()
    g = BELGraph(name='test', version='0.0.0')
    c = 0
    for ixn in m.session.query(ChemGeneIxn).limit(limit):
        res = add_chemical_gene_interaction(g, ixn)

        if res:
            c += 1

        if size is not None and c > size:
            break

    return g
