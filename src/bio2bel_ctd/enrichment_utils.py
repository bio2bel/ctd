# -*- coding: utf-8 -*-

import logging

from pybel.constants import DECREASES, INCREASES, REGULATES, VARIANTS
from pybel.dsl import (
    abundance as abundance_dsl, activity, complex_abundance as complex_abundance_dsl, fragment, gene as gene_dsl, pmod,
    protein as protein_dsl, reaction, rna as rna_dsl, translocation,
)

log = logging.getLogger('bio2bel_ctd')


def get_dsl_chemical(ixn):
    """Returns a PyBEL DSL object for the chemical from the interaction

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: pybel.dsl.abundance
    """
    return abundance_dsl(
        namespace='MESH',
        name=str(ixn.chemical.chemical_name),
        identifier=str(ixn.chemical.chemical_id)
    )


def get_dsl_gene(ixn):
    """Returns a PyBEL DSL object for the gene from the interaction

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: pybel.dsl.gene or pybel.dsl.rna or pybel.dsl.protein
    """
    kwargs = dict(
        namespace='ENTREZ',
        name=str(ixn.gene.gene_symbol),
        identifier=str(ixn.gene.gene_id)
    )

    forms = list(ixn.gene_forms)

    # do checking here too?
    form = forms[0].gene_form

    if form == 'gene':
        return gene_dsl(**kwargs)

    if form == 'mRNA':
        return rna_dsl(**kwargs)

    if form == 'protein':
        return protein_dsl(**kwargs)

    raise ValueError('unknown form: {}'.format(form))


def get_citation(ixn):
    """Extracts the first PubMed identifier for this interaction

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: str
    """
    pubmed_models = list(ixn.pubmed_ids)
    pubmed_model = pubmed_models[0]
    return str(pubmed_model.pubmed_id)


def _ixn_is_changes_entity(ixn, interaction_action, gene_form):
    """
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :param str interaction_action:
    :param str gene_form:
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


def _ixn_is_changes_mrna(ixn, interaction_action):
    """
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :param str interaction_action:
    :rtype: bool
    """
    return _ixn_is_changes_entity(ixn, interaction_action, 'mRNA')


def _ixn_is_changes_protein(ixn, interaction_action):
    """
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :param str interaction_action:
    :rtype: bool
    """
    return _ixn_is_changes_entity(ixn, interaction_action, 'protein')


def ixn_is_increases_mrna(ixn):
    """Checks if the interaction results in the increase of the expression of the mRNA of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_mrna(ixn, 'increases^expression')


def ixn_is_decreases_mrna(ixn):
    """Checks if the interaction results in the decrease of the expression of the mRNA of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_mrna(ixn, 'decreases^expression')


def ixn_is_regulates_mrna(ixn):
    """Checks if the interaction results in the regulation of the expression of the mRNA of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_mrna(ixn, 'affects^expression')


def ixn_is_increases_protein(ixn):
    """Checks if the interaction results in the increase of the expression of the protein of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_protein(ixn, 'increases^expression')


def ixn_is_decreases_protein(ixn):
    """Checks if the interaction results in the decrease of the expression of the protein of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_protein(ixn, 'decreases^expression')


def ixn_is_increases_activity(ixn):
    """Checks if the interaction results in the decrease of the activity of the protein of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_protein(ixn, 'increases^activity')


def ixn_is_decreases_activity(ixn):
    """Checks if the interaction results in the decrease of the activity of the protein of the gene

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool
    """
    return _ixn_is_changes_protein(ixn, 'decreases^activity')


def ixn_is_increases_phosphorylation(ixn):
    return _ixn_is_changes_protein(ixn, 'increases^phosphorylation')


def ixn_is_decreases_phosphorylation(ixn):
    return _ixn_is_changes_protein(ixn, 'decreases^phosphorylation')


def ixn_is_increases_hydroxylation(ixn):
    return _ixn_is_changes_protein(ixn, 'increases^hydroxylation')


def ixn_is_decreases_hydroxylation(ixn):
    return _ixn_is_changes_protein(ixn, 'decreases^hydroxylation')


def _add_ixn_causal_expression(graph, ixn, edge_type):
    """Adds an interaction to the graph after

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :param str edge_type: Either :data:`pybel.constants.INCREASES` or :data:`pybel.constants.DECREASES`
    :return: The hash of the added edge
    :rtype: str
    """
    gene = get_dsl_gene(ixn)
    chemical = get_dsl_chemical(ixn)
    citation = get_citation(ixn)

    return graph.add_qualified_edge(
        chemical,
        gene,
        edge_type,
        evidence=ixn.interaction,
        citation=citation,
        annotations={
            'Species': str(ixn.organism_id)
        }
    )


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


def _add_ixn_causal_actvity(graph, ixn, edge_type):
    """Adds an interaction to the graph after

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :param str edge_type: Either :data:`pybel.constants.INCREASES` or :data:`pybel.constants.DECREASES`
    :return: The hash of the added edge
    :rtype: str
    """
    gene = get_dsl_gene(ixn)
    chemical = get_dsl_chemical(ixn)
    citation = get_citation(ixn)

    return graph.add_qualified_edge(
        chemical,
        gene,
        edge_type,
        evidence=ixn.interaction,
        citation=citation,
        annotations={
            'Species': str(ixn.organism_id)
        },
        object_modifier=activity(),
    )


def add_ixn_increases_activity(graph, ixn):
    """Adds an interaction that represents the chemical increasing the activity of a protein

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_causal_actvity(graph, ixn, INCREASES)


def add_ixn_decreases_activity(graph, ixn):
    """Adds an interaction that represents the chemical decreasing the activity of a protein

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    return _add_ixn_causal_actvity(graph, ixn, DECREASES)


def _add_ixn_changes_pmod(graph, ixn, edge_type, pmod_name):
    chemical = get_dsl_chemical(ixn)

    protein = get_dsl_gene(ixn)
    protein[VARIANTS] = [pmod(name=pmod_name)]

    citation = get_citation(ixn)

    return graph.add_qualified_edge(
        chemical,
        protein,
        edge_type,
        evidence=ixn.interaction,
        citation=citation,
        annotations={
            'Species': str(ixn.organism_id)
        },
    )


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


def add_ixn_binding(graph, ixn):
    """Adds an interaction that represents the chemical binding to a protein

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    chemical = get_dsl_chemical(ixn)
    protein = get_dsl_gene(ixn)
    chemical_protein_complex = complex_abundance_dsl(members=[
        chemical,
        protein
    ])

    citation = get_citation(ixn)

    return graph.add_qualified_edge(
        chemical,
        chemical_protein_complex,
        INCREASES,
        evidence=ixn.interaction,
        citation=citation,
        annotations={
            'Species': str(ixn.organism_id)
        },
    )


def ixn_is_decreasing_expression_and_activity(ixn):
    """

    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :rtype: bool

    Chemical: 103D5R
    Gene: hypoxia inducible factor 1 alpha subunit
    Interaction: 103D5R results in decreased expression of and results in decreased activity of HIF1A protein
    Action: decreases^activity
    Action: decreases^expression
    Gene Form: protein
    """
    raise NotImplementedError  # TODO


def ixn_is_binding_and_increasing(ixn):
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


def add_ixn_affect_localization(graph, ixn):
    """

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    chemical = get_dsl_chemical(ixn)
    protein = get_dsl_gene(ixn)
    citation = get_citation(ixn)

    return graph.add_qualified_edge(
        chemical,
        protein,
        REGULATES,
        evidence=ixn.interaction,
        citation=citation,
        annotations={
            'Species': str(ixn.organism_id)
        },
        object_modifier=translocation({}, {})
    )


def add_ixn_increases_cleavage(graph, ixn):
    """

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    chemical = get_dsl_chemical(ixn)
    protein = get_dsl_gene(ixn)

    cleaved_protein = get_dsl_gene(ixn)
    cleaved_protein[VARIANTS] = [fragment()]

    citation = get_citation(ixn)

    cleavage_reaction = reaction(
        reactants=[protein],
        products=[cleaved_protein]
    )

    return graph.add_qualified_edge(
        chemical,
        cleavage_reaction,
        INCREASES,
        evidence=ixn.interaction,
        citation=citation,
        annotations={
            'Species': str(ixn.organism_id)
        },
    )


def add_ixn_increases_chemical_synthesis(graph, ixn):
    """

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    :return: The hash of the added edge
    :rtype: str
    """
    chemical = get_dsl_chemical(ixn)
    protein = get_dsl_gene(ixn)
    citation = get_citation(ixn)

    return graph.add_qualified_edge(
        protein,
        chemical,
        INCREASES,
        evidence=ixn.interaction,
        citation=citation,
        annotations={
            'Species': str(ixn.organism_id)
        },
    )


def add_chemical_gene_interaction(graph, ixn):
    """Adds a chemical-gene interaction to the BEL graph

    :param pybel.BELGraph graph: A BEL graph
    :param pyctd.manager.models.ChemGeneIxn ixn: A chemical-gene interaction
    """
    if ixn_is_increases_mrna(ixn) or ixn_is_increases_protein(ixn):
        return add_ixn_increases_expression(graph, ixn)

    if ixn_is_decreases_mrna(ixn) or ixn_is_decreases_protein(ixn):
        return add_ixn_decreases_expression(graph, ixn)

    if ixn_is_regulates_mrna(ixn):
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
