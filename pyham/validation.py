"""Check a species tree and an OrthoXML file against each other before analysing them with pyHam.

Most problems users hit when loading data come from a species tree that doesn't match the OrthoXML
file: species named differently, levels missing from the tree, or a group nesting that contradicts
the tree. :func:`validate_input` checks for all of these in one go and returns a
:class:`ValidationReport` listing every problem found, instead of failing at the first one.

The checks run in four stages:

1. **Species tree**: can it be read, does every node have a name, are the names unique?
2. **OrthoXML**: duplicate species or gene ids, references to undefined genes, genes used twice.
3. **Tree vs. OrthoXML**: is every OrthoXML species a leaf of the tree, and is every level the
   OrthoXML claims for its HOGs (``TaxRange``/``taxid``) a node of the tree?
4. **Trial build**: if nothing above is an error, every family is built exactly as :obj:`pyham.Ham`
   would, but a family that fails is recorded and skipped instead of stopping the whole parse.

The same checks are available from the command line::

    pyham-validate species_tree.nwk hogs.orthoxml
"""

import argparse
import collections
import difflib
import gzip
import io
import os
import re
import sys
import textwrap
import types
from dataclasses import dataclass, field
from typing import List, Optional
from xml.etree.ElementTree import XMLParser, iterparse, ParseError

__all__ = ["validate_input", "ValidationReport", "ValidationIssue", "InputValidationError"]

ERROR = "error"
WARNING = "warning"
NOTE = "note"
_SEVERITY_ORDER = {ERROR: 0, WARNING: 1, NOTE: 2}

_NS = "{http://orthoXML.org/2011/}"
_PHYLOXML_TAG_ATTR = {"clade_name": "name", "taxonomy_scientific_name": "scientific_name", "taxonomy_code": "code"}
_SPECIES_CODE_RE = re.compile(r"^[A-Z][A-Z0-9]{2,4}$")
_MAX_FUZZY_LOOKUPS = 50


@dataclass
class ValidationIssue:
    """One problem, or noteworthy fact, found in the inputs.

    Attributes:
        | severity (:obj:`str`): ``'error'`` (pyHam will fail, or silently build something wrong),
          ``'warning'`` (pyHam works, but probably not as intended) or ``'note'`` (worth knowing).
        | code (:obj:`str`): stable identifier of the kind of issue, e.g. ``'species-not-in-tree'``.
        | summary (:obj:`str`): one-sentence description.
        | items (:obj:`list` of :obj:`str`): every affected entity, one description each. The printed
          report shows only the first few; this list is always complete.
        | hint (:obj:`str`): suggested fix, when one can be inferred from the data.
        | data (:obj:`dict`): machine-readable details, specific to each code.
    """
    severity: str
    code: str
    summary: str
    items: List[str] = field(default_factory=list)
    hint: Optional[str] = None
    data: dict = field(default_factory=dict)


class ValidationReport:
    """Result of :func:`validate_input`.

    Print it for a human-readable report, or inspect :attr:`issues` (or :attr:`errors`,
    :attr:`warnings`, :attr:`notes`) programmatically.

    Attributes:
        | inputs (:obj:`dict`): short description of each input that was checked.
        | issues (:obj:`list` of :obj:`ValidationIssue`): everything found, in the order it was found.
        | stats (:obj:`dict`): counts gathered along the way (species, genes, families, ...).
        | trial_build (:obj:`str`): outcome of the trial build, or why it was skipped. None if not applicable.
    """

    def __init__(self):
        self.inputs = {}
        self.issues = []
        self.stats = {}
        self.trial_build = "not run"

    def add(self, severity, code, summary, items=None, hint=None, data=None):
        issue = ValidationIssue(severity, code, summary, list(items or []), hint, dict(data or {}))
        self.issues.append(issue)
        return issue

    @property
    def errors(self):
        return [i for i in self.issues if i.severity == ERROR]

    @property
    def warnings(self):
        return [i for i in self.issues if i.severity == WARNING]

    @property
    def notes(self):
        return [i for i in self.issues if i.severity == NOTE]

    @property
    def ok(self):
        """True if no errors were found (warnings and notes are allowed)."""
        return not self.errors

    def raise_for_errors(self):
        """Raise :obj:`InputValidationError` if the report contains any error."""
        if self.errors:
            raise InputValidationError(self)

    def format(self, max_items=10):
        """Render the report as text, showing at most `max_items` affected entities per issue."""
        counts = [_plural(len(self.errors), "error"), _plural(len(self.warnings), "warning"),
                  _plural(len(self.notes), "note")]
        lines = ["pyHam input check: " + ", ".join(counts)]
        width = max((len(k) for k in self.inputs), default=0)
        for label, desc in self.inputs.items():
            lines.append(f"  {label + ':':<{width + 1}} {desc}")

        for issue in sorted(self.issues, key=lambda i: _SEVERITY_ORDER[i.severity]):
            lines.append("")
            lines.append(f"{issue.severity.upper()} [{issue.code}] {issue.summary}")
            shown = issue.items if max_items is None else issue.items[:max_items]
            for item in shown:
                first, *rest = item.split("\n")
                lines.append(f"    - {first}")
                if rest:
                    lines.extend(f"      {line}" for line in textwrap.dedent("\n".join(rest)).split("\n"))
            if len(shown) < len(issue.items):
                lines.append(f"    ... and {len(issue.items) - len(shown)} more")
            if issue.hint:
                lines.append(f"  hint: {issue.hint}")

        if not self.issues:
            lines += ["", "No problems found."]
        if self.trial_build is not None:
            lines += ["", f"Trial build: {self.trial_build}"]
        return "\n".join(lines)

    def __str__(self):
        return self.format()

    def __repr__(self):
        return (f"<ValidationReport: {len(self.errors)} errors, {len(self.warnings)} warnings, "
                f"{len(self.notes)} notes>")


class InputValidationError(KeyError):
    """Raised when the inputs given to pyHam contain errors. The full :obj:`ValidationReport` is
    available as the `report` attribute.

    Subclasses :obj:`KeyError` because that is what pyHam raised for these problems before.
    """

    def __init__(self, report, footer=None):
        message = report.format()
        if footer:
            message += "\n\n" + footer
        super().__init__(message)
        self.report = report

    def __str__(self):
        return self.args[0]


_IRREGULAR_PLURALS = {"leaf": "leaves", "tree leaf": "tree leaves", "family": "families", "species": "species"}


def _plural(n, word):
    return f"{n:,} {word if n == 1 else _IRREGULAR_PLURALS.get(word, word + 's')}"


def _counted(n, word, verb_singular, verb_plural):
    """e.g. _counted(2, 'leaf', 'has', 'have') -> '2 leaves have'"""
    return f"{_plural(n, word)} {verb_singular if n == 1 else verb_plural}"


def _src_label(src, is_string):
    if is_string:
        return "<string>"
    return str(src)


# ----------------------------------------------------------------------------------------------
# Names in the species tree
# ----------------------------------------------------------------------------------------------

@dataclass
class _TreeNames:
    """Names available in a species tree, under the chosen settings and under the alternatives."""
    leaves: list
    internal: list
    # internal node labels as written in the file (newick only), even when use_internal_name=False
    raw_internal_labels: Optional[set] = None
    # phyloxml only: {tag: set of values} for every name tag, to suggest a better choice
    alt_leaf_names: dict = field(default_factory=dict)
    alt_internal_names: dict = field(default_factory=dict)
    leaf_tag: Optional[str] = None
    internal_tag: Optional[str] = None

    @classmethod
    def from_taxonomy(cls, taxonomy):
        leaves, internal = [], []
        alt_leaf = collections.defaultdict(set)
        alt_internal = collections.defaultdict(set)
        for node in taxonomy.tree.traverse():
            (leaves if node.is_leaf else internal).append(node.name)
            for tag, attr in _PHYLOXML_TAG_ATTR.items():
                value = node.props.get(attr) if attr != "name" else None
                if value:
                    (alt_leaf if node.is_leaf else alt_internal)[tag].add(value)
        return cls(leaves, internal, alt_leaf_names=dict(alt_leaf), alt_internal_names=dict(alt_internal))


def _normalize_name(name):
    return " ".join(name.strip().strip("'\"").replace("_", " ").lower().split())


def _describe_name_difference(a, b):
    if a.replace("_", " ") == b.replace("_", " "):
        return "underscores vs spaces"
    if a.lower() == b.lower():
        return "letter case"
    if a.strip().strip("'\"").strip() == b.strip().strip("'\"").strip():
        return "quotes or surrounding whitespace"
    return "underscores/spaces, letter case or quotes"


def _species_issue(report, unmapped, n_species, tree_names, all_species=None):
    """Add a 'species-not-in-tree' error for the (name, reason) pairs in `unmapped`."""
    leaves = [n for n in tree_names.leaves if n]
    internal = {n for n in tree_names.internal if n}
    by_norm = collections.defaultdict(list)
    for leaf in leaves:
        by_norm[_normalize_name(leaf)].append(leaf)

    items, suggestions = [], {}
    differences = collections.Counter()
    unexplained = fuzzy_lookups = 0
    for name, _reason in unmapped:
        if name in internal:
            detail = "is an internal node of the tree, not a leaf"
            unexplained += 1
        elif by_norm.get(_normalize_name(name)):
            match = by_norm[_normalize_name(name)][0]
            difference = _describe_name_difference(name, match)
            differences[difference] += 1
            suggestions[name] = match
            detail = f'the tree has "{match}" ({difference})'
        else:
            unexplained += 1
            detail = None
            if fuzzy_lookups < _MAX_FUZZY_LOOKUPS:
                fuzzy_lookups += 1
                close = difflib.get_close_matches(name, leaves, n=1, cutoff=0.8)
                if close:
                    suggestions[name] = close[0]
                    detail = f'closest tree leaf: "{close[0]}"'
                else:
                    detail = "no similar leaf in the tree"
        items.append(f'"{name}": {detail}' if detail else f'"{name}"')

    n_matched = n_species - len(unmapped)
    hints = []
    if unmapped and not unexplained:
        kinds = " and ".join(differences)
        hints.append(f"All of these differ from a tree leaf only in {kinds}; make the names identical in "
                     f"both files.")

    unmapped_names = [name for name, _ in unmapped]
    if unexplained and leaves:
        xml_codes = sum(bool(_SPECIES_CODE_RE.match(n)) for n in unmapped_names) / len(unmapped_names)
        tree_codes = sum(bool(_SPECIES_CODE_RE.match(n)) for n in leaves) / len(leaves)
        if xml_codes >= 0.8 and tree_codes <= 0.2:
            hints.append("The OrthoXML names species by short codes (like \"HUMAN\"), but the tree leaves are not "
                         "codes.")
        elif xml_codes <= 0.2 and tree_codes >= 0.8:
            hints.append("The tree leaves are short species codes (like \"HUMAN\"), but the OrthoXML uses other "
                         "names.")

    if all_species is not None:
        species_set = set(all_species)
        for tag, values in tree_names.alt_leaf_names.items():
            n_alt = len(species_set & values)
            if tag != tree_names.leaf_tag and n_alt > n_matched:
                hints.append(f"phyloxml_leaf_name_tag='{tag}' would match {n_alt:,} of {n_species:,} species "
                             f"(currently {n_matched:,}).")

    if not hints:
        if n_matched == 0:
            hints.append("None of the OrthoXML species are in the tree. Is this the species tree the HOGs were "
                         "inferred with?")
        else:
            hints.append("Every OrthoXML species must be a leaf of the species tree, with exactly the same name. "
                         "Use the species tree the HOGs were inferred with.")

    verb = "is not a leaf" if len(unmapped) == 1 else "are not leaves"
    report.add(ERROR, "species-not-in-tree",
               f"{len(unmapped):,} of {n_species:,} OrthoXML species {verb} of the species tree.",
               items, " ".join(hints), {"species": unmapped_names, "suggestions": suggestions})


def species_not_in_tree_error(unmapped, n_species, taxonomy):
    """Build the :obj:`InputValidationError` that :obj:`pyham.Ham` raises when OrthoXML species are
    missing from the species tree."""
    report = ValidationReport()
    report.trial_build = None
    _species_issue(report, unmapped, n_species, _TreeNames.from_taxonomy(taxonomy))
    return InputValidationError(report, footer="Run pyham.validate_input() with the same inputs to check "
                                               "everything else as well.")


# ----------------------------------------------------------------------------------------------
# Stage 1: species tree
# ----------------------------------------------------------------------------------------------

def _detect_tree_format(tree, tree_format):
    if tree_format is not None:
        return tree_format
    if isinstance(tree, (str, os.PathLike)) and os.path.isfile(tree):
        name = str(tree).lower()
        if name.endswith(".gz"):
            name = name[:-3]
        return "phyloxml" if name.endswith((".xml", ".phyloxml")) else "newick"
    return "newick_string"


def _read_newick(tree, tree_format, report):
    from .taxonomy import create_tree_from_newick
    try:
        # use_internal_name=True keeps the labels as written; they're replaced later if needed
        return create_tree_from_newick(tree, tree_format=tree_format, use_internal_name=True)
    except Exception as e:
        report.add(ERROR, "tree-unreadable", f"The species tree can't be read as newick: {e}",
                   hint="Names that contain spaces, quotes, parentheses, commas, colons or semicolons must be "
                        "quoted in newick, e.g. 'Aster yellows witches''-broom phytoplasma' (an apostrophe "
                        "inside single quotes is written twice).")
        return None


def _read_phyloxml(tree, report):
    from .parsers import PhyloXMLToETE
    factory = PhyloXMLToETE()
    parser = XMLParser(target=factory)
    open_ = gzip.open if str(tree).endswith(".gz") else open
    try:
        with open_(tree, "rt") as fh:
            for line in fh:
                parser.feed(line)
        parser.close()
    except Exception as e:
        report.add(ERROR, "tree-unreadable", f"The species tree can't be read as phyloXML: {e}")
        return None
    if factory.root is None:
        report.add(ERROR, "tree-unreadable", "The phyloXML file contains no <clade> elements.")
        return None
    return factory.root


def _name_field(node, tag):
    """Value of a phyloxml name field; the clade name is kept aside as node.name gets overwritten."""
    return node.props.get("_clade_name") if tag == "clade_name" else node.props.get(_PHYLOXML_TAG_ATTR[tag])


def _node_label(node):
    """The node's name or, failing that, any other name field it has."""
    if node.name:
        return node.name
    for tag in _PHYLOXML_TAG_ATTR:
        if _name_field(node, tag):
            return f"{tag}={_name_field(node, tag)!r}"
    return "?"


def _describe_node(node, max_names=3):
    names = [_node_label(leaf) for leaf in node.leaves()]
    shown = ", ".join(names[:max_names]) + (", ..." if len(names) > max_names else "")
    return f"{len(names)} leaves: {shown}" if len(names) > 1 else f"leaf {shown}"


def _describe_unnamed_leaf(node, alternatives):
    values = [f"{tag}={_name_field(node, tag)!r}" for tag in alternatives if _name_field(node, tag)]
    if values:
        return "leaf with " + ", ".join(values)
    return f"unnamed leaf next to {_describe_node(node.up)}" if node.up else "unnamed leaf"


def _check_tree(tree, tree_format, use_internal_name, leaf_tag, internal_tag, report):
    """Stage 1. Returns (Taxonomy or None, _TreeNames or None)."""
    from .taxonomy import Taxonomy

    is_phyloxml = tree_format == "phyloxml"
    root = _read_phyloxml(tree, report) if is_phyloxml else _read_newick(tree, tree_format, report)
    if root is None:
        return None, None

    names = _TreeNames([], [], leaf_tag=leaf_tag if is_phyloxml else None,
                       internal_tag=internal_tag if is_phyloxml else None)
    if is_phyloxml:
        alt_leaf, alt_internal = collections.defaultdict(set), collections.defaultdict(set)
        for node in root.traverse():
            if node.name:
                node.add_prop("_clade_name", node.name)
            for tag, attr in _PHYLOXML_TAG_ATTR.items():
                value = node.props.get(attr)
                if value:
                    (alt_leaf if node.is_leaf else alt_internal)[tag].add(value)
        names.alt_leaf_names, names.alt_internal_names = dict(alt_leaf), dict(alt_internal)
        leaf_attr, internal_attr = _PHYLOXML_TAG_ATTR[leaf_tag], _PHYLOXML_TAG_ATTR[internal_tag]
        for node in root.traverse("postorder"):
            if node.is_leaf:
                node.name = node.props.get(leaf_attr) or None
            elif use_internal_name:
                node.name = node.props.get(internal_attr) or None
    else:
        names.raw_internal_labels = {n.name for n in root.traverse() if not n.is_leaf and n.name}

    if not use_internal_name:
        for node in root.traverse("postorder"):
            if not node.is_leaf:
                node.name = "/".join(child.name or "?" for child in node.get_children())

    leaves = list(root.leaves())
    internal = [n for n in root.traverse() if not n.is_leaf]
    names.leaves = [n.name for n in leaves]
    names.internal = [n.name for n in internal]
    report.stats.update(tree_leaves=len(leaves), tree_internal_nodes=len(internal))

    def tag_hint(alternatives, current, n_nodes, option, noun):
        better = [tag for tag, values in alternatives.items() if tag != current and len(values) >= n_nodes]
        return f"Every {noun} has a value for {option}='{better[0]}'." if better else ""

    unnamed_leaves = [n for n in leaves if not n.name]
    if unnamed_leaves:
        hint = None
        if is_phyloxml:
            hint = (f"pyHam names leaves from phyloxml_leaf_name_tag='{leaf_tag}'. " +
                    tag_hint(names.alt_leaf_names, leaf_tag, len(leaves), "phyloxml_leaf_name_tag", "leaf")).strip()
        report.add(ERROR, "unnamed-leaves", f"{_counted(len(unnamed_leaves), 'leaf', 'has', 'have')} no name in "
                   f"the species tree.",
                   [_describe_unnamed_leaf(n, names.alt_leaf_names) for n in unnamed_leaves], hint)

    if use_internal_name:
        # newick roots may be unnamed; phyloxml roots need a name like every other node
        unnamed_internal = [n for n in internal if not n.name and (is_phyloxml or not n.is_root)]
        if unnamed_internal:
            hint = "Set use_internal_name=False to let pyHam name internal nodes after their children."
            if is_phyloxml:
                hint = (tag_hint(names.alt_internal_names, internal_tag, len(internal),
                                 "phyloxml_internal_name_tag", "internal node") + " " + hint).strip()
            report.add(ERROR, "unnamed-internal-nodes",
                       f"{_counted(len(unnamed_internal), 'internal node', 'has', 'have')} no name in the species "
                       f"tree, but use_internal_name=True needs one for every internal node.",
                       [f"node above {_describe_node(n)}" for n in unnamed_internal], hint)

    leaf_counts = collections.Counter(n for n in names.leaves if n)
    dup_leaves = {n: c for n, c in leaf_counts.items() if c > 1}
    if dup_leaves:
        report.add(ERROR, "duplicate-leaf-names", f"{_counted(len(dup_leaves), 'leaf name', 'occurs', 'occur')} "
                   f"more than once in the species tree.", [f'"{n}" ({c} times)' for n, c in dup_leaves.items()],
                   "Each species must appear exactly once in the species tree.")

    if use_internal_name:
        internal_counts = collections.Counter(n for n in names.internal if n)
        dup_internal = {n: c for n, c in internal_counts.items() if c > 1}
        if dup_internal:
            numeric = sum(_is_number(n) for n in internal_counts) / max(len(internal_counts), 1)
            if numeric >= 0.5:
                hint = ("Most internal node labels are numbers, probably support values rather than taxon "
                        "names. Set use_internal_name=False to let pyHam name internal nodes after their children.")
            else:
                hint = ("pyHam identifies ancestral levels by name, so internal names must be unique. Rename "
                        "them, or set use_internal_name=False to name internal nodes after their children.")
            report.add(ERROR, "duplicate-internal-names",
                       f"{_counted(len(dup_internal), 'internal node name', 'occurs', 'occur')} more than once in the "
                       f"species tree.",
                       [f'"{n}" ({c} times)' for n, c in dup_internal.items()], hint)

        shared = sorted(set(leaf_counts) & {n for n in names.internal if n})
        if shared:
            report.add(WARNING, "leaf-and-internal-name", f"{_counted(len(shared), 'name', 'is', 'are')} used for both "
                       f"a leaf and an internal node of the species tree, so looking up a level by name is ambiguous.",
                       [f'"{n}"' for n in shared])

    if report.errors:
        return None, names
    try:
        return Taxonomy.from_ete_tree(root), names
    except Exception as e:
        report.add(ERROR, "tree-invalid", f"pyHam can't use this species tree: {e}")
        return None, names


def _is_number(text):
    try:
        float(text)
        return True
    except (TypeError, ValueError):
        return False


# ----------------------------------------------------------------------------------------------
# Stages 2 and 3: OrthoXML content, and OrthoXML vs. tree
# ----------------------------------------------------------------------------------------------

class _LevelResolver:
    """Decides whether a HOG's claimed level is a node of the tree, the same way the parser does."""

    def __init__(self, taxonomy, tree_names):
        self.taxonomy = taxonomy
        self.names = {n for n in tree_names.leaves + tree_names.internal if n} if tree_names else set()
        self._cache = {}

    def resolves(self, key, value):
        if (key, value) not in self._cache:
            self._cache[(key, value)] = self._resolves(key, value)
        return self._cache[(key, value)]

    def _resolves(self, key, value):
        if self.taxonomy is None:
            return key != "taxonId" and value in self.names
        if key == "taxonId":
            try:
                return len(list(self.taxonomy.nodes_by_attr(id=int(value)))) == 1
            except ValueError:
                return False
        try:
            self.taxonomy.get_node_by_name(value)
            return True
        except KeyError:
            return False


def _scan_orthoxml(open_orthoxml, resolver, report):
    """Stage 2 and the level part of stage 3, in one streaming pass. Returns the species names in
    file order and the unresolved levels ({label: [count, first family]}), or (None, None) if the
    file can't be parsed."""
    species, gene_species = [], {}
    duplicate_genes = collections.defaultdict(list)
    references = collections.defaultdict(list)
    unknown_refs = []
    unresolved_levels = collections.OrderedDict()
    n_families = n_groups = n_claimed = 0
    current_species = family = None
    frames = []

    try:
        with open_orthoxml() as fh:
            for event, elem in iterparse(fh, events=("start", "end")):
                tag = elem.tag
                if event == "start":
                    if tag == _NS + "species":
                        current_species = elem.get("name")
                        species.append(current_species)
                    elif tag == _NS + "gene":
                        gene_id = elem.get("id")
                        if gene_id in gene_species:
                            duplicate_genes[gene_id].append(current_species)
                        else:
                            gene_species[gene_id] = current_species
                    elif tag in (_NS + "orthologGroup", _NS + "paralogGroup"):
                        if not frames:
                            n_families += 1
                            family = elem.get("id") or f"#{n_families} (without id)"
                        claims = {}
                        if tag == _NS + "orthologGroup":
                            n_groups += 1
                            if elem.get("taxonId") is not None:
                                claims["taxonId"] = elem.get("taxonId")
                        frames.append((tag, claims))
                    elif tag == _NS + "property" and frames and frames[-1][0] == _NS + "orthologGroup":
                        if elem.get("name") in ("TaxRange", "taxid"):
                            frames[-1][1].setdefault(elem.get("name"), elem.get("value"))
                    elif tag == _NS + "geneRef" and frames:
                        gene_id = elem.get("id")
                        if gene_id not in gene_species:
                            unknown_refs.append((gene_id, family))
                        else:
                            references[gene_id].append(family)
                else:
                    if tag in (_NS + "orthologGroup", _NS + "paralogGroup"):
                        group_tag, claims = frames.pop()
                        if group_tag == _NS + "orthologGroup" and claims:
                            n_claimed += 1
                            if not any(resolver.resolves(k, v) for k, v in claims.items()):
                                label = claims.get("TaxRange") or (
                                    f"taxid={claims['taxid']}" if "taxid" in claims else f"taxonId={claims['taxonId']}")
                                entry = unresolved_levels.setdefault(label, [0, family])
                                entry[0] += 1
                    elem.clear()
    except (ParseError, OSError) as e:
        report.add(ERROR, "orthoxml-unreadable", f"The OrthoXML file can't be read: {e}")
        return None, None

    report.stats.update(species=len(species), genes=len(gene_species), families=n_families, hogs=n_groups,
                        hogs_with_level=n_claimed)

    dup_species = {n: c for n, c in collections.Counter(species).items() if c > 1}
    if dup_species:
        report.add(ERROR, "duplicate-species", f"{_counted(len(dup_species), 'species', 'is', 'are')} defined more "
                   f"than once in the OrthoXML.", [f'"{n}" ({c} times)' for n, c in dup_species.items()],
                   "Merge the genes of each species into a single <species> element.")
    if duplicate_genes:
        report.add(ERROR, "duplicate-gene-ids", f"{_counted(len(duplicate_genes), 'gene id', 'is', 'are')} defined "
                   f"more than once in the OrthoXML.",
                   [f'gene id {g} (in {gene_species[g]}, and again in {", ".join(map(str, sp))})'
                    for g, sp in duplicate_genes.items()],
                   "Gene ids must be unique across the whole file.")
    if unknown_refs:
        report.add(ERROR, "unknown-gene-reference", f"{_counted(len(unknown_refs), 'geneRef', 'points', 'point')} "
                   f"to a gene id that isn't defined in the <species> section.",
                   [f"geneRef {g} in family {f}" for g, f in unknown_refs])
    reused = {g: fams for g, fams in references.items() if len(fams) > 1}
    if reused:
        report.add(ERROR, "gene-referenced-twice", f"{_counted(len(reused), 'gene', 'is', 'are')} referenced more "
                   f"than once in the groups. pyHam would place such a gene in several HOGs and keep only the last "
                   f"as its parent.",
                   [f"gene {g} ({gene_species[g]}) in families {', '.join(dict.fromkeys(fams))}"
                    for g, fams in reused.items()],
                   "Each gene must belong to at most one group.")
    singletons = len(gene_species) - len(references)
    if singletons and n_families:
        report.add(NOTE, "singleton-genes", f"{singletons:,} of {len(gene_species):,} genes are not in any group. "
                   f"pyHam keeps them as genes without a HOG, so HOG-based tools don't apply to them.")

    report.stats["levels_not_in_tree"] = sum(c for c, _ in unresolved_levels.values())
    if n_groups and not n_claimed:
        report.add(NOTE, "no-level-annotations", "No HOG states its taxonomic level (TaxRange/taxid), so pyHam "
                   "infers each HOG's level as the common ancestor of its member species. The nesting of the "
                   "groups must then agree with the species tree, e.g. come from gene trees reconciled with it.")
    return species, unresolved_levels


def _level_issue(report, unresolved, n_claimed, tree_names, use_internal_name, tree_format):
    if not unresolved:
        return
    n_hogs = sum(count for count, _ in unresolved.values())
    items = [f'"{label}": {_plural(count, "HOG")}, e.g. in family {family}'
             for label, (count, family) in sorted(unresolved.items(), key=lambda kv: -kv[1][0])]
    labels = set(unresolved)

    hint = None
    if tree_format != "phyloxml" and not use_internal_name and tree_names.raw_internal_labels is not None:
        n_match = len(labels & tree_names.raw_internal_labels)
        if n_match and n_match >= len(labels) / 2:
            hint = (f"The tree's own internal node names match {n_match:,} of these levels. Set "
                    f"use_internal_name=True so pyHam uses them.")
        elif not tree_names.raw_internal_labels:
            hint = ("The species tree has no internal node names, so pyHam can't use the levels stated in the "
                    "OrthoXML. Use a tree with named internal nodes (e.g. the one the HOGs were inferred with) and "
                    "use_internal_name=True.")
    elif tree_format == "phyloxml":
        best = max(((len(labels & values), tag) for tag, values in tree_names.alt_internal_names.items()),
                   default=(0, None))
        if best[0] and best[0] >= len(labels) / 2 and (best[1] != tree_names.internal_tag or not use_internal_name):
            hint = (f"phyloxml_internal_name_tag='{best[1]}' with use_internal_name=True would match "
                    f"{best[0]:,} of these levels.")
    if hint is None:
        hint = ("These levels are not in the species tree. If it is a pruned or different version of the tree the "
                "HOGs were inferred with, pyHam infers these HOGs' levels from their member species instead, which "
                "can differ from the OrthoXML or make families fail to build. Prefer the original species tree.")

    report.add(WARNING, "level-not-in-tree",
               f"{n_hogs:,} of {n_claimed:,} HOGs state a taxonomic level that isn't a node of the species tree "
               f"({_plural(len(unresolved), 'distinct level')}). pyHam infers their level from their member species "
               f"instead.", items, hint, {"levels": {k: v[0] for k, v in unresolved.items()}})


def _check_species_mapping(report, species, taxonomy, tree_names):
    if taxonomy is not None:
        unmapped, used = [], collections.defaultdict(list)
        for name in dict.fromkeys(species):
            try:
                used[taxonomy.get_extant_taxa_by_name(name)].append(name)
            except KeyError as e:
                unmapped.append((name, e.args[0] if e.args else str(e)))
        shared = {node: names for node, names in used.items() if len(names) > 1}
        if shared:
            report.add(ERROR, "species-share-leaf", f"{_counted(len(shared), 'tree leaf', 'is', 'are')} matched by "
                       f"more than one OrthoXML species.",
                       [f'leaf "{node.name}": {", ".join(names)}' for node, names in shared.items()])
    else:
        leaves = set(tree_names.leaves)
        unmapped = [(name, "") for name in dict.fromkeys(species) if name not in leaves]
    if unmapped:
        _species_issue(report, unmapped, len(set(species)), tree_names, all_species=species)
    n_matched = len(set(species)) - len(unmapped)
    unused = len(tree_names.leaves) - n_matched
    if unused > 0 and n_matched > 0:
        report.add(NOTE, "tree-leaves-without-species", f"{unused:,} of {len(tree_names.leaves):,} tree leaves have "
                   f"no species in the OrthoXML. This is fine: pyHam only uses the part of the tree it needs.")


# ----------------------------------------------------------------------------------------------
# Stage 4: trial build
# ----------------------------------------------------------------------------------------------

def _make_tolerant_parser_class():
    from .parsers import OrthoXMLParser
    from .abstractgene import TaxonomicConflictError

    class TolerantOrthoXMLParser(OrthoXMLParser):
        """OrthoXMLParser that records a family that fails to build, skips the rest of it, and
        continues with the next family. The objects it builds are only good for validation."""

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.family_failures = []

        def start(self, tag, attrib):
            depth = len(self.hog_stack)
            try:
                super().start(tag, attrib)
            except (InputValidationError, TaxonomicConflictError):
                raise
            except Exception as e:
                if tag == _NS + "orthologGroup" and len(self.hog_stack) == depth:
                    self.hog_stack.append(0)  # keep the stack balanced for the matching end()
                self._record_failure(e)

        def end(self, tag):
            try:
                super().end(tag)
            except (InputValidationError, TaxonomicConflictError):
                raise
            except Exception as e:
                self._record_failure(e)

        def _record_failure(self, exc):
            self.family_failures.append((self.current_family_id, exc))
            self.paralog_stack = []
            self.in_paralogGroup = None
            self.paralogyNode = None
            if self.hog_stack:
                self.skip_this_hog = True
            else:
                self.current_family_id = None

    return TolerantOrthoXMLParser


def _trial_build(open_orthoxml, taxonomy, report):
    from .parsers import sniff_id_schema, _format_one_conflict

    with open_orthoxml(text=True) as fh:
        id_schema = sniff_id_schema(fh)
    target = _make_tolerant_parser_class()(types.SimpleNamespace(taxonomy=taxonomy), id_schema=id_schema)
    parser = XMLParser(target=target)
    try:
        with open_orthoxml(text=True) as fh:
            for line in fh:
                parser.feed(line)
    except Exception as e:
        report.add(ERROR, "trial-build-crashed", f"Building the HOGs stopped outside of any family: "
                   f"{type(e).__name__}: {e}")
        report.trial_build = "failed"
        return

    n_families = report.stats.get("families", 0)
    failures = target.family_failures
    if failures:
        items, messages = [], {}
        for family, exc in failures:
            message = " ".join(str(exc).split())
            messages[family] = f"{type(exc).__name__}: {message}"
            items.append(f"family {family}: {type(exc).__name__}: {_shorten(message, 600)}")
        hint = None
        if any("contradicts the species tree" in m or "older than itself" in m for m in messages.values()):
            hint = ("The group nesting of these families doesn't fit the species tree. If the OrthoXML was "
                    "derived from gene trees, reconcile them with this species tree first.")
        report.add(ERROR, "family-build-failed", f"{len(failures):,} of {_plural(n_families, 'family')} can't be "
                   f"built by pyHam.", items, hint, {"failures": messages})

    conflicts = target.taxonomic_conflicts
    if conflicts:
        n_families_hit = len({c.family_id for c in conflicts})
        report.add(ERROR, "taxonomic-conflicts", f"{len(conflicts):,} HOGs in {n_families_hit:,} families state a "
                   f"level that conflicts with the species tree.",
                   [_format_one_conflict(c) for c in conflicts],
                   "pyHam raises TaxonomicConflictError for these. The OrthoXML and the species tree disagree on "
                   "where these species belong; use the species tree the HOGs were inferred with.",
                   {"conflicts": conflicts})

    collapsed = target.collapsed_same_level_hogs
    if collapsed:
        report.add(NOTE, "same-level-hogs-collapsed", f"{len(collapsed):,} HOGs sit at the same level as their "
                   f"parent and are merged into it, as is usual for 'augmented' HOGs.",
                   [str(h) for h in collapsed])

    n_failed = len({f for f, _ in failures} | {c.family_id for c in conflicts})
    report.stats["families_built"] = n_families - n_failed
    if n_failed:
        report.trial_build = f"{n_families - n_failed:,} of {_plural(n_families, 'family')} built without errors"
    else:
        report.trial_build = f"ok, all {_plural(n_families, 'family')} built"


def _shorten(text, width):
    return text if len(text) <= width else text[:width - 3].rstrip() + "..."


# ----------------------------------------------------------------------------------------------
# Entry points
# ----------------------------------------------------------------------------------------------

def validate_input(tree_file, hog_file, tree_format=None, use_internal_name=False,
                   phyloxml_leaf_name_tag="taxonomy_scientific_name",
                   phyloxml_internal_name_tag="taxonomy_scientific_name",
                   orthoXML_as_string=False, trial_build=True):
    """Check a species tree and an OrthoXML file for problems that would stop pyHam, or make it
    silently build something different from what the files say.

    Takes the same arguments as :obj:`pyham.Ham`, so the check reflects what Ham would do with them.
    Unlike Ham, it doesn't stop at the first problem: it reports everything it finds.

    Args:
        | tree_file (:obj:`str`): path to the species tree file, or a newick string.
        | hog_file (:obj:`str`): path to the OrthoXML file (may be gzipped), or its content if
          `orthoXML_as_string` is True.
        | tree_format (:obj:`str`, optional): 'newick', 'newick_string' or 'phyloxml'. Defaults to
          guessing from the file extension.
        | use_internal_name (:obj:`bool`, optional): as for :obj:`pyham.Ham`. Defaults to False.
        | phyloxml_leaf_name_tag, phyloxml_internal_name_tag (:obj:`str`, optional): as for :obj:`pyham.Ham`.
        | orthoXML_as_string (:obj:`bool`, optional): `hog_file` is the OrthoXML content. Defaults to False.
        | trial_build (:obj:`bool`, optional): also build every family, as Ham would, and report the
          ones that fail. This takes about as long as loading the data with Ham. Defaults to True.

    Returns:
        :obj:`ValidationReport`: print it for a readable summary, or inspect its `issues`.
    """
    if tree_file is None or hog_file is None:
        raise TypeError("validate_input needs both a species tree and an OrthoXML file.")
    for tag in (phyloxml_leaf_name_tag, phyloxml_internal_name_tag):
        if tag not in _PHYLOXML_TAG_ATTR:
            raise ValueError(f"{tag!r} is not a phyloxml name tag. Options: {', '.join(_PHYLOXML_TAG_ATTR)}.")

    report = ValidationReport()
    tree_format = _detect_tree_format(tree_file, tree_format)

    def open_orthoxml(text=False):
        if orthoXML_as_string:
            return io.StringIO(hog_file) if text else io.BytesIO(hog_file.encode("utf-8"))
        open_ = gzip.open if str(hog_file).endswith(".gz") else open
        return open_(hog_file, "rt" if text else "rb")

    taxonomy, tree_names = _check_tree(tree_file, tree_format, use_internal_name, phyloxml_leaf_name_tag,
                                       phyloxml_internal_name_tag, report)
    tree_desc = f"{_src_label(tree_file, tree_format == 'newick_string')} ({tree_format.split('_')[0]}"
    if tree_names is not None:
        naming = ("internal names from the file" if use_internal_name
                  else "internal nodes named after their children")
        tree_desc += (f", {_plural(report.stats['tree_leaves'], 'leaf')}, "
                      f"{_plural(report.stats['tree_internal_nodes'], 'internal node')}, {naming}")
    report.inputs["species tree"] = tree_desc + ")"

    resolver = _LevelResolver(taxonomy, tree_names)
    species, unresolved_levels = _scan_orthoxml(open_orthoxml, resolver, report)
    if species is not None:
        report.inputs["OrthoXML"] = (f"{_src_label(hog_file, orthoXML_as_string)} "
                                     f"({_plural(report.stats['species'], 'species')}, "
                                     f"{_plural(report.stats['genes'], 'gene')}, "
                                     f"{_plural(report.stats['families'], 'family')})")
        if tree_names is not None:
            _check_species_mapping(report, species, taxonomy, tree_names)
            _level_issue(report, unresolved_levels, report.stats["hogs_with_level"], tree_names,
                         use_internal_name, tree_format)
    else:
        report.inputs["OrthoXML"] = _src_label(hog_file, orthoXML_as_string)

    if not trial_build:
        report.trial_build = "not run (trial_build=False)"
    elif report.errors or taxonomy is None or species is None:
        report.trial_build = "skipped, fix the errors above first"
    else:
        _trial_build(open_orthoxml, taxonomy, report)
    return report


def main(argv=None):
    """Command line entry point: ``pyham-validate TREE ORTHOXML [options]``."""
    p = argparse.ArgumentParser(
        prog="pyham-validate",
        description="Check a species tree and an OrthoXML file for problems before analysing them with pyHam. "
                    "Exits with status 1 if any error is found.")
    p.add_argument("tree", help="species tree file (newick or phyloXML)")
    p.add_argument("orthoxml", help="OrthoXML file (may be gzipped)")
    p.add_argument("--tree-format", choices=["newick", "phyloxml"],
                   help="format of the tree file (default: guessed from the extension)")
    p.add_argument("--use-internal-name", action="store_true",
                   help="use the internal node names of the tree, as Ham(use_internal_name=True)")
    p.add_argument("--phyloxml-leaf-name-tag", choices=list(_PHYLOXML_TAG_ATTR), default="taxonomy_scientific_name")
    p.add_argument("--phyloxml-internal-name-tag", choices=list(_PHYLOXML_TAG_ATTR),
                   default="taxonomy_scientific_name")
    p.add_argument("--no-trial-build", action="store_true", help="skip building the families (faster)")
    p.add_argument("--max-items", type=int, default=10, help="entries shown per issue (default: 10; 0 for all)")
    args = p.parse_args(argv)

    report = validate_input(args.tree, args.orthoxml, tree_format=args.tree_format,
                            use_internal_name=args.use_internal_name,
                            phyloxml_leaf_name_tag=args.phyloxml_leaf_name_tag,
                            phyloxml_internal_name_tag=args.phyloxml_internal_name_tag,
                            trial_build=not args.no_trial_build)
    print(report.format(max_items=args.max_items or None))
    return 1 if report.errors else 0


if __name__ == "__main__":
    sys.exit(main())
