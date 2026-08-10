"""Detection and synthesis helpers for the HOG id "scheme" used by an OrthoXML file.

Three schemes are recognised:

- LOFT: after a duplication, paralogous branches get a ``.<num><letter>`` suffix
  appended (e.g. ``1074943`` -> ``1074943.2a`` -> ``1074943.2a.7b``). Non-duplication
  levels keep the same id as they move up/down the tree, so two different taxonomic
  levels can end up sharing an id.
- LOFT_TAXID: LOFT plus a trailing ``_<taxid>`` tied to the specific taxonomic level,
  which restores per-level uniqueness (e.g. ``HOG:0402997.2b.1a_47424``).
- GENERIC: no inferable structure; ids are opaque and used as-is.
"""

import re

LOFT = "LOFT"
LOFT_TAXID = "LOFT_TAXID"
GENERIC = "GENERIC"
SCHEMES = {LOFT, LOFT_TAXID, GENERIC}

# <fam><subhog>[_<taxid>], e.g. "HOG:0402997.2b.1a_47424" or "HOG:0402997" or "3".
# Kept identical to the grammar Ham.get_hog_by_id used to inline, for backward compatibility.
HOG_ID_RE = re.compile(r"^(?P<fam>(HOG:[A-Z]?)?(\d+))(?P<subhog>[0-9a-z.]*)(_(?P<taxid>-?\d+))?$")

# Strict LOFT dot-chain shape: one or more ".<digits><letters>" groups.
_LOFT_CHAIN_RE = re.compile(r"^(\.\d+[a-zA-Z]+)+$")

_CHAIN_SHAPE_THRESHOLD = 0.8
_TAXID_SUFFIX_THRESHOLD = 0.5


def detect_id_scheme(sample_ids):
    """Classify a batch of raw HOG id strings into one of LOFT / LOFT_TAXID / GENERIC.

    Ids without a subhog part (e.g. a bare top-level id like "HOG:0402997") are
    ignored: they say nothing about the nested-level convention used by a file.
    """
    matches = [HOG_ID_RE.match(s) for s in sample_ids]
    matches = [m for m in matches if m]

    with_subhog = [m for m in matches if m.group('subhog')]
    if not with_subhog:
        return GENERIC

    chain_shaped = [m for m in with_subhog if _LOFT_CHAIN_RE.match(m.group('subhog'))]
    if len(chain_shaped) / len(with_subhog) < _CHAIN_SHAPE_THRESHOLD:
        return GENERIC

    with_taxid = [m for m in chain_shaped if m.group('taxid') is not None]
    if len(with_taxid) / len(chain_shaped) > _TAXID_SUFFIX_THRESHOLD:
        return LOFT_TAXID
    return LOFT


def _taxid_for_taxon(taxon):
    """Best-effort numeric/string id for an ete4 taxonomy node.

    Newick-derived trees carry no id props at all; phyloxml-derived trees store it
    under 'taxon_id'; orthoxml-embedded taxonomies store it under 'id'. The `_taxid`
    suffix must stay numeric to remain parseable by `HOG_ID_RE`/`Ham.get_hog_by_id`, so
    this returns None (rather than e.g. the node's name) when no numeric id is found.
    """
    for key in ('id', 'taxon_id'):
        value = taxon.props.get(key)
        if value is None:
            continue
        try:
            return int(value)
        except (TypeError, ValueError):
            continue
    return None


def make_missing_hog_id(scheme, anchor_id, taxon):
    """Build an id for a HOG synthesized for a taxonomic level absent from the XML.

    Reuses the LOFT dot-chain of `anchor_id` (a neighboring real HOG's id) and, under
    LOFT_TAXID, swaps in the taxid of `taxon` so each synthesized level gets a
    distinct, resolvable id. Falls back to returning `anchor_id` unchanged (today's
    behavior) under LOFT/GENERIC, if `anchor_id` doesn't parse, or if no numeric taxid
    can be found for `taxon` (e.g. a newick tree with no embedded taxon ids) -- a
    non-unique id is preferable to a synthesized one that can't be looked back up.
    """
    if scheme != LOFT_TAXID or anchor_id is None:
        return anchor_id

    m = HOG_ID_RE.match(anchor_id)
    if not m:
        return anchor_id

    taxid = _taxid_for_taxon(taxon)
    if taxid is None:
        return anchor_id

    base = m.group('fam') + m.group('subhog')
    return f"{base}_{taxid}"
