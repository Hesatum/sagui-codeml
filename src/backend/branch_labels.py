"""Branch labels for the Branch and Branch-site models: the Newick that codeml
receives and a plain description of which branches each label covers.

The tree window writes every labelled branch explicitly with '#N' (never '$N'),
so what codeml tests is exactly what the window shows. A label on a node covers
the branch leading to that node; a label on a whole clade is written on its stem
and on every branch inside it, like codeml's '$N' (PAML manual, "clade label").
"""

import io
import re
from typing import Callable, Dict, List, Optional

try:
    from Bio import Phylo
except ImportError:          # the window and the run check for Biopython themselves
    Phylo = None

STEM, CLADE = 'stem', 'clade'
_LABEL_RE = re.compile(r'\s*[#$](\d+)\s*$')


def parent_map(tree) -> Dict[object, object]:
    return {child: parent for parent in tree.find_clades() for child in parent.clades}


def effective_labels(tree, marks: Dict[object, tuple]) -> Dict[object, int]:
    """Label of every branch (keyed by the clade the branch leads to) from the marks
    {clade: (number, STEM | CLADE)}. A clade mark covers its stem and the branches
    inside it; a mark nearer the tips wins, as in codeml. The root has no branch."""
    parents = parent_map(tree)
    labels = {}
    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        if clade in marks:
            labels[clade] = marks[clade][0]
            continue
        node = parents.get(clade)
        while node is not None:
            mark = marks.get(node)
            if mark and mark[1] == CLADE:
                labels[clade] = mark[0]
                break
            node = parents.get(node)
    return labels


def _clean_name(name: Optional[str]) -> str:
    return _LABEL_RE.sub('', name or '').strip()


def to_newick(tree, labels: Dict[object, int]) -> str:
    """Topology with tip names and '#N' on each labelled branch; no branch lengths or
    support values (codeml takes the labelled tree without lengths)."""
    def write(clade) -> str:
        if clade.is_terminal():
            text = _clean_name(clade.name)
        else:
            text = '(' + ','.join(write(c) for c in clade.clades) + ')'
        if clade in labels:
            text += f" #{labels[clade]}"
        return text
    return write(tree.root) + ';'


def read_labels(newick: str):
    """(tree, {clade: number}) from a Newick with '#N' and '$N' labels, applying the
    PAML rules: '$N' covers the whole clade, '#N' beats '$N', labels nearer the tips
    win. None if the tree cannot be read."""
    if Phylo is None or not newick or '(' not in newick:
        return None
    text = newick[newick.find('('):]
    text = re.sub(r':\s*-?[0-9.]+(?:[eE][-+]?\d+)?', '', text)
    # Bio.Phylo keeps ') #1' as the node name only without the space
    text = re.sub(r'\)\s+([#$]\d+)', r')\1', text)
    text = re.sub(r'(?<=[^\s(),;])\s+([#$]\d+)', r'\1', text)
    try:
        tree = Phylo.read(io.StringIO(text), 'newick')
    except Exception:
        return None
    marks = {}
    for clade in tree.find_clades():
        m = re.search(r'([#$])(\d+)\s*$', clade.name or '')
        if m:
            marks[clade] = (int(m.group(2)), CLADE if m.group(1) == '$' else STEM)
    # '#' beats '$' on the same node: a node written '#N' is never a clade mark here
    return tree, effective_labels(tree, marks)


def _tips(clade, limit: int = 3) -> str:
    names = [_clean_name(t.name) or '?' for t in clade.get_terminals()]
    if len(names) <= limit:
        return ', '.join(names)
    return ', '.join(names[:limit]) + f" +{len(names) - limit}"


_TEXT = {
    'en': {
        'tip': "the branch to {tip}",
        'stem': "the ancestral branch of {tips}",
        'clade': "the whole clade {tips}, ancestral branch included",
        'inner': "some of the branches of the clade {tips}",
        'none': "no labelled branches",
        'line': "#{num}: {parts} ({n} branch{es} in total)",
    },
    'pt': {
        'tip': "o ramo terminal de {tip}",
        'stem': "o ramo ancestral de {tips}",
        'clade': "o clado inteiro {tips}, com o ramo ancestral",
        'inner': "parte dos ramos do clado {tips}",
        'none': "nenhum ramo marcado",
        'line': "#{num}: {parts} ({n} ramo{es} no total)",
    },
}


def describe(tree, labels: Dict[object, int], lang: str = 'en') -> List[str]:
    """One line per label number saying which branches it covers, e.g.
    '#1: the ancestral branch of Hsa, Ptr, Ppa (1 branch in total)'."""
    t = _TEXT.get(lang, _TEXT['en'])
    if not labels:
        return [t['none']]
    parents = parent_map(tree)
    lines = []
    for num in sorted(set(labels.values())):
        covered = [c for c, n in labels.items() if n == num]
        tops = [c for c in covered if labels.get(parents.get(c)) != num]
        parts = []
        for top in tops:
            if top.is_terminal():
                parts.append(t['tip'].format(tip=_clean_name(top.name) or '?'))
                continue
            inside = [c for c in top.find_clades() if c is not top]
            n_in = sum(1 for c in inside if labels.get(c) == num)
            key = 'stem' if n_in == 0 else ('clade' if n_in == len(inside) else 'inner')
            parts.append(t[key].format(tips=_tips(top)))
        n = len(covered)
        es = '' if n == 1 else ('es' if lang != 'pt' else 's')
        lines.append(t['line'].format(num=num, parts='; '.join(parts), n=n, es=es))
    return lines


def describe_newick(newick: str, lang: str = 'en') -> List[str]:
    """describe() for a Newick string; [] if it cannot be read."""
    parsed = read_labels(newick)
    if parsed is None:
        return []
    return describe(*parsed, lang=lang)
