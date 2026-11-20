"""Branch labels from the tree window: what codeml receives is what the window shows.

Before, a click on an internal node wrote #1 on the tips of the clade only, so the
ancestral branch and the inner branches stayed background while the window coloured
the whole clade."""
import io
from types import SimpleNamespace

import pytest
from Bio import Phylo

from src.backend import branch_labels as bl
from src.backend.methods_text import build_methods_text

TREE = "((Hsa,(Ptr,Ppa)),(Ggo,(Ppy,(Mmu,Cja))));"


def _tree(text=TREE):
    return Phylo.read(io.StringIO(text), 'newick')


def _newick(marks_by_names, text=TREE):
    """marks given as {(tip, tip): (num, scope)} on the common ancestor, or {tip: …}."""
    t = _tree(text)
    marks = {}
    for names, mark in marks_by_names.items():
        clade = t.common_ancestor(*names) if isinstance(names, tuple) else t.find_any(name=names)
        marks[clade] = mark
    return bl.to_newick(t, bl.effective_labels(t, marks))


def test_first_click_labels_the_ancestral_branch_only():
    assert _newick({('Hsa', 'Ppa'): (1, bl.STEM)}) == "((Hsa,(Ptr,Ppa)) #1,(Ggo,(Ppy,(Mmu,Cja))));"


def test_second_click_labels_the_whole_clade_with_its_stem():
    assert _newick({('Hsa', 'Ppa'): (1, bl.CLADE)}) == \
        "((Hsa #1,(Ptr #1,Ppa #1) #1) #1,(Ggo,(Ppy,(Mmu,Cja))));"


def test_a_tip_labels_its_terminal_branch():
    assert _newick({'Ggo': (2, bl.STEM)}) == "((Hsa,(Ptr,Ppa)),(Ggo #2,(Ppy,(Mmu,Cja))));"


def test_label_nearer_the_tips_wins():
    assert _newick({('Hsa', 'Ppa'): (1, bl.CLADE), 'Ptr': (2, bl.STEM)}) == \
        "((Hsa #1,(Ptr #2,Ppa #1) #1) #1,(Ggo,(Ppy,(Mmu,Cja))));"


def test_written_tree_reads_back_to_the_same_labels():
    nwk = _newick({('Hsa', 'Ppa'): (1, bl.CLADE), ('Ppy', 'Cja'): (2, bl.STEM)})
    tree, labels = bl.read_labels(nwk)
    assert bl.to_newick(tree, labels) == nwk


def test_dollar_label_follows_the_paml_manual():
    """PAML manual, "clade label": the two trees of each pair are equivalent."""
    pairs = [("(((rabbit, rat) $1, human), goat_cow, marsupial);",
              "(((rabbit #1, rat #1) #1, human), goat_cow, marsupial);"),
             ("((((rabbit, rat) $2, human #3), goat_cow) $1, marsupial);",
              "((((rabbit #2, rat #2) #2, human #3) #1, goat_cow #1) #1, marsupial);")]
    for dollar, hashes in pairs:
        a, b = bl.read_labels(dollar), bl.read_labels(hashes)
        assert bl.to_newick(*a) == bl.to_newick(*b)


def test_descriptions_say_which_branches():
    t = _tree()
    clade = t.common_ancestor('Hsa', 'Ppa')
    stem = bl.describe(t, bl.effective_labels(t, {clade: (1, bl.STEM)}))
    whole = bl.describe(t, bl.effective_labels(t, {clade: (1, bl.CLADE)}))
    assert stem == ["#1: the ancestral branch of Hsa, Ptr, Ppa (1 branch in total)"]
    assert whole == ["#1: the whole clade Hsa, Ptr, Ppa, ancestral branch included (5 branches in total)"]
    assert bl.describe(t, {}, lang='pt') == ["nenhum ramo marcado"]


def test_methods_text_names_the_labelled_branches():
    text = build_methods_text(version="0.4.0", codeml_version="4.10.10", models=['Branch-site', 'Branch-site_null'],
                              ctl={'CodonFreq': 2, 'kappa': 2, 'cleandata': 1}, omega0=0.5, pruned=True,
                              family_sizes={}, n_genes=1,
                              foreground={'Branch-site': ["#1: the ancestral branch of Hsa, Ptr, Ppa (1 branch in total)"]})
    assert "Branch labels (Branch-site model): #1, the ancestral branch of Hsa, Ptr, Ppa" in text


def test_pruning_keeps_labels_on_internal_branches():
    from src.backend.codeml_backend import CodemlBatchAnalysis
    tree, _ = bl.read_labels("((Hsa #1,(Ptr #1,Ppa #1) #1) #1,(Ggo,(Ppy,(Mmu,Cja))));")
    CodemlBatchAnalysis._prune_keep_labels(tree, {'Hsa', 'Ptr', 'Ggo', 'Ppy', 'Mmu', 'Cja'})
    out = io.StringIO()
    Phylo.write(tree, out, 'newick')
    assert bl.to_newick(*bl.read_labels(out.getvalue())) == \
        "((Hsa #1,Ptr #1) #1,(Ggo,(Ppy,(Mmu,Cja))));"


@pytest.fixture
def window():
    """TreeLabelWindow._on_pick on a stand-in object (no display needed)."""
    gui = pytest.importorskip('src.gui.main_gui')
    W = gui.TreeLabelWindow
    log = []
    w = SimpleNamespace(tree=_tree(), marks={}, mode='branchsite',
                        parent=SimpleNamespace(append_log=log.append),
                        _draw_tree=lambda: None, _refresh_legend=lambda: None)
    w._get_clade_name = lambda c: W._get_clade_name(w, c)
    w._labels = lambda: W._labels(w)
    w.click = lambda clade: W._on_pick(w, SimpleNamespace(artist=SimpleNamespace(clade_obj=clade)))
    return w


def test_clicks_cycle_stem_clade_off(window):
    clade = window.tree.common_ancestor('Hsa', 'Ppa')
    newick = lambda: bl.to_newick(window.tree, window._labels())  # noqa: E731
    window.click(clade)
    assert newick() == "((Hsa,(Ptr,Ppa)) #1,(Ggo,(Ppy,(Mmu,Cja))));"
    window.click(clade)
    assert newick() == "((Hsa #1,(Ptr #1,Ppa #1) #1) #1,(Ggo,(Ppy,(Mmu,Cja))));"
    window.click(clade)
    assert newick() == TREE
    window.click(window.tree.root)            # the root has no branch
    assert window.marks == {}
