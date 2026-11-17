"""Bugs found by the validation for the manuscript (validation/debate/sintese/S5_sintese.md,
section "Bugs do programa"), one test each."""
import itertools
import json
import platform
import re
from pathlib import Path

import pytest

from src.backend import lrt_stats, messages
from src.backend.codeml_backend import CodemlBatchAnalysis
from src.backend.methods_text import build_methods_text
from src.backend.sites_parser import SitesParser
from tests.test_runner_failures import FAKE, _app, fake_codeml  # noqa: F401 (fixture)

MLC = Path(__file__).resolve().parent / 'data' / 'mlc'

needs_fake = pytest.mark.skipif(platform.system() == 'Windows',
                                reason="the fake codeml uses a shebang (Linux/macOS)")


# 1. a BEB site whose amino acid in the first sequence is a gap ---------------------

def test_beb_site_with_gap_in_first_sequence_is_kept():
    """lysin, cleandata = 0: codeml lists '134 -  0.960*' (gap in the 1st sequence)."""
    sites = SitesParser.parse_sites_from_file(MLC / 'lysin_clean0_M2a.txt', 'BEB')
    assert len(sites) == 33
    row = sites[sites['position'] == 134].iloc[0]
    assert row['amino_acid'] == '-' and row['pr_w_gt_1'] == pytest.approx(0.960)
    assert row['significance'] == '*'


def test_branch_site_beb_site_with_gap_is_kept(tmp_path):
    f = tmp_path / 'bs.txt'
    f.write_text("Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005)\n"
                 "Positive sites for foreground lineages Prob(w>1):\n"
                 "    12 K 0.970*\n    40 - 0.991**\n    41 ? 0.960*\n\n")
    sites = SitesParser.parse_sites_from_file(f, 'BEB')
    assert list(sites['position']) == [12, 40, 41]
    assert list(sites['amino_acid']) == ['K', '-', '?']


# 2. 2Δℓ < 0 is re-optimised, then marked as a failure ------------------------------

NEG_FAKE = FAKE.replace(
    "elif mode == 'idle':",
    "elif mode == 'neg':\n"
    "    w0 = float(re.search(r'^\\s*omega\\s*=\\s*([\\d.]+)', ctl, re.M).group(1))\n"
    "    lnl = -4100.0\n"
    "    if sys.argv[1].endswith('_M8.ctl'):   # M8 misses the optimum unless it starts at ω 2.5\n"
    "        good = w0 >= 2.5 and not os.environ.get('FAKE_NEG_ALWAYS')\n"
    "        lnl = -4099.0 if good else -4110.0\n"
    "    with open(out, 'w') as f:\n"
    "        f.write('CODONML (in paml version 9.9fake)\\nns =  10  ls = 300\\n')\n"
    "        f.write('lnL(ntime: 17  np: 22):  %.6f      +0.000000\\n' % lnl)\n"
    "elif mode == 'idle':")


@pytest.fixture
def neg_codeml(tmp_path):
    import stat
    import sys
    path = tmp_path / 'fake_codeml_neg'
    path.write_text(NEG_FAKE.format(python=sys.executable))
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


@needs_fake
def test_negative_lrt_is_reoptimised(tmp_path, neg_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'neg')
    app = _app(tmp_path, neg_codeml, models=('M8a', 'M8'))
    app.config['beta_starts'] = 1
    app.run_batch_analysis()
    assert app.results['gene']['M8']['lnL'] == -4099.0
    assert '-4099.000000' in (tmp_path / 'out' / 'M8' / 'gene_M8_results.txt').read_text()
    assert app._lrt_pvalues[('M8a', 'M8')]['gene'] < 1


@needs_fake
def test_negative_lrt_that_persists_is_a_failure_not_p_1(tmp_path, neg_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'neg')
    monkeypatch.setenv('FAKE_NEG_ALWAYS', '1')
    app = _app(tmp_path, neg_codeml, models=('M8a', 'M8'))
    app.config['beta_starts'] = 1
    app.run_batch_analysis()
    out = tmp_path / 'out'
    assert 'gene' not in app._lrt_pvalues[('M8a', 'M8')]       # out of p, q and BH
    lrt = (out / 'LRT_results.txt').read_text()
    assert 'FAILED' in lrt and 'p-value = 1' not in lrt
    header, row = (out / 'analysis_summary.tsv').read_text().splitlines()[:2]
    cols = dict(zip(header.split('\t'), row.split('\t')))
    assert cols['p_M8a_vs_M8'] == 'NA' and cols['q_M8a_vs_M8'] == 'NA'
    assert 'M8 vs M8a: optimisation failure' in (out / 'genes_status.tsv').read_text()
    # the rebuilt LRT file applies the same rule
    CodemlBatchAnalysis._regenerate_lrt_results(out)
    assert 'FAILED' in (out / 'LRT_results.txt').read_text()


def test_optimization_failure_threshold():
    assert lrt_stats.optimization_failed(-28.28)
    assert not lrt_stats.optimization_failed(-0.001)    # rounding, set to 0
    assert not lrt_stats.optimization_failed(3.0)


# 3. gene_verdict follows METHODS.md ----------------------------------------------

def _methods_rule(sig, s2, s8):
    """METHODS.md, table "Evidence"."""
    a, b, c = sig[('M1a', 'M2a')], sig[('M7', 'M8')], sig[('M8a', 'M8')]
    if (c and s8 > 0) or (a and s2 > 0):
        return 'positive'
    if b and c is False and not a:
        return 'neutral'
    if b and c is None and not a:
        return 'no_m8a'
    if a or b or c:
        return 'weak'
    return 'none'


def test_gene_verdict_matches_methods_in_every_combination():
    from src.gui.summary_tab import gene_verdict, test_label
    diverging = []
    for a, b, c in itertools.product([True, False, None], repeat=3):
        for s2, s8 in itertools.product([0, 2], repeat=2):
            sig = {('M1a', 'M2a'): a, ('M7', 'M8'): b, ('M8a', 'M8'): c}
            sig_by_pair = {k: v for k, v in sig.items() if v is not None}
            rows = [(test_label(pair), v, '', '', '', s2 if pair[1] == 'M2a' else s8, 3.0, 0.1)
                    for pair, v in sig.items() if v is not None]
            got = gene_verdict(sig_by_pair, rows)[0]
            if got != _methods_rule(sig, s2, s8):
                diverging.append(((a, b, c), (s2, s8), got))
    assert diverging == []


# 4. run_config.json records the values the run used ------------------------------

@needs_fake
def test_run_config_has_the_codonfreq_used_and_the_warm_start(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml, models=('M0', 'M8'))
    app.config.update({'CodonFreq': 7, 'warm_start_m0': True, 'beta_starts': 1})
    tree = (Path(__file__).resolve().parent / 'data' / 'gene_example.nwk').read_text().strip()
    monkeypatch.setattr(app, '_extract_fitted_tree', lambda path: tree)
    monkeypatch.setattr(app, '_extract_kappa', lambda path: 2.5)
    app.run_batch_analysis()
    assert 'CodonFreq = 7' in (tmp_path / 'out' / 'M8' / 'gene_M8.ctl').read_text()
    rc = json.loads((tmp_path / 'out' / 'run_config.json').read_text())
    for m in ('M0', 'M8'):
        assert rc['ctl_parameters_by_model'][m]['CodonFreq'] == 7
    start = rc['starts_by_gene']['gene']['M8']
    assert start['kappa'] == 2.5 and start['fix_blength'] == 1
    assert start['omega0'] in CodemlBatchAnalysis._WARM_START_OMEGA_TRIALS
    assert rc['starts_by_gene']['gene']['M0']['fix_blength'] == 0


# 5. summary ω of M2a/M8 -------------------------------------------------------------

@pytest.mark.parametrize('name, expected', [
    ('hiv4_M2a.txt', 1.9777), ('hiv4_M8.txt', 1.9785),    # 4 taxa: NEB lines were read as branches
    ('ds0_M2a.txt', 999.0), ('ds0_M8.txt', 998.99), ('ds0_M0.txt', 999.0),   # dS = 0
])
def test_summary_omega_small_trees_and_ds_zero(name, expected):
    assert SitesParser.extract_omega_robust(MLC / name) == pytest.approx(expected, abs=0.01)


# 6. columns removed by cleandata are reported ---------------------------------------

@needs_fake
def test_cleandata_removed_columns_are_reported(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'ok')
    app = _app(tmp_path, fake_codeml, models=('M0',))
    fas = tmp_path / 'in' / 'gene.fasta'
    lines = fas.read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l and not l.startswith('>'))
    lines[i] = '---NNN' + lines[i][6:]           # codons 1 and 2 of the first sequence
    fas.write_text('\n'.join(lines) + '\n')
    app.run_batch_analysis()
    out = tmp_path / 'out'
    assert re.search(r'cleandata.*2 of 300 codon', (out / 'batch_analysis_log.txt').read_text())
    assert '2 of 300 codon' in (out / 'genes_status.tsv').read_text()
    assert '2 of 300 codon columns' in (out / 'methods_text.txt').read_text()


def test_methods_text_reports_cleandata_columns():
    txt = build_methods_text(version="1.0", codeml_version="4.10.10", models=['M0'],
                             ctl={'CodonFreq': 2, 'kappa': 2, 'cleandata': 1}, omega0=0.5,
                             pruned=True, family_sizes={}, n_genes=2,
                             cleandata_removed={'a': (13, 135), 'b': (0, 100)})
    assert "13 of 235 codon columns (5.5%) were removed" in txt


# 7. several starts for M7, M8 and M8a ------------------------------------------------

@needs_fake
def test_beta_models_run_several_starts_and_keep_the_best(tmp_path, fake_codeml, monkeypatch):
    monkeypatch.setenv('FAKE_CODEML_MODE', 'omega')     # best lnL when the start ω is 1
    app = _app(tmp_path, fake_codeml, models=('M1a', 'M7'))
    app.config['beta_starts'] = 4
    calls = []
    real = app._run_single_analysis
    monkeypatch.setattr(app, '_run_single_analysis',
                        lambda *a, **k: calls.append(k.get('model_name')) or real(*a, **k))
    app.run_batch_analysis()
    assert calls.count('M7') == 4 and calls.count('M1a') == 1
    assert app.results['gene']['M7']['lnL'] == -4100.0
    rc = json.loads((tmp_path / 'out' / 'run_config.json').read_text())
    assert rc['options']['beta_starts'] == 4
    assert 'M7 was fitted from 4 starts' in (tmp_path / 'out' / 'methods_text.txt').read_text()


def test_default_number_of_beta_starts():
    assert CodemlBatchAnalysis.DEFAULT_BETA_STARTS >= 5
