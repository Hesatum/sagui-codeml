# SAGUI methods

What SAGUI does to the data, for writing a methods section or reviewing a
manuscript that used it. Installation, usage and the list of output files are in
the [README](README.md). Every run records its own settings in `run_config.json`,
in the `.ctl` of each gene and model, and in `LRT_results.txt`. This describes
version 0.4.0; cite the version and commit recorded in `run_config.json`.

## Input

Alignments: one file per gene (`.fasta`, `.fas`, `.fa`, `.fna`, `.phy`, `.phylip`).
PHYLIP may be strict (10-column names) or relaxed, sequential or interleaved. If a
gene has both a FASTA and a PHYLIP file, the FASTA file is used and the other is
reported as ignored.

Tree: Newick, rooted or unrooted, with or without branch lengths and with or
without an `N 1` header. A per-gene tree (in the tree folder or the alignments
folder) replaces the general tree for that gene. It is paired with the gene by file
name, without extensions and the names IQ-TREE and RAxML give (`GENE.nwk`,
`GENE.fasta.treefile`, `GENE.contree`, `RAxML_bestTree.GENE`, `GENE.raxml.bestTree`)
and ignoring case; near names are reported, never paired. The pairs are listed in
`run_config.json` (`per_gene_trees`). Branch labels for the Branch and Branch-site
models (`#1`, `#2` …) are set in the "Label branches" window, or come from a per-gene
tree that has them (`Name#1` or `Name #1`); site models get the tree without labels.

codeml is searched in this order: `--codeml` or the setting in the window, the
`SAGUI_CODEML` variable, `bin/codeml(.exe)`, then `PATH`.

## Checks before running

| Finding | Consequence |
|---|---|
| fewer than 3 sequences, unequal lengths, length not a multiple of 3, repeated name | the gene does not run and is reported as failed |
| internal stop codon (not in the last codon) | the codon is masked (codeml treats the column as missing data) if "Ignore stop codons" is on, `--ignore-stop-codons` is given, or "Continue anyway" is chosen in the data check; otherwise the gene fails with the sequence and codon position. Masked codons are listed in `genes_status.tsv` and `methods_text.txt` |
| sequence missing from the tree | with automatic pruning (default) the sequence is left out of that gene and the closest tree name is suggested; without pruning codeml fails |
| tree taxon missing from the alignment | with automatic pruning it is removed from that gene's tree |

`--strict` makes the command-line mode stop if any problem is found.

## What codeml receives

For each gene and model, codeml gets a FASTA copy of the alignment (only the
sequences present in the tree) and the pruned tree. For the site models (M0, M1a,
M2a, M7, M8, M8a) a rooted tree is unrooted, as PAML requires.

Every parameter is written to the `.ctl` file. Leaving one out changes results:
without `ncatG`, codeml uses 4 beta categories when M7 or M8 runs alone and 10 when
`NSsites = 7 8` runs together (checked in PAML 4.9j and 4.10.10).

| Parameter | Default | Note |
|---|---|---|
| `seqtype` | 1 | codons |
| `CodonFreq` | 2 (F3x4) | 0 Fequal, 1 F1x4, 2 F3x4, 3 F61, 4 F1x4MG, 5 F3x4MG, 6 FMutSel0, 7 FMutSel |
| `estFreq` | 0 | observed frequencies |
| `icode` | 0 | universal code |
| `cleandata` | 1 | drop columns with a gap, ambiguity or stop codon in any sequence |
| `kappa` / `fix_kappa` | 2 / 0 | κ estimated from 2 |
| `omega` / `fix_omega` | 0.5 / 0 | ω estimated from 0.5 (M8a and the branch-site null: 1 / 1) |
| `ncatG` | 10 | beta categories (M7, M8, M8a) |
| `fix_alpha` / `alpha` / `Malpha` | 1 / 0 / 0 | no gamma rate variation |
| `clock` | 0 | no clock |
| `fix_blength` | 0 | branch lengths estimated from scratch (1 with `--warm-start-m0`) |
| `method` | 0 | simultaneous optimization |
| `getSE`, `RateAncestor` | 0, 0 | |
| `Small_Diff` | 0.5×10⁻⁶ | convergence criterion |

| Model | model | NSsites | Note |
|---|---|---|---|
| M0 | 0 | 0 | |
| M1a | 0 | 1 | ω₁ = 1 imposed by codeml |
| M2a | 0 | 2 | |
| M7 | 0 | 7 | beta(p, q) |
| M8 | 0 | 8 | beta plus a free-ω class |
| M8a | 0 | 8 | `fix_omega = 1`, `omega = 1` (Swanson et al. 2003) |
| Branch | 2 | 0 | labelled tree, branch lengths removed |
| Branch-site and null | 2 | 2 | null with `fix_omega = 1`, `omega = 1` |

Some starting values are set by codeml itself (`codeml.c`, PAML 4.10.10):

- In M8, an initial ω below 1 (the default 0.5 included) is replaced by a random
  value between 2 and 3, and codeml prints `initial w for M8:NSbetaw>1 reset.`
- In M7, M8 and M8a, the beta parameters start from random values (p between 0.2
  and 1.2, q between 1 and 2; p₀ starts at 0.9). Two runs of the same `.ctl` can
  therefore stop at different optima. SAGUI runs each of these models
  `--beta-starts` times (5 by default; "Multiple starts" below).
- With `fix_blength = 0` codeml ignores the input branch lengths. The tree file
  written for codeml keeps the lengths it was given and adds none.

## Execution

codeml runs with standard input closed, so a "Press Enter" prompt (for example
after a stop codon) does not block it.

Each run has a time limit, the same in the window and on the command line. By
default it is 5 times the expected time for that model and gene size (taxa ×
codons), and at least 30 minutes; the calibration is in
[docs/timing_benchmark.md](docs/timing_benchmark.md). A run is also stopped if
codeml, counting its child processes, uses no CPU for `--idle-timeout` seconds
(300 by default).

A gene fails if any requested model exits with an error, writes no output, reports
no lnL, or is stopped by a limit. The partial output of a failed run is kept as
`*_results_FAILED.txt` and does not enter any test.

### Multiple starts (M7, M8, M8a)

M7, M8 and M8a are each run `--beta-starts` times (default 5; 1 gives a single
run), and the run with the highest lnL is kept with its files. Each start is a
separate codeml run, so codeml draws new random starting values for the beta
parameters; the starts after the first also change the initial ω (0.2, 1.0, 2.5,
0.05, 1.5), which M8 uses. Starts are not fully independent: with HIV *env* (13 taxa),
a single run of M7 reached the best lnL in about 1 run in 5 to 1 in 10 and changing
the initial κ or ω did not help, and 3 of 6 runs with 5 starts reached it. The lnL of
every start is listed in the log and in `run_config.json` (`starts_by_gene`). If the
first start fails (error, time limit), the model is not run again.

When the output folder already holds a finished run whose `.ctl`, alignment and
tree are byte-for-byte the ones about to be used, codeml is not run again and that
output is reused; the log says so. `--rerun-all` runs everything again.

## Likelihood ratio tests

2Δℓ = 2(lnL_alternative − lnL_null). The alternative contains the null, so at
their optima 2Δℓ ≥ 0. A value between −0.01 and 0 is rounding of the printed lnL
and is set to 0 (p = 1). A value below −0.01 means the alternative stopped below
the optimum of the null: the alternative is run again from five other initial ω
(0.2, 1.0, 2.5, 0.05, 1.5, each with new random starting values) and kept if its
lnL is higher. If 2Δℓ is still below −0.01, the test is an optimization failure:
`LRT_results.txt` marks it FAILED, it has no p- or q-value in `analysis_summary.tsv`,
it is left out of the Benjamini-Hochberg correction of that test, and
`genes_status.tsv` notes it. It is never reported as p = 1.

| Test | df | Null distribution |
|---|---|---|
| M1a vs M2a | 2 | χ²₂ |
| M7 vs M8 | 2 | χ²₂ |
| M8a vs M8 | 1 | χ²₁ (*) |
| M0 vs M1a | 1 | χ²₁ (*) |
| M0 vs Branch | number of foreground groups | χ²_df |
| Branch-site null vs Branch-site | 1 | χ²₁ (*) |

(*) The null lies on the boundary of the parameter space (ω = 1 fixed, or p₀ = 1 for
M0 vs M1a), so the asymptotic distribution is a 50:50 mixture of χ²₀ and χ²₁ (Self
& Liang 1987). SAGUI uses χ²₁, which is conservative, as the PAML manual
recommends for the branch-site test; the mixture p-value is printed in
`LRT_results.txt` for reference.

p-values come from the χ² survival function (`scipy.stats.chi2.sf`), so very
small values such as 10⁻²² are reported as such instead of 0.

Multiple testing is corrected with Benjamini-Hochberg (`scipy.stats.false_discovery_control`)
within each pair of models, across the genes of one run. Genes split over several
runs form separate families, so a methods section should say how genes were
grouped.

M8 can beat M7 only because some sites are neutral (ω = 1), which the M7 beta
cannot fit; M8a vs M8 does not have this problem. The results panel warns when
M7 vs M8 is significant and M8a vs M8 is not, and the log warns when M8 runs
without M8a.

## Sites under selection

Sites come from the Bayes Empirical Bayes section (BEB; Yang, Wong & Nielsen
2005) of M2a, M8 and Branch-site, or from NEB when there is no BEB. `*` marks
Pr(ω>1) ≥ 0.95 and `**` ≥ 0.99, as in codeml.

With `cleandata = 1` codeml numbers sites after removing columns. SAGUI
applies the same rule (a column is removed if any sequence has something other
than A, C, G, T there, or a stop codon), checks the result against the number of
sites codeml reports (`ls` in its output), and writes `GENE_MODEL_sitemap.json`.
Tables show both the alignment position and the codeml position. If the check
fails, only the codeml position is shown, with a warning. The number of codon
columns removed is written to the log, to `genes_status.tsv`, to `methods_text.txt`
(total over the run, and the range per gene) and to `run_config.json`
(`cleandata_removed_codons`).

A site is listed whatever the residue of the first sequence at that position
(with `cleandata = 0` it can be `-`, `?` or `*`).

"ω and p₁ of the positive class" are the last class of the "MLEs of dN/dS (w) for
site classes" table of M2a or M8. The gene's mean ω is not used as a criterion.
The ω column of `analysis_summary.tsv` is the dN/dS of the "dN & dS for each branch"
table (for M1a, M2a, M7, M8 and M8a, the mean over site classes); 999 is codeml's
value when dS = 0 and is kept.

## Options for many genes

`--skip-beb` stops M2a and M8 when codeml starts BEB. lnL, np and ω, which the LRT
uses, are already written by then; the site table falls back to NEB.

`--two-pass` runs every gene with `--skip-beb`, applies the LRT and BH correction,
then reruns with BEB only the genes with q below `--sig-threshold` (0.05). Outputs
go to `OUT/pass1_screen/` and `OUT/pass2_beb/`.

Warm start (`--warm-start-m0`, or "Start from M0" in the advanced settings) is off by
default: every model is fitted from scratch, also when M0 is selected. When it is on,
M0 is fitted first for each gene (also when it is not selected) and its κ and branch
lengths are the starting values of the other models (`fix_blength = 1`); each model is
fitted from initial ω 0.2, 1.0 and 2.5 (M7, M8 and M8a: more starts up to
`--beta-starts`) and the best lnL is kept. The κ, `fix_blength` and initial ω of the
kept start of each gene and model are in `run_config.json` (`starts_by_gene`). Results are not
guaranteed to match a fit from scratch. `methods_text.txt` says when it was used.

## Evidence column in the results panel

Each gene gets one entry in "Evidence (all tests)" from its positive-selection tests (q < 0.05):

| Evidence | When |
|---|---|
| Supported | M8 vs M8a significant with at least one M8 BEB site at Pr(ω>1) ≥ 0.95, or M2a vs M1a significant with at least one M2a BEB site at Pr(ω>1) ≥ 0.95 |
| Not confirmed by M8a | M8 vs M7 significant but M8 vs M8a not (and M2a vs M1a not) |
| Possible (no M8a) | M8 vs M7 significant and M8a not run (and M2a vs M1a not) |
| Weak signal | a test significant, but none of the above (M8 sites do not count through M8 vs M7) |
| Not detected | no test significant; this is not evidence that the gene is free of selection |

The rows are checked from the top. A test marked as an optimization failure counts
as not run.

"Copy significant genes" copies the genes marked "Supported" for the selected
site test, or the genes significant in the selected branch or branch-site test.

## Citing

Each run writes `methods_text.txt`: a methods paragraph filled in with that run's
SAGUI version and commit, codeml version, models, parameters, tests with df,
and the number of genes in each BH correction. The commit is also shown under
About, by `sagui-cli.sh --version`, and in `run_config.json` (`sagui_commit`).

## References

- Yang Z (2007) PAML 4: phylogenetic analysis by maximum likelihood. *Mol Biol Evol* 24:1586–1591.
- Yang Z, Nielsen R, Goldman N, Pedersen A-MK (2000) Codon-substitution models for
  heterogeneous selection pressure at amino acid sites. *Genetics* 155:431–449. (M7, M8)
- Swanson WJ, Nielsen R, Yang Q (2003) Pervasive adaptive evolution in mammalian
  fertilization proteins. *Mol Biol Evol* 20:18–20. (M8a)
- Wong WSW, Yang Z, Goldman N, Nielsen R (2004) Accuracy and power of statistical
  methods for detecting adaptive evolution in protein coding sequences. *Genetics*
  168:1041–1051. (M1a, M2a, M8a)
- Yang Z, Wong WSW, Nielsen R (2005) Bayes empirical Bayes inference of amino acid
  sites under positive selection. *Mol Biol Evol* 22:1107–1118. (BEB)
- Zhang J, Nielsen R, Yang Z (2005) Evaluation of an improved branch-site likelihood
  method for detecting positive selection at the molecular level. *Mol Biol Evol*
  22:2472–2479. (branch-site)
- Self SG, Liang K-Y (1987) *J Am Stat Assoc* 82:605–610.
- Benjamini Y, Hochberg Y (1995) *J R Stat Soc B* 57:289–300.
