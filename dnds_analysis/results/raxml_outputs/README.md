# CESA dN/dS Analysis in Bryophytes using PAML

This repository contains data, scripts, and output files for estimating rates of molecular evolution in *Physcomitrium patens* and related bryophyte CESA gene sequences using the codeml program from the PAML package.

---

## Project Overview

Cellulose synthase (CESA) genes are essential for plant cell wall biosynthesis. This project investigates the evolutionary pressures acting on CesA genes across bryophyte lineages by identifying codon sites under positive selection.

---

##  Analysis Workflow

1. **Sequence Retrieval**

   - CESA coding sequences from *P. patens* and related bryophyte species were identified via BLASTX and curated for stop codons and frame shifts.

2. **Alignment**

   - Codon-based nucleotide alignments generated using [TranslatorX](http://translatorx.co.uk/), guided by peptide alignments.
   - Regions with excessive gaps were manually removed in **AliView**, preserving codon structure and equal sequence length.
   - The final aligned FASTA file was converted to **PHYLIP format (`.phy`)**, as required by **codeml**.

3. **Phylogenetic Tree Construction**

   - Trees were built with **RAxML** using the aligned amino acid sequences.

4. **dN/dS Analysis**

   - Conducted using **codeml** (PAML) with site models:
     
     - M1a (Nearly Neutral) vs M2a (Positive Selection)
     - M7 (Beta) vs M8 (Beta + ω > 1)

5. **Model Comparison**

   - Likelihood Ratio Tests (LRTs) were used to compare models and identify statistically significant positive selection.

---

### Dependencies

- [PAML 4.9](http://abacus.gene.ucl.ac.uk/software/paml.html) — tested with codeml
- Add `codeml` to `scripts/` or ensure it's in your system PATH

---

## Requirements

- PAML v4.9 or later
- RAxML
- TranslatorX
- R (for visualization/mapping scripts)

---

## Example codeml Command

```bash
codeml codeml_configs/codeml_M8.ctl
```

---

## Results


---


## Citation

If using this workflow or data, please cite:

- Yang Z. (2007). PAML 4: Phylogenetic Analysis by Maximum Likelihood.
- TranslatorX: Abascal et al. (2010)
- Our publication when available



