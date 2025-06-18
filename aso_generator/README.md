# ASO Generator (`aso_generator/` folder)

This folder provides a web-based interface for generating and evaluating **Antisense Oligonucleotide (ASO) candidates** targeting a specific gene or transcript. It allows researchers to input a gene symbol or Ensembl transcript ID, customize design parameters, and receive optimized ASO sequences with key performance metrics and visualizations.

---

## Key Features

- **Gene or Transcript Input**: Accepts both gene names (e.g., `TP53`) and Ensembl transcript IDs (e.g., `ENST00000357033`) via the Ensembl REST API.
- **Target Region Selection**: Enables targeting either the entire transcript or a specific exon.
- **Custom Design Parameters**:
  - ASO length (window size)
  - GC content filtering
  - Maximum homopolymer length
  - ASO type (`SSO`, `GAPMER`)
- **Melting Temperature Calculation**: Calculates Tm using Biopython’s Nearest-Neighbor model.
- **Off-Target Screening**: Uses Bowtie2 to align candidate ASOs to a genomic fragment and estimate off-target potential.
- **Visualization Outputs**:
  - GC% distribution
  - Tm distribution
  - ASO length boxplot
  - Filtering pipeline summary
  - ASO positions on the targeted exon
- **Export**: Downloadable CSV file with full ASO candidate data.

---

## File Overview

- `app.py`  
  The main Streamlit web application. Handles user input, runs the ASO generation pipeline, and renders the interactive interface.

- `aso_generator.py`  
  Core backend class responsible for:
  - Sequence extraction from Ensembl
  - Exon parsing and indexing
  - Fragment generation with user-defined filters
  - Tm and GC% calculation
  - Bowtie2-based off-target alignment
  - Final output formatting and visualization

- `requirements.txt`
  Lists all necessary dependencies to run the ASO Generator, including:
  - streamlit
  - requests
  - pandas
  - matplotlib
  - seaborn
  - biopython  

  ---

## How It Works

1. **Sequence Retrieval**  
 Retrieves the cDNA or genomic sequence from Ensembl using the Ensembl REST API.

2. **Exon Identification**  
 Expands transcript annotations to identify exon positions and coordinates.

3. **Candidate Generation**  
 Slides windows of size 16–25nt over the selected region, filters by GC% and homopolymer rules, and generates ASOs.

4. **ASO Evaluation**  
 - Calculates melting temperature (`Tm`)
 - Aligns each ASO to a local genomic fragment using **Bowtie2** (optional if genome is downloaded)

5. **Output and Visualization**  
 Outputs candidates as a table and generates visuals to assess GC%, Tm, length, exon mapping, and filtering steps.

---

## Live Demo

You can try the tool online via Hugging Face Spaces:  
🔗 [https://huggingface.co/spaces/jpatsiukova/ASO-generator](https://huggingface.co/spaces/jpatsiukova/ASO-generator)

The dashboard was developed in **Python using Streamlit**, and deployed using Hugging Face’s hosting infrastructure. It integrates live querying of Ensembl for gene sequences and local alignment using Bowtie2.

---

## Usage

To run locally:

```bash
pip install -r requirements.txt
streamlit run app.py
```

Ensure you have `bowtie2` installed and in your `PATH` for the off-target alignment step.

---

## Planned Extensions
- Add support for multi-exon targeting  
- Integrate more chemical modification options  
- Link with ML models for predictive ASO scoring  
- Export ready-to-order oligo files  

---

## Related Modules
- [`/modelling/`](../modelling): ML models for ASO efficacy prediction  
- [`/database/`](../database): Scripts for mining ASO sequences and metadata  
- [`/data/`](../data): Streamlit dashboard visualizing final ASO dataset  

---

This tool is designed to support research in RNA-targeting therapeutics and facilitate rapid ASO design and screening.
