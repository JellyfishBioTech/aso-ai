# ASO-AI: A Platform for Antisense Oligonucleotide Discovery and Design

ASO-AI is an end-to-end computational framework for discovering and evaluating antisense oligonucleotides (ASOs) for therapeutic applications. The platform integrates data mining, machine learning, and interactive tools to support both exploratory analysis and gene-specific ASO design.

Antisense oligonucleotides (ASOs) are short, synthetic RNA-like molecules that bind to specific mRNA targets to modify gene expression. They hold significant promise for treating genetic disorders, cancer, and other complex diseases. This project was created to streamline the ASO development process and support researchers working in RNA-targeting therapeutics.

---

## Project Overview

The project includes three core components:

### 1. Data Collection & Integration

Data was mined from multiple structured and unstructured sources, including:
- [eSkip-Finder](https://www.eskip-finder.org/)
- LNA gapmer studies
- siRNA-Mod database
- PubMed/PMC full texts
- Lens.org and Google Patents

Extracted data includes ASO sequences, target genes, modifications, delivery methods, cell models, and experimental efficacy metrics. Large language models (LLMs) were used to extract structured metadata from textual sources such as abstracts, full articles, and patent descriptions.

All data was cleaned, standardized, and merged into a unified dataset in CSV format.

---

### 2. Machine Learning Modelling

With the consolidated dataset, predictive models were developed to identify sequence and contextual features that influence ASO efficacy. Techniques applied include:
- Preprocessing and normalization
- Classical ML models (Random Forest, XGBoost)
- Feature importance and explainability (e.g., SHAP)
- Evaluation and optimization of regression performance

Modelling notebooks are located in the [`modelling/`](modelling) folder.

---

### 3. ASO Generator Tool

A web-based generator allows users to input either a gene name or Ensembl transcript ID and receive optimized ASO candidate sequences. The app automatically:
- Queries Ensembl to fetch mRNA sequence and exon boundaries
- Filters candidates based on window size, GC content, and homopolymer rules
- Calculates melting temperature (Tm)
- Runs Bowtie2 to assess off-target alignment
- Outputs a table of ranked ASO candidates with genomic coordinates and quality metrics
- Generates interactive plots, including GC content distribution, melting temperature histogram, and exon-level mapping of ASO positions

The generator is built with Python and Streamlit and can be launched locally or deployed remotely.
Related ASO design scripts and interface code are available in the [`aso_generator/`](aso_generator) folder.

---

## Interactive Dashboards

An exploratory dashboard was developed to help visualize and filter the final ASO dataset across key parameters such as:
- Efficiency vs. concentration
- Target genes and modifications
- Distribution by source (journal vs patent)

The full dashboard code is available in the [`data/`](data) folder.

---

## Repository Structure

```
aso-ai/
│
├── database/             # Scripts for data mining from Lens.org, PubMed, Google Patents
│
├── data/                 # Streamlit dashboard code for dataset exploration
│
├── modelling/            # ML notebooks for training efficacy prediction models
│
├── aso_generator/        # ASO design tool: backend + Streamlit UI
│
├── .env.example          # Environment variable template
├── requirements.txt      # Common dependencies
└── README.md             # Project overview (this file)
```

---

## Requirements & Setup

- Python 3.8+
- Streamlit
- Biopython
- Bowtie2 (for off-target alignment)

To install dependencies and launch apps:
```bash
pip install -r requirements.txt
streamlit run aso_generator/app.py
```

## Planned Extensions

- Add support for multi-exon and intron-spanning ASOs  
- Integrate RNA secondary structure accessibility into design  
- Link generator to predictive ML models  
- Enable export of ready-to-order oligonucleotide sequences  

---

For feedback or collaboration, feel free to contact via GitHub.

---
