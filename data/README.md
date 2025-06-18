# ASO Data Dashboard (`data/` folder)

This folder contains the source code for an **interactive dashboard** that visualizes a comprehensive dataset of antisense oligonucleotides (ASOs). The dashboard was built using **Python and Streamlit**, and is currently deployed on **Hugging Face Spaces** for public access and demonstration.

---

## Live Dashboard

You can explore the live version of the dashboard here:

**[https://huggingface.co/spaces/jpatsiukova/ASO-dashboard](https://huggingface.co/spaces/jpatsiukova/ASO-dashboard)**

> Hosted on Hugging Face Spaces, this deployment allows users to interactively filter and visualize the ASO dataset without any local installation.

---

## Purpose

The dashboard enables researchers and developers to **interactively explore the ASO dataset** collected from scientific articles, patents, and clinical sources. It helps to:

- Identify trends in ASO efficiency across species and gene targets
- Compare experimental data by modification type, delivery method, and concentration
- Visualize the structure and scale of the underlying dataset
- Facilitate hypothesis generation for further modeling

---

## Dataset Overview

The visualized dataset combines entries from:
- **eSkip-Finder**
- **Gapmer**
- **siRNA-Mod**
- Additional extracted data from **PubMed**, **Lens.org**, and **Google Patents**

Each ASO record contains metadata such as:
- Target gene and exon/intron region
- Organism and cell model
- ASO type and chemical modifications
- Concentration and inhibition efficiency
- Experimental context (e.g., delivery method, source type)

---

## Files

- `app.py` – Main Streamlit app script for rendering the dashboard
- `requirements.txt` – Python dependencies needed to run the app locally

---

## How to Run Locally

To run the dashboard on your machine:

1. Clone the repository and navigate to the `data/` directory:

   ```bash
   git clone https://github.com/<your_username>/aso-ai.git
   cd aso-ai/data
   pip install -r requirements.txt
   streamlit run app.py
   ```
