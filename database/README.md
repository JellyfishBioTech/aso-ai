# database/

This folder contains all scripts and notebooks used for **data collection, extraction, validation, and integration** of antisense oligonucleotide (ASO)-related information from various public resources. It represents the first stage of the project: **creating a unified structured dataset** to be used for exploratory analysis, modeling, and ASO design.

---

## Folder Overview

### `google_bigquery/`
Contains scripts to extract ASO-relevant patents from the [Google Patents Public Datasets](https://console.cloud.google.com/marketplace/product/google_patents_public_datasets/patents-public-data) via BigQuery.

- `bigquery_google_setup.py` — Google Cloud authentication and BigQuery client setup.
- `google_patent.ipynb` — Sample notebook for querying patent metadata and full-text related to ASOs.

---

### `lens_org/`
Gathers patent data from [Lens.org](https://www.lens.org/) and uses OpenAI’s language models to extract structured ASO-related information.

- `lens_org_api_query_openai.ipynb` — Queries Lens API, retrieves abstract/claim/description text, and extracts features such as gene target, sequence, efficiency using GPT.
- `README.md` — Describes Lens API usage and prompt structure.

---

### `pubmed/`
Extracts full-text and metadata from PubMed Central (PMC) articles mentioning ASOs using the NCBI Entrez API.

- `pubmed_full_texts_extract_data.py` — Parses abstracts and body text from PMC XML responses using BeautifulSoup.
- `pubmed_validating_data.py` — Validates and standardizes extracted PubMed data.
- `test.ipynb` — Optional/test notebook for prototyping.
- `README.md` — Description of the PubMed extraction pipeline.

---

### Other Files

- `patent.py` — Integration script that processes and merges ASO-relevant data from patent databases.
- `.env.example` — Template for storing API keys securely.
- `requirements.txt` — List of dependencies to recreate the environment.
- `README.md` — This file.

---

## Output

The output from all subfolders is merged into a **single structured CSV file** containing:
- Target gene and exon/intron location
- Oligo sequence and chemical modification
- Efficiency data (e.g., inhibition %, IC50)
- Source type (journal, patent, clinical study)
- Delivery method, species, and cell line
- Link to original data

This dataset is used in downstream stages of the project (e.g., ML modeling, dashboarding, and ASO generation).

---

## API Access Setup

Some scripts require API access:

- `Lens.org` — `LENS_API_KEY`
- `OpenAI` — `OPENAI_API_KEY`
- `NCBI Entrez` — `ENTREZ_EMAIL`

Store your credentials in a `.env` file based on the `.env.example` provided.

```bash
LENS_API_KEY=your_lens_api_key
OPENAI_API_KEY=your_openai_key
ENTREZ_EMAIL=your_email@example.com
