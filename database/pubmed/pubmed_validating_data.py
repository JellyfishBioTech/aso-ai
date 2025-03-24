import pandas as pd
import requests
import time
from bs4 import BeautifulSoup
from Bio import Entrez
from tqdm import tqdm

Entrez.email = "julia.patsyukova@gmail.com"

df = pd.read_excel("/Users/julia_patsiukova/Downloads/output_table_FINAL.xlsx")
df['PMID'] = df['id'].str.extract(r'PMID\s*(\d+)')
pmid_list = df['PMID'].dropna().unique().tolist()

not_found_pmids = []

def get_pmcid_from_pmid(pmid):
    try:
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        record = Entrez.read(handle)
        for link in record[0].get('LinkSetDb', []):
            if link['LinkName'] == 'pubmed_pmc':
                return link['Link'][0]['Id']
    except:
        pass
    return None

def fetch_full_article(pmcid):
    url = f"https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:{pmcid}&metadataPrefix=pmc"
    try:
        response = requests.get(url)
        if response.status_code != 200:
            return None
        soup = BeautifulSoup(response.content, "xml")

        title = soup.find("article-title")
        abstract = soup.find("abstract")
        authors = soup.find_all("contrib", {"contrib-type": "author"})
        date = soup.find("pub-date")
        body = soup.find("body")

        title_text = title.get_text(separator=" ") if title else ""
        abstract_text = abstract.get_text(separator=" ") if abstract else ""
        body_text = body.get_text(separator="\n") if body else ""
        author_list = [a.find("surname").text + " " + a.find("given-names").text
                       for a in authors if a.find("surname") and a.find("given-names")]
        pub_date = ""
        if date:
            year = date.find("year").text if date.find("year") else ""
            month = date.find("month").text if date.find("month") else ""
            day = date.find("day").text if date.find("day") else ""
            pub_date = f"{year}-{month}-{day}".strip("-")

        return {
            "PMCID": f"PMC{pmcid}",
            "Author": ", ".join(author_list),
            "Name of the article": title_text,
            "Abstract": abstract_text,
            "Date": pub_date,
            "Full text": body_text,
        }
    except:
        return None

results = []
seen_pmcs = set()

for pmid in tqdm(pmid_list):
    pmcid = get_pmcid_from_pmid(pmid)
    if not pmcid or pmcid in seen_pmcs:
        not_found_pmids.append(pmid)
        continue
    article = fetch_full_article(pmcid)
    if article:
        results.append(article)
        seen_pmcs.add(pmcid)
    else:
        not_found_pmids.append(pmid)
    time.sleep(0.5)


df_results = pd.DataFrame(results)