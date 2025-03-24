from Bio import Entrez

Entrez.email = "julia.patsyukova@gmail.com"

queries = [
    "antisense oligonucleotides",
    "ASO"
]

MAX_RESULTS = 100

results = []

def search_pmc(query):
    handle = Entrez.esearch(db="pmc", term=query, retmax=MAX_RESULTS)
    record = Entrez.read(handle)
    return record["IdList"]

def fetch_full_article(pmcid):
    url = f"https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:{pmcid}&metadataPrefix=pmc"
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

seen_ids = set()
for query in queries:
    print(f"Searching for: {query}")
    pmc_ids = search_pmc(query)
    for pmcid in tqdm(pmc_ids):
        if pmcid in seen_ids:
            continue  # skip duplicates
        article = fetch_full_article(pmcid)
        if article:
            results.append(article)
            seen_ids.add(pmcid)

df = pd.DataFrame(results)