import os
import time
import pandas as pd
from Bio import Entrez

# Configure Entrez
Entrez.email = os.environ.get('ENTREZ_EMAIL', 'your_email@example.com')  # Replace with your email or set as env var
TOOL_NAME = 'BioTextClassifier'

# Search parameters
SEARCH_TERM = 'neurological disorders OR behavioral disorders'
MAX_COUNT = 1000
BATCH_SIZE = 100


def fetch_pubmed_ids(term, max_count):
    handle = Entrez.esearch(db='pubmed', term=term, retmax=max_count, tool=TOOL_NAME, email=Entrez.email)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']  # type: ignore


def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='xml', tool=TOOL_NAME, email=Entrez.email)
    records = Entrez.read(handle)
    handle.close()
    return records['PubmedArticle']


def parse_article(article):
    medline = article['MedlineCitation']
    article_data = medline['Article']
    pmid = str(medline['PMID'])
    title = article_data.get('ArticleTitle', '')
    abstract = ''
    if 'Abstract' in article_data and 'AbstractText' in article_data['Abstract']:
        abstract = ' '.join(article_data['Abstract']['AbstractText'])
    journal = article_data['Journal']['Title'] if 'Journal' in article_data else ''
    year = ''
    if 'Journal' in article_data and 'JournalIssue' in article_data['Journal']:
        year = article_data['Journal']['JournalIssue'].get('PubDate', {}).get('Year', '')
    authors = []
    if 'AuthorList' in article_data:
        for author in article_data['AuthorList']:
            if 'LastName' in author and 'ForeName' in author:
                authors.append(f"{author['ForeName']} {author['LastName']}")
    return {
        'PMID': pmid,
        'Title': title,
        'Abstract': abstract,
        'Journal': journal,
        'Year': year,
        'Authors': '; '.join(authors)
    }


def main():
    print(f"Fetching up to {MAX_COUNT} PubMed abstracts for: {SEARCH_TERM}")
    ids = fetch_pubmed_ids(SEARCH_TERM, MAX_COUNT)
    print(f"Found {len(ids)} articles.")
    all_articles = []
    for start in range(0, len(ids), BATCH_SIZE):
        batch_ids = ids[start:start+BATCH_SIZE]
        print(f"Fetching details for records {start+1} to {start+len(batch_ids)}...")
        articles = fetch_details(batch_ids)
        for article in articles:
            all_articles.append(parse_article(article))
        time.sleep(0.5)  # Be polite to NCBI servers
    df = pd.DataFrame(all_articles)
    df.to_csv('data_ingestion/pubmed_abstracts.csv', index=False)
    print(f"Saved {len(df)} abstracts to data_ingestion/pubmed_abstracts.csv")


if __name__ == '__main__':
    main() 