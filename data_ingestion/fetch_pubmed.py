import os
import time
import pandas as pd
from Bio import Entrez

# Configure Entrez
Entrez.email = os.environ.get('ENTREZ_EMAIL', 'your_email@example.com')  # Replace with your email or set as env var
TOOL_NAME = 'BioTextClassifier'

# Search parameters
SEARCH_TERM = (
    '((neurological disorder OR neurological disease OR behavioral disorder OR behavioral disease) '
    'AND ("Phase II" OR "efficacy of" OR "therapeutic impact" OR "clinical response" OR '
    '"treatment outcome" OR "drug therapy" OR "side effects" OR "mechanism of action")) '
    'AND ("Therapeutic Uses"[MeSH Terms] OR "Drug Therapy"[MeSH Terms] OR "Treatment Outcome"[MeSH Terms]) '
    'AND (journal:"Nature Reviews Neuroscience" OR journal:"Neuron" OR journal:"Annual Review of Neuroscience" '
    'OR journal:"Trends in Neurosciences" OR journal:"Nature Neuroscience" OR journal:"Brain" '
    'OR journal:"Molecular Psychiatry" OR journal:"Biological Psychiatry" OR journal:"Progress in Neurobiology" '
    'OR journal:"Neuropsychopharmacology")'
)
MAX_COUNT = 2000
BATCH_SIZE = 100
SUPPLEMENT_MARGIN = 50  # Extra articles to fetch in case of NaN values


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


def fetch_with_supplementation(term, target_count):
    """Fetch articles with supplementation to handle NaN abstracts"""
    print(f"Fetching up to {target_count + SUPPLEMENT_MARGIN} PubMed abstracts for: {term}")
    
    # Fetch initial batch with extra margin
    ids = fetch_pubmed_ids(term, target_count + SUPPLEMENT_MARGIN)
    print(f"Found {len(ids)} articles.")
    
    all_articles = []
    valid_articles = 0
    
    for start in range(0, len(ids), BATCH_SIZE):
        batch_ids = ids[start:start+BATCH_SIZE]
        print(f"Fetching details for records {start+1} to {start+len(batch_ids)}...")
        articles = fetch_details(batch_ids)
        
        for article in articles:
            parsed_article = parse_article(article)
            all_articles.append(parsed_article)
            
            # Count valid articles (those with abstracts)
            if parsed_article['Abstract'] and parsed_article['Abstract'].strip():
                valid_articles += 1
                
                # If we have enough valid articles, we can stop
                if valid_articles >= target_count:
                    print(f"Reached target of {target_count} valid articles. Stopping fetch.")
                    break
        
        time.sleep(0.5)  # Be polite to NCBI servers
        
        # Check if we have enough valid articles
        if valid_articles >= target_count:
            break
    
    # Filter to keep only articles with valid abstracts
    valid_articles_list = [article for article in all_articles if article['Abstract'] and article['Abstract'].strip()]
    
    print(f"Fetched {len(all_articles)} total articles")
    print(f"Found {len(valid_articles_list)} articles with valid abstracts")
    
    # If we still don't have enough, try to fetch more
    if len(valid_articles_list) < target_count:
        print(f"Warning: Only found {len(valid_articles_list)} valid articles out of {target_count} requested")
    
    return valid_articles_list[:target_count]


def main():
    # Fetch articles with supplementation
    all_articles = fetch_with_supplementation(SEARCH_TERM, MAX_COUNT)
    
    df = pd.DataFrame(all_articles)
    df.to_csv('data_ingestion/pubmed_abstracts.csv', index=False)
    print(f"Saved {len(df)} abstracts to data_ingestion/pubmed_abstracts.csv")
    
    # Verify no NaN abstracts
    nan_count = df['Abstract'].isna().sum()
    if nan_count > 0:
        print(f"Warning: Still found {nan_count} articles with NaN abstracts")
    else:
        print("Success: All articles have valid abstracts!")


if __name__ == '__main__':
    main() 