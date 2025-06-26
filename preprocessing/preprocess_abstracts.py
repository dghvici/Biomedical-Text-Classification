import pandas as pd
import spacy
from spacy.lang.en.stop_words import STOP_WORDS

# Load spaCy English model
nlp = spacy.load('en_core_web_sm')

def clean_text(text):
    doc = nlp(text)
    tokens = [
        token.lemma_.lower() for token in doc
        if not token.is_stop and not token.is_punct and not token.is_space
    ]
    return ' '.join(tokens)

def main():
    df = pd.read_csv('data_ingestion/pubmed_abstracts.csv')
    df['Cleaned_Abstract'] = df['Abstract'].astype(str).apply(clean_text)
    df.to_csv('preprocessing/cleaned_abstracts.csv', index=False)
    print('Saved cleaned abstracts to preprocessing/cleaned_abstracts.csv')

if __name__ == '__main__':
    main() 