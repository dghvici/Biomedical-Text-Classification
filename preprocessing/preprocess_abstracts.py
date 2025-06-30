import pandas as pd
import spacy
from spacy.lang.en.stop_words import STOP_WORDS

# Load spaCy English model
nlp = spacy.load('en_core_web_sm')

def clean_text(text):
    # Handle empty or None text
    if not text or pd.isna(text) or text.strip() == '':
        return ''
    
    doc = nlp(text)
    tokens = [
        token.lemma_.lower() for token in doc
        if not token.is_stop and not token.is_punct and not token.is_space
    ]
    return ' '.join(tokens)

def main():
    df = pd.read_csv('data_ingestion/pubmed_abstracts.csv')
    
    # Check for NaN or empty abstracts before processing
    nan_count = df['Abstract'].isna().sum()
    empty_count = (df['Abstract'] == '').sum()
    
    print(f"Processing {len(df)} articles...")
    print(f"Found {nan_count} NaN abstracts and {empty_count} empty abstracts")
    
    # Filter out articles without abstracts
    df_valid = df[df['Abstract'].notna() & (df['Abstract'] != '')].copy()
    print(f"Proceeding with {len(df_valid)} articles with valid abstracts")
    
    # Clean the abstracts
    df_valid['Cleaned_Abstract'] = df_valid['Abstract'].apply(clean_text)  # type: ignore
    
    # Save the cleaned data
    df_valid.to_csv('preprocessing/cleaned_abstracts.csv', index=False)
    print(f'Saved {len(df_valid)} cleaned abstracts to preprocessing/cleaned_abstracts.csv')
    
    # Verify no NaN values in cleaned data
    final_nan_count = df_valid['Cleaned_Abstract'].isna().sum()  # type: ignore
    if final_nan_count == 0:
        print("✅ SUCCESS: No NaN values in cleaned abstracts!")
    else:
        print(f"⚠️  WARNING: {final_nan_count} NaN values found in cleaned abstracts")

if __name__ == '__main__':
    main() 