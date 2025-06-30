import pandas as pd
import spacy
import re
from collections import Counter


nlp = spacy.load('en_ner_bc5cdr_md')

def extract_entities(text):
    """
    Extract DISEASE and CHEMICAL entities using scispaCy's en_ner_bc5cdr_md model.
    Returns a dict with lists for each entity type.
    """
    doc = nlp(text)
    entities = {'DISEASE': [], 'CHEMICAL': []}
    for ent in doc.ents:
        if ent.label_ in entities:
            entities[ent.label_].append(ent.text)
    return entities

def extract_research_type(text):
    """
    Improved keyword-based extraction for research type.
    Uses regex and expanded keyword list.
    """
    research_types = {
        'clinical_trial': [r'clinical trial', r'randomi[sz]ed', r'\brct\b', r'phase [i1]{1,3}|phase iv'],
        'review': [r'review', r'meta[- ]?analysis', r'systematic review'],
        'case_study': [r'case study', r'case report'],
        'observational': [r'observational', r'cohort', r'cross[- ]sectional']
    }
    text_lower = text.lower()
    found_types = []
    for research_type, patterns in research_types.items():
        for pattern in patterns:
            if re.search(pattern, text_lower):
                found_types.append(research_type)
                break
    return found_types

def extract_trial_phase(text):
    """
    Improved regex for clinical trial phases.
    """
    phase_patterns = {
        'phase_1': r'phase\s*[i1]\b',
        'phase_2': r'phase\s*[i1]{2}\b',
        'phase_3': r'phase\s*[i1]{3}\b',
        'phase_4': r'phase\s*iv\b|phase\s*4\b'
    }
    text_lower = text.lower()
    found_phases = []
    for phase, pattern in phase_patterns.items():
        if re.search(pattern, text_lower):
            found_phases.append(phase)
    return found_phases

def main():
    # Load cleaned abstracts
    df = pd.read_csv('preprocessing/cleaned_abstracts.csv')
    print("Extracting categories from abstracts...")
    categories = []
    for idx, row in df.iterrows():
        abstract = row['Abstract']  # Use original abstract for NER
        # Extract all named entities
        entities = extract_entities(abstract)
        # Extract research type
        research_types = extract_research_type(abstract)
        # Extract trial phases
        trial_phases = extract_trial_phase(abstract)
        # Combine all categories
        all_categories = entities['DISEASE'] + entities['CHEMICAL'] + research_types + trial_phases
        categories.append({
            'PMID': row['PMID'],
            'Categories': '; '.join(all_categories),
            'Research_Type': '; '.join(research_types),
            'Trial_Phase': '; '.join(trial_phases),
            'Diseases': '; '.join(entities['DISEASE']),
            'Chemicals': '; '.join(entities['CHEMICAL'])
        })
    # Create categories DataFrame
    categories_df = pd.DataFrame(categories)
    # Merge with original data
    result_df = pd.merge(df, categories_df, on='PMID', how='left')
    # Save results
    result_df.to_csv('preprocessing/abstracts_with_categories.csv', index=False)
    # Print summary statistics
    print(f"Processed {len(result_df)} abstracts")
    print(f"Found {len([c for c in result_df['Categories'] if c])} abstracts with categories")
    # Show most common entities
    all_entities = []
    for ents in result_df['Categories'].dropna():
        all_entities.extend(ents.split('; '))
    if all_entities:
        print("\nTop 10 most common named entities:")
        for entity, count in Counter(all_entities).most_common(10):
            print(f"  {entity}: {count}")
    print("\nSaved results to preprocessing/abstracts_with_categories.csv")

if __name__ == '__main__':
    main() 