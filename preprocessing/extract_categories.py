import pandas as pd
import spacy
import re
from collections import Counter

# Load spaCy English model
nlp = spacy.load('en_core_web_sm')

def extract_entities(text):
    """Extract medical entities and research-related terms from text"""
    doc = nlp(text)
    
    # Extract named entities
    entities = {
        'DISEASE': [],
        'CONDITION': [],
        'TREATMENT': [],
        'DRUG': [],
        'PROCEDURE': []
    }
    
    for ent in doc.ents:
        if ent.label_ in entities:
            entities[ent.label_].append(ent.text.lower())
    
    return entities

def extract_research_type(text):
    """Extract research type indicators"""
    research_types = {
        'clinical_trial': ['clinical trial', 'randomized', 'rct', 'phase'],
        'review': ['review', 'meta-analysis', 'systematic review'],
        'case_study': ['case study', 'case report'],
        'observational': ['observational', 'cohort', 'cross-sectional']
    }
    
    text_lower = text.lower()
    found_types = []
    
    for research_type, keywords in research_types.items():
        if any(keyword in text_lower for keyword in keywords):
            found_types.append(research_type)
    
    return found_types

def extract_trial_phase(text):
    """Extract clinical trial phases"""
    phase_patterns = {
        'phase_1': r'phase\s*[i1]',
        'phase_2': r'phase\s*[i1]{2}',
        'phase_3': r'phase\s*[i1]{3}',
        'phase_4': r'phase\s*[iv4]'
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
    
    # Extract categories
    print("Extracting categories from abstracts...")
    
    categories = []
    for idx, row in df.iterrows():
        abstract = row['Abstract']  # Use original abstract for better entity recognition
        
        # Skip rows with NaN abstracts
        if abstract != abstract or abstract is None:  # NaN check: NaN != NaN
            categories.append({
                'PMID': row['PMID'],
                'Categories': '',
                'Research_Type': '',
                'Trial_Phase': '',
                'Diseases': '',
                'Treatments': ''
            })
            continue
        
        # Extract entities
        entities = extract_entities(abstract)
        
        # Extract research type
        research_types = extract_research_type(abstract)
        
        # Extract trial phases
        trial_phases = extract_trial_phase(abstract)
        
        # Combine all categories
        all_categories = []
        for entity_type, entity_list in entities.items():
            all_categories.extend(entity_list)
        all_categories.extend(research_types)
        all_categories.extend(trial_phases)
        
        categories.append({
            'PMID': row['PMID'],
            'Categories': '; '.join(all_categories),
            'Research_Type': '; '.join(research_types),
            'Trial_Phase': '; '.join(trial_phases),
            'Diseases': '; '.join(entities['DISEASE'] + entities['CONDITION']),
            'Treatments': '; '.join(entities['TREATMENT'] + entities['DRUG'])
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
    
    # Show most common categories
    all_categories = []
    for cats in result_df['Categories'].dropna():
        all_categories.extend(cats.split('; '))
    
    if all_categories:
        print("\nTop 10 most common categories:")
        for category, count in Counter(all_categories).most_common(10):
            print(f"  {category}: {count}")
    
    print("\nSaved results to preprocessing/abstracts_with_categories.csv")

if __name__ == '__main__':
    main() 