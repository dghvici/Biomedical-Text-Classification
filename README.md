# Biomedical Text Classification Pipeline

This project implements a scalable ML pipeline for classifying biomedical research abstracts (e.g., from PubMed) into categories such as treatment type, disease area, or trial phase. The pipeline is designed for cloud-native deployment (Azure-ready).

## Project Structure
- `data_ingestion/`: Scripts to fetch and preprocess data from PubMed
- `preprocessing/`: Text cleaning and feature engineering
- `models/`: Model training and evaluation (TensorFlow + scikit-learn)
- `api/`: Model serving (FastAPI/Flask)
- `docker/`: Dockerfiles and deployment scripts
- `k8s/`: Kubernetes manifests

## Getting Started
1. Install dependencies: `pip install -r requirements.txt`
2. Run data ingestion: `python data_ingestion/fetch_pubmed.py`

## Cloud Deployment
- Designed for Azure ML and AKS, but can be adapted to AWS or GCP.

## Disclaimer
This project uses data from NCBI PubMed via the Entrez API. Please review the [NCBI Disclaimer and Copyright notice](https://www.ncbi.nlm.nih.gov/About/disclaimer.html) before using or redistributing any data. Abstracts in PubMed may be protected by copyright. All users are expected to adhere to the terms and conditions asserted by the copyright holder. This project is for research and educational purposes only.

## Project Logic & Data Flow

- **Entity Extraction (NER):**
  - Use the original abstracts for named entity recognition (NER) and category extraction.
  - Cleaning (lemmatization, stopword removal) can remove important context and harm NER performance.
- **Modeling/Classification:**
  - Use the cleaned abstracts (lowercased, lemmatized, stopwords removed) for feature engineering and as input to ML models.
- **Conclusion:**
  - **Keep both versions:**
    - Use the original abstracts for NER/category extraction.
    - Use the cleaned abstracts for modeling/classification.

This approach ensures optimal performance for both entity extraction and downstream machine learning tasks. 