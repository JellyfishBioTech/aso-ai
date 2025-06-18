# ASO Modelling (`modelling/` folder)

This folder contains machine learning notebooks focused on modeling and predicting the efficacy of antisense oligonucleotides (ASOs) using structured data collected from scientific publications, patents, and databases.

---

## Purpose

The goal of this part of the project was to analyze the performance of ASOs and explore factors affecting their biological activity. The notebooks cover preprocessing, feature engineering, model training, and evaluation for predicting ASO efficiency and off-target behavior.

---

## Contents

- **`Comprehensive_ASO_Analysis_Complete.ipynb`**  
  A full exploratory and modeling pipeline notebook including:
  - Data cleaning and preprocessing
  - Feature engineering (modifications, GC%, melting temperature, etc.)
  - Exploratory Data Analysis (EDA) with visualizations
  - Machine learning models (Random Forest, Logistic Regression, XGBoost)
  - Model evaluation and interpretability

- **`ASO_ML_Analysis_Notebook_Fixed.ipynb`**  
  A streamlined and optimized version of the modeling notebook. Focuses on:
  - Performance comparison of models
  - Parameter tuning
  - Handling categorical and continuous features
  - Export-ready outputs for further integration

---

## Modeling Objectives

- Predict **ASO efficiency** (e.g., % inhibition or IC50)
- Explore the impact of:
  - ASO length and GC content
  - Chemical modifications (e.g., LNA, PMO)
  - Target gene region (exon/intron)
  - Organism and delivery method

---

## Techniques Used

- Random Forest Regressor & Classifier
- XGBoost
- Logistic Regression
- GridSearchCV for hyperparameter tuning
- SHAP (optional) for model interpretation
- Data encoding and imputation
- Train/test splitting and cross-validation

---

## Requirements

Notebooks require the following Python libraries:

```bash
pandas
numpy
scikit-learn
xgboost
matplotlib
seaborn
pip install -r requirements.txt
```
