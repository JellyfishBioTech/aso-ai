import os
from google.cloud import bigquery
import db_dtypes

# Step 1: Install required dependencies (Run these manually in the terminal if needed)
# Uncomment and run these in a separate terminal if installation issues persist.
# !pip install google-cloud-bigquery
# !pip install --no-cache-dir --only-binary ":all:" --upgrade pyarrow
# !pip install --upgrade db-dtypes

# If using a Mac with an M1/M2 chip, make sure required system dependencies are installed.
# Run these in the terminal, not in Python:
# arch -arm64 brew install cmake ninja llvm boost rapidjson

# Step 2: Set Google Cloud credentials
GOOGLE_CREDENTIALS_PATH = "your_path_to_the_creds.json"
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = GOOGLE_CREDENTIALS_PATH

# Step 3: Initialize BigQuery client
client = bigquery.Client()

# Step 4: Verify authentication by printing the active project
print(f"Authenticated with project: {client.project}")

# Step 5: Run a sample BigQuery query (Retrieves 10 rows from the patents-public-data dataset)
query = """
SELECT
  *
FROM
  `patents-public-data.patents.publications`
LIMIT 10
"""

try:
    df = client.query(query).to_dataframe()
    print("Query executed successfully. Here are the first few rows:")
    print(df.head())
except Exception as e:
    print(f"Error executing query: {e}")
