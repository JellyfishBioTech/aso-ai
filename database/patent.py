import os

# Set env vars for patent_client BEFORE importing patent_client
# TODO: Get API keys from .env
os.environ["PATENT_CLIENT_EPO_API_KEY"] = "PUT_YOUR_PATENT_CLIENT_EPO_API_KEY_HERE"
os.environ["PATENT_CLIENT_EPO_API_SECRET"] = "PUT_YOUR_PATENT_CLIENT_EPO_SECRET_HERE"

# Import and initialize patent_client
from patent_client import USApplication, Inpadoc

# USApplication works without API key, but we need key for USPTO Open Data Portal
# TODO: Get API key for USPTO Open Data Portal
google_apps = USApplication.objects.filter(first_named_applicant="Google LLC")
patent = google_apps[0]

# Print the first Google patent
print(google_apps[0].appl_id, google_apps[0].patent_title)

# Print Inpadoc data for a patent
pub = Inpadoc.objects.get("EP3082535A1")
print(pub.biblio.title)
