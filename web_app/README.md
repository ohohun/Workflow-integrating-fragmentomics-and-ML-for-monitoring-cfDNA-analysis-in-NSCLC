# cfDNA Monitor
Version: Final V.8.2

## Project Description
cfDNA Monitor is a production-grade web application designed for monitoring longitudinal circulating tumor DNA (cfDNA) in Non-Small Cell Lung Cancer (NSCLC) patients. The platform provides clinical-style analytics, including tumor fraction distribution, treatment response assessment via machine learning, and fragmentomics profiles.

## Main Technologies
- **Backend**: FastAPI (Python)
- **Frontend**: HTML5, CSS3 (Vanilla), JavaScript (ES6+)
- **Data Source**: Google Cloud Storage (GCS)
- **Deployment**: Google Cloud Run (Serverless)

## How to Run Locally
1. Install dependencies: `pip install -r requirements.txt`
2. Set environment variables (optional): `USE_LOCAL_DATA=1` for local CSV testing.
3. Start the server: `uvicorn backend.main:app --reload`
4. Access via: `http://localhost:8000`

## How to Deploy
Deploy to Cloud Run using the provided Dockerfile:
`gcloud run deploy cfdna-monitor --source . --region asia-southeast1`

## Demo Access
- **Login Credentials**:
  - Username: `admin`
  - Password: `admin1234`
- **Live URL**: [https://cfdna-monitor-118512591233.asia-southeast1.run.app/](https://cfdna-monitor-118512591233.asia-southeast1.run.app/)

## Note
“Research use only. Not for diagnostic use.”
