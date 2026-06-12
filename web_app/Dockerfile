# Use slim python image for efficiency
FROM python:3.9-slim

# Prevent python from generating .pyc files
ENV PYTHONDONTWRITEBYTECODE=1
# Ensure stdout/stderr are sent directly to terminal (Cloud Run logging)
ENV PYTHONUNBUFFERED=1
# Force local runtime mode to use bundled data
ENV USE_LOCAL_DATA=1

WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy backend and frontend code
COPY backend ./backend
COPY frontend ./frontend

# Expose the default Cloud Run port
EXPOSE 8080

# Start command: Point to the backend.main module
# Cloud Run automatically sets the PORT environment variable
CMD ["sh", "-c", "uvicorn backend.main:app --host 0.0.0.0 --port ${PORT:-8080}"]
