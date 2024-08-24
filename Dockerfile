FROM python:3.11-slim

WORKDIR /for_testing

COPY requirements.txt requirements.txt

RUN apt-get update && \
    apt-get install -y libpq-dev build-essential python3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt


COPY . .

CMD ["pytest", "testing/test.py"]
