FROM python:3.11-slim

WORKDIR /for_testing

COPY requirements.txt requirements.txt

RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["pytest", "testing/test.py"]