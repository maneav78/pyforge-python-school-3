import boto3
import json
import os
from dotenv import load_dotenv

load_dotenv(".env")

session = boto3.Session(
    aws_access_key_id=os.getenv('AWS_ACCESS_KEY_ID'),
    aws_secret_access_key=os.getenv('AWS_SECRET_ACCESS_KEY'),
    region_name=os.getenv('AWS_REGION')  
)


lambda_client = session.client('lambda')

event = {
"name": "Student" 
}


response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event) 
)


response_payload = json.loads(response['Payload'].read())

print("Response: ", response_payload)
