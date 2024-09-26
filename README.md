# AWS EC2 Deployment Guide

## Prerequisites

1. **AWS EC2 Instance**: Ensure you have an EC2 instance set up.
2. **GitHub Repository**: Your project must be hosted in a GitHub repository.
3. **Docker**: Your application should be containerized using Docker.
4. **GitHub Secrets**: You will need to store sensitive credentials in GitHub Secrets for secure access. The required secrets are:

- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`
- `SSH_HOST`
- `SSH_USERNAME`
- `SSH_PRIVATE_KEY`

## Steps for Deployment

### 1. Configure AWS Credentials

I have used the `aws-actions/configure-aws-credentials@v3` action to configure AWS credentials within GitHub Actions. This allows secure access to AWS resources from GitHub.

```yaml
- name: Configuration
   uses: aws-actions/configure-aws-credentials@v3
   with:
      aws-access-key-id: ${{ secrets.AWS_ACCESS_KEY_ID }}
      aws-secret-access-key: ${{ secrets.AWS_SECRET_ACCESS_KEY }}
      aws-region: us-east-1
```

### 2. SSH into EC2 Instance and Transfer Files

GitHub Actions will establish an SSH connection to your EC2 instance and transfer the necessary files using SCP (Secure Copy).

```yaml
- name: Copy files to EC2    
   run: |
      scp -o StrictHostKeyChecking=no -r * ${{ secrets.SSH_USERNAME }}@${{ secrets.SSH_HOST }}:/home/${{ secrets.SSH_USERNAME }}/app

```

### 3. Deploy the Project on EC2

The deployment is handled by the `deploy.sh` script. This script prepares the EC2 instance by setting up the environment and running Docker Compose to start the application.
