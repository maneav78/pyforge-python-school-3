#!/bin/sh
set -e

echo "Starting deploy.sh..."

# Directly check for Docker in /usr/bin/docker
if ! [ -x "/usr/bin/docker" ]; then
    echo "Docker not found in /usr/bin. Installing Docker…"
    curl -fsSL https://get.docker.com -o install-docker.sh
    sudo sh install-docker.sh
    sudo usermod -aG docker $USER
    newgrp docker
else
    echo "Docker is already installed in /usr/bin. Skipping installation."
fi

cd /home/ubuntu/app || { echo "App directory not found. Exiting."; exit 1; }

sudo sed -i '/DB_URL/d' /etc/environment
sudo sed -i '/DB_USER/d' /etc/environment
sudo sed -i '/DB_PASSWORD/d' /etc/environment
sudo sed -i '/DB_NAME/d' /etc/environment
sudo sed -i '/MAIN_DATABASE_URL/d' /etc/environment
sudo sed -i '/TEST_DB_NAME/d' /etc/environment
sudo sed -i '/TEST_DATABASE_URL/d' /etc/environment

echo "Setting environment variables..."
echo "DB_URL='${DB_URL}'" | sudo tee -a /etc/environment
echo "DB_USER='${DB_USER}'" | sudo tee -a /etc/environment
echo "DB_PASSWORD='${DB_PASSWORD}'" | sudo tee -a /etc/environment
echo "DB_NAME='${DB_NAME}'" | sudo tee -a /etc/environment
echo "MAIN_DATABASE_URL='${MAIN_DATABASE_URL}'" | sudo tee -a /etc/environment
echo "TEST_DB_NAME='${TEST_DB_NAME}'" | sudo tee -a /etc/environment
echo "TEST_DATABASE_URL='${TEST_DATABASE_URL}'" | sudo tee -a /etc/environment

. /etc/environment

if [ ! -f docker-compose.yaml ]; then
    echo "docker-compose.yml not found in the current directory. Exiting."
    exit 1
fi

running_containers=$(/usr/bin/docker ps -q)

if [ -n "$running_containers" ]; then
    echo "Running containers found. Restarting with --force-recreate…"
    /usr/bin/docker compose down || { echo "Failed to stop containers. Exiting."; exit 1; }
    /usr/bin/docker compose up --build --force-recreate -d || { echo "Failed to restart containers. Exiting."; exit 1; }
else
    echo "No running containers found. Starting containers…"
    /usr/bin/docker compose up --build -d || { echo "Failed to start containers. Exiting."; exit 1; }
fi

echo "Waiting for container to initialize…"
sleep 10

/usr/bin/docker ps
echo "Deployment script completed."

