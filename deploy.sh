#!/bin/bash
sudo apt-get update
sudo apt-get install docker.io -y

sudo curl -L "https://github.com/docker/compose/releases/download/v2.21.0/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose

docker compose up --build -d

