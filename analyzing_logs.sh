container_name="pyforge-python-school-3_web1"
sleep 10

container_hash=$(docker ps -aqf "name=$container_name")

logs=$(docker logs $container_hash)
echo "$logs" 
echo "$logs" | grep "10 passed" || exit 1

