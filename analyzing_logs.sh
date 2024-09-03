container_name="pyforge-python-school-3-web1"
sleep 7

container_hash=$(docker ps -aqf "name=$container_name")

logs=$(docker logs $container_hash)
echo "$logs" | grep "9 passed" || exit 1

