# Example script how an image was built

cd /Users/iosubi/Documents/docker/env3
ls
# [out] Dockerfile	env3.yml
docker build -t iraiosub/nf-riboseq .
docker push iraiosub/nf-riboseq