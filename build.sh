set -ex

# SET THE FOLLOWING VARIABLES
# docker hub username
USERNAME=mattjvincent
# image name
IMAGE=cc_proteins

docker build -t $USERNAME/$IMAGE:latest .



