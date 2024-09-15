# Docker alternative

This docker file creates a new ubuntu image and installs super-focus and mmseqs2. Then there is a directory called `/superfocus` where you can mount your data.

This image is particularly good for singularity!

# Building this image

### Cleaning up before you start
> If you have been using docker recently, you might want to do a complete nuke
> of the docker directory and start again. This is the easy way
> and also updates docker!

```
sudo -H bash
apt remove -y docker.io && rm -rf /var/lib/docker && apt install -y docker.io
```

### Remove any existing docker images

```
docker images
docker images | awk '{print $3}' | grep -v IM | xargs docker rmi -f {}
```

### build a new image
docker build -t superfocus . > build.out 2> build.err

### tag it and upload to docker quay
```
docker login -u $USER
TAG=37c06e3a9f8a     ## ADD THE BUILD MD5SUM here
docker tag superfocus $USER/superfocus:v1.6.0_${TAG}

docker push $USER/superfocus:v1.6.0_${TAG}
```

#  Running the image

You can run the image to see what's there

```
docker run -i -t $USER/superfocus:v1.6.0_${TAG}  /bin/bash
```



# Running Superfocus with singularity

We can use this dockerfile and convert it  into a `.sif` image.

The current tag is 5abbef95a6ca


# Load the singularity module and create a new image


NOTE: There is (currently) an issue with singularity 4.1.0. If it persists, here is how to load a previous image:

```
module load pawseyenv/2023.08
module load singularity/3.11.4-slurm
```

Otherwise, this will make a superfocus `.sif`

```
export TAG=5abbef95a6ca
module load singularity/4.1.0-slurm
mkdir sif tmp
# remember that singularity needs the full path (not a relative path, like tmp)
export SINGULARITY_TMPDIR=$PWD/tmp/
singularity pull --dir sif docker://linsalrob/superfocus:v1.6.0_${TAG}
```

Now, you need to link the directories so that everything will work.

Assuming you have

```
/home/edwa0468/Projects/superfocus
/home/edwa0468/Projects/superfocus/fasta/


/home/edwa0468/Projects/superfocus_db/db/
/home/edwa0468/Projects/superfocus_db/db/database_PKs.txt
/home/edwa0468/Projects/superfocus_db/db/static/mmseqs2
```


Then this should work in a slurm job

```
singularity exec --bind $PWD/:/superfocus,/home/edwa0468/Projects/superfocus_db:/superfocus_db sif/superfocus_v1.6.0_5abbef95a6ca.sif \
        superfocus -q /superfocus/fasta -t 16 -dir /superfocus/superfocus -a mmseqs2 -db DB_95 -b /superfocus_db 
```
