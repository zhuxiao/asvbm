## Docker
Clone this repository and build a Docker image as follows.

-------------------
If the current directory contains a clone of the asvbm repository, asvbm can be run in Docker as follows:
```
$ docker build -t asvbm:v1.2.0 .	
$ docker images
REPOSITORY        TAG       IMAGE ID       CREATED         SIZE
asvbm             v1.2.0    22ab0947a2ba   1 second ago    <...>
$ docker run -it asvbm:v1.2.0
OR
$ docker run -it --name xxx -v 'pwd':/data_test asvbm:v1.2.0 bash
OR
$ docker run -it --name xxx -v 'pwd':/data_test asvbm:v1.2.0 ./asvbm -m 50000 -T "tool1;tool2;tool3" /data_test/user_sv.vcf /data_test/user_sv1.vcf /data_test/user_sv2.vcf /data_test/benchmark_sv.vcf /data_test/ref.fa -o /data_test/output
```

The -v argument mounts the current directory as /data_test in the Docker image. The output should also appear in the current directory.

The default Docker image is based on Ubuntu. 

Clone this repository and build a Docker image as follows. A pre-built docker image can be found here: https://hub.docker.com/xzhu2/xxx. It can be obtained by running:

```
$ docker pull xzhu2/xxx
```

The complete list of dependencies / packages to install beforehand can be found in the Dockerfile.