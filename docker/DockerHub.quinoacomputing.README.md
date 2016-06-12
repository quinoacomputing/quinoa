### What is Quinoa?
[Quinoa](https://github.com/quinoacomputing/quinoa) is a set of computational tools that enables research and numerical analysis in fluid dynamics. At this time it is a test-bed to experiment with various algorithms using fully asynchronous runtime systems.

### Purpose
This image contains all executables, configured for a single computer, and thus intended to quickly try Quinoa on a multi-core workstation. To run on clusters of networked compute nodes you must build from source, see the [README](https://github.com/quinoacomputing/quinoa/blob/master/README.md) for build instructions.

### Usage
1. Run the [one of the containers](https://hub.docker.com/r/quinoacomputing/quinoa/tags) on your local machine
```
docker run -ti quinoacomputing/quinoa:alpine
```
2. Run Quinoa executables inside the container, e.g.,
```
charmrun +p4 /usr/local/bin/unittest -v
```

### Organization
- All Quinoa executables are in `/usr/local/bin/`.
- Inside the container there is a minimal install of [Alpine Linux](http://www.alpinelinux.org).
- Install [Alpine packages](https://pkgs.alpinelinux.org) using [apk](https://wiki.alpinelinux.org/wiki/Alpine_Linux_package_management).
