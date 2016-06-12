### **What is Quinoa?**
[Quinoa](https://github.com/quinoacomputing/quinoa) is a set of computational tools that enables research and numerical analysis in fluid dynamics. At this time it is a test-bed to experiment with various algorithms using fully asynchronous runtime systems.

### **Purpose**
The images in this repository are configured for a single computer, and thus intended to quickly try the executables ([:alpine](https://hub.docker.com/r/quinoacomputing/quinoa/tags/)) or setup a complete development environment ([:debian](https://hub.docker.com/r/quinoacomputing/quinoa/tags/)) on a multi-core workstation. For production runs on clusters of networked compute nodes you should build from source, see the [README](https://github.com/quinoacomputing/quinoa/blob/master/README.md) or one of the [docker files](https://github.com/quinoacomputing/quinoa/tree/master/docker) for instructions.

### **Usage**
1. Run the [one of the containers](https://hub.docker.com/r/quinoacomputing/quinoa/tags) on your local machine
```
docker run -ti quinoacomputing/quinoa:alpine
```
2. Run Quinoa executables inside the container, e.g.,
```
charmrun +p4 /usr/local/bin/unittest -v
```

### **Organization**

#### [:alpine](https://hub.docker.com/r/quinoacomputing/quinoa/tags/)
- _Purpose:_ **Containerized release**
- _Operating system:_ [Alpine Linux](http://www.alpinelinux.org), install additional [packages](https://pkgs.alpinelinux.org) using [apk](https://wiki.alpinelinux.org/wiki/Alpine_Linux_package_management)
- Minimalistic container with only static executables (installed in `/usr/local/bin/`) + OpenMPI

#### [:debian](https://hub.docker.com/r/quinoacomputing/quinoa/tags/)
- _Purpose:_ **Simple development environment**
- _Operating system:_ [Debian Linux](https://www.debian.org), install additional [packages](https://packages.debian.org/testing/) using [apt-get](https://www.debian.org/doc/manuals/debian-faq/ch-pkgtools.en.html)
- System-wide GNU compilers + OpenMPI, source, third-party libraries, documentation, build objects
