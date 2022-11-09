# ePSproc + JupyterLab builds with Docker

Dockerfiles for ePSproc, based on [Jupyter Docker Stacks images](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html), are in `/docker`.


## Quick build container only

### Basic container + ePSproc install.

```
docker build -t jupyterlab_epsproc .
docker run -p 8888:8888 jupyterlab_epsproc
```

where the port mapping is `host:container`.


### Basic container + ePSproc + PEMtk install.

Includes PEMtk install, with libmsym build for symmetrized harmonics support.

```
docker build -t jupyterlab_epsproc_pemtk -f Dockerfile-PEMtk .
docker run -p 8888:8888 jupyterlab_epsproc_pemtk
```

where the port mapping is `host:container`.


## Other options

### Terminal in container

Use `exec -it <container> bash` to attach to a running container, or `run -it <container> bash` to spin one up, and connect to the terminal.

E.g. for named container as above: `docker exec -it jupyterlab_epsproc bash`


### Build with Compose (includes some extra options)

#### Compose + local install (for development)

See `docker-compose.local.yml` for options, including port and storage mapping.

The default case pulls source image, installs dependencies, and mounts local folders for source code and notebooks - see the compose file to set the mount paths.

```
docker-compose -f docker-compose.local.yml build
docker-compose -f docker-compose.local.yml up
```

Once up, run `installlocal.sh` in the container to install the local code in editable mode.

A one-liner for this (assuming a `github` directory mounted with local code): `docker exec jupyterlab-ePSproc-dev /home/jovyan/github/ePSproc/docker/localinstall.sh`, or run via a terminal in the container, e.g. `docker exec -it jupyterlab-ePSproc-dev bash` to connect.

Note the script defaults to `~/github` for the install, pass a different path if required, e.g. `docker exec jupyterlab-ePSproc-dev /home/jovyan/github/ePSproc/docker/localinstall.sh /path/to/repos`.



#### Compose + ePSproc install

See `docker-compose.ePSproc.yml` for options, including port and storage mapping.

```
docker-compose -f docker-compose.ePSproc.yml build
```

NOTE: this compose version is currently not working.


## Notes

For more details & options, see:

- [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/running.html).
- [Dockerfile reference](https://docs.docker.com/engine/reference/builder/).
- [Docker Compose reference](https://docs.docker.com/compose/compose-file/compose-file-v3/).

For use with full JupyterHub deployment, see https://github.com/phockett/jupyterhub-docker


## TODO

- Further testing.
- Different install types. (Here using `pip install git+git://github.com/phockett/ePSproc@dev`).
- Deploy to Docker Cloud/Jupyter Stacks repos.
