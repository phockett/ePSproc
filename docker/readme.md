# ePSproc + JupyterLab builds with Docker

Dockerfiles for ePSproc, based on [Jupyter Docker Stacks images](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html), are in `/docker`.


## Quick build container only

### Basic container + ePSproc install.

```
docker build -t jupyterlab_epsproc .
docker run -p 8888:8888 jupyterlab_epsproc
```

where the port mapping is `host:container`.


### Basic container + ePSproc +PEMtk install.

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


### Quick build with Compose (includes some extra options)

See `docker-compose.ePSproc.yml` for options, including port and storage mapping.

```
docker-compose -f docker-compose.ePSproc.yml build
```

NOTE: the compose version is currently not working.


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
