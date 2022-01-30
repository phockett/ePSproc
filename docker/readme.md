# ePSproc + JupyterLab builds with Docker

Dockerfiles for ePSproc, based on [Jupyter Docker Stacks images](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html), are in `/docker`.


Quick build container only:


```
docker build -t jupyterlab_epsproc .
docker run -p 8888:8888 jupyterlab_epsproc
```

where the port mapping is `host:container`.



Quick build with Compose (includes some extra options):


```
docker-compose -f docker-compose.ePSproc.yml build
```

NOTE: the compose version is currently not working.


For more details & options, see:

- [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/running.html).
- [Dockerfile reference](https://docs.docker.com/engine/reference/builder/).
- [Docker Compose reference](https://docs.docker.com/compose/compose-file/compose-file-v3/).

For use with full JupyterHub deployment, see https://github.com/phockett/jupyterhub-docker


## TODO

- Further testing.
- Different install types. (Here using `pip install git+git://github.com/phockett/ePSproc@dev`).
- Deploy to Docker Cloud/Jupyter Stacks repos.
