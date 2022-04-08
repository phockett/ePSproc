# Jupyterlab with ePSproc install.
# Based on OpenDreamKit JupyterHub distro, https://github.com/defeo/jupyterhub-docker/
# Uses Jupyter Stacks scipy notebook as base container, https://github.com/jupyter/docker-stacks/tree/master/scipy-notebook

FROM jupyter/scipy-notebook


# USER $NB_UID


# Add conda env hook
COPY ./conda-activate.sh /usr/local/bin/before-notebook.d/

#*** Add ePSproc + dependencies
# Note some may become/are optional, see also requirements.txt
RUN conda install --quiet --yes -c conda-forge spherical_functions
RUN conda install --quiet --yes xarray=0.17
RUN conda install --quiet --yes seaborn=0.9
RUN conda install --quiet --yes holoviews

# Note optional @branch here
# RUN pip install git+git://github.com/phockett/ePSproc@dev
RUN pip install git+https://github.com/phockett/ePSproc@dev
