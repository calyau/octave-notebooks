FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install Octave and dependencies
RUN apt-get update && apt-get install -y \
    octave \
    octave-signal \
    octave-image \
    octave-statistics \
    octave-control \
    octave-io \
    gnuplot-qt \
    python3 \
    python3-pip \
    python3-venv \
    jupyter-notebook \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install the Octave Jupyter kernel
RUN pip3 install --break-system-packages metakernel==0.30.2 octave_kernel==0.35.1

# Create workspace
WORKDIR /notebooks

# Jupyter configuration
RUN jupyter notebook --generate-config && \
    echo "c.NotebookApp.ip = '0.0.0.0'" >> /root/.jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.port = 8888" >> /root/.jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.open_browser = False" >> /root/.jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.allow_root = True" >> /root/.jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.token = ''" >> /root/.jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.password = ''" >> /root/.jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.notebook_dir = '/notebooks'" >> /root/.jupyter/jupyter_notebook_config.py

# Configure Octave for inline plotting
RUN mkdir -p /root/.config/octave && \
    echo "graphics_toolkit('gnuplot');" >> /root/.octaverc && \
    echo "setenv('GNUTERM', 'dumb');" >> /root/.octaverc

EXPOSE 8888

CMD ["jupyter", "notebook"]
