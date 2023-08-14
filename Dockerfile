# Use a specific Debian Buster image with a slim tag to ensure stability
FROM debian:buster-slim

# Set the default shell to bash
SHELL ["/bin/bash", "-c"]

# Install required packages in a single RUN instruction
RUN mkdir -p /usr/share/man/man1 && \
    apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        build-essential \
        gnupg \
        libfftw3-dev \
        default-jdk \
        curl \
        python3 \
        python3-dev \
        python3-pip \
        # Additional packages for Google Cloud SDK
        apt-transport-https \
        ca-certificates \
        software-properties-common && \
    # Add Google Cloud SDK repository
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install google-cloud-sdk -y && \
    # Symlink python3 to python
    ln -s /usr/bin/python3 /usr/bin/python && \
    # Upgrade pip, setuptools, and wheel
    python -m pip install --upgrade pip setuptools wheel --no-cache-dir && \
    # Install required Python packages
    python -m pip install 'scanpy[leiden]' tangram-sc pytorch --no-cache-dir && \
    # Clean up to reduce image size
    apt-get -qq autoremove -y && \
    apt-get -qq clean && \
    rm -rf /var/lib/apt/lists/*

# Set the entrypoint or CMD as needed for your specific application

