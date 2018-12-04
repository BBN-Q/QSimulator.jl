FROM library/julia:1.0.2

RUN apt-get update && apt-get install -yq bzip2

# The ~/.julia from the host is mounted here
ENV JULIA_DEPOT_PATH /opt/.julia
RUN mkdir -p /opt/.julia

# The project source is mounted here
WORKDIR /project-mount

# COPY Manifest.toml .
# COPY Project.toml .
# COPY scripts ./scripts

# RUN /qsimulator/scripts/build.sh

# COPY . .

# RUN julia /qsimulator/scripts/setup.jl /qsimulator

CMD ["julia"]
