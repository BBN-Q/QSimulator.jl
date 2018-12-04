FROM library/julia:1.0.2

RUN apt-get update && apt-get install -yq bzip2

WORKDIR /qsimulator

COPY Manifest.toml .
COPY Project.toml .
COPY scripts ./scripts

RUN /qsimulator/scripts/build.sh

COPY . .

CMD ["./runtests.sh"]
