FROM library/julia:1.0.2

RUN apt-get update && apt-get install -yq bzip2

WORKDIR /workspace

COPY . .

CMD ["/workspace/ci/runtests.sh"]
