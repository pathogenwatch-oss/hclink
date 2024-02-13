FROM python:3.11-slim

ARG VERSION
ARG API_KEY

#RUN apt update && \
#    apt install -y --no-install-recommends pip3 && \
#    rm -rf /var/lib/apt/lists/*
#
RUN --mount=type=bind,source=requirements.txt,target=/tmp/requirements.txt \
    pip install --requirement /tmp/requirements.txt && \
    pip cache purge && \
    rm -rf ~/.cache/pip

RUN mkdir /hclink

WORKDIR /hclink

COPY hclink.py /hclink/

COPY lib /hclink/lib

COPY LICENSE /hclink/LICENSE

ENV API_KEY=${API_KEY}
ENV VERSION=${VERSION}

RUN python hclink.py build ${VERSION} ${API_KEY} --clean

ENTRYPOINT ["python", "hclink.py", "assign", "-"]

