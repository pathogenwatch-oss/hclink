FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS build

# Include `--secret id=api,src=api_key.txt` in the docker build command to mount the API key.

COPY src uv.lock pyproject.toml .env schemes.json LICENSE README.md /hclink/

WORKDIR /hclink

RUN uv venv && uv build --wheel && mkdir /build && mv LICENSE README.md schemes.json dist/*.whl /build/

FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS hclink

# SPECIES = senterica or ecoli
ARG SPECIES
ARG VERSION

COPY --from=build /build /hclink

WORKDIR /hclink

RUN uv venv
RUN . .venv/bin/activate && uv pip install hclink-"${VERSION}"-py3-none-any.whl

ENV VERSION=${VERSION}
ENV SPECIES=${SPECIES}

RUN --mount=type=secret,id=api,env=API_KEY \
    . .venv/bin/activate && \
    hclink build ${VERSION} "${API_KEY}" -s ${SPECIES} --clean

ENTRYPOINT ["/bin/bash", "-c", "source .venv/bin/activate && hclink assign -"]
