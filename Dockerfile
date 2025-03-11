FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

# Include `--secret id=api,src=api_key.txt` in the docker build command to mount the API key.

# SPECIES = senterica or ecoli
ARG SPECIES
ARG VERSION

RUN mkdir /hclink

COPY src uv.lock pyproject.toml .env schemes.json LICENSE README.md /hclink/

WORKDIR /hclink

RUN uv sync --frozen

ENV VERSION=${VERSION}
ENV SPECIES=${SPECIES}

RUN --mount=type=secret,id=api,env=API_KEY uv run hclink build ${VERSION} "${API_KEY}" -s ${SPECIES} --clean

ENTRYPOINT ["uv", "run", "hclink", "assign", "-"]

