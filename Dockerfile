FROM rust:1.80-bookworm AS build

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    cmake \
    libhts-dev \
    nodejs \
    npm \
    pkg-config \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

RUN cmake -S . -B build -DMITO_BUILD_TESTS=ON \
  && cmake --build build -j2 \
  && ctest --test-dir build --output-on-failure

RUN cargo build --release --workspace

RUN npm ci \
  && npm --workspace visualization-lib run build \
  && npm --workspace web run build

FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    libhts3 \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=build /src/target/release/mito-cli /usr/local/bin/mito-cli
COPY --from=build /src/target/release/mito-server /usr/local/bin/mito-server
COPY --from=build /src/core/data /app/core/data
COPY --from=build /src/web/dist /app/web/dist

ENV MITO_CLINICAL_ANNOTATIONS=/app/core/data/clinical_annotations.tsv
ENV MITO_SERVER_ADDR=0.0.0.0:8080
EXPOSE 8080
CMD ["mito-server"]
