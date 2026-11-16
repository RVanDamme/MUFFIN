#!/usr/bin/env bash
# Installs the pinned toolchain for the MUFFIN nf-core rebuild into ~/.local/bin.
# Safe to re-run. Needs Java 17+ and Python 3.10+ already present.
set -euo pipefail

NXF_VER_PIN="${NXF_VER_PIN:-26.10.6}"
NFCORE_VER_PIN="${NFCORE_VER_PIN:-4.2.2}"
BIN="${HOME}/.local/bin"
mkdir -p "${BIN}"
export PATH="${BIN}:${PATH}"

java -version 2>&1 | grep -q 'version "\(1[7-9]\|[2-9][0-9]\)' \
    || { echo "Java 17+ is required (see the install-nextflow skill)"; exit 1; }

if ! nextflow -version 2>/dev/null | grep -q "${NXF_VER_PIN}"; then
    curl -fsSL https://get.nextflow.io -o "${BIN}/nextflow"
    chmod +x "${BIN}/nextflow"
    NXF_VER="${NXF_VER_PIN}" nextflow -version >/dev/null
fi

if ! nf-core --version 2>/dev/null | grep -q "${NFCORE_VER_PIN}"; then
    python3 -m pip install --user --quiet "nf-core==${NFCORE_VER_PIN}" \
        || python3 -m pip install --user --quiet --break-system-packages "nf-core==${NFCORE_VER_PIN}"
fi

if ! command -v nf-test >/dev/null; then
    (cd "${BIN}" && curl -fsSL https://get.nf-test.com | bash >/dev/null) \
        || echo "WARN: nf-test could not be downloaded here (GitHub release assets blocked?). Install it manually: https://www.nf-test.com"
fi

echo "--- toolchain ---"
NXF_VER="${NXF_VER_PIN}" nextflow -version 2>&1 | grep -i 'version'
nf-core --version 2>&1 | tail -1
nf-test version 2>/dev/null | grep -i 'nf-test' | head -1 || echo "nf-test: NOT INSTALLED"
for engine in docker apptainer singularity; do
    command -v "${engine}" >/dev/null && echo "container engine found: ${engine}"
done
docker info >/dev/null 2>&1 || echo "NOTE: no running Docker daemon here; real (non-stub) runs need a machine with Docker or Apptainer."
echo "Add to your shell: export PATH=${BIN}:\$PATH NXF_VER=${NXF_VER_PIN}"
