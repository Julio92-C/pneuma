#!/usr/bin/env bash
# Toolchain for PNEUMA development: Nextflow, nf-core tools, nf-test, pytest.
set -euo pipefail
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/ && sudo chmod +x /usr/local/bin/nextflow
python3 -m pip install --user --upgrade "nf-core>=3.3" pytest
curl -fsSL https://get.nf-test.com | bash
sudo mv nf-test /usr/local/bin/
nextflow -version && nf-core --version && nf-test version
