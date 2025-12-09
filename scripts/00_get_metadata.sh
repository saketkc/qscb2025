#!/bin/bash
set -eox pipefail
pip install git+https://github.com/saketkc/pysradb@develop -q
pysradb metadata GSE192740 --detailed --saveto data/metadata/gex_GSE192740_metadata.tsv
pysradb metadata SRP352824 --detailed --saveto data/metadata/fastq_GSE192740_metadata.tsv
