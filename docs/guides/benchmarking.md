# Benchmarking with MucOneUp

[MucOneUp](https://github.com/berntpopp/muconeup) is a companion tool that generates simulated MUC1 VNTR haplotypes and PacBio HiFi amplicon reads with known ground truth, for benchmarking and validation of open-pacmuci.

---

## Prerequisites

| Tool | Version | Purpose | Install |
|------|---------|---------|---------|
| MucOneUp | >= 0.44.0 | VNTR simulation + amplicon reads | `pip install muc-one-up` |
| pbsim3 | latest | PacBio HiFi read simulation | `conda install -c bioconda pbsim3` |
| minimap2 | >= 2.28 | Read alignment | `conda install -c bioconda minimap2` |
| samtools | >= 1.21 | BAM processing | `conda install -c bioconda samtools` |

!!! note "MucOneUp config.json"
    MucOneUp requires a `config.json` file defining repeat sequences and probabilities. This file is included in the [MucOneUp repository](https://github.com/berntpopp/muconeup).

---

## Generating Test Data

open-pacmuci includes a test data generation script that creates 10 simulated samples covering a range of mutation types and allele configurations.

### Run the Generator

```bash
# Ensure pbsim3 is on PATH
export PATH="/path/to/conda/envs/env_pacbio/bin:$PATH"

# Generate all test samples
python scripts/generate_testdata.py
```

### Test Samples Generated

| Sample | Allele 1 | Allele 2 | Mutation | Purpose |
|--------|----------|----------|----------|---------|
| `sample_dupc_60_80` | 60 repeats | 80 repeats | dupC at repeat 25 | Common mutation, well-separated alleles |
| `sample_dupa_60_80` | 60 repeats | 80 repeats | dupA at repeat 25 | Alternative single-base mutation |
| `sample_insg_60_80` | 60 repeats | 80 repeats | insG at repeat 25 | G insertion mutation |
| `sample_dupcccc_60_80` | 60 repeats | 80 repeats | insCCCC at repeat 25 | Multi-base insertion |
| `sample_del_60_80` | 60 repeats | 80 repeats | del18_31 at repeat 25 | Large deletion |
| `sample_normal_60_80` | 60 repeats | 80 repeats | None | Negative control |
| `sample_homozygous_60_60` | 60 repeats | 60 repeats | dupC at repeat 25 | Same-length alleles |
| `sample_asymmetric_25_140` | 25 repeats | 140 repeats | dupC at repeat 10 | Extreme length asymmetry |
| `sample_short_25_30` | 25 repeats | 30 repeats | dupC at repeat 10 | Very close alleles |
| `sample_long_120_140` | 120 repeats | 140 repeats | dupC at repeat 50 | Long alleles |

All samples are generated at **200x template coverage** with realistic PCR length bias.

**Output location:** `tests/data/generated/`

---

## Running Integration Tests

Integration tests exercise the full pipeline against the generated test data.

```bash
# Ensure external tools are on PATH
export PATH="/path/to/conda/envs/env_pacbio/bin:$PATH"

# Run integration tests
uv run pytest tests/integration/ -v --no-cov
```

---

## Validation Results

Tested against MucOneUp-simulated PacBio HiFi amplicon data:

| Metric | Result |
|--------|--------|
| VNTR classification (given correct sequence) | 20/20 haplotypes perfect (100% exact match) |
| Allele length detection (gap >= 5 repeats) | Exact match |
| Allele length detection (gap 3-4 repeats) | Within +/- 2 repeats |
| Full pipeline with Clair3 (dupC sample) | dupC detected at correct position, 100% confidence |
| Close allele pairs (gap 3-12 repeats) | 12/12 resolved by indel-valley splitting |

### Known Limitations

!!! warning "Edge cases in allele detection"
    These limitations affect allele detection from BAM files, not classification:

    - **`sample_asymmetric_25_140`**: The 140-repeat allele is undetectable (only ~2 reads due to extreme PCR bias favoring the 25-repeat allele)
    - **`sample_short_25_30`**: Alleles 5 repeats apart may merge into one cluster when the gap between them contains no zero-count contigs

These are fundamental limitations of the PCR amplicon approach, not software bugs.

---

## Custom Benchmarking

### Generate Specific Test Cases

Use MucOneUp directly to create custom scenarios:

```bash
# Generate a diploid haplotype with specific lengths and mutation
muconeup --config config.json simulate \
  --out-base custom_test \
  --out-dir my_output/ \
  --num-haplotypes 2 \
  --fixed-lengths 45 \
  --fixed-lengths 65 \
  --mutation-name dupC \
  --mutation-targets 1,20 \
  --output-structure \
  --seed 42

# Simulate PacBio HiFi amplicon reads
muconeup --config config.json reads amplicon \
  custom_test.001.simulated.fa \
  --coverage 200 \
  --out-base custom_reads

# Run open-pacmuci on the simulated data
open-pacmuci run \
  --input custom_reads.amplicon.bam \
  --output-dir custom_results/ \
  --clair3-model /path/to/models/hifi

# Compare to ground truth
jq '.mutations' custom_test.001.simulation_stats.json
cat custom_results/repeats.txt
```

---

## Next Steps

- **[Known Mutations](../reference/mutations.md)** -- Full mutation catalog
- **[CLI Reference](../reference/cli.md)** -- All command options
- **[Core Concepts](../getting-started/concepts.md)** -- Pipeline architecture
