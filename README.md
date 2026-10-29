# Lab 4 — Database-Driven Patient Feature Table

## Overview

The objective of this lab is to build a reproducible, patient-level
feature table from the supplied TCGA Lung Adenocarcinoma PanCancer Atlas
SQLite database. The final table contains **566 distinct patients (566
rows and 11 columns)** and no patients were dropped from the source
tables. This report describes the methods and results for tasks 1-5.

| Gene | Non-Silent mutation records | Patients with at least one recorded mutation |
|:---|:---|:---|
| TP53 | 305 | 287 |
| KRAS | 172 | 168 |
| EGFR | 85 | 69 |
| STK11 | 81 | 74 |
| KEAP1 | 103 | 102 |

There are 746 qualifying mutation records in these five genes,
distributed over 700 nonzero patient–gene combinations. Mutation records
and patients differ in that a patient can have multiple records in one
gene or records in multiple genes.

## Tools, Input, and Execution

Tools: Python standard library (`sqlite3`, `csv`, `collections`,
`pathlib`, `json`, `hashlib`, `argparse`); `sqlite3` for joins and
aggregation. Python dictionaries for converting gene counts from rows
into columns. Built and tested with **Python 3.14.7** and **SQLite
3.53.4**.

``` bash
# Check python and sqlite versions to ensure reproducibility
python3 --version  # version 3.14.7
sqlite3 --version  # version 3.53.4
```

From the project directory:

``` bash
# Run `build_features.py` and `verify_reproducibility`.
python3 scripts/build_features.py
python3 scripts/verify_reproducibility.py
```

Secure Hash Algorithm (SHA-256) is a hash function that converts any
data input into a unique 256-bit string. `SHA256SUMS.txt` leverages this
function so that relative paths entered on the command line follow the
shell's working directory.The default database and output paths are
resolved relative to the script's project directory, not the shell's
working directory. The script expects the existing four-table database
and 566 patients. It opens the input and confirms its SHA-256 is
unchanged after running. It does not select an arbitrary sample if
future data have multiple samples, which would require a documented
clinical selection policy.

**SHA-256 Input** (256-bit string assigned to `cohort.db` database):

``` text
6da3f91401918dda543aced9fa6f411dcecf204552cebfc6e54ab7d8f387e8e7
```

The SHA-256 file contains the analysis database input. Running the build
overwrites generated files in the specified output directory but does
not change the database.

## Task 1 — Inspected Schema and Grain

Schema was obtained from `sqlite_master`. Exact column types,
primary-key flags, and declared foreign keys are saved in
`outputs/schema.json`.

| Table     | Rows    | Grain                                       | Primary key |
|:----------|:--------|:--------------------------------------------|:------------|
| patients  | 566     | One patient                                 | patient_id  |
| samples   | 566     | One tumor sample                            | sample_id   |
| mutations | 200,789 | One recorded somatic mutation in one sample | mutation_id |
| genes     | 17,756  | One stored gene-symbol entry                | gene_symbol |

### Columns

| Table | Column | SQLite type | Meaning |
|----|----|----|----|
| patients | patient_id | TEXT | De-identified TCGA patient barcode |
| patients | sex | TEXT | Recorded sex |
| patients | age_at_diagnosis | REAL | Age in years |
| patients | stage | TEXT | Recorded pathologic tumor stage |
| patients | vital_status | TEXT | Recorded survival status |
| patients | followup_months | REAL | Source OS_MONTHS survival/follow-up time in months |
| samples | sample_id | TEXT | Tumor sample barcode |
| samples | patient_id | TEXT | Required patient reference |
| samples | tissue_source_site | TEXT | Collecting institution/source site |
| samples | tmb_nonsynonymous | REAL | Source nonsynonymous tumor mutation burden |
| samples | aneuploidy_score | REAL | Source aneuploidy score |
| samples | msi_mantis_score | REAL | Source MANTIS microsatellite-instability score |
| mutations | mutation_id | INTEGER | Mutation-record identifier |
| mutations | sample_id | TEXT | Required sample reference |
| mutations | gene_symbol | TEXT | Required gene reference |
| mutations | variant_classification | TEXT | Functional classification, including Silent |
| mutations | variant_type | TEXT | Variant type, e.g. SNP or DEL |
| mutations | protein_change | TEXT | Source protein-change annotation |
| mutations | chromosome | TEXT | Chromosome label |
| mutations | start_pos | INTEGER | GRCh37/hg19, 1-based start position |
| mutations | end_pos | INTEGER | GRCh37/hg19, 1-based end position |
| mutations | t_ref_count | INTEGER | Tumor reference-allele read count |
| mutations | t_alt_count | INTEGER | Tumor alternate-allele read count |
| genes | gene_symbol | TEXT | Stored gene symbol |
| genes | entrez_gene_id | INTEGER | Entrez Gene identifier |
| genes | chromosome | TEXT | Gene chromosome annotation |

### Relationships

- `samples.patient_id → patients.patient_id`: a many-to-one in the
  declared schema. The observed database has exactly one sample for each
  of its 566 patients. `samples.patient_id` is not declared unique, so
  one-to-one must be checked rather than assumed.
- `mutations.sample_id → samples.sample_id`: many mutations can belong
  to a sample. Four samples have zero mutation rows.
- `mutations.gene_symbol → genes.gene_symbol`: many mutations can
  reference one stored gene entry.

`PRAGMA integrity_check` returned `ok` and `PRAGMA foreign_key_check`
returned no issues. An exact mutation-to-gene join matched all 200,789
mutation rows. 2,431 mutation rows have a symbol differing from
`UPPER(TRIM(gene_symbol))` and each symbol also exists as stored keys in
`genes`.

Normalized gene keys have no collisions in this database. None of those
differences affects the five target genes: their stored symbols are
already canonical. Normalization is applied to mutation gene symbols
during aggregation without changing raw data. No gene-annotation join is
needed for the feature table, which avoids unnecessary complexity in the
join.

## Task 2 — SQL queries and expected versus observed counts

The runnable queries are in `sql/`. Results and the executed join audit
are included in the `outputs/` folder.

### Query A: patient–gene mutation counts

`sql/patient_gene_counts.sql` creates a five-gene list, aggregates
eligible mutation rows to sample–gene grain, and then left-joins those
counts to every patient–gene combination. `COUNT(mutation_id)` counts
records; `COALESCE` converts an absent aggregate to zero. Summing the
sample-level counts allows mutation counts to be defined across samples,
although the feature build deliberately requires at most one sample for
sample-level clinical fields.

P**rotein-altering is defined here as any classification in this
supplied coding/filtered database other than `Silent`**. Included
classes are Missense_Mutation, Nonsense_Mutation, Nonstop_Mutation,
Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins,
Splice_Site, and Translation_Start_Site.

Splice effects are treated as having protein-altering potential. There
are no NULL classifications in this input and the script would fail if
NULL or undefined classes appear. No records are removed on the basis of
a missing protein-change string, and no undocumented de-duplication is
performed.

Expected result: **566 × 5 = 2,830 rows**, one per patient and target
gene. Actual result: **2,830**. All five genes appear for every patient,
including patients with zero qualifying records.

### Query B: patients per tissue source site and stage

`sql/site_stage_counts.sql` left-joins patients to samples and uses
`COUNT(DISTINCT p.patient_id)` per observed source-site/stage
combination. It does not join mutations, which would inflate counts.
Unknown stages remain in NULL groups, exported as `NA`. Tissue source
site identifies the collecting institution, not an anatomic tumor site.

Before running the grouped SQL, the script independently enumerates
observed site/stage pairs from the raw patient and sample tables in
Python. Expected result: **142 observed combinations**. Actual result:
**142 rows**; their patient counts sum to **566**. Unobserved site/stage
combinations are not generated as zero rows. The 142-row result is
smaller than 566 because its grain is a group rather than a patient.
This is expected and does not indicate that patients have been excluded.

### Join counts and grain checks

| Step | Rows before | Expected after | Actual after | Result grain |
|:---|:---|:---|:---|:---|
| patients LEFT JOIN samples | 566 | 566 | 566 | patient, after verifying cardinality |
| patients CROSS JOIN five genes | 566 | 2,830 | 2,830 | patient × gene |
| patient–gene grid LEFT JOIN samples | 2,830 | 2,830 | 2,830 | patient × gene × sample |
| LEFT JOIN sample–gene mutation aggregate | 2,830 | 2,830 | 2,830 | patient × gene × sample |
| Query A final GROUP BY | 2,830 | 2,830 | 2,830 | patient × gene |
| Query B GROUP BY | 566 | 142 | 142 | observed site × stage |
| Python reshape and combination | 566 | 566 | 566 | patient |

The mutation aggregation reduces 746 eligible raw records to 700
observed sample–gene groups before joining; this is deliberate
aggregation. The aggregate has at most one matching row per sample–gene
key, so the subsequent left join cannot multiply the grid. Full executed
counts, including the raw-record inspection join, are in
`outputs/join_audit.tsv`.

Optional-side match conditions are in `ON`. The `WHERE` inside the
independent mutation aggregation defines eligible records before the
left join; it never filters the patient grid. Placing
`m.gene_symbol='KRAS'` and the non-Silent condition in an outer `WHERE`
would return only 168 patients instead of 566. A deliberately tested
`COUNT(*)` on the left join gives 1 for a sample with no KRAS records
because it counts the retained placeholder row, while counting the
mutation key or using the pre-aggregated count avoids that error.

## Task 3 — Build the Feature Table with SQL and Python

1.  Inspect schema, integrity, classifications, and sample cardinality.
2.  Execute `patient_clinical.sql` to retain all patients and attach
    sample attributes.
3.  Execute `patient_gene_counts.sql` to obtain the complete 2,830-row
    long table.
4.  Reshape in Python: use a dictionary keyed by patient ID with five
    gene counts. Combine each patient's SQL-derived clinical row with
    that patient's dictionary. This is Python long-to-wide reshaping,
    with no SQL pivot.
5.  Verify unique IDs, preserved patient set, nonnegative integer gene
    counts, original clinical values, and independent raw-record counts.
6.  Export deterministic, sorted TSVs with an explicit `NA` token for
    missing values.

`outputs/features.tsv` columns, in order:

``` text
patient_id | sex | age_at_diagnosis | stage | tissue_source_site | tmb_nonsynonymous | TP53 | KRAS | EGFR | STK11 | KEAP1
```

The first 6 fields come from the clinical SQL, while the 5 gene fields
are counts after reshaping in Python. The TSV has a header plus 566 rows
of data.

## Task 4 — Missing data and verification

| Feature or gap | Observed missing/no-record count | Decision |
|:---|:---|:---|
| sex | 52 patients | Keep patient; export SQL NULL as NA |
| age_at_diagnosis | 71 patients | Keep patient; export SQL NULL as NA |
| stage | 54 patients | Keep patient and missing-stage group; export NA |
| tissue_source_site | 0 patients | Preserve source values; export NA if ever absent |
| tmb_nonsynonymous | 0 patients | Preserve source values, including zero |
| Sample absent | 0 patients | Left join would retain patient; sample attributes would be NA |
| All mutation rows absent | 4 samples | Retain patient and use zero recorded counts for all five genes |
| Qualifying mutation absent for a gene | 2,130 of 2,830 patient–gene pairs | Use zero recorded counts |

Clinical values were not imputed, and no complete-case filtering was
used. Missing sex, age, and stage can overlap within the same patient,
so those counts must not be added to estimate excluded patients.

**Zero is a count of records, not a claim of confirmed wild-type
status.** It means that the supplied database contains no qualifying
mutation record for that patient–gene pair. `NA` would instead mean that
a count is undefined or unavailable under a different measurement-based
definition. Those four are listed in `outputs/validation.json` and their
zero counts should not be interpreted that there were no mutations. No
separate clinical measurements are replaced with zero. The source TMB
values are copied rather than recomputed or used to infer coverage.

### Raw Record Checks

Raw mutation records were inspected, then compared to the final feature
rows.

| Patient | Raw evidence | Expected and observed feature result |
|----|----|----|
| TCGA-05-4382 | EGFR mutation IDs 2635 (p.R222L) and 2636 (p.E545Q); TP53 ID 1770 (p.V272G); all Missense_Mutation | TP53=1, KRAS=0, EGFR=2, STK11=0, KEAP1=0 |
| TCGA-17-Z023 | TP53 ID 21335 is Nonsense_Mutation (p.E198\*); ID 21336 is Silent (p.V197=) | TP53=1, confirming Silent exclusion |
| TCGA-64-1678 | Sample TCGA-64-1678-01 exists; mutation table contains no rows for it | All five gene counts=0; patient retained with female, age 70, stage IIA |
| TCGA-17-Z000 | Missing sex, age, and stage in patients; sample still exists | All three clinical fields remain NA; patient retained |

Raw target-gene rows and corresponding feature rows are in
`outputs/patient_check_mutations.tsv` and
`outputs/patient_check_features.tsv`. A patient with no matching
target-gene records correctly has no row in the raw mutation extract but
still appears in the feature check.

The script additionally traverses all raw mutation rows in Python, maps
samples to patients without the feature SQL, and compares every one of
the **2,830** patient–gene values. It also checks all **142** site–stage
counts independently and verifies every clinical field against raw-table
dictionaries. No patients are lost at the 566-row clinical join or final
combination. The separate reproduction check rebuilds in a temporary
directory and compares TSV/schema bytes.

## Task 5 — Data Dictionary

`data_dictionary.tsv` contains one row per final column with its
meaning, units, source table/column, and missing-value policy. Age is in
years and counts reflect mutation records per patient. TMB uses the
conventional mutations-per-megabase interpretation. The SQLite schema
and supplied clinical metadata do not encode a unit, so this unit
interpretation is documented rather than asserted as a constraint.
Tissue source site is a categorical collecting-site name. Stage is
categorical. IA, IB, IIA, IIB, IIIA, and IIIB are preserved rather than
collapsed.

|  |  |  |  |  |
|:---|:---|:---|:---|:---|
| **column_name** | **meaning** | **units** | **source_table** | **missing_value_policy** |
| patient_id | De-identified TCGA patient barcode; unique feature-table key | identifier | patients.patient_id | No missing values allowed |
| sex | Recorded sex; female or male | category | patients.sex | SQL NULL exported as NA; not imputed |
| age_at_diagnosis | Age at diagnosis copied without rounding | years | patients.age_at_diagnosis | SQL NULL exported as NA; not imputed |
| stage | Recorded AJCC pathologic tumor stage; original substages retained | category | patients.stage | SQL NULL exported as NA; not imputed |
| tissue_source_site | Institution/source site collecting tissue and clinical metadata; not anatomic tumor location | category | samples.tissue_source_site | SQL NULL exported as NA; no missing values observed |
| tmb_nonsynonymous | Source-reported nonsynonymous tumor mutation burden; not recalculated from five genes | mutations per megabase $mut/Mb$; conventional TMB interpretation, not encoded in SQLite schema | samples.tmb_nonsynonymous | SQL NULL exported as NA; recorded zero retained |
| TP53 | Recorded non-Silent TP53 mutation rows across the patient's samples | mutation records per patient | mutations via samples.patient_id | 0 when no qualifying row recorded; no NA |
| KRAS | Recorded non-Silent KRAS mutation rows across the patient's samples | mutation records per patient | mutations via samples.patient_id | 0 when no qualifying row recorded; no NA |
| EGFR | Recorded non-Silent EGFR mutation rows across the patient's samples | mutation records per patient | mutations via samples.patient_id | 0 when no qualifying row recorded; no NA |

## Structure Notes for Reproducibility

``` text
HDS-applied-computing-lab4/
  README.md
  AI_USAGE.md
  data_dictionary.tsv
  requirements.txt
  SHA256SUMS.txt
  data/raw/lab4-cbioportal-luad/
    cohort.db
    SOURCE.md
  sql/
    patient_clinical.sql
    patient_gene_counts.sql
    site_stage_counts.sql
  scripts/
    build_features.py
    verify_reproducibility.py
  outputs/
    features.tsv
    patient_gene_counts.tsv
    site_stage_counts.tsv
    join_audit.tsv
    patient_check_mutations.tsv
    patient_check_features.tsv
    schema.json
    validation.json
```

`AI_USAGE.md` includes actual prompts used for assistance and debugging
observed errors. `verify_reproducibility.py` checks assigned file
strings or hashes before rebuilding.

## Sources

- [cBioPortal
  study](https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018):
  Lung Adenocarcinoma (TCGA, PanCancer Atlas), study ID
  `luad_tcga_pan_can_atlas_2018`
- [cBioPortal clinical sample
  metadata](https://github.com/cBioPortal/datahub/blob/master/public/luad_tcga_pan_can_atlas_2018/data_clinical_sample.txt):
  tissue source site describes the collecting institution. Source
  metadata were inspected without replacing the supplied database
