# Lab 4 - AI Usage

AI Assistant Used: ChatGPT (GPT-6.1 Sol model) / Codex

## Prompt 1 — Draft SQL while preserving the cohort

### **Prompt:**

> *"I have SQLite tables patients(patient_id), samples(sample_id,
> patient_id), and mutation (mutation_id, sample_id, gene_symbol,
> variant_classification). Write a SQLite query with one row for every
> patient and each of TP53, KRAS, EGFR, STK11, and KEAP1. Count mutation
> records excluding Silent, including patients with zero. Explain the
> grain of each step and predict the result row count for 566 patients.
> Do not pivot in SQL."*

**Accepted:** Generate the complete patient × five-gene grid, aggregate
eligible mutations to sample–gene grain, and left-join the aggregate.
Use `COALESCE` for absent recorded counts. The independent aggregation's
`WHERE` defines mutation eligibility and optional-side conditions of the
cohort-preserving joins remain in `ON`.

**Rejected:** Start from the mutation table or use an outer `WHERE` to
require a matching mutation. That excludes patients without matching
records; An additional suggestion involved aggregation with SQL only as
the final pivot, which was rejected since the task requires Python or R
reshaping.

**Verification:**

- Expected and observed long result: 566 × 5 = 2,830 rows.
- Every patient has exactly five distinct gene entries (final feature
  key is unique).
- Expected and observed final wide table: 566 patients and 11 columns.
- Independently iterating raw mutation records in Python agrees with all
  2,830 SQL-derived counts.
- Raw-record inspection of TCGA-05-4382: EGFR IDs 2635 and 2636 give
  EGFR=2; TP53 ID 1770 gives TP53=1. Both values match the feature row.

## Prompt 2 — LEFT JOIN silently loses zero-count patients

### **Prompt:**

> *"This query runs but gives 168 patients when I expected 566:*
>
> ``` sql
> SELECT p.patient_id, COUNT(m.mutation_id)
> FROM patients p 
> LEFT JOIN samples s ON s.patient_id=p.patient_id 
> LEFT JOIN mutations m 
>   ON m.sample_id=s.sample_id 
> WHERE m.gene_symbol='KRAS' AND m.variant_classification<>'Silent' 
> GROUP BY p.patient_id; 
> ```
>
> *Explain why the row count changes. Show a corrected query and explain
> what to count so patients with no KRAS mutation get zero rather than
> one."*

**Accepted:** Illustrated by the delivered approach: keep the patient
table on the preserved side and move conditions on optional matches into
`ON`, or independently pre-aggregate eligible mutation records before
left-joining them. Count a non-NULL mutation key, not the placeholder
row.

A simple corrected diagnostic query is:

``` sql
SELECT p.patient_id, COUNT(m.mutation_id) AS kras_count
FROM patients p
LEFT JOIN samples s ON s.patient_id=p.patient_id
LEFT JOIN mutations m
  ON m.sample_id=s.sample_id
 AND UPPER(TRIM(m.gene_symbol))='KRAS'
 AND m.variant_classification<>'Silent'
GROUP BY p.patient_id;
```

**Rejected:** `DISTINCT` on patient rows without fixing eligibility
logic; an inner join; filling every missing value in the entire feature
table with zero. None preserves the intended meaning of both clinical NA
and mutation zero.

**Verification:** the bad filter returns 168 distinct patients, while
the final pipeline preserves all 566. A deliberately tested `COUNT(*)`
for TCGA-64-1678 produces 1 despite no KRAS mutation row and reported a
count of 0. There are exactly four samples with no mutation records, and
their patients remain in the final table.

## Prompt 3 — Gene name mismatch

### **Prompt:**

> *"Some mutation gene symbols contain lowercase letters or surrounding
> spaces. How can I check whether exact matches are losing records,
> normalize the key without editing raw data, and avoid creating a
> complex many-to-many join? Please distinguish checking foreign-key
> validity from gene-name standardization."*

**Accepted:** inspect exact and normalized keys and use
`UPPER(TRIM(gene_symbol))` inside the mutation aggregation. Do not
mutate the source or join on a normalized key without checking its
uniqueness.

**Rejected:** all 2,431 noncanonical mutation symbols are unmatched
foreign keys. In this file they have exact counterparts in `genes`, so
the foreign-key check passes and the exact mutation-to-gene join matches
all 200,789 rows. There are no normalized gene-key collisions. None of
the noncanonical symbols belongs to the five target genes; normalization
therefore does not change the five-gene totals for this input.

**Verification:** schema inspection, `PRAGMA foreign_key_check`, exact
join count, normalized-key duplicate check, and independent Python
normalization during all 2,830 reference comparisons. Normalization
addresses case/space only; no synonym/alias remapping is claimed.

## Prompt 4 — Count patients rather than mutation rows

### **Prompt:**

> *"I am attempting to determine patients per tissue source site and
> stage. A sample can have many mutation rows. Which tables are needed
> and how to retain patients with NULL stage. Before querying, how can I
> predict the expected number of site-stage groups without assuming
> every site has every stage?"*

**Accepted:** Query patients LEFT JOIN samples, group by site and stage,
and count distinct patient IDs. Do not join mutations. Independently
enumerate observed category pairs in Python before executing the grouped
SQL.

**Rejected:** Joining mutations to derive the patient denominator;
dropping NULL stages; generating a full site × stage Cartesian product
and describing its size as the expected observed-group count.

**Verification:** Predicted 142 observed groups, matching the observed
output of 142. Group counts sum to 566. Independent Python grouping
agrees with every group. Washington University / missing stage contains
52 patients. One-sample-per-patient was verified, making group
membership exclusive in this input.

## Prompt 5 — Missing versus zero, and Silent exclusion

### **Prompt:**

> *"My clinical data have missing sex, age, and stage. Four samples have
> no mutation rows, and most patients have no record for a specified
> target gene. How might these missing data be handled without dropping
> patients. How can I verify Silent mutations are excluded while
> clinical NULLs remain missing?"*

**Accepted:** Retain clinical NULLs as NA, use zero for the defined
count of qualifying records, and explicitly list the four samples with
no mutation rows and state the coverage limitation. Preserve source TMB,
including zeros.

**Verification:**

- Before and after clinical joins and final feature combination: 566
  patients.
- Final missing clinical counts: sex=52, age=71, stage=54; checks
  confirm that source clinical values were unchanged.
- TCGA-17-Z023 has TP53 ID 21335 (Nonsense_Mutation, p.E198\*) and ID
  21336 (Silent, p.V197=). Its TP53 feature equals 1, not 2.
- TCGA-17-Z000 retains NA sex, age, and stage while remaining in the
  final table.
- TCGA-64-1678 retains its feature row and five zero recorded counts,
  with its sample listed among the four lacking mutation records.

## Reproducibility Evidence

The build script records exact counts in `outputs/join_audit.tsv`,
schema in `outputs/schema.json`, and cohort-level validation in
`outputs/validation.json`. The patient-check TSVs expose the raw records
behind the reviewed examples. `scripts/verify_reproducibility.py`
validates packaged checksums and rebuilds into a temporary directory,
comparing every output TSV and schema byte-for-byte and comparing
validation content except runtime versions.
