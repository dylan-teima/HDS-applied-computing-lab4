#!/usr/bin/env python3
"""Read an existing SQLite database; query in SQL, reshape in Python, verify."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import platform
import sqlite3
import sys

ROOT = Path(__file__).resolve().parents[1]
GENES = ['TP53', 'KRAS', 'EGFR', 'STK11', 'KEAP1']
CLASSES = {'Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins',
           'Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation',
           'Splice_Site','Translation_Start_Site','Silent'}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_tsv(path, columns, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(columns)
        for row in rows:
            w.writerow(['NA' if v is None else v for v in row])

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--db', type=Path, default=ROOT/'data/raw/lab4-cbioportal-luad/cohort.db')
    ap.add_argument('--out', type=Path, default=ROOT/'outputs')
    args = ap.parse_args()
    db = args.db.resolve()
    require(db.is_file(), f'Database not found: {db}; supply the existing cohort.db, do not rebuild it.')
    digest = sha(db)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db.as_uri()+'?mode=ro', uri=True)
    con.row_factory = sqlite3.Row
    def rows(q): return [dict(r) for r in con.execute(q)]
    def scalar(q): return con.execute(q).fetchone()[0]
    def query(name): return rows((ROOT/'sql'/name).read_text())

    require(scalar('PRAGMA integrity_check') == 'ok', 'Database integrity check failed')
    require(not rows('PRAGMA foreign_key_check'), 'Foreign-key orphans found')
    schema = {}
    table_counts = {}
    for t in ['patients','samples','mutations','genes']:
        table_counts[t] = scalar(f'SELECT COUNT(*) FROM {t}')
        schema[t] = {'columns': rows(f'PRAGMA table_info({t})'),
                     'foreign_keys': rows(f'PRAGMA foreign_key_list({t})'),
                     'ddl': scalar(f"SELECT sql FROM sqlite_master WHERE name='{t}'")}
    (out/'schema.json').write_text(json.dumps(schema, indent=2)+'\n')
    n = table_counts['patients']
    require(n == 566, f'Expected this assignment database to contain 566 patients, found {n}')
    samples = rows('SELECT * FROM samples ORDER BY sample_id')
    patients = rows('SELECT * FROM patients ORDER BY patient_id')
    sample_counts = Counter(s['patient_id'] for s in samples)
    require(max(sample_counts.values(), default=0) <= 1,
            'Multiple samples per patient: define a sample-selection policy before building clinical features.')
    classes = rows('SELECT variant_classification, COUNT(*) AS n FROM mutations GROUP BY 1')
    require(all(r['variant_classification'] in CLASSES for r in classes),
            'Unknown/NULL mutation classification: review protein-altering definition')
    require(scalar('SELECT COUNT(*) FROM patients WHERE patient_id IS NULL') == 0, 'Missing patient key')
    clinical = query('patient_clinical.sql')
    site_expected = len({(s['tissue_source_site'], p['stage'])
                         for p in patients for s in samples if s['patient_id']==p['patient_id']}
                        | {(None,p['stage']) for p in patients if not sample_counts[p['patient_id']]})
    site = query('site_stage_counts.sql')
    long = query('patient_gene_counts.sql')
    require(len(clinical)==n and len({r['patient_id'] for r in clinical})==n, 'Clinical join changed grain')
    require(len(long)==n*5, 'Patient-gene grid incomplete')
    require(len(site)==site_expected, 'Unexpected site-stage grouping count')
    require(sum(r['patient_count'] for r in site)==n, 'Site-stage totals do not reconcile')

    # Python reshaping: no SQL pivot and no imputation of clinical fields.
    wide = defaultdict(dict)
    for r in long:
        require(r['gene_symbol'] not in wide[r['patient_id']], 'Duplicate patient-gene key')
        wide[r['patient_id']][r['gene_symbol']] = r['mutation_count']
    features = [{**r, **{g:wide[r['patient_id']][g] for g in GENES}} for r in clinical]
    require({r['patient_id'] for r in features} == {p['patient_id'] for p in patients}, 'Lost patient IDs')
    for r in features:
        require(all(isinstance(r[g],int) and r[g]>=0 for g in GENES), 'Invalid gene count')

    # Independent reference: iterate raw tables, join via Python dictionaries,
    # and apply the recorded classification rule without reusing SQL outputs.
    sample_patient = {s['sample_id']:s['patient_id'] for s in samples}
    independent = Counter()
    raw_sample_n = Counter()
    altered_symbols = 0
    for m in con.execute('SELECT sample_id,gene_symbol,variant_classification FROM mutations'):
        raw_sample_n[m['sample_id']] += 1
        g = m['gene_symbol'].strip().upper()
        altered_symbols += g != m['gene_symbol']
        if g in GENES and m['variant_classification'] != 'Silent':
            independent[(sample_patient[m['sample_id']],g)] += 1
    require(all(r[g]==independent[(r['patient_id'],g)] for r in features for g in GENES),
            'SQL/Python independent mutation counts disagree')
    lookup = {s['patient_id']:s for s in samples}
    patient_map = {p['patient_id']:p for p in patients}
    for r in features:
        p = patient_map[r['patient_id']]
        for k in ['sex','age_at_diagnosis','stage']:
            require(r[k]==p[k], f'Clinical value changed: {k}')
        for k in ['tissue_source_site','tmb_nonsynonymous']:
            require(r[k]==lookup.get(r['patient_id'],{}).get(k), f'Sample value changed: {k}')
    independent_site = Counter((lookup.get(p['patient_id'],{}).get('tissue_source_site'),p['stage']) for p in patients)
    require(independent_site==Counter({(r['tissue_source_site'],r['stage']):r['patient_count'] for r in site}),
            'Independent site-stage counts disagree')
    no_mut = [s for s in samples if raw_sample_n[s['sample_id']]==0]
    manual_id = next(r['patient_id'] for r in features if any(r[g]>1 for g in GENES))
    missing_id = next(p['patient_id'] for p in patients if p['stage'] is None)
    silent_sample = scalar("SELECT sample_id FROM mutations WHERE gene_symbol IN ('TP53','KRAS','EGFR','STK11','KEAP1') AND variant_classification='Silent' ORDER BY mutation_id LIMIT 1")
    check_ids = list(dict.fromkeys([manual_id, no_mut[0]['patient_id'], missing_id, sample_patient[silent_sample]]))
    placeholders = ','.join('?' for _ in check_ids)
    # Raw rows for inspectable patient checks; intentionally not aggregated.
    raw = [dict(r) for r in con.execute(f'''SELECT s.patient_id,m.* FROM samples s
        JOIN mutations m ON m.sample_id=s.sample_id
        WHERE s.patient_id IN ({placeholders})
          AND UPPER(TRIM(m.gene_symbol)) IN ('TP53','KRAS','EGFR','STK11','KEAP1')
        ORDER BY s.patient_id,m.gene_symbol,m.mutation_id''',check_ids)]
    write_tsv(out/'patient_check_mutations.tsv', ['patient_id']+list(schema['mutations']['columns'][i]['name'] for i in range(len(schema['mutations']['columns']))),
              [[r[k] for k in ['patient_id']+[c['name'] for c in schema['mutations']['columns']]] for r in raw])
    write_tsv(out/'patient_check_features.tsv', list(features[0]), [[r[k] for k in features[0]] for r in features if r['patient_id'] in check_ids])
    write_tsv(out/'features.tsv', list(features[0]), [[r[k] for k in features[0]] for r in features])
    write_tsv(out/'patient_gene_counts.tsv', ['patient_id','gene_symbol','mutation_count'], [[r[k] for k in ['patient_id','gene_symbol','mutation_count']] for r in long])
    write_tsv(out/'site_stage_counts.tsv', ['tissue_source_site','stage','patient_count'], [[r[k] for k in ['tissue_source_site','stage','patient_count']] for r in site])

    audit=[]
    def record(step, before, expected, got, grain):
        require(expected==got, f'Unexpected join count at {step}')
        audit.append([step,before,expected,got,grain])
    record('patients LEFT JOIN samples', n, sum(max(1,sample_counts[p['patient_id']]) for p in patients),len(clinical),'patient (at most one sample verified)')
    record('patients CROSS JOIN five target genes', n,n*5,scalar("SELECT COUNT(*) FROM patients CROSS JOIN (SELECT 'TP53' UNION ALL SELECT 'KRAS' UNION ALL SELECT 'EGFR' UNION ALL SELECT 'STK11' UNION ALL SELECT 'KEAP1')"),'patient x gene')
    record('patient-gene grid LEFT JOIN samples',n*5,len(clinical)*5,scalar('SELECT COUNT(*)*5 FROM patients p LEFT JOIN samples s ON s.patient_id=p.patient_id'),'patient x gene x sample')
    eligible = scalar("SELECT COUNT(*) FROM mutations WHERE variant_classification <> 'Silent' AND UPPER(TRIM(gene_symbol)) IN ('TP53','KRAS','EGFR','STK11','KEAP1')")
    observed = scalar("SELECT COUNT(*) FROM (SELECT sample_id,UPPER(TRIM(gene_symbol)) FROM mutations WHERE variant_classification <> 'Silent' AND UPPER(TRIM(gene_symbol)) IN ('TP53','KRAS','EGFR','STK11','KEAP1') GROUP BY 1,2)")
    record('LEFT JOIN sample-gene aggregate (before final GROUP BY)',n*5,n*5,scalar("""WITH mc AS (SELECT sample_id,UPPER(TRIM(gene_symbol)) g,COUNT(*) n FROM mutations WHERE variant_classification <> 'Silent' AND UPPER(TRIM(gene_symbol)) IN ('TP53','KRAS','EGFR','STK11','KEAP1') GROUP BY 1,2), tg(g) AS (VALUES ('TP53'),('KRAS'),('EGFR'),('STK11'),('KEAP1')) SELECT COUNT(*) FROM patients p CROSS JOIN tg LEFT JOIN samples s ON s.patient_id=p.patient_id LEFT JOIN mc ON mc.sample_id=s.sample_id AND mc.g=tg.g"""),'patient x gene x sample')
    record('Raw-record patient check JOIN',len(check_ids),sum(1 for m in con.execute('SELECT sample_id,gene_symbol FROM mutations') if sample_patient[m['sample_id']] in check_ids and m['gene_symbol'].strip().upper() in GENES),len(raw),'selected mutation record')
    record('Task 2 patient-gene result',n*5,n*5,len(long),'patient x gene')
    record('Task 2 site-stage result',len(clinical),site_expected,len(site),'observed site x stage')
    record('Python wide reshape and combine',n,n,len(features),'patient')
    write_tsv(out/'join_audit.tsv',['step','rows_before','expected_rows_after','actual_rows_after','grain_after'],audit)
    summary = {
        'table_counts':table_counts,'db_sha256':digest,
        'missing_features':{k:sum(r[k] is None for r in features) for k in features[0]},
        'genes':{g:{'mutation_count':sum(r[g] for r in features), 'patients_with_recorded_mutation':sum(r[g]>0 for r in features)} for g in GENES},
        'mutation_classes':classes,'noncanonical_gene_symbol_rows':altered_symbols,
        'eligible_mutation_rows':eligible,'observed_sample_gene_groups':observed,
        'samples_with_no_mutations':[s['sample_id'] for s in no_mut],
        'patient_checks':check_ids,'site_stage_groups':len(site),
        'checks':{'integrity':'ok','foreign_key_violations':0,'independent_patient_gene_comparisons':n*5,
                  'independent_site_stage_comparisons':len(site),'patient_ids_preserved':True},
        'runtime':{'python':platform.python_version(),'sqlite':sqlite3.sqlite_version},
    }
    require(sha(db)==digest, 'Input database changed')
    (out/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    con.close()
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
