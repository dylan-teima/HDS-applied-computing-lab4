-- Output grain: one patient x target gene, including zero counts.
-- WHERE clause included inside this independent aggregation to define eligible mutations;
-- it does not filter patients from the later LEFT JOIN.
WITH target_genes(gene_symbol) AS (
    VALUES ('TP53'), ('KRAS'), ('EGFR'), ('STK11'), ('KEAP1')
), mutation_counts AS (
    SELECT sample_id, UPPER(TRIM(gene_symbol)) AS gene_symbol,
           COUNT(mutation_id) AS mutation_count
    FROM mutations
    WHERE variant_classification <> 'Silent'
      AND UPPER(TRIM(gene_symbol)) IN ('TP53','KRAS','EGFR','STK11','KEAP1')
    GROUP BY sample_id, UPPER(TRIM(gene_symbol))
)
SELECT p.patient_id, g.gene_symbol,
       COALESCE(SUM(mc.mutation_count), 0) AS mutation_count
FROM patients p
CROSS JOIN target_genes g
LEFT JOIN samples s ON s.patient_id = p.patient_id
LEFT JOIN mutation_counts mc
  ON mc.sample_id = s.sample_id AND mc.gene_symbol = g.gene_symbol
GROUP BY p.patient_id, g.gene_symbol
ORDER BY p.patient_id, g.gene_symbol;
