-- Output grain: one observed tissue-source-site x stage combination.
-- NULL stage is retained as its own GROUP BY bucket.
SELECT s.tissue_source_site, p.stage,
       COUNT(DISTINCT p.patient_id) AS patient_count
FROM patients p
LEFT JOIN samples s ON s.patient_id = p.patient_id
GROUP BY s.tissue_source_site, p.stage
ORDER BY s.tissue_source_site, p.stage;
