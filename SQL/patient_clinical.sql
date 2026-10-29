-- One row per patient only after checking <= 1 sample per patient.
SELECT p.patient_id, p.sex, p.age_at_diagnosis, p.stage,
       s.tissue_source_site, s.tmb_nonsynonymous
FROM patients p
LEFT JOIN samples s ON s.patient_id = p.patient_id
ORDER BY p.patient_id;
