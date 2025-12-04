# Benchmark Approach Workflows

## Corrected Multi-Platform Workflows

All approaches assume **multiple genotyping platforms/batches** as input.

---

## APPROACH A: Traditional - Thorough QC Before (Per-Platform)

**Philosophy**: QC each platform thoroughly, then merge, then impute.
**This is what most GWAS do.**

```
Platform 1 ──┬──► Thorough QC ──┬──► INTERSECT ──► Ref Align ──► Impute ──► R² filter ──► Basic QC
Platform 2 ──┤                  │
Platform 3 ──┘                  │
                                └── (QC per-platform BEFORE merge)

STEPS:
1. Per-platform thorough QC:
   - Call rate 95%
   - MAF 1%
   - HWE SKIPPED (or optional)
   - Het ±3 SD
   - Relatedness

2. INTERSECT variants across platforms

3. QC on intersected data:
   - Re-check call rate after merge
   - Final relatedness check

4. Reference alignment (Rayner script)

5. Liftover to hg38

6. Impute (TOPMed/AllOfUs/Michigan)

7. Post-imputation:
   - R² > 0.3 filter
   - Basic call rate check
```

---

## APPROACH B: Southam 2011 - Minimal QC Before, Thorough After

**Philosophy**: Don't filter variants before imputation (hurts imputation quality).
**Based on**: Southam et al. 2011 EJHG

```
Platform 1 ──┬──► Minimal QC ──┬──► INTERSECT ──► Ref Align ──► Impute ──► R² filter ──► THOROUGH QC
Platform 2 ──┤   (call rate   │                                              │
Platform 3 ──┘    only 95%)   │                                              └── MAF, Het, Relatedness
                              │
                              └── (Minimal QC - just call rate)

STEPS:
1. Per-platform minimal QC:
   - Call rate 95% only
   - NO MAF filter
   - NO HWE filter

2. INTERSECT variants across platforms

3. Reference alignment (Rayner script)

4. Liftover to hg38

5. Impute

6. Post-imputation THOROUGH QC:
   - R² > 0.3 filter
   - Call rate 95%
   - MAF 1%
   - Het ±3 SD
   - Relatedness
```

---

## APPROACH C: Intersect First, Then QC (Common Practice)

**Philosophy**: Get common variants first, then QC on combined data.
**Ensures identical variant sets across platforms.**

```
Platform 1 ──┬──► INTERSECT FIRST ──► Thorough QC ──► Ref Align ──► Impute ──► R² filter ──► Basic QC
Platform 2 ──┤        │
Platform 3 ──┘        └── (Find common variants BEFORE any QC)

STEPS:
1. INTERSECT variants across platforms FIRST

2. Merge samples on common variants

3. Thorough QC on merged data:
   - Call rate 95%
   - MAF 1%
   - HWE SKIPPED
   - Het ±3 SD
   - Relatedness

4. Reference alignment (Rayner script)

5. Liftover to hg38

6. Impute

7. Post-imputation:
   - R² > 0.3 filter
   - Basic call rate check
```

---

## APPROACH D: Intersect First, QC After Imputation (Verma/Charon)

**Philosophy**: Intersect for consistency, but delay QC until after imputation.
**Based on**: Verma 2014 workflow + Charon 2021 QC timing findings

```
Platform 1 ──┬──► INTERSECT FIRST ──► Minimal QC ──► Ref Align ──► Impute ──► R² filter ──► THOROUGH QC
Platform 2 ──┤        │
Platform 3 ──┘        └── (Intersect first, minimal pre-imputation QC)

STEPS:
1. INTERSECT variants across platforms

2. Merge samples on common variants

3. Minimal QC:
   - Call rate 95% only

4. Reference alignment (Rayner script)

5. Liftover to hg38

6. Impute

7. Post-imputation THOROUGH QC:
   - R² > 0.3 filter
   - Call rate 95%
   - MAF 1%
   - Het ±3 SD
   - Relatedness
```

---

## APPROACH E: Our Pipeline - 1-Step

**Philosophy**: Union within platform, impute separately, then intersect.
**Key innovations**: MagicalRsq-X, skip HWE, union merge preserves rare variants.

```
Platform 1 ──► Union batches ──► Minimal QC ──► Ref Align ──► Impute ──┐
Platform 2 ──► Union batches ──► Minimal QC ──► Ref Align ──► Impute ──┼──► INTERSECT ──► MagicalRsq-X ──► Thorough QC
Platform 3 ──► Union batches ──► Minimal QC ──► Ref Align ──► Impute ──┘
              (UNION within)    (call rate)                (separate)     (after impute)   (ancestry-   (final)
                                                                                           calibrated)

STEPS:
1. Per-platform: UNION merge batches (keep all variants)

2. Per-platform: Minimal QC (call rate 95%)

3. Per-platform: Reference alignment

4. Per-platform: Impute SEPARATELY

5. INTERSECT imputed data across platforms

6. MagicalRsq-X filter (ancestry-calibrated, replaces R²)

7. Thorough QC:
   - Call rate 95%
   - MAF 1%
   - HWE SKIPPED
   - Het ±3 SD
   - Relatedness (GENESIS PCRelate)
```

---

## APPROACH F: Our Pipeline - 2-Step (TRUE APPROACH)

**Philosophy**: Two imputation passes to recover variants lost in merge.
**Key**: MagicalRsq-X after EACH imputation step.

```
Platform 1 ──► Union ──► Min QC ──► Ref Align ──► IMPUTE 1 ──► MagicalRsq-X ──┐
Platform 2 ──► Union ──► Min QC ──► Ref Align ──► IMPUTE 1 ──► MagicalRsq-X ──┼──► INTERSECT ──► IMPUTE 2 ──► MagicalRsq-X ──► Thorough QC
Platform 3 ──► Union ──► Min QC ──► Ref Align ──► IMPUTE 1 ──► MagicalRsq-X ──┘
                                                    │                               (merge)      (re-impute   (ancestry-     (final)
                                                    │                                             fills gaps)  calibrated)
                                                    └── MagicalRsq-X after 1st impute

STEPS:
1. Per-platform: UNION merge batches

2. Per-platform: Minimal QC (call rate 95%)

3. Per-platform: Reference alignment

4. Per-platform: IMPUTE (pass 1)

5. Per-platform: MagicalRsq-X filter (ancestry-calibrated)

6. INTERSECT across platforms

7. RE-IMPUTE (pass 2) - fills gaps from intersect merge

8. MagicalRsq-X filter (ancestry-calibrated)

9. Thorough QC:
   - Call rate 95%
   - MAF 1%
   - HWE SKIPPED
   - Het ±3 SD
   - Relatedness (GENESIS PCRelate)
```

---

## Summary Comparison Table

| Approach | Pre-Impute QC | Merge Strategy | QC Timing | Post-Impute Filter | Re-imputation |
|----------|---------------|----------------|-----------|-------------------|---------------|
| **A** | Thorough (per-platform) | Intersect | Before + Basic After | R² > 0.3 | No |
| **B** | Minimal (call rate) | Intersect | After only (thorough) | R² > 0.3 | No |
| **C** | Thorough (after intersect) | Intersect First | Before + Basic After | R² > 0.3 | No |
| **D** | Minimal (call rate) | Intersect First | After only (thorough) | R² > 0.3 | No |
| **E** | Minimal (call rate) | Union→Impute→Intersect | After only (thorough) | MagicalRsq-X | No |
| **F** | Minimal (call rate) | Union→Impute→Intersect→Impute | After only (thorough) | MagicalRsq-X ×2 | Yes |

---

## Standardized Settings (All Approaches)

| Setting | Value | Notes |
|---------|-------|-------|
| Call rate | 95% (geno=0.05, mind=0.05) | Standard |
| MAF | 1% when applied | Standard |
| HWE | **SKIPPED** | For mixed cohorts |
| Relatedness | KING > 0.125 | Same for all |
| Reference alignment | Rayner script | All approaches |
| R² threshold | 0.3 | Traditional (A-D) |
| MagicalRsq-X | 0.3 (calibrated) | Our pipeline (E-F) |
