# Benchmark Approach Workflows

## Key Differentiators

| Factor | Options | Approaches |
|--------|---------|------------|
| **Merge timing** | BEFORE vs AFTER imputation | A,B = before; C,D,E,F = after |
| **QC timing** | BEFORE vs AFTER imputation | A,C = before; B,D,E,F = after |
| **Quality filter** | R² vs MagicalRsq-X | A,B,C,D = R²; E,F = MagicalRsq-X |

## Summary Table

| Approach | Merge Timing | QC Timing | Filter | Key Test |
|----------|--------------|-----------|--------|----------|
| **A** | BEFORE impute | BEFORE impute (thorough) | R² > 0.3 | Traditional baseline |
| **B** | BEFORE impute | AFTER impute (thorough) | R² > 0.3 | Southam 2011 - QC timing |
| **C** | **AFTER impute** | **BEFORE impute** (thorough) | R² > 0.3 | **QC before, merge after** |
| **D** | AFTER impute | AFTER impute | R² > 0.3 | Same as ours but R² filter |
| **E** | AFTER impute | AFTER impute | MagicalRsq-X | Our 1-step |
| **F** | AFTER impute | AFTER impute | MagicalRsq-X ×2 | Our 2-step |

---

## APPROACH A: Traditional - Merge Before, QC Before

**Philosophy**: QC each platform thoroughly, then merge, then impute.
**This is the standard GWAS workflow.**

```
Platform 1 ──► Thorough QC ──┐
Platform 2 ──► Thorough QC ──┼──► INTERSECT/MERGE ──► Ref Align ──► Impute ──► R² > 0.3 ──► Basic QC
Platform 3 ──► Thorough QC ──┘
              (per-platform)      (BEFORE impute)                    (traditional)
```

**Steps:**
1. Per-platform thorough QC (95% call rate, MAF, het, relatedness)
2. INTERSECT variants across platforms
3. Merge samples
4. Reference alignment (Rayner script)
5. Liftover → Impute
6. R² > 0.3 filter
7. Basic post-imputation QC (call rate check)

---

## APPROACH B: Southam 2011 - Merge Before, QC After

**Philosophy**: Don't filter variants before imputation (hurts imputation quality).
**Based on**: Southam et al. 2011 EJHG

```
Platform 1 ──► Minimal QC ──┐
Platform 2 ──► Minimal QC ──┼──► INTERSECT/MERGE ──► Ref Align ──► Impute ──► R² > 0.3 ──► Thorough QC
Platform 3 ──► Minimal QC ──┘
              (call rate     (BEFORE impute)                    (traditional)   (MAF, het,
               only 95%)                                                         relatedness)
```

**Steps:**
1. Per-platform minimal QC (95% call rate ONLY)
2. INTERSECT variants across platforms
3. Merge samples
4. Reference alignment
5. Liftover → Impute
6. R² > 0.3 filter
7. **Thorough** post-imputation QC (MAF, het, relatedness)

---

## APPROACH C: Thorough QC Before → Merge After Imputation

**Philosophy**: Apply thorough QC on each platform BEFORE imputation, impute separately, then merge.
**KEY COMPARISON**: Same as D but with thorough QC BEFORE imputation instead of after.

```
Platform 1 ──► THOROUGH QC ──► Ref Align ──► Impute ──┐
Platform 2 ──► THOROUGH QC ──► Ref Align ──► Impute ──┼──► INTERSECT/MERGE ──► R² > 0.3 ──► Light QC
Platform 3 ──► THOROUGH QC ──► Ref Align ──► Impute ──┘
              (before impute)              (separately)    (AFTER impute)    (traditional)   (call rate)
```

**Steps:**
1. Per-platform THOROUGH QC (95% call rate, het, relatedness)
2. Per-platform reference alignment
3. Per-platform imputation (SEPARATELY)
4. INTERSECT/MERGE across platforms (AFTER imputation)
5. R² > 0.3 filter (traditional)
6. Light post-imputation QC (call rate only - thorough was done before)

**C vs D comparison**: Does QC timing matter? (before vs after imputation, both merge after)

---

## APPROACH D: Merge AFTER Imputation (R² filter)

**Philosophy**: Impute each platform separately, then merge. Uses traditional R² filter.
**KEY COMPARISON**: Same workflow as E/F but with R² instead of MagicalRsq-X.

```
Platform 1 ──► Minimal QC ──► Ref Align ──► Impute ──┐
Platform 2 ──► Minimal QC ──► Ref Align ──► Impute ──┼──► INTERSECT/MERGE ──► R² > 0.3 ──► Thorough QC
Platform 3 ──► Minimal QC ──► Ref Align ──► Impute ──┘
              (call rate)                (separately)    (AFTER impute)    (traditional)
```

**Steps:**
1. Per-platform minimal QC (95% call rate)
2. Per-platform reference alignment
3. Per-platform imputation (SEPARATELY)
4. INTERSECT/MERGE across platforms (AFTER imputation)
5. R² > 0.3 filter (traditional)
6. Thorough post-imputation QC

**D vs E/F shows the value of MagicalRsq-X over traditional R².**

---

## APPROACH E: Our Pipeline - 1-Step (MagicalRsq-X)

**Philosophy**: Merge after imputation with ancestry-calibrated quality filter.

```
Platform 1 ──► Minimal QC ──► Ref Align ──► Impute ──┐
Platform 2 ──► Minimal QC ──► Ref Align ──► Impute ──┼──► INTERSECT/MERGE ──► MagicalRsq-X ──► Thorough QC
Platform 3 ──► Minimal QC ──► Ref Align ──► Impute ──┘
              (call rate)                (separately)    (AFTER impute)    (ancestry-
                                                                           calibrated)
```

**Steps:**
1. Per-platform minimal QC (95% call rate)
2. Per-platform reference alignment
3. Per-platform imputation (SEPARATELY)
4. INTERSECT/MERGE across platforms (AFTER imputation)
5. **MagicalRsq-X** filter (ancestry-calibrated)
6. Thorough post-imputation QC (HWE skipped for admixed)

---

## APPROACH F: Our Pipeline - 2-Step (MagicalRsq-X ×2)

**Philosophy**: Two imputation passes with MagicalRsq-X after each.
**Key**: Re-imputation fills gaps from merge.

```
Platform 1 ──► Min QC ──► Ref Align ──► IMPUTE 1 ──► MagicalRsq-X ──┐
Platform 2 ──► Min QC ──► Ref Align ──► IMPUTE 1 ──► MagicalRsq-X ──┼──► INTERSECT ──► Call rate ──► IMPUTE 2 ──► MagicalRsq-X ──► Thorough QC
Platform 3 ──► Min QC ──► Ref Align ──► IMPUTE 1 ──► MagicalRsq-X ──┘
                                        (1st pass)   (filter 1)      (merge)      (adjust)    (2nd pass)   (filter 2)     (final)
```

**Steps:**
1. Per-platform minimal QC (95% call rate)
2. Per-platform reference alignment
3. Per-platform IMPUTE (pass 1)
4. **MagicalRsq-X** filter (pass 1)
5. INTERSECT/MERGE across platforms
6. Call rate adjustment before 2nd imputation
7. RE-IMPUTE (pass 2) - fills gaps from merge
8. **MagicalRsq-X** filter (pass 2)
9. Thorough post-imputation QC

---

## Key Comparisons

| Comparison | Question |
|------------|----------|
| **A vs B** | Does QC timing matter when merging BEFORE impute? |
| **A/B vs C/D** | Does merge timing matter? (before vs after impute) |
| **C vs D** | Does QC timing matter when merging AFTER impute? **(KEY)** |
| **D vs E** | Does MagicalRsq-X improve over R²? **(KEY)** |
| **E vs F** | Does 2-step re-imputation help? |

### Merge Timing Groups
- **Merge BEFORE** (A, B): Traditional approach - intersect variants first
- **Merge AFTER** (C, D, E, F): Our approach - impute separately, then merge

---

## Standardized Settings (All Approaches)

| Setting | Value | Notes |
|---------|-------|-------|
| Call rate | 95% (geno=0.05, mind=0.05) | Standard |
| MAF | 1% when applied | Standard |
| HWE | **SKIPPED** | For mixed cohorts (option to enable) |
| Relatedness | KING > 0.125 | Same for all |
| Reference alignment | Rayner script | All approaches |
| R² threshold | 0.3 | Traditional (A-D) |
| MagicalRsq-X | 0.3 (calibrated) | Our pipeline (E-F) |

---

## Output Tracks

Each approach produces two output tracks:

1. **GWAS track**: MAF > 1% variants for standard GWAS
2. **RVAS track**: All variants (including rare) for rare variant analysis
