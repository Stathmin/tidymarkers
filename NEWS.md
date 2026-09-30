# kaspexample (kasp-example)

Демонстрационный пайплайн подбора KASP-маркеров для *Arabidopsis thaliana*
на чистом R. Здесь — краткая история значимых изменений; подробный инвентарь и
очередь работ — в `ROADMAP.md`.

## 0.3.0

### Дизайн
- Движок InDel: аллельные шаблоны и скан стыка (`design_as_primers_indel()`,
  `allele_templates()`); размещение выбирает физика (`DECIPHER`), не ручное правило.
- Substitution-блоки / MNP: вариант — строка REF/ALT, единая геометрия через
  `as_anchor()` / `block_in_seg()`.
- `design_as_primers_indel()` и SNP-путь унифицированы в `assemble_kasp_assay()`
  (пер-аллельный 3'-футпринт, ампликон — худший из двух аллелей).

### Специфичность
- `specificity_seed()`: точный путь `PDict`/`vcountPDict`, мисматч-скан с
  `BiocParallel`-хуком (`BPPARAM=`).
- Батч: `design_kasp_markers()` с вердиктом на локус и параллелью по локусам.

### Термодинамика / QC
- `terminal_dg()` — 3'-концевая стабильность (GC-clamp, `OligoArrayAux`).
- `wobble_mismatch()` — флаг слабого G:T / A:C неспаривания (симметрично по REF/ALT).

### Ввод / вывод
- `normalise_allele()` — IUPAC-коды сворачиваются к `ACGTN`, а не отбрасываются.
- `kasp_mix_recipe()` — tidy-рецепт смеси (AS:AS:COM = 2:2:5) и реакции.

## 0.2.0
- `kasp_design_checks()` и `KASP_LIMITS`; пороги LGC/EasyKASP.
- `in_silico_pcr()` (`matchProbePair`), SSR-скрин, `strand_gc_hint()`.
- Пакетный ввод VCF (`vroom`) и FASTA с нормализацией имён.

## 0.1.0
- Каркас `renv`, дизайн AS/COM через `DECIPHER`/`TmCalculator`/`OligoArrayAux`,
  библиотечная графика (`ggmsa`, `gggenes`, `ggrepel`, `patchwork`), экспорт IDT/LGC.
