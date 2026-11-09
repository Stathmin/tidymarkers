# tidymarkers — roadmap (открытая работа)

> Внутренний документ разработки: в собранный пакет не входит (`.Rbuildignore`).
> Пользовательская документация — `README.md` и виньетки.

Только незакрытое; история изменений — `NEWS.md`.

Обозначения: `[ ]` очередь · `[~]` частично · `[x]` сделано · `[—]` проверено
и **отложено** (с причиной) · `[!]` проверено и **опровергнуто**.

**Формат открытого пункта:** **правка** (файл и функция) · **проверка** (тест, id
в `dev/check-criteria.R`, артефакт) · **блокер** (что не даёт начать и чем
снимается) · **следствие** (что меняется, если пункт не делать). Пункт без
одного из четырёх остаётся в секции как описание, а не в очереди.

**О провенансе.** DOI/название — это провенанс, а не доказательство. Число
попадает в `inst/extdata/`, только если оно прочитано из скачанного файла и у
него есть `data-raw/<paper>.R`, который утверждает его против итогов статьи и
падает при расхождении. Где такого утверждения построить нельзя — файл не
пишется, причина записывается.

## Очередь

Открытая работа ведётся по секциям ниже, каждая с чекбоксом; `Y`/`YP` — тоже
очередь (пункты `[ ]`/`[~]` рабочие, не идеи). История — `NEWS.md`, провенанс
источников — `inst/extdata/science_ledger.csv`.

## Приоритет

Наука раньше архитектуры: `Z7` (намеренные вторые мисматричи) → `Z8`
(калибровка); `ZA6` шаги 1–4 — контракты данных, шаги 5–8 — после `Z8`.
Дисциплина «рефакторинг не сдвигает вердикт»: шаг принимается только с
неизменными вердиктами на золотом наборе.

## Решения

Все решения `Q1`/`Q2` приняты; каждое — `*.decision.md` в `data-raw/proposed/`
(статус, правка, тест). `.proposal.md` — запись доказательств, после решения не
правится. Пока решение не принято, артефакт говорит «данных нет», а не «данные
против».

| решение                  | итог                                                            |
| ------------------------ | --------------------------------------------------------------- |
| `arms_mismatch_position` | −3 по умолчанию; −2 — fallback для слабого терминального класса |
| `rpa_anchor_rule`        | якорь 4 nt; веса ранжируют, не предсказывают частоту            |
| `crispr_seed_rule`       | правило позиции; величина сигнала `unknown`                     |
| `hrm_resolvability`      | `dTm >= 0.5 °C` для LightCycler 96 и среднего по классу         |
| `kasp_z8_scope`          | конверсии не объединяются; `assay_x_sample` — не ассей          |

## Открытая работа — сводка

### §U. Другие PCR-маркеры и внешние приоры

- [x] **CAPS/dCAPS и ARMS/AS-PCR.** `R/marker_caps.R`, `R/marker_arms.R`;
      «смена нуклеазы» = изошизомеры (`restriction_enzymes()` отдаёт одну
      строку на пару (сайт, `cut`), остальные — в `$isoforms`);
      `design_arms()` падает на отрицательном `stab_ats`. `test-caps.R`,
      `test-arms.R`. Остаток — откат на запасную нуклеазу, отдельная правка.
- [x] **SNaPshot / single-base extension.** 19 критериев `family_sbe_*`
      (8 `primary_literature`, 10 `vendor_guideline`, 1 `rule_of_thumb`) в
      `criterion_spec()`; `SNAPSHOT_EXT_TM` → `c(50, 60)`,
      `SNAPSHOT_EXT_LENGTHS` → `17:25`; провенанс —
      `inst/extdata/snapshot_multiplex_design_rules.csv` (7 open-access статей).
      `test-snapshot*`. Провенанс — это rung, а не калибровка.
- [x] **Провенанс порогов + термодинамическая сверка.** Лестница
      `CRITERION_EVIDENCE` (`ZW4.1`) и `inst/extdata/criterion_provenance.csv`
      (212 строк: `id` → `source` + `reason_class` + причина; 41 названа, 171 —
      `package assumption` с явным классом). `lobato2018` (RPA ×4) влит.
      Внешний Tm-якорь — литературный пробел (`ZW3.4`), не пробел провенанса.

### §V. Модели считывания (readout simulation)

- [~] **HRM — плавление; dPCR — пуассоновское разбиение; CRISPR-Dx —
  коллатеральное расщепление; ARMS/CAPS/indel/SSR — виртуальный гель;
  LAMP/RPA — изотермическая кинетика.** `marker_readout()`
  (`R/marker_readout.R`) возвращает `simulated = TRUE` только для dPCR
  (`poisson_lambda()`, `poisson_positive_fraction()`, `dpcr_partition_table()`;
  `test-readout-dpcr.R` воспроизводит λ Dong и разведение MPS). У остальных
  четырёх химий — `FALSE` с причиной: HRM есть набор из 28 ампликонов, но нет
  Mg/соли (**partial**); LAMP есть Tp × температура, но `copies_per_reaction`
  = `NA` (**partial**, `Q1` не закрыт); CRISPR-Dx — Nalefski 2024 печатает
  kcat/KM (`readout_crispr_nalefski2024_pass2.conditions.csv`, 109 строк), но не
  ряд «концентрация → доля расщеплённого»; виртуального геля нет (bp→мм с
  % геля, вольтажом и временем).
    - **правка:** код начинается только когда есть набор «условия → исход»,
      который модель обязана воспроизвести; иначе `simulated = TRUE` — ложь.
    - **проверка:** тот же набор, что читает `test-readout-physics.R`.
    - **блокер:** `Q1` — наборы «условия → исход» в открытом доступе; отсутствие
      подтверждено поиском (`v-readout-conditions.md`), это не доказательство
      несуществования.
    - **следствие:** `simulated = FALSE` — честная граница, а не забытый код.

### §W. Повторяющийся аудит

- [x] DRY, дрейф документов, styler/lintr/`R CMD check --as-cran`, «дерево ==
      код», guard «single-base only» (критерий `family_single_base`), общие
      фикстуры, `in_silico_pcr()` с IUPAC/mismatch. Аудит живёт в
      `dev/check.R --quick` (38 s).

### §X. Слой плана

- [x] **Двухпроходная дизайн-проводка.** `two_pass_schemes()`/`two_pass()` в
      `R/marker_plan.R`; второй проход у snp, indel, multi, lax, duplex, caps,
      arms, hrm, indel_size, ssr, sanger. Меняет число кандидатов, не вердикт
      (`ZW0`). `dev/test.R plan` — 139 pass.

### Прочее

- [x] **CRAN-готовность: зависимости и документация сборки.** `DESCRIPTION`
      (`Authors@R` с ORCID, `Imports`/`Suggests` под реальные вызовы), `@param`
      для всех аргументов `\usage`, `.Rbuildignore`, фикстуры; `R CMD check
--as-cran` — 0 ERROR, 0 WARNING, 1 NOTE (New submission).
- [x] **`ZV`. CRAN-остаток, проход 1.** `.github/workflows/R-CMD-check.yaml` —
      пять ног (Ubuntu release/oldrel-1/devel, macOS release, Windows release),
      `error-on: warning`; `cran-comments.md`; `revdep/README.md`;
      `CITATION.cff` + «How to cite»; `inst/WORDLIST` (527) + шаг 4 гейта;
      тарбол 2.3 MB при лимите 5 MB.
- [x] **CRAN-остаток, проход 2.** `vignettes/index.Rmd` + ссылка в
      `docs/INDEX.md`; `https://www.vniisb.ru` не отвечает (`000`), поэтому в
      `README.md` — только-HTTP `http://www.vniisb.ru/ru/news/` (`200`);
      `dev/prose-tics.R` (семь уровней LLM-прозы), `--check` гейтит `filler`,
      `comment`, `verbosity`; `Suggests` перепроверены — все используются.
- [x] **Языковая чистота, проход 1.** `dev/lang_purity.md` — разбор уровней и
      очередь; `comment_verbose` (6, плотные пояснения ограничений),
      `prose_verbose` (0; четыре длинных абзаца в `markers-gel.Rmd`/
      `markers-kasp.Rmd` разбиты), оба гейтятся.
- [x] **Языковая чистота, проход 2.** `meta` → 0, `hedge` → 4 (все
      технические), `punctuation` → 0, `negation` (397) прочитан целиком,
      оставлен; гейт `encoding` (невалидный UTF-8, U+FFFD, мойбейк, zero-width,
      управляющие — 0); 14 WARN «Could not resolve link» исправлены причиной.
- [ ] **Windows/R-hub: сборка OligoArrayAux; бейджи/CI.**
      правка: `data-raw/install-oligoarrayaux.R` (unix-only) +
      `.github/workflows/R-CMD-check.yaml`. проверка: сборка на чистой машине;
      `dev/install-deps.R` ставит всё без ручных шагов. блокер: нет локального
      Windows. следствие: в CI движок недоступен → `capabilities()` = `unknown`
      (честно), но ни один внешний тест не исполняется.
- [x] **Скорость.** `dev/test.R --full` двухфазный: 126 файлов параллельно
      (`TESTTHAT_CPUS` = 8), затем 4 SOCK-файла серийно; список SOCK-файлов
      защищён grep-guard. Замер 2026-11-06: **8m 10s**, 131 файл, 36894 pass,
      0 fail (`dev/logs/latest-summary.md`).

## Z. Внешний научный аудит (окт. 2026) — код, требования и физика

Каждый пункт закрывается кодом, тестом и (где есть данные) калибровкой.
Приоритеты: **P0** — неверная биология/вердикт; **P1** — специфичность,
популяция, калибровка; **P2** — расширение методов.

### Z0. Фактические ошибки физики/биологии [P0]

- [x] **`Z0.0` Cas12a collateral — ssDNase, репортёр ssDNA.**
      `tests/testthat/test-wording.R` запрещает `dsDNase` и dsDNA-коллатераль
      по виньеткам; виньетки уже были корректны.
- [x] **`Z0.1` Cas13 — не «PAM-free» без оговорок.** `crispr_guide_candidates()`
      несёт `pfs`/`pfs_side` (Cas13a 3′ H, Cas13b 5′ D), Cas12a — `NA`; длина
      спейсера — параметр; `test-crispr*` 181 pass. `crisprBase` PFS не
      описывает — поле своё.
- [ ] **`Z0.2` dPCR Tm-шкала — калибровочный долг, не норма.** Стандарт:
      праймеры 58–62 °C, проба +8–10 °C; текущий default 44–50 °C.
    - правка: **окно не менять**, пока `ZW3.4` не объяснён (иначе подгонка под
      смещение); `@param tm_range` в `R/marker_dpcr.R` называет default
      калибровочным долгом.
    - проверка: `test-dpcr.R` — default `c(44, 50)` пиннится и назван долгом.
    - блокер: `ZW3.4` — смещение разложено: `Schildkraut` −16.9 °C, NN-консистентный
      `SantaLucia1998-2` −9.6 °C при 25 нМ, −5.4 °C при 1 мкМ тотальной нити
      (`dnac=`), остаток ≈−5.4 — метка/развёртка. Внешнего измеренного якоря на
      KASP/dPCR-условиях нет, поэтому окно не двигается. Перенос координаты на
      SL2 (с пере-якорением окон) — отдельный пункт `ZW3.10`.
    - следствие: смена окна «на глаз» замаскирует биас и сдвинет вердикты всех
      dPCR-ассays без данных.
- [x] **`Z0.3` LAMP F1c/B1c > F2/B2.** `lamp_fuse()` фильтрует `t1 > t2`; FIP 3′
      несёт ALT, а не REF.
- [x] **`Z0.4` KASP: «Taq не продлевает 3′ мисматч» — кинетика, не бинарность.**
      `marker_mismatch_severity()` (Kwok 1990) + строка «AS 3' mismatch base
      pair» с классом; `test-kasp-mismatch-type.R`, `test-kasp-checks.R` (70 pass).

### Z1–Z3. Вердикт, пороги, условия реакции [P0–P1] — закрыто

`Evaluation` — один конструктор; `NA` не проходит отбор; каждый критерий —
data-object (`R/marker_criteria.R`); `reaction_conditions()` доходит до call
sites. Детали — `NEWS.md`.

### Z4. Движки и capabilities [P1]

- [x] **`Z4.1` Нет движка — `unknown`, не тихий `pass`.** `lamp_geometry()`,
      `design_snapshot_marker()`, `pair_structure_ok()`/`rpa_structure_ok()`
      дают `NA` + причину без OligoArrayAux; `test-capabilities.R`,
      `test-lamp.R`, `test-snapshot.R` (493 pass).
- [x] **`Z4.2` Ecosystem reuse, not reimplementation.** Кандидаты: `crisprVerse`,
      `DECIPHER`, `rprimer`, `TmCalculator` (27 NN, salt/DMSO) как
      Tm-cross-check; `PolyMarker` — benchmark. Внешние Primer3/MFEprimer/
      NUPACK/ViennaRNA — только через `capabilities()`. Потребителя в коде пока
      нет (вторая зависимость без вызова = мёртвый код).

### Z5. Mismatch-tolerant specificity [P1]

- [x] **`Z5.0` `in_silico_pcr()`** ищет IUPAC-праймеры, `max_mismatch`, склеивает
      реальные продукты; degenerate покрыт тестом (`"ACGS"`).
- [x] **Политика и вердикт:** `specificity_policy()`, `specificity_verdict()`
      (seed/exact/BLAST/FM) с capability-флагом и `unknown`.
- [x] **`Z5.1` Position-aware off-target.** `primer_hit_table()` →
      `mismatch_profile()` + `plus_to_primer()`; 3′-терминальный промах штраф
      **12**, серединный — **1**; `test-specificity-orientation.R` (36 pass).
- [ ] **`Z5.2` Family-specific mismatch models.** RPA (несколько мисматчей),
      CRISPR (seed) и KASP (3′-якорь) сведены к одной политике.
    - правка: `specificity_policy()` — модель на семейство.
    - проверка: `test-specificity.R` по семействам.
    - блокер: нет пер-семейной предсказуемой величины. CRISPR — гейт по числу
      мисматчей/длине спейсера (3 мисматча в 23-nt → 11.0 % WT; 28-nt не
      различает одиночный), не fold-сигнал. RPA — только бинарная вероятность
      0.59–0.984. LAMP — бинарно + Tt. HRM — dTm локус-специфичен. CAPS — оси нет.
    - следствие: `specificity_verdict()` отвечает одной шкалой на разной физике.

### Z6. Популяционный слой [P1]

- [x] `marker_attach_haplotypes()` + `marker_population_checks()`; CAPS/CRISPR
      предупреждают без популяции и при вариантах в сайте/guide-PAM; нет
      частоты = `unknown`.
- [x] **`Z6.1` Гаплотип/частота во все схемы.** CAPS site, CRISPR guide/PAM,
      LAMP/RPA primer sites, SSR null alleles; приёмка — инвариант `ZW0`.
- [x] **`Z6.2` Полиплоидия — отдельная ось, не `pop_freq`.** `$dosage` — поле
      первого класса (`marker_dosage()`, `marker_dosage_gate()`); строка
      `allele dosage` в гейте; `test-polyploid-dosage.R` (204 pass),
      `test-polyploid-scope.R`. Движок — `Z10.2`.

### Z7. Намеренные вторые мисматричи [P1] — первый научный пункт после `ZW`

Приоры собраны как проверяемые предсказания: `inst/extdata/z7_priors.csv`
(`data-raw/z7-priors.R`, `test-z7-priors.R`) — восемь предсказаний с дословной
цитатой и `calibrated = FALSE`. Противоречия не усреднены: Ye 2001 (−2),
Liu 2012 (3-е от 3′-конца; offset −2), Lefever 2019 (4-е; offset −3) хранятся
все три.

Блокер измерен: `z7_engine_probe.csv` — против **матчевой** матрицы потеря
3′-концевого основания стоит 0.3–0.7 % эффективности DECIPHER; вся
дискриминация в **off**-аллеле (G:A 0.72, C:A 0.66, A:G и A:A ≈ 0). Порядок Kwok
воспроизводится, но это проверка **пробы**, не калибровка.

- [x] **`Z7.1` Параметр deliberate mismatch.** Позиция решена (`Q2`): −3 по
      умолчанию, −2 для слабого терминального класса; `inner$ref`/`inner$alt`
      несут `stab_at`, `weak_terminal`, `stab_policy`; `design_as_primers()`
      отвергает отрицательный `stab_ats`. `test-arms-mismatch-position.R` (43).
- [x] **`Z7.2` Совместимость якоря по химии.** `z7_priors.csv` несёт
      `anchor_chemistry` (Kwok — Taq-кинетика; Lefever — real-time AS-PCR, dCq;
      Liu — AS-PCR end-point; Ye — tetra-primer ARMS) и `transfer_assumption`
      («перенос на KASP не измерен»). `test-z7-priors.R` (49 pass).
- [x] **`Z7.3` Внешний якорь определён.** Liu 2012 (1686 праймеров; 3-е
      основание взвешенно 38.6 против 2-го 29.8 и 4-го 28.4, пересчитывается из
      таблицы) + Lefever + Kwok. Нет набора «SNP → сработало/нет» для
      KASP-переходов — поэтому `calibrated = FALSE`; закроется через `Z8.2`.
- [x] **`Z7.5` Добавленный мисматч не ослабляет другие гейты.** Мисматч —
      одним основанием внутри праймера; набор гейтов по id/стадии/уровню до и
      после совпадает; `ZW0 M2` держится и с мисматчем.
- [x] **`Z7.6` Физика, которой в модели нет.** Накопление продукта ослабляет
      дискриминацию. `marker_artefact()` пишет `limitations` (всегда массив,
      пустой без намеренного мисматча); поле не в `artefact_id()`.
      `test-report.R`, `test-plan.R`.

### Z8. Калибровка и benchmark [P1] — после `Z7`

Провалы собраны, вероятности — нет. `inst/extdata/conversion_rates.csv` несёт
долю конверсии **вместе с провалами** (`n_failed`, `n_fail_*`, `decomposes`).
Полностью разложены **3** строки: Makhoul S3, Yamagata Exp. 3 (24 → 18),
Grewal 2022 (363 = 698 − 335). Плюс Ayalew (TaqMan 39/50, KASP 49/50,
rhAmp 49/50; причины качественные), CAPS watermelon (532 → 271), soybean
EST-CAPS, SNP2CAPS, HRM TP53, AS-LAMP (с цензурированием), ddPCR JAK2
(LoB 0.002 % / LoD 0.005 % / LoQ 0.01 %). Выход пока **порядковый**.

- [ ] **`Z8.1` Отрицательные примеры — достаточно ли их?** Сейчас нет
      (3 разложенные строки); `yardstick`-метрики — только после этого.
    - правка: `conversion_rates.csv` + `data-raw/<paper>.R` с `outcome = failed`
      и разложенными провалами.
    - проверка: `test-z8-outcomes.R` — число `fully_decomposed` названо в отчёте.
    - блокер: **новые статьи** (открытого доступа мало найдено).
    - следствие: любая вероятность конверсии сегодня — одно наблюдение.
- [x] **`Z8.2` Harness/`targets`.** `dev/targets/_targets.R` (`panel_outcomes()`)
      пишет `inst/extdata/z8_gate_outcomes.csv`; 36 строк, `tar_make()`
      байт-в-байт. Все 36 — `worked` (0 `failed`), поэтому `z8_gate_rate()`
      (`dev/z8-harness.R`) даёт `rate = NA`, а не «1.0».
- [x] **`Z8.3` Не объединять доли конверсии** между культурами/исследованиями.
      `conversion_rates.csv` несёт `conversion_definition`
      (`design`/`manufacture`/`work`/`usable`) и `scope`; guard в
      `test-z8-outcomes.R` (114 pass).
- [x] **`Z8.4` Разделение по источнику, не по строке.** Ключ `(source, row)`;
      `z8_split_by_source()` в `dev/z8-harness.R`; ни один `source` не
      попадает на обе стороны.
- [x] **`Z8.5` Analytical performance: отдельные поля.** `R/marker_performance.R`
      (`ANALYTICAL_QUANTITIES`: LoB, LoD, LoQ, MAF, false-positive rate,
      partition count); `scheme_def_dpcr()` объявляет шесть; `test-report.R` —
      три определения не сливаются, `LoB < LoD <= LoQ`.
- [x] **`Z8.6` dPCR:** 4-инструментный набор разбиений (λ совпадает в пределах
      0.005) + панель из 7 ассays; Пуассон проверяется аналитически.

### Z9. Пробелы по методам [P1–P2]

- [x] **`Z9.1` LAMP: paired ALT set.** `paired` (FIP 3′ на REF; `status` =
      `complete`/`ref_not_placed`/`not_applicable`); `test-lamp.R`.
- [x] **`Z9.2` HRM: порог `dTm` измеряем.** `hrm_resolvability()` — 0.5 °C по
      среднему класса (LightCycler 96), `marginal` не повышается до `pass`;
      `family_hrm_resolvable` в каждой HRM-оценке. Остаток (melt-домен,
      heteroduplex) — `§V`.
- [x] **`Z9.3` TaqMan MGB/LNA.** `probe_modification = c("none", "mgb", "lna")`;
      `mgb` сдвигает окно по Kutyavin 2000 (`rule_of_thumb`,
      `taqman_kutyavin2000_mgb_tm.csv`), `lna` отклоняется. `test-taqman-mgb.R`.
- [x] **`Z9.4` Sanger:** объявленный предел чтения вместо самодельной модели
      трассы; ручной deconvolver удалён; ICE/TIDE/Poly Peak — внешние.
- [x] **`Z9.5` Sanger: гетерозиготные индели.** `$readout` несёт
      `SANGER_DETECTION` (`status = "not_modelled"`, `floor = 0.15`,
      `mixed = 0.20`, `calibrated = FALSE`) + внешние deconvolver'ы;
      `test-sanger.R`.
- [x] **`Z9.6` CRISPR: несколько спейсеров.** `crispr_guide_candidates()`
      ранжирует спейсеры; кинетика — `simulated = FALSE` до набора данных.
- [x] **`Z9.7` dPCR:** разбиение и LoB/LoD/LoQ — `Z8.6`/`Z8.5`; Tm-окно —
      долг `Z0.2`.

### ZW. Поперечные инварианты и провенанс [P0–P1]

#### ZW0. Метаморфические тесты

- [x] **Обмен REF/ALT** (и FAM/HEX): вердикт совпадает (`test-invariants.R`).
- [x] **Сдвиг координат** ничего не меняет (критерии по целому контигу
      исключены по имени).
- [x] **Порядок/перестановка входов** не меняет вердикт.
- [x] **`NA` ≠ `pass`** — ни в одной схеме.
- [x] **Реверс-комплемент локуса.** Модель локуса несёт `strand` (в
      `target_cols()`); `marker_design_at()` один раз реверс-комплементирует
      минус-страндовый локус до схемы. `ZW0 M1` в `test-invariants.R` проверяет,
      что плюс- и минус-запись одного локуса дают один вердикт.
- [x] **`ZW0.1` Где инвариант намеренно не держится.** `criterion_spec()` несёт
      `asymmetry_reason` (9 строк); `test-invariants.R` требует непустую
      причину; 7 названы в `docs/decision-audit.md`.

#### ZW1. Каждый гейт должен быть способен провалиться

- [x] **`ZW1.1`** `test-falsifiable.R` + `dev/check-criteria.R`: у каждого
      criterion-id есть валящая и проходящая фикстура; 137/137.
- [x] **`ZW1.2`** Новый criterion-id без фикстуры **блокирует PR**
      (`dev/check-criteria.R` в `--quick`).
- [x] **`ZW1.3`** Аудит критериев идёт по умолчанию.

#### ZW2. Штамп версий движков

- [x] `marker_engine_versions()`/`marker_engine_stamp()`/`..._field()` для семи
      движков; недоступный движок — `NA`, а не выпадение.
- [x] **`ZW2.2`** `engine_stamp` в каждой `Evaluation` (`marker_with_stamp()`) и
      в манифесте артефакта.
- [x] **`ZW2.3`** `marker_replay()` сверяет штамп до вердиктов
      (`artefact_engine_match()`), предупреждает с именами сдвинувшихся
      движков; манифест без штампа — `NA`.

#### ZW3. Независимый Tm-оракул как тест

- [x] **`ZW3.1`** `helper-tm-tolerance.R` + `test-tm-oracle.R` (GC-края,
      гомополимеры).
- [x] **`ZW3.2`** Primer3 `oligotm()` как третий оракул; `skip` без него.
- [x] **`ZW3.3`** Допуск — по **паре моделей**, не константа: SantaLucia↔Allawi
      max ΔTm 0.46 °C (допуск 0.5), Sugimoto +1.5…+7.3, Breslauer от 9.5.
- [x] **`ZW3.4` Измеренный Tm — смещение разложено, не скорректировано.**
      `tm_mm2020_s1_measured.csv` (4096 19-меров, 50 mM NaCl). `Schildkraut`
      читает −16.9 °C; NN-консистентный `SantaLucia1998-2` −9.6 °C при
      дефолтных 25 нМ каждой нити; при условии статьи (1 мкМ тотальной нити,
      `dnac=500`) −5.4 °C; остаток ≈−5.4 — флуоресцентная метка / первая
      производная. `calibration_stance` = «recorded, not corrected» сохранён;
      `primer_tm(dnac=)` добавлен. ZW3.5: `SantaLucia1998-2` воспроизводит
      печатную Primer3-панель сои до 0.08 °C — это референсная модель, но окна
      KASP/dPCR выписаны на координате `Schildkraut` и не двигаются без
      пере-якорения (отдельный пункт).
- [x] **`ZW3.5`** Реальные праймеры: 85 олиго соевой панели, max |ΔTm| 0.079 °C
      (r 0.9995); `salt_method = "SantaLucia1998-2"` воспроизводит напечатанное
      до 0.12 °C.
- [x] **`ZW3.6` Off-target соевой панели — по геному.** `genome_hits` считается
      из GCF_000004515.6 (обе цепи, 20 хромосом), печатные числа сохранены как
      `genome_hits_printed`; все 30 чисел воспроизводятся.
- [x] **`ZW3.7` 23 Sanger-ампликона — продукт из генома.** `in_silico_pcr()`
      подтверждает ровно один точный продукт; координаты совпадают во всех 23.
      Исправлен `exact_pattern_hits()` (фантомные хиты 1–2 п.н.).
- [x] **`ZW3.8` Качественные гейты панели — движком.** 13 KASP-наборов × 27 =
      351 строка (`soybean_worked_panel_kasp_gates.csv`); off-target пройден
      всеми 13; провалы (COM Tm 11/13, дискриминация 5, ампликон > 100 п.н. 4,
      GC 4, гомополимер 4, SSR 3, инвертированный повтор 3) — находки.
- [x] **`ZW3.9` Фиксированные праймеры и аннотация — общий интерфейс.**
      `assemble_primer_assay()` (`R/marker_assemble.R`), `marker_features()`
      (`R/marker_annotate.R`); 13 наборов пересобраны и нарисованы.
- [ ] **`ZW3.10` Пере-якорение Tm-координаты (flip на `SantaLucia1998-2`).**
      Дефолт `TM_SALT_METHOD` (`R/marker_core.R`) — `"Schildkraut2010"`; окна
      KASP (`KASP_LIMITS$as_tm` 55–65, `com_tm` 58–63) и dPCR `tm_range` (44–50)
      выписаны на этой координате. `ZW3.5`: NN-консистентный
      `SantaLucia1998-2` — пара к `DNA_NN_SantaLucia_2004` (Primer3/oligotm,
      Biopython), воспроизводит печатную сою до 0.08 °C и лежит в основе `ZW3.4`.
      правка: сменить дефолт **одновременно** с пере-якорением окон (вывести из
      SL2-шкалы, не «подогнать под биас»), пере-пин `tm_tolerances()`
      (`helper-tm-tolerance.R` — таблица снята на Schildkraut) и фикстур
      `test-taqman.R` (центральность), `test-selection.R`,
      `test-plot-methods.R`; обновить ожидания биаса `test-tm-mm2020.R` /
      `test-tm-oracle.R`.
      проверка: `dev/test.R --full`; `ZW0` — ни один вердикт не сдвигается молча
      (любой сдвиг фиксируется в `NEWS.md`).
      блокер: нет внешнего измеренного якоря на KASP/dPCR-условиях; пере-якорение
      переносит координату на SL2, а не калибрует (`fitted_here` остаётся 0).
      следствие: без flip 13 тестов держат старую координату; наивный flip без
      пере-якорения роняет `snp_fixture()` в doom на стадии pair (AS 55.3–63.8,
      COM ≤ 63, `COM ≥ AS`), то есть молча сдвигает дизайн.

#### ZW4. Уровень доказательности вместо `calibrated`-флага

- [x] **Лестница.** `CRITERION_EVIDENCE` — фактор
      `rule_of_thumb` < `vendor_guideline` < `primary_literature` <
      `fitted_here`; `criterion_spec()` несёт `evidence`.
- [x] **Миграция без потери.** Все 191 строка были `rule_of_thumb`;
      `calibrated` — deprecated-аргумент, не удалён.
- [x] **Вердикт несёт уровень.** `evaluation_evidence()` — самая слабая ступень
      среди прошедших проверок (`NA`, если ни одна).
- [x] **`ZW4.5`** `marker_decision_text(evidence = TRUE)` добавляет ступень к
      каждой проверке; по умолчанию вывод прежний.
- [x] **`ZW4.6` Источники названы построчно.** `criterion_provenance.csv`
      (212 строк) даёт `id` → `source` + `reason_class` + причина;
      `test-criteria.R` сверяет с реестром. **Итог:** `rule_of_thumb` **71**,
      `primary_literature` **128**, `vendor_guideline` **13**,
      `package assumption` **63**, `fitted_here` **0**.
    - **Tier-разделение (фазы A и B):** tier-строка = `target` (печатная точка
      отсчёта, её называет источник) + политика агрегации (`limit` =
      «maximise (per-tier bottleneck min)», `CRITERION_TIER_POLICY_EVIDENCE =
"rule_of_thumb"`, объявлена один раз в коде). Ступень несёт `target`,
      политика остаётся `rule_of_thumb` честно. Так поднято 36 из 54 tier-строк;
      18 остаются `rule_of_thumb` (Kwok-ранжирование, симулятор, метрика пакета
      без печатного порога).
    - **Фазы C–E:** `R/marker_bncc.R` — политика BNCC (regret [0,1], полосы
      A/B/C, взвешенный Чебышёв, фронт Парето) как опция
      `select_with_ladder(bncc = )`; лексикографический порядок по умолчанию,
      `BNCC_CALIBRATION` инертен (`fitted = FALSE`). Фаза E подняла две
      CAPS-строки (`liu2016_caps`, `vcf2caps2021`). Остальные 53 не-tier строки
      честно `rule_of_thumb`: источник ранжирует, а не печатает порог, либо
      якорь без читаемого вложения в библиотеке (операторское действие, ADR-0037).
    - **Фаза D (калибровка → `fitted_here`) заблокирована:** нет мокрых исходов;
      протокол §6 и хук `BNCC_CALIBRATION` готовы и инертны.

### Z10. Структура, плоидность, аналитика [P2]

- [x] **`Z10.1` Amplicon/target folding** целиком, не праймера
      (`R/marker_thermo.R`); движок через `capabilities()`.
- [x] **`Z10.2` Dosage/copy-number.** Данные:
      `polyploid_jaimes2024_sugarcane_dosage.csv` (42/76 размечены),
      `polyploid_cotton2023_kasp_screen.csv` (8/11 — отрицательный, в средние не
      входит). Ось в `marker_population_checks()`, не `pop_freq`;
      `test-polyploid-dosage.R` (51 pass). Движок `fitPoly` — Suggests, только
      с реализацией (`ZA4`).
- [x] **`Z10.3` GRanges/VCF-backed interchange** (`R/marker_io.R`).

### Z11. Документация [P2]

- [x] **`Z11.1`** «Which numbers are calibrated» в
      `vignettes/reading-results.Rmd` — уровни `rule_of_thumb` /
      `primary_literature` / `vendor_guideline` / `calibrated`.
- [x] **`Z11.2`** Cas12a/Cas13/dPCR/LAMP/KASP формулировки исправлены;
      `test-wording.R` по регуляркам запрещённого.

## ZA. Целевые интерфейсы и швы (codebase-design)

План перестройки **модулей и швов**, не список функций. Принцип: табличные
данные — `tibble`/`vctrs`; поведение — S3-generic; граф — `tidygraph`; рисование
— грамматика `ggplot2`; внешний движок — adapter с `capabilities()`. Не тащить
`R6`/`data.table`.

**Состояние:** шаги 1–4 частично. `Evaluation`, `MarkerContext`/
`ReactionConditions`/`Engines`, `criterion_spec()` есть; `Scheme` — адаптер
(`marker_scheme_adapter()`, `scheme_designer()`). Шаги 5–8 заперты на `Z8`.
Наука (`Z`, `Z7`) выше всего этого.

### ZA1. Целевой путь вызова

```r
context <- marker_context(
  genome, variants = NULL, population = NULL,
  conditions = reaction_conditions("pcr"), engines = marker_engines()
)
scheme <- marker_scheme("caps")
m  <- scheme$design(snp, context)
ev <- scheme$evaluate(m, context)
tr <- scheme$tracks(m, context)
dc <- scheme$decision()
ro <- scheme$readout(m, context)
```

Сегодня `marker_scheme("caps")` возвращает строку реестра, а не объект с
`$design()`; вызов идёт через `design_caps_marker(genome, snp)`.

### ZA2. Модули и швы

| Модуль                           | Tidy/ggplot-форма                                    | Заменяет                                    | Зависимость                       | Состояние                                             |
| -------------------------------- | ---------------------------------------------------- | ------------------------------------------- | --------------------------------- | ----------------------------------------------------- |
| `MarkerContext`                  | classed list, `rlang::list2`, `withr`                | ручной `...`-плумбинг                       | существующие                      | есть (`marker_run_context()`)                         |
| `Scheme`                         | S3 + `tibble`-реестр                                 | `design_<scheme>_marker` + `switch(scheme)` | нет                               | есть (`marker_scheme_adapter()`, `scheme_designer()`) |
| `Criterion`/`CriterionSet`       | `tibble` + list-column evaluator, `vctrs::vec_rbind` | дубли порогов в design/check/DAG            | `vctrs`                           | есть (`criterion_spec()`, 191 строка)                 |
| `Evaluation`                     | `tibble` + `forcats` status                          | три evaluator-пути                          | `vctrs`, `forcats`                | есть                                                  |
| `Population`                     | `tibble` + `GenomicRanges` adapters                  | мутация `marker_attach_haplotypes()`        | `fitPoly` позже (Suggests)        | частично                                              |
| `Engines`                        | S3-adapter + `capabilities()`                        | глобальный `PCR_COND`, молчаливые skip      | существующие + внешние            | есть (`marker_engines()`, `engine_has()`; `Z4.1`)     |
| `DecisionSchema`/`DecisionState` | `tidygraph` + Graphviz                               | ручная раскладка и статусы                  | `ggraph` (Suggests)               | нет (заперто `Z8`)                                    |
| `TrackSpec`/`Figure`             | `ggplot2` + `ggenes`/`patchwork`                     | стиль внутри схем                           | `scales`                          | частично                                              |
| `Readout`/`Calibration`          | model-result + `yardstick`                           | строки «not modeled»                        | `yardstick`, `targets` (Suggests) | нет (заперто `Z8`)                                    |
| YAML-технология                  | `yaml`                                               | ручная сериализация                         | `yaml`                            | нет (`Y1`)                                            |
| Отчёты/таблицы                   | `knitr::kable`/`gt`                                  | ручной текст                                | опционально Suggests              | частично                                              |

#### ZA2.1–ZA2.3 — сделано (коммит `33ebedf`)

- **`capabilities()` спрашивают.** Один шов `engine_has(capability, engines)`
  (`R/marker_core.R`); необъявленная/`"unknown"` — `FALSE`, никогда не проход;
  `test-capabilities.R` (10 pass). Остаток: 13 литералов
  `has_oligoarrayaux()` в шести файлах — переносить вместе с `Z0.1`/`Z5.1`.
- **`marker_generic_checks()`** перенесена из `R/marker_plot_methods.R` в
  `R/marker_checks.R`.
- **Семь стадий селекции** — одна запись на стадию: `marker_builtin_stages()`
  (`R/marker_stages.R`); `criterion_spec()` и `marker_decision_stages()` читают
  её; членство перечисляется, не выводится из семейства (`duplex` семейства
  `multiplex` и `asq` тоже проходят стадии `as`/`pair`).

### ZA3. Home-rolled → проверенный импорт

`spec_bind()` → `purrr::list_rbind()`/`vctrs::vec_rbind()`; `format_check_value()`
→ `scales::percent()`/`number()`; progress → `cli::cli_progress_*()`;
cache key → `rlang::hash()`; IUPAC → `Biostrings::matchPattern(fixed = FALSE)`,
`PDict`; интервалы → `IRanges`/`GenomicRanges`; граф → `tidygraph` + Graphviz
`dot`, `ggraph` для R-view; test doubles → `testthat::local_mocked_bindings()`;
таблицы в доках → `knitr::kable()`/`gt`; метрики → `yardstick`; resampling →
`rsample`; pipeline → `targets`.

### ZA4. Политика зависимостей

- **Imports при прямом использовании:** `vctrs`, `scales`, `yaml`.
- **Suggests:** `yardstick`, `rsample`, `targets`, `fitPoly`, `crisprBase`,
  `crisprDesign`, `crisprScore`, `crisprBowtie`, `crisprBwa`, `vdiffr`, `gt`.
- **Уже есть:** `Biostrings`, `IRanges`, `GenomicRanges`, `Rsamtools`,
  `BiocParallel`, `DECIPHER`, `TmCalculator`, `rBLAST`, `Rsubread`, `rprimer`,
  `dplyr`, `tidyr`, `purrr`, `forcats`, `tibble`, `ggplot2`, `ggenes`, `ggforce`,
  `ggrepel`, `patchwork`, `tidygraph`, `cli`, `withr`, `rlang`, `glue`,
  `processx`.
- **Не добавлять без второго adapter:** `R6`, `data.table`, `DiagrammeR`,
  `openPrimeR` (сначала maintenance). Тяжёлые `crispr*`/`fitPoly` — только при
  реальной реализации.
- [x] **`crisprBase` — источник данных, не обёртка.** Есть `SpCas9`, `SaCas9`,
      `SpG`, `AsCas12a`, `enAsCas12a`, `MAD7`, `CasRx` (Cas13d), `Csm`; Cas13 —
      только Cas13d, PFS не поле. Реестр `cas` им не заменяется.

### ZA5. Что не делать tidy/ggplot

`Scheme`/`Engines` — S3-generic, не data-frame; `DecisionSchema` — графовые
данные + state, рисование отдельным adapter; `Readout` — typed result
(`simulated`, `lod`, `maf`, `stderr`), не ggplot. Не вводить seam для одного
adapter (deletion test).

### ZA6. Порядок миграции

**До калибровки — шаги 1–4 (контракты данных):**

1. [x] `vctrs`-typed `Criterion`/`Evaluation`. `criterion_spec()` возвращает
       `marker_criteria` на `vctrs` (`vctrs` в Imports); вторая ось убрана —
       `marker_check_classes()` выводится из трёх `KASP_CHECK_*`. Колонка `strand`
       в модели локуса добавлена вместе с нормализацией в `marker_design_at()`
       (см. `ZW0`), поэтому она больше не объявляет ориентацию, которой дизайн
       не следует.
2. [x] `MarkerContext` + `ReactionConditions` + `Engines` (без тяжёлых deps).
3. [x] `Scheme` S3-adapter; схемы мигрируют по одной.
       `scheme_designer(scheme)` — один резолвер для встроенных и обнаруженных
       дизайнеров (им пользуются `marker_design_at()` и `artefact_designer()`);
       `marker_scheme_adapter()` даёт `design/evaluate/tracks/decision/readout`
       для всех 20 схем. `primer` — без локусного дизайнера намеренно (C13);
       `duplex`/`asq` резолвятся (`design_duplex_marker`/`design_asq_marker`).
4. [x] `Population` в context; `marker_attach_haplotypes()` — compatibility
       adapter; инварианты `ZW0` изолированы на шве.

**После калибровки (Z7/Z8) — шаги 5–8, не раньше:**

5. [ ] `DecisionSchema`/`DecisionState`; `tidygraph` + Graphviz. Сейчас DAG
       раскладывается вручную и работает; смена движка графа не меняет вердикт.
6. [ ] `TrackSpec`/`Figure`; `scales`; стиль только в renderer (косметика).
7. [ ] `Readout` + `yardstick`/`targets` — **только после `Z8`**: пока есть лишь
       положительные примеры, `yardstick`-метрики будут мерить смещение.
8. [ ] `yaml` export; затем optional `fitPoly`/`crispr*`.

Правило перехода: шаги 5–6 не начинаются, пока нет хотя бы одной внешней
фикстуры на `Z7` и отрицательных примеров в `Z8` (или явного решения отложить
`Z8`, зафиксированного здесь).

### ZA7. Acceptance / test surface

Вызовы в тестах идут тем же швом; contract-test на
`design/evaluate/tracks/decision/readout`; `NA` никогда не `pass`, `unknown` не
успех; у каждого criterion один bounds-object и `source`/`calibrated`; fake
thermo/specificity engines; golden-тесты — criterion-id, evaluation-table,
track-rows, не пиксели; после шага — `graphify update .` и `dev/check.R`. Граф —
ориентир, не доказательство: `stage_gates`/`stage_relax`/`stage_tiers` —
локальные векторы, граф о них молчит.

### ZA8. Docs/artifacts graph и freshness

Graphify индексирует код и doc-headings, но **не** файлы документов и **не**
артефакты; sidecar-merge провалился (`graph.json` — dict, `merge-graphs` ждёт
node-link list). Поэтому docs/artifacts — через manifest + generated Markdown.

- [x] **`docs/INDEX.md`** — единый вход, каждая ссылка проверяется.
- [x] **`docs/artifacts.yml`** — источник истины: `path`, `kind`, `generator`,
      `inputs`, `input_hash`, `output_hash`, `docs`, `code`, `status`.
- [x] **`dev/check-docs.R`** — нет orphan-документа, нет битой ссылки, каждый
      `include_graphics()` в manifest.
- [x] **`dev/check-artifacts.R`** — нет orphan-файла, все пути существуют,
      хеши совпадают, иначе **stale**.
- [x] **`docs/generated/artifacts.md`** — heading label = путь артефакта.
- [x] **Gate.** `dev/check.R` вызывает docs/artifacts checks; graphify —
      opt-in, не CI-зависимость.
- [x] **`ZA8.1` `targets` как держатель freshness.** `dev/targets/_targets.R` +
      `dev/targets.R` (`--make`/`--outdated`/`--check`): `engine` → `figures` →
      `manifest`. Хеши не переписаны на `tar_meta()` намеренно.
- [x] **`ZA8.2` Graphify health.** `dev/graphify-health.R --write` пишет
      `docs/generated/graphify-health.md`; `--check` падает только на регрессиях.
      Регенерация без `--no-cluster` обязательна. «Бог-узлы»: `tibble::tibble()`
      (171), `cli::cli_abort()` (94), `glue::glue()` (59), из нашего —
      `plot_primers()` (43). Это ограничение ориентации, не доказательство.
- [x] **`ZA8.3` Имена сообществ и `supersedes`.** `dev/relabel-bonsai.R` (424
      сообщества `named`, `echo-honest` → 0); `dev/handrolled.yml` +
      `dev/graphify-handrolled.R` (8011 узлов / 10688 рёбер / 2 `supersedes`).
- **Состояние:** `docs/artifacts.yml` хеширует весь `R/`, поэтому любая правка
  `R/` делает 20 SVG stale: `dev/render-gallery.R`, затем
  `dev/artifacts.R --write`. PNG git-ignored; SVG коммитятся.

## ZP. Параллельный rollout (coordination)

Правило: **ветки изолируют историю, не рабочее дерево.** Один clone = один
checkout. Все волны слиты, worktree один (`main`); новые ветки — от `main`.

### ZP0. Base snapshot [x]

Snapshot (`wip/roadmap-base`, `2785097`), Wave 1 (`e5a581a`), Wave 0 (`5e31a1a`)
слиты; `main` не трогали; ветка удалена после слияния.

### ZP1. Волны и ownership

- **Wave 0** (`5e31a1a`, последовательно): `Evaluation`/status,
  `MarkerContext`/`ReactionConditions`/`Engines`, KASP-plot нейтральный статус,
  единый `marker_evaluation_table()`. Один владелец shared-файлов.
- **Lean test loop:** `dev/test.R` — 1414 pass, 1m33s; `--scheme crispr` ~1m10s;
  `--full` — ≈8 мин (см. `dev/logs/latest-summary.md`); `check.R --quick` — 38 s,
  `--no-render` — 44 s. Тайминги пишутся в `dev/logs/`.
- **Wave 1** (`e5a581a`, disjoint paths): docs-artifacts, crispr-physics,
  sanger-physics, lamp-physics — по одной ветке на владельца.
- **Wave 1b** (слито): `feat/z2-criterion` (реестр, уровни без дрейфа),
  `feat/z3-conditions`, `feat/z5-specificity`, `feat/z9-sanger`.
- **Wave 2a** (слито `87433a8`, `0368d3b`): z3-callsites (`scheme_conditions()`,
  LAMP 65 °C, RPA 42 °C), gates-general (`MARKER_FAMILY_CHECKS`), degenerate-pcr
  (позиции мисматчей в 5′→3′ ориентации праймера). Integration: lean 1676 pass,
  затронутые схемы 2267 pass, 0 fail.
- **Wave 2b** (слито): `haplotype-track` — частота и гаплотип во всех схемах;
  нет частоты — `NA`, не вес 1.
- **Следующая волна — по одной ветке из `main`, один владелец файла.**

| Ветка                | Owns                                                                                | Пункты                                                                                               | Блокер                            |
| -------------------- | ----------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- | --------------------------------- |
| `feat/q1-withdraw`   | `R/marker_taqman.R`, `R/marker_readout.R`, `R/marker_polyploid.R`, их тесты         | `Q1`                                                                                                 | нет                               |
| `feat/kasp-mismatch` | `R/marker_checks.R`, `tests/testthat/test-kasp-checks.R`                            | `Z0.4`                                                                                               | нет                               |
| `feat/crispr-pfs`    | `R/marker_crispr.R`, `tests/testthat/test-crispr.R`                                 | `Z0.1`                                                                                               | нет                               |
| `feat/z8-negative`   | `inst/extdata/conversion_rates.csv`, `data-raw/`, `R/marker_evaluation*.R`          | `Z8.1`–`Z8.4`                                                                                        | новые статьи с `outcome = failed` |
| `feat/plots-style`   | `R/marker_plot_methods.R`, `R/marker_plots.R`, `tests/testthat/test-plot-methods.R` | `ZA8`/`ZW2.2`-доводка; оставшиеся call sites `dimer_thermo()`/`primer_tm()` (`:155`, `:350`, `:757`) | нет                               |
| `feat/za-scheme`     | `R/marker_schemes.R`, `R/marker_criteria.R`                                         | `ZA6` шаг 1 (единая ось id), шаг 3                                                                   | `ZA6` шаг 3 — после `Z0.4`, `ZA8` |
| `feat/zp5-tests`     | `tests/testthat/test-specificity.R`, `test-oligo-dg.R`, `test-dimer.R`              | `ZP5`                                                                                                | нет                               |
| vignette-ветки       | по одному файлу `vignettes/*.Rmd`                                                   | —                                                                                                    | последними; `Q1` не блокер        |

### ZP2. Merge order

Одна ветка из таблицы за раз, база — `main`; после каждого merge —
`Rscript dev/check.R --no-render`, полный gate — в integration worktree;
конфликты решать там, не в агентских. Vignette-ветки — последними. Удалять
ветку после merge и `git worktree list`. **Не squash:** golden-тесты и
`docs/artifacts.yml` привязаны к коммитам, сквош делает 20 SVG stale.

### ZP3. Agent brief

```text
Worktree: /abs/path/tidymarkers-<slug>
Branch: feat/<slug> (already checked out; do NOT checkout elsewhere)
Owns: <paths>
Forbidden: <paths>
Env (every shell call):
  PATH=/opt/R-4.6.1/bin:<clone>/.tools/bin:$PATH
  UNAFOLDDAT=<clone>/.tools/share/oligoarrayaux
  R_LIBS_USER=<clone>/.tmp-rt/lib
  library(tidymarkers) — stale 0.4.0; use pkgload::load_all(".")
Tests: Rscript dev/test.R <filter>; Rscript dev/check.R --no-render
Never run the full suite while iterating (8+ min): test_dir(filter=...) or dev/test.R.
Commit only on feat/<slug>; verify git branch --show-current before commit.
Commit message via `git commit -F file` — never -m with backticks or $().
```

### ZP4. Anti-patterns

N веток в одном clone без worktree; worktree от stale `HEAD` при грязном дереве;
два агента в одной виньетке/helper; параллельная dPCR-калибровка до Tm-движка;
merge без gate; удаление worktree до merge; новые worktree от
`wip/roadmap-base` (ветки нет); squash-merge.

### ZP5. Здоровье тестов

- `--full` был 10 → ≈17 мин после волны Z3/Z5. **Цель ≤ 10 мин достигнута:**
  двухфазный прогон (параллель + серийные SOCK-файлы) даёт **8m 10s**, 131 файл,
  36894 pass, 0 fail (`dev/logs/latest-summary.md`, 2026-11-06).
- Внешние сверки (BLAST+, Rsubread, rprimer, SOCK) не ослаблены; контроль — по
  `dev/logs/latest-summary.md` (топ медленных файлов), не по ощущениям.

## Y. Следующий слой: конструктор технологий, цели по области, специфичность

Весь раздел — `calibration_evidence = "rule_of_thumb"`: пороги не калиброваны
(`Z8`).

- [x] **Y0.** `tidykasp` → `tidymarkers`; `kasp_*` → `marker_*` (0.4.0).
- [x] **Y1.** Конструктор технологий: `marker_technology()` +
      `marker_filter_blocks()`; блок несёт `applies_to`, несовместимая смесь —
      ошибка, собранная технология — `rule_of_thumb`. Тест `test-technology.R`.
- [x] **Y2.** Цели по области: `marker_plan(region=)` + `region_targets()`;
      область расширяется только в варианты внутри блока (вне блока — ошибка).
      проверка: `test-plan.R` (область и haplotype→signal).
- [~] **Y3.** Специфичность как политика технологии: `marker_technology(specificity=)`
  и `homeolog_axes()`/`subgenome` есть; порог по гексаплоиду (A/B/D) не
  калиброван. блокер: `Z5.2`/`Z8`. следствие: порог чужой физики.
- [x] **Y4.** Ранжирование и Парето-фронт: `marker_pareto()` (критерий отбора) +
      `marker_informativeness()`; `linkage_risk` — `rule_of_thumb` (`Z8`).
      Тесты `test-pareto.R`, `test-informativeness.R`.
- [x] **Y5.** Инпуты: `read_variants_vcf()` индексно, `depth` (`Y5.3`); края
      закрыты — пустой/без данных → 0 строк, `.`/`*`/символьные отброшены,
      мультиаллельный расщеплён, `#CHROM` без 1000-строчного лимита, заголовок
      без CHROM/POS/REF/ALT — читаемая ошибка. проверка: `test-io.R`.
- [ ] **Y6.** Валидация и обратная задача. блокер: нет набора опубликованных
      маркеров для предсказательной проверки. следствие: валидация не начата.
- [x] **Y7.** Артефакт «почему этот маркер» — слой объяснения поверх `§X`.
      `marker_explain()`: проверки как упорядоченные причины (fail/unknown
      первыми) + `why.csv` в артефакте. проверка: `test-plan.R` (артефакт).

## YP. Распространение: контейнер, MCP, артефакт, обратная связь

- [x] **YP1.** Контейнер: `docker/Dockerfile` на pinned `rocker/r-ver:4.6.1`,
      `docker/entrypoint.sh`, `.dockerignore`; ставит только `Imports` через
      `pak` (бинарники PPM), `R CMD INSTALL`, `ENTRYPOINT` — stdio MCP.
      проверка: образ собран (`docker build -f docker/Dockerfile`, ~5.5 мин) и
      smoke в контейнере — `initialize`/`tools/call` вернули Tm 58.6 и KASP-микс.
- [x] **YP2.** MCP-слой: `inst/mcp/tools.R` — реестр 38 инструментов,
      покрывающий рабочий процесс виньеток: вход (`genome_info`, `read_variants`,
      `mask_positions`, `fasta`-пути), классификация (`variant_type`,
      `marker_site_class`), термодинамика/QC (`primer_tm`, `primer_thermo`,
      `oligo_dg`, `terminal_dg`, `dimer_thermo`, `hairpin_dg_thermo`,
      `primer_set_compatibility`, `fold_amplicon`, `pcr_efficiency`,
      `wobble_mismatch`), специфичность (`in_silico_pcr`, `find_caps_enzymes`,
      `screen_snp_targets`), проектирование (`design_marker`, `design_scheme`,
      `design_markers` — по `genome_path`/`variants_path`), решение
      (`marker_criteria`, `marker_decision_text`, `rescue_strategy`,
      `engine_status`, `engine_capabilities`, `design_chemistries`) и конфигурация схем.
      Схемные аргументы несут `enum` из `marker_schemes()`/`rescue_strategy()`
      (проектирование — 18: без `primer` и `duplex`). Не входят
      сознательно: графики и объекты `marker_plan()`/`marker_solve()` (не JSON),
      писатели (`marker_export`/`write_marker_results` — не пишут вне монтирования)
      и Ensembl-запросы (сеть). `inst/mcp/server.R` — stdio JSON-RPC 2.0.
      проверка: `test-mcp.R` вызывает каждый tool (448 pass) и валидирует
      каждый `inputSchema` (`required` ⊆ `properties`, у каждого своя `type`).
- [x] **YP3.** `marker_artefact()`/`marker_replay()` — манифест, провенанс,
      переигровка вердикта.
- [x] **YP4.** Риски агентского вызова: `verdict`/`relaxed`/`rejected` — поля
      структуры, не текст; README — что монтируется и какая сеть.
      проверка: `test-plan.R` пиннит типизированные колонки `outcome`/`verdict`/
      `relaxed`/`violated`/`risk`; `README.md` §Container and MCP server +
      `docker/README.md` — mount read-only `/data`, сеть не нужна на рантайме.
- [x] **YP5.** `marker_record_outcome()` дописывает исход в `outcomes.csv`,
      привязанный к `id`; `fails` без `notes` — ошибка.
