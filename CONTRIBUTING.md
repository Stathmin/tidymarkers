# Contributing

Спасибо! Это учебный пакет, но правила как в обычном R-пакете.

## Быстрый цикл разработки

```r
# окружение: Imports + dev/test-инструменты, без тяжёлых Suggests
# (быстрее, чем renv::restore(): без исторического 262-пакетного dev-лока)
# Rscript dev/install-deps.R                  # зависимости + сам пакет
# Rscript dev/install-deps.R --with-suggests  # полный набор, для R CMD check

# загрузить пакет из исходников, без установки
pkgload::load_all()

# быстрые тесты; --full добавляет внешние сверки (BLAST/Rsubread/аннотации/SOCK)
# Rscript dev/test.R --full
testthat::test_local()

# документация из roxygen -> NAMESPACE + man/
roxygen2::roxygenise()

# полная проверка
#   R CMD build . && R CMD check --no-manual tidymarkers_*.tar.gz
```

## Стиль

- **tidyverse-стайл**: `|>` (base pipe), `` ` ``, `snake_case`, `snake_case()` для
  глаголов; никаких `*_dfr()` — `purrr::list_rbind()`.
- Внутри `dplyr`/`ggplot2` используем местоимение `.data$col` (не голые имена).
- Комментарии — по делу: объясняют **почему**, не «что».
- Форматирование: `styler::style_pkg()`; линт: `lintr::lint_package()` (конфиг `.lintr`).

## Документация

- Каждая экспортируемая функция: `@param` на все аргументы, `@return`, `@family`
  (группа попадает в «See also» и в pkgdown), `@examples` (тяжёлые — `\dontrun{}`).
- Новое поведение → новый тест в `tests/testthat/`.
- Значимое изменение → запись в `NEWS.md` и в `ROADMAP.md`.

## Модель для правок

`tidymarkers` — обычный каталог R-пакета (`DESCRIPTION`, `R/`, `man/`, `vignettes/`,
`tests/`). CI в пакет не входит; при выделении в отдельный репозиторий настройте его
в корне репозитория.
