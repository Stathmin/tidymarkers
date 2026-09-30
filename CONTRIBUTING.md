# Contributing

Спасибо! Это учебный пакет, но правила как в обычном R-пакете.

## Быстрый цикл разработки

```r
# загрузить пакет из исходников, без установки
pkgload::load_all()

# тесты (наследуют namespace, внутренние функции доступны без :::)
testthat::test_local()

# документация из roxygen -> NAMESPACE + man/
roxygen2::roxygenise()

# полная проверка
#   R CMD build . && R CMD check --no-manual kaspexample_*.tar.gz
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

`kasp-example` — обычный каталог пакета (`DESCRIPTION`, `R/`, `man/`, `vignettes/`,
`tests/`). Файл `.github/workflows/R-CMD-check.yaml` лежит внутри пакета как шаблон:
при выделении пакета в отдельный репозиторий перенесите его в корень репозитория.
