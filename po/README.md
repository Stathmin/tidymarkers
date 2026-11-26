# Локализация (i18n) — сообщения

Статус: **все пользовательские сообщения локализуются** через GNU gettext; в каталоге
`ru` 262 сообщения (199 переведено, 8 fuzzy, 55 новых без перевода — см. `po/R-ru.po`).
Справка (`.Rd`) намеренно не переводится (см. `ROADMAP.md` §Q).

## Как устроено

- Строки-шаблоны обёрнуты в `gettext()`; домен `R-tidymarkers` R подставляет по namespace.
  `rlang::abort/warn/inform` и `cli::cli_*` potools сам **не** распознаёт — их строки
  обёрнуты вручную.
- `abort_msg()` (`R/marker_i18n.R`) — единственный путь для пользовательских abort'ов:
  переводит сообщение и каждый bullet **один раз** через `gettext()` и передаёт в
  `cli::cli_abort()`. Для извлекателя он зарегистрирован как custom translation
  function: `potools::po_extract(custom_translation_functions = list(
  R = "abort_msg:...\\class,call,parent"))` — все аргументы, кроме `class`, `call`,
  `parent`, считаются переводимыми строками. Без `R =` в списке potools аргумент
  игнорирует и msgid'ы `abort_msg()` в шаблон не попадают.
- `po/R-tidymarkers.pot` — шаблон (`potools::po_extract()`).
- `po/R-ru.po` — перевод (`potools::po_update()`).
- `inst/po/ru/LC_MESSAGES/R-tidymarkers.mo` — скомпилированный каталог
  (`potools::po_compile()`); именно он ставится и виден R.

## Цикл (нужны system gettext-tools + potools)

```r
# код -> po/R-tidymarkers.pot; `R =` обязателен, иначе abort_msg() не виден
potools::po_extract(custom_translation_functions = list(
  R = "abort_msg:...\\class,call,parent"
))
potools::po_update()    # pot -> po/R-ru.po   (merge, сохраняя перевод)
potools::po_compile()   # po  -> inst/po/ru/LC_MESSAGES/R-tidymarkers.mo
```

Затем `R CMD INSTALL .` и:

```r
Sys.setenv(LANGUAGE = "ru"); library(tidymarkers)
```

## Грабли (проверено на этом пакете)

1. **`gettext()` оборачивает ШАБЛОН, а не результат.** С `glue` — только
   `glue::glue(gettext("..."))`, **не** `gettext(glue::glue("..."))`: иначе `glue`
   подставит значения до перевода, и lookup промахнётся (тихо останется английский).
2. **Одна строка = одно сообщение.** `gettext("a", "b")` рантайм склеивает, а
   извлекатель даёт два `msgid` → промах. Весь шаблон — одним литералом.
3. **Никаких `%` в сообщениях.** `tools::checkPoFile()` читает `%` как форматную
   директиву и ругается на несовпадение; переформулировать (мы убрали «100%»).
4. **`R CMD INSTALL` не компилирует `po/`.** `.mo` кладёт `po_compile()` в `inst/po/`,
   откуда он и ставится; без этого `LANGUAGE=ru` молчит.
5. **Ярлыки `print()`** (`<marker_assay>`, `AS-ref`, `sequence(s):`, …) намеренно
   оставлены на английском — это фиксированные диагностические ярлыки/терминология.
6. Числа — `ngettext()`; `Plural-Forms` для `ru` (3 формы) прописан в `R-ru.po`.
