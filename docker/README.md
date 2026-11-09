# Контейнер: tidymarkers + MCP-сервер

Pinned-образ на `rocker/r-ver:4.6.1` — той же R, на которую нацелен lockfile.
Ставятся только `Imports` (через `pak` и бинарники Posit Package Manager);
тяжёлые `Suggests` (аннотация, `rBLAST`, `Rsubread`, `rprimer`, `crisprVerse`) в
образ не входят, и `capabilities()` честно отдаёт `unknown` для того, что им
нужно. Внешние движки (`OligoArrayAux`, `BLAST+`) в образ не ставятся:
`oligo_dg()`/`terminal_dg()` возвращают `NA` с предупреждением.

## Сборка

```sh
docker build -t tidymarkers:0.4.0 .
```

## Запуск

Контейнер stdio-only: сервер читает JSON-RPC 2.0 со stdin и пишет ответ в
stdout. Порт не открывается.

```sh
printf '%s\n' '{"jsonrpc":"2.0","id":1,"method":"tools/list"}' \
  | docker run --rm -i tidymarkers:0.4.0
```

## Что монтируется и какая сеть

- **Данные.** Геном/аннотация/кэш монтируются **только для чтения** в `/data`
  (`-v /path/to/data:/data:ro`). Ничего в образ не пишется. Дизайн-инструменты
  принимают путь монтированного FASTA/VCF: `genome_path` (`/data/genome.fa`),
  `variants_path` (`/data/variants.vcf`).
- **Сеть.** На рантайме не нужна: сервер не делает запросов. Сборка качает
  пакеты из PPM/CRAN/Bioconductor. Демо-данные (`ath_genome()`) качаются в
  per-user кэш при первом вызове — если сеть не нужна, монтируйте кэш или не
  вызывайте демо.
- **Инструменты.** Движки вроде `OligoArrayAux` сознательно не в образе; при
  необходимости смонтируйте `.tools` и передайте `oligoarrayaux_dir`
  (`/data/.tools`, префикс с `bin/` или сам `bin/`) — движок ищется только на
  время вызова.

## MCP-клиент

```json
{
    "mcpServers": {
        "tidymarkers": {
            "command": "docker",
            "args": [
                "run",
                "--rm",
                "-i",
                "-v",
                "/path/to/data:/data:ro",
                "tidymarkers:0.4.0"
            ]
        }
    }
}
```

Инструменты — `inst/mcp/tools.R` (реестр) и `inst/mcp/server.R` (цикл stdio).
Добавляя инструмент, добавьте его вызов в `tests/testthat/test-mcp.R`: тест
проверяет, что **каждый** зарегистрированный инструмент вызывается.
