#!/bin/sh
# Start the stdio MCP server from the installed package.  Arguments after the
# image name are ignored; the server takes its input on stdin.
set -eu
server="$(Rscript -e 'cat(system.file("mcp/server.R", package = "tidymarkers"))')"
exec Rscript "$server"
