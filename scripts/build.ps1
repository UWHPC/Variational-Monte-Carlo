cmake --build build
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
cmake -S . -B build
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
