ctest --test-dir build --output-on-failure
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
