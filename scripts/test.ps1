param([string]$BuildType = "Debug")

cmake -S . -B build-tests -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=$BuildType `
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
cmake --build build-tests --target vmc_tests
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
ctest --test-dir build-tests --output-on-failure
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }