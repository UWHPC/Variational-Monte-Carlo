cmake -S . -B build -DBUILD_TESTING=OFF
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
cmake --build build --target vmc
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
