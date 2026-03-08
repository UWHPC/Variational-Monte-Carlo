param([string]$BuildType = "Release")

cmake -S . -B build -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=$BuildType
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
cmake --build build --target vmc
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
