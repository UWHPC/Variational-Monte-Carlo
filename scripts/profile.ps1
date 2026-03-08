param(
    [string]$BuildType = "RelWithDebInfo",
    [string]$BuildDir = "build-prof"
)

cmake -S . -B $BuildDir -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=$BuildType -DPROFILE_MODE=ON
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
cmake --build $BuildDir --target vmc
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& "./$BuildDir/vmc"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
