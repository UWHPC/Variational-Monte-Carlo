param(
    [string]$BuildDir = "build-prof"
)

cmake -S . -B $BuildDir -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
cmake --build $BuildDir --target vmc
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Profiler build ready: ./$BuildDir/vmc"
