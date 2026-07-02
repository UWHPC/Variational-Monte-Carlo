param(
    [string]$BuildDir = "",
    [switch]$Cuda
)

$vcvars = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat'
if (-not (Test-Path $vcvars)) {
    Write-Error "vcvars64.bat not found at '$vcvars'. Install the MSVC Build Tools, or edit this script if it's installed elsewhere."
    exit 1
}

if ($BuildDir -eq "") {
    $BuildDir = if ($Cuda) { "build-prof-cuda" } else { "build-prof" }
}
$cudaFlag = if ($Cuda) { "ON" } else { "OFF" }

cmd /c "`"$vcvars`" && cmake -S . -B $BuildDir -G Ninja -DVMC_ENABLE_CUDA=$cudaFlag -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo && cmake --build $BuildDir --target vmc"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host "Profiler build ready: ./$BuildDir/vmc"
