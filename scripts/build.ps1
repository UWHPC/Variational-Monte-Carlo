param(
    [string]$BuildType = "Release",
    [switch]$Cuda
)

$vcvars = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat'
if (-not (Test-Path $vcvars)) {
    Write-Error "vcvars64.bat not found at '$vcvars'. Install the MSVC Build Tools, or edit this script if it's installed elsewhere."
    exit 1
}

$buildDir = if ($Cuda) { "build-cuda" } else { "build" }
$cudaFlag = if ($Cuda) { "ON" } else { "OFF" }

cmd /c "`"$vcvars`" && cmake -S . -B $buildDir -G Ninja -DVMC_ENABLE_CUDA=$cudaFlag -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=$BuildType && cmake --build $buildDir --target vmc"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
