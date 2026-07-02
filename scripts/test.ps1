param([string]$BuildType = "Debug")

$vcvars = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat'
if (-not (Test-Path $vcvars)) {
    Write-Error "vcvars64.bat not found at '$vcvars'. Install the MSVC Build Tools, or edit this script if it's installed elsewhere."
    exit 1
}

cmd /c "`"$vcvars`" && cmake -S . -B build-tests -G Ninja -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_EXPORT_COMPILE_COMMANDS=ON && cmake --build build-tests --target vmc_tests && ctest --test-dir build-tests --output-on-failure"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
