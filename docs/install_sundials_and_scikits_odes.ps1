# 使用 CMake 构建 Sundials 并安装 scikits.odes
# 请确保你已安装 Visual Studio（含 C++ Build Tools）和 CMake

$ErrorActionPreference = "Stop"

# 项目路径设置
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$SundialsVer = "7.3.0"
$SundialsZip = "sundials-$SundialsVer.tar.gz"
$SundialsUrl = "https://github.com/LLNL/sundials/releases/download/v$SundialsVer/$SundialsZip"
$SundialsSrcDir = "$ProjectRoot\sundials-src"
$SundialsBuildDir = "$ProjectRoot\sundials-build"
$SundialsInstallDir = "$ProjectRoot\sundials-install"
$VenvPath = "$ProjectRoot\.venv"

# 下载源码包
if (-not (Test-Path "$ProjectRoot\$SundialsZip")) {
    Write-Host "📥 正在下载 Sundials v$SundialsVer ..."
    Invoke-WebRequest -Uri $SundialsUrl -OutFile "$ProjectRoot\$SundialsZip"
}

# 解压
if (-not (Test-Path $SundialsSrcDir)) {
    Write-Host "📦 正在解压 Sundials..."
    tar -xzf "$ProjectRoot\$SundialsZip" -C $ProjectRoot
    Rename-Item -Path "$ProjectRoot\sundials-$SundialsVer" -NewName "sundials-src"
}

# 创建构建目录
if (-not (Test-Path $SundialsBuildDir)) {
    New-Item -ItemType Directory -Path $SundialsBuildDir | Out-Null
}

# 配置 CMake
Write-Host "⚙️ 使用 CMake 配置 Sundials..."
Push-Location $SundialsBuildDir
cmake $SundialsSrcDir `
    -DCMAKE_INSTALL_PREFIX=$SundialsInstallDir `
    -DCMAKE_BUILD_TYPE=Release `
    -DBUILD_SHARED_LIBS=ON `
    -DSUNDIALS_BUILD_STATIC_LIBS=OFF `
    -DEXAMPLES_ENABLE=OFF `
    -DSUNDIALS_INDEX_TYPE=int32 `
    -DSUNDIALS_PRECISION=double `
    -DSUNDIALS_BUILD_ARKODE=ON `
    -DSUNDIALS_BUILD_CVODE=ON `
    -DSUNDIALS_BUILD_IDA=ON `
    -DSUNDIALS_BUILD_KINSOL=ON `
    -DSUNDIALS_BUILD_TESTS=OFF
Pop-Location

# 编译并安装
Write-Host "🔨 开始构建 Sundials..."
cmake --build $SundialsBuildDir --config Release --target install

# 激活虚拟环境
Write-Host "🐍 激活虚拟环境"
& "$VenvPath\Scripts\Activate.ps1"

# 设置编译路径
$env:SUNDIALS_INCLUDE_DIR = "$SundialsInstallDir\include"
$env:SUNDIALS_LIB_DIR = "$SundialsInstallDir\lib"

# 安装 scikits.odes
Write-Host "📦 开始安装 scikits.odes ..."
pip install scikits.odes

Write-Host "`n✅ 完成！scikits.odes 已在虚拟环境中安装成功。"
