@echo off
chcp 65001 >nul
title RNA 分析网页 · 一键安装

echo 检查 WSL（Windows 的 Linux 子系统）...
wsl --cd "%CD%" bash -lc "echo ok" >nul 2>&1
if errorlevel 1 goto NO_WSL
echo WSL 已就绪。

REM 按物理内存自动配置 WSL 内存（分析需要大内存），上限 30GB
for /f "delims=" %%i in ('powershell -NoProfile -Command "$m=[math]::Floor((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB*0.75); if($m -gt 30){$m=30}; Write-Output $m"') do set RAMGB=%%i
REM 只在没有内存设置时追加；先备份原 .wslconfig 到 %USERPROFILE%\.wslconfig.bak，
REM 并用 UTF8 写回（SEC-07 修复：旧的 ascii 编码会把文件里的中文/其他非 ASCII
REM 字符变成 '?'，直接破坏已有的 WSL 配置）。
REM BUG-08 修复：swap 只在没有已有 swap 设置时才追加，不再无条件追加导致重复键。
for /f "usebackq delims=" %%w in (`powershell -NoProfile -Command "$p=Join-Path $env:USERPROFILE '.wslconfig'; $l=@(); if(Test-Path $p){$l=Get-Content $p}; if($l -match '^\s*memory\s*='){ Write-Output 'WSLCFG_KEEP' } else { if(Test-Path $p){ Copy-Item $p ($p + '.bak') -Force }; if(-not ($l -match '^\[wsl2\]')){ $l+='[wsl2]' }; $l+=('memory=%RAMGB%GB'); if(-not ($l -match '^\s*swap\s*=')){ $l+=('swap=16GB') }; Set-Content -Encoding UTF8 -Path $p -Value $l; Write-Output 'WSLCFG_CHANGED' }"`) do set WSLCFG=%%w
if "%WSLCFG%"=="WSLCFG_CHANGED" goto WSLCFG_CHANGED
echo 沿用现有 WSL 内存设置
goto INSTALL_BEGIN

:WSLCFG_CHANGED
echo 已更新 WSL 内存设置（原配置已备份为 %USERPROFILE%\.wslconfig.bak），重启 WSL 使内存生效...
wsl --shutdown
ping -n 3 127.0.0.1 >nul

:INSTALL_BEGIN
echo.
echo ============ 开始安装（全程自动，约 30-60 分钟，请勿关窗口） ============
wsl --cd "%CD%" bash -lc "bash setup_env.sh"
if errorlevel 1 goto INSTALL_FAIL
echo.
echo ==================== 安装结束 ====================
echo.
echo 接下来：
echo   1. 双击『打开分析网页.bat』启动网页
echo   2. 浏览器自动打开 http://localhost:8501
echo   3. 详细用法看『使用说明.md』
pause
exit /b 0

:INSTALL_FAIL
echo.
echo [安装失败] 安装失败，请把窗口内容截图发给懂技术的人。
echo.
pause
exit /b 1

:NO_WSL
echo.
echo [未检测到 WSL] 这台电脑还没有 Linux 子系统，需要先装一次：
echo.
echo   1. 右键点"开始"菜单 → 选"终端(管理员)" 或 "PowerShell(管理员)"
echo   2. 输入下面这行，回车：
echo.
echo      wsl --install -d Ubuntu
echo.
echo   3. 按提示重启电脑
echo   4. 重启后再双击本文件继续安装
echo.
pause
exit /b 1
