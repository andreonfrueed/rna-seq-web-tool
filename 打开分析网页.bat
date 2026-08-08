@echo off
chcp 65001 >nul
title RNA 分析网页

curl -f -s -o nul http://localhost:8501 && goto OPEN

echo 正在启动分析网页，首次约 20-30 秒，请稍候...
REM SEC-08 修复：日志改放工作目录并收紧权限（600），不再写全局可读的 /tmp
wsl bash -lc "mkdir -p ~/rna_web_workspace/logs && cd ~/rna_web_app && (umask 077; nohup bash run.sh > ~/rna_web_workspace/logs/web.log 2>&1 &) ; sleep 5"

set N=0
:WAIT
ping -n 4 127.0.0.1 >nul
curl -f -s -o nul http://localhost:8501 && goto OPEN
set /a N+=1
if %N% GEQ 25 goto FAIL
goto WAIT

:FAIL
echo.
echo [启动失败] 等了约 2 分钟网页还没起来。
echo 排查方法：在 WSL 终端里运行  cat ~/rna_web_workspace/logs/web.log  查看报错。
pause
exit /b 1

:OPEN
start http://localhost:8501
echo 已打开 http://localhost:8501
echo 提示：关掉本窗口不会影响分析运行，分析在后台独立跑。
ping -n 5 127.0.0.1 >nul
