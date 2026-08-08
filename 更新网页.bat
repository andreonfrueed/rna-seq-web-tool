@echo off
chcp 65001 >nul
title RNA 分析网页 · 更新并重启

echo 正在停止旧服务、同步最新网页文件并重启，约 15 秒...
REM SEC-06 修复：
REM 1. 停止旧进程时优先用 pid 文件精确终止（旧写法 pkill -f 'streamlit[r]un'
REM    会误杀同机其他人的 streamlit 进程），pid 不存在才退回 pkill；
REM 2. 同步前先备份旧版网页到 ~/rna_web_app.bak（旧写法 rm -rf 后若复制失败，
REM    网页代码全丢，只能重装）；复制成功才删备份，失败自动回滚。
wsl --cd "%CD%" bash -lc "
set -e
if [ -f ~/rna_web_app/.web.pid ]; then
  PID=\$(cat ~/rna_web_app/.web.pid)
  kill \$PID 2>/dev/null || true
  rm -f ~/rna_web_app/.web.pid
  sleep 1
fi
pkill -f 'streamlit[r]un' 2>/dev/null || true
sleep 1
rm -rf ~/rna_web_app.bak
if [ -d ~/rna_web_app ]; then mv ~/rna_web_app ~/rna_web_app.bak; fi
if cp -r app ~/rna_web_app; then
  rm -rf ~/rna_web_app.bak
else
  echo '复制失败，回滚到旧版本...'
  rm -rf ~/rna_web_app
  mv ~/rna_web_app.bak ~/rna_web_app
  exit 1
fi
cd ~/rna_web_app
mkdir -p ~/rna_web_workspace/logs
(umask 077; nohup bash run.sh > ~/rna_web_workspace/logs/web.log 2>&1 &)
sleep 8
"
if errorlevel 1 goto FAIL

curl -f -s -o nul http://localhost:8501
if errorlevel 1 goto FAIL

start http://localhost:8501
echo 更新完成，网页已打开。
ping -n 4 127.0.0.1 >nul
exit /b 0

:FAIL
echo.
echo [更新失败] 请把窗口内容截图发给懂技术的人，或双击『一键安装.bat』重新安装。
pause
exit /b 1
