@echo off
title RNA 分析网页 · 一键环境测试

echo 正在复制迷你测试样例并运行环境自检（约 1 分钟）...
wsl --cd "%CD%" bash -lc "mkdir -p ~/rna_web_workspace/uploads && cp -f example/fastq/*.fq.gz ~/rna_web_workspace/uploads/ && cd app && ~/miniconda3/envs/pyseqrna/bin/python selfcheck.py"
if errorlevel 1 goto FAIL
echo.
echo ============ 环境测试通过 ============
echo 迷你样例已放进上传目录。
echo 打开网页 → 数据与参考文件 → 点"扫描已有文件"，即可看到 CON_1/CON_2/TREAT_1/TREAT_2。
echo 详细说明见 example\\README.md。
pause
exit /b 0

:FAIL
echo.
echo [测试未通过] 请把窗口内容截图发给懂技术的人。
echo 如果提示找不到 miniconda3，说明还没执行过一键安装，请先双击 一键安装.bat。
pause
exit /b 1
