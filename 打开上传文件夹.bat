@echo off
chcp 65001 >nul
title 打开上传文件夹
wsl bash -lc "mkdir -p ~/rna_web_workspace/uploads && explorer.exe $(wslpath -w ~/rna_web_workspace/uploads) 2>/dev/null; true"
echo 已打开上传文件夹（大文件直接拖进去，然后在网页里点"扫描已有文件"）。
ping -n 3 127.0.0.1 >nul
