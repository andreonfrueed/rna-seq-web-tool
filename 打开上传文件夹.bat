@echo off
chcp 65001 >nul
title 打开上传文件夹

REM 注意：下面写死的 ~/rna_web_workspace 必须与 web_config.yaml 的 workspace_dir 保持一致
REM （仓库里是 app/web_config.yaml，部署后是 WSL 侧 ~/rna_web_app/web_config.yaml），
REM 且注意环境变量 RNA_WEB_WORKSPACE 会覆盖该配置。将来改工作区目录要一起改，否则又会跳回系统「文档」目录。
for /f "usebackq delims=" %%P in (`wsl bash -lc "mkdir -p ~/rna_web_workspace/uploads && wslpath -w ~/rna_web_workspace/uploads"`) do set "UPLOAD_WIN=%%P"

if not defined UPLOAD_WIN (
    echo 打开失败：无法访问 WSL 工作区，请先运行 一键环境测试.bat
    pause
    exit /b 1
)

explorer "%UPLOAD_WIN%"
echo 已打开上传文件夹（大文件直接拖进去，然后在网页里点"扫描已有文件"）。
ping -n 3 127.0.0.1 >nul
