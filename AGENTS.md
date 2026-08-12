# AGENTS.md

本文件为 AI 助手（Claude Code、Codex、Hermes 等）在本仓库工作时的配置说明。

## Agent skills

### Issue tracker

Issues and specs for this repo live as GitHub issues. Use the `gh` CLI for all operations. See `docs/agents/issue-tracker.md`.

### Triage labels

The skills speak in terms of five canonical triage roles, mapped to this repo's label strings: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `CONTEXT.md` + `docs/adr/` at the repo root, created lazily by `/domain-modeling` when terms or decisions get resolved. See `docs/agents/domain.md`.
