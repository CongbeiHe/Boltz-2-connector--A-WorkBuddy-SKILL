#!/usr/bin/env bash
# boltz2-structure-prediction — 一键安装
# 用法:
#   ./install.sh            # 装 Boltz CLI + 登录 + 装进 ~/.workbuddy/skills/
#   ./install.sh --skill    # 只装 skill（跳过 CLI 安装）
#   ./install.sh --cli      # 只装 CLI + 登录
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SKILL_DIR="$HOME/.workbuddy/skills/boltz2-structure-prediction"

install_cli() {
  echo "== 安装 boltz-api CLI =="
  if command -v boltz-api >/dev/null 2>&1; then
    echo "已安装: $(command -v boltz-api)"
  else
    curl -fsSL https://install.boltz.bio/boltz-api/install.sh | sh
  fi
  echo
  echo "== 登录（浏览器设备码流程）=="
  echo "提示: 沙箱/受限环境下务必加 BOLTZ_API_TEST_DISABLE_KEYRING=1"
  BOLTZ_API_TEST_DISABLE_KEYRING=1 boltz-api auth login
}

install_skill() {
  echo "== 安装 skill -> $SKILL_DIR =="
  mkdir -p "$HOME/.workbuddy/skills"
  rm -rf "$SKILL_DIR"
  cp -R "$HERE" "$SKILL_DIR"
  echo "已安装。新开一个 WorkBuddy 会话即可使用（说\"用 Boltz 跑结构预测\"即可触发）。"
}

case "${1:-all}" in
  --cli)  install_cli ;;
  --skill) install_skill ;;
  all)    install_cli; install_skill ;;
  *) echo "用法: ./install.sh [--cli|--skill]"; exit 1 ;;
esac
echo "完成。"
