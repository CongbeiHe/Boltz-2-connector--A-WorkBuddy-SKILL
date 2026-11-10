# 踩坑清单（认证 / 会话 / 沙箱）

## 1. macOS 沙箱下凭据写入失败（最高优先级）

**现象**：登录流程走到最后报 `refresh_token_store_failed`；`~/Library/Application Support/boltz-api/` 下只留下 `credentials.json.tmp-*` 孤儿文件；所有后续命令返回 `auth_required`，且 refresh token 丢失不可恢复。

**原因**：沙箱拦截 keychain 写入与 `Application Support/boltz-api/` 目录写入。

**解决**：
```bash
export BOLTZ_API_TEST_DISABLE_KEYRING=1   # 强制文件后端凭据
# 并在沙箱外执行（WorkBuddy 中即允许 Bash 工具的沙箱豁免）
```
孤儿 tmp 文件可手工改名为 `credentials.json` 救回一次（有效但 refreshable:false，过期前够用）。

**铁律：所有 boltz-api 调用都带上这个环境变量。** 本 skill 的 `boltz_pipeline.py` 已内置。

## 2. 会话约 1 小时过期

- `auth validate` 返回 `valid:true` 不代表 10 分钟后还 true。
- `refreshable:false` 的会话过期后无法刷新，只能重走 `auth login --device-code`。
- 设备码流程**有时自动批准，有时需要人工在浏览器点确认**（同一账号同一天内两种都出现过）。pipeline 脚本会打印 `verification_uri_complete` 与 user code，等人工确认，超时（默认 180 s）则提示手动执行。

## 3. 后台长任务与 token 轮换冲突

**现象**：一个 shell 里的后台 `run`（>2 min）跑着，另一个 shell 的凭据失效；或后台 `run` 中途 `auth_required` 失败退出。

**原因**：refresh token 是轮换式的——CLI 刷新后写回文件，多 shell 并发时互相踩；且刷新写入可能再次被沙箱拦。

**解决**：
- 长任务尽量**单 shell、脱沙箱**执行；
- 远端任务不受本地掉线影响——`run` 失败后：
  ```bash
  boltz-api predictions:structure-and-binding list      # 找最新 sab_pred_XXX
  boltz-api download-results --id sab_pred_XXX --name <name> --root-dir <dir>
  ```
  即可补下已完成的结果。`boltz_pipeline.py` 内置了这个恢复逻辑（从 run 目录的 `.boltz-run.json`/`run.json` 或远端 list 找 id）。

## 4. 约束与指标陷阱

- `force:false` 软约束实测几乎不起作用（5 条全超标 2–13 Å）。要几何必须 `force:true`。
- `force:true` 的势场参与扩散采样 → **ipTM 被"带约束"偏置**，引用时必须注明。
- `binding_metrics.binding_confidence` 与 ipTM 可以相差 2–3 个数量级（ipTM 0.918 vs confidence 0.001），当前不要用它做结论。
- 单 sample 的约束命中可能是 seed 侥幸：同 seed 复现一致，换 seed 距离可能翻倍。**下结论前先跑 5 sample 看中位数**。
- 5 样本距离分布极窄（σ≈1 Å）却全部违反约束 → 是锚点残基几何与模型偏好构象系统性不相容，force 也压不动（实测 4 Å 硬约束 vs 实测 9–19 Å 纹丝不动）。
- `--model boltz-2.2` 无法用于实际提交（见 api-schema.md）。
