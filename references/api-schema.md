# Boltz API — 实测输入 Schema 与 CLI 参考

> 来源：2026-11-10 对 boltz-api 0.43.0 的服务端校验报错逐字段反推（官方文档只描述 UI，无 JSON schema）。所有字段名均为服务端严格校验后的实测有效值。

## CLI 子命令

```
boltz-api auth login|validate|status|whoami
boltz-api predictions:structure-and-binding estimate-cost|start|run|list|retrieve
boltz-api download-results --id <pred_id> --name <name> --root-dir <dir>
```

- 输入传法：`--input @json://相对或绝对路径.json`
- 版本查询：`boltz-api cli version`；`--debug` 可看校验细节。

## 顶层字段（严格校验，未知字段直接拒绝）

| 字段 | 类型 | 说明 |
|---|---|---|
| `entities` | array | 必填，至少一个 |
| `binding` | object | 蛋白-蛋白/配体-蛋白结合建模 |
| `num_samples` | int | 采样数；每 sample 单独计费 |
| `constraints` | array | 距离接触约束 |
| `bonds` | array | 共价键（环化、交联等） |
| `model_options` | object | seed / recycling_steps / sampling_steps / step_scale |

## entities 元素

```json
{"type": "protein", "chain_ids": ["A"], "value": "MKTAYIA..."}
```

- `type`：实测有效值 `protein`（`proteinChain`/`PROTEIN` 报错）
- `chain_ids`：字符串数组；`value`：序列。允许 `modifications`。
- 实体声明顺序 = token 顺序（PAE 矩阵按此排列，每残基 1 token）。

## binding

```json
{"type": "protein_protein_binding", "binder_chain_ids": ["B"]}
{"type": "ligand_protein_binding", "binder_chain_id": "B"}
```

注意是 `binder_chain_ids`（蛋白，复数数组）/ `binder_chain_id`（配体，单数），二者不可混用。没有 target 字段——目标是其余所有链。

## constraints（contact 类型）

```json
{
  "type": "contact",
  "token1": {"type": "polymer_contact", "chain_id": "B", "residue_index": 1},
  "token2": {"type": "polymer_contact", "chain_id": "A", "residue_index": 25},
  "max_distance_angstrom": 4.0,
  "force": true
}
```

- **`residue_index` 是 0 基**；CIF 里残基号是 1 基，写输入时必须减 1。
- `force:false`（软约束）实测几乎无效——5/5 条全部大幅超标；要几何就用 `force:true`（势场在扩散采样中起作用，会偏置 ipTM）。
- token 的 `type` 实测有效值：`polymer_contact`（链间/链内接触约束用这个）。

## bonds（共价键，用于环肽等）

```json
{
  "atom1": {"type": "polymer_atom", "chain_id": "B", "residue_index": 6, "atom_name": "C"},
  "atom2": {"type": "polymer_atom", "chain_id": "B", "residue_index": 0, "atom_name": "N"}
}
```

- `residue_index` 同样 0 基；`atom_name` 是 PDB 命名（C/N/CA/...）。
- 环肽首尾酰胺环化：末位残基的 `C` ↔ 首位残基的 `N`。实测出结构后键长 1.24–1.58 Å（<1.8 Å 判为成键）。

## model_options

```json
{"seed": 42, "recycling_steps": 3, "sampling_steps": 200}
```

同 seed 可逐原子复现（已验证）。

## 成本（2026-11 实测）

| 规模 | 单 sample |
|---|---|
| 合成小双链（33+26 aa） | $0.025 |
| 461+187 aa 二元复合物 | $0.04 |
| 461+7+187 aa 三元 + 3 约束 | $0.04 |

`num_samples:5` 相应 ×5。运行时长：小体系 <1 min；655 token / 5 sample 约 26 min。

## 输出文件

`<root-dir>/<name>/outputs/files/prediction/`：

- `<pred_id>_predicted.cif`：best sample 结构
- `sample_N_predicted_structure.cif`：各 sample 结构
- `sample_N_pae.npz`：各 sample PAE 矩阵（npy 数组，zip 容器）
- `metrics.json`：`all_sample_results[N].metrics{ptm,iptm,complex_plddt,...}`、`best_sample.metrics`、`binding_metrics.binding_confidence`
- CIF 的 B-factor 列 = pLDDT（0–1 标度）。

## 已知服务端限制

- `--model boltz-2.2`：`estimate-cost` 接受，`start`/`run` 报 `Unsupported parameter: body.model`。实际可跑的只有默认 `boltz-2.1`。
- MSA 服务端自动生成，输入里无需（也无法）提供。
- `retrieve --id X --transform status -r` 可轮询状态。
