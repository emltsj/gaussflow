# GaussFlow

**面向计算化学的 Gaussian 16 自动化平台——确定性输入、科学验收、结构化结果。**

GaussFlow 把分子（SMILES、批量 `.gjf`，或可选的 AI 生成结构）变成完整的 Gaussian 16 工作流：3D 结构生成、受控写输入、排队计算、解析输出、科学验收、Li⁺ 结合能（含 BSSE）、氧化还原电位、前线轨道、热化学，以及 Excel/CSV 导出。

本仓库只包含**通用、可复用的预设**，不包含任何论文专用的生产计算链。

[English README](README.md)

---

## 为什么用 GaussFlow？

| 问题 | 做法 |
|---|---|
| 手写 `.gjf` 容易出错 | 确定性输入生成 + 路由白名单 + Gaussian 语法自动纠正 |
| Normal termination 不等于结果可用 | 检查 SCF、几何优化、虚频、原子/电荷/自旋、BSSE 解析 |
| 长任务/软件崩溃会丢结果 | 每个 job/study 实时写入 fsync 的 JSONL 日志；崩溃后可恢复 |
| 加 F 的位置不明确 | RDKit 逐 H 枚举 + 对称去重 + CH3/CH2/OH 位点标签 |
| 不同体系需要不同计算链 | 通用预设 + 可选的 AI 生成计算链（经路由白名单校验） |
| 批量对比麻烦 | 每个体系一个 batch，统一输出结合能/氧化还原/轨道/热化学表格 |

---

## 功能概览

- **三种分子输入**
  - SMILES（RDKit ETKDGv3 生成 3D + 力场清理）
  - 一个或多个 `.gjf` / `.com` 文件，自动解析路由、电荷、多重度
  - 可选的 AI 生成（DeepSeek；**需要自备 API Key**）
- **计算预设**
  - `opt_freq` — B3LYP/6-31G(d) 气相 opt+freq（默认，通用起点）
  - `standard` — opt+freq → 气相 SP → 溶剂相 SP，B3LYP/6-31G(d)，Water PCM
  - 可选的 AI 生成自定义计算链，校验后保存为 `ai_*` 方案
- **Li⁺ 结合能与 BSSE**
  - 在 O/F 位点放置 Li⁺，计算配合物、单体、Li⁺ 参考
  - 气相层级 `counterpoise=2` BSSE 校正
- **氧化还原**
  - cation/anion 各自 opt+freq + SP；相对电位 `E = −ΔG_sol(eV) + offset`
- **F 取代位点枚举**
  - 逐个 H 替换；canonical SMILES 对称去重
  - 标签如 `C1（CH3 端基）`、`C2（CH2）`、`O3（OH）`、`C1×2 + O3`
  - 支持 F=0–3 系列，每个结构一个 `.gjf`
- **科学验收**
  - Normal termination、Error/QPErr、SCF 收敛、优化收敛、虚频、原子/电荷/多重度、
    最终阶段能量、BSSE 解析
- **结果与导出**
  - Li⁺ 结合能（未校正/BSSE 校正）、氧化还原电位、HOMO/LUMO/Gap、热化学
  - Excel、CSV、`analysis.json`、Markdown 报告、XYZ 结构
- **崩溃恢复**
  - 每个 job 和完成的 study 实时写入 JSONL 日志
  - 从日志恢复数据库中缺失的已完成 study
  - 重新打开 GUI 后可继续被中断的排队任务
- **两种入口**
  - 桌面 GUI：`python GaussFlow.py` / `run_gui.cmd`
  - 命令行/Agent：`python gflow.py <command> --json`

---

## 环境要求

- Windows 10/11 64 位
- **Gaussian 16**（商业软件，需自行安装并持有许可）
- **GaussView 6** — 可选，仅用于人工查看结构
- **Shermo** — 可选，热化学后处理
- Python 3.10+（实测 3.14）用于源码运行
- Python 依赖：`numpy`、`scipy`、`rdkit`、`cclib`、`openpyxl`、`customtkinter`、`Pillow`、`psutil`

GaussFlow **不隶属于** Gaussian, Inc. 或 GaussView；请自行遵守其许可协议。

---

## 安装（源码）

```bash
git clone https://github.com/emltsj/gaussflow.git
cd gaussflow

python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt

python gflow.py doctor --json
python gflow.py detect --json
```

`detect` 会搜索环境变量、PATH、常见安装目录、注册表和运行中的进程。如果找不到 Gaussian，可在设置页手动指定。

启动 GUI：

```bash
python GaussFlow.py
# 或：run_gui.cmd
```

打包 Windows 单文件 exe（可选）：

```bash
build_exe.cmd
# 产物：dist\GaussFlow.exe
```

---

## 快速开始

### 1. 查看预设

```bash
python gflow.py recipes --json
```

### 2. 准备分子

```bash
python gflow.py prepare --smiles "CCO" --name ethanol --route "B3LYP/6-31G(d)" --workflow opt_freq --json
```

这一步只生成结构、路由预览和 study id，不提交计算。

### 3. 提交与轮询

```bash
python gflow.py submit --study st_xxx --json
python gflow.py status --study st_xxx --json
python gflow.py result --study st_xxx --json
```

### 4. 创建完整计算体系

```bash
python gflow.py recipe --smiles "O" --name water --recipe standard --li --redox --submit --json
python gflow.py analysis --batch bt_xxx --json
python gflow.py export-system --batch bt_xxx --json
```

### 5. 批量导入 gjf

```bash
python gflow.py gjf-batch --files "a.gjf;b.gjf" --submit --json
python gflow.py gjf-batch --dir "D:\inputs\gjf" --submit --json
```

### 6. F 取代系列

```bash
# 只列出对称不等价的单 F 位点（不写结构）
python gflow.py fluorinate --smiles "CCO" --n-f 1 --describe --json

# 生成 F=2 的全部对称不等价结构 + gjf
python gflow.py fluorinate --smiles "CCO" --n-f 2 --json
```

---

## AI 功能（可选）

AI 功能是**可选**的。仓库中**不包含任何 API Key**。如需使用，请自备 DeepSeek API Key：

- GUI：**设置 → DeepSeek 大模型** → 填入 Key → **测试连接** → **保存设置**
- 命令行/环境变量：`DEEPSEEK_API_KEY`
- 或写入本地 `config.json` 的 `llm.api_key`

Key 只保存在本机 `config.json`，该文件已被 **git 忽略**。请勿提交 API Key。

AI 可以：

- 用自然语言生成一个或多个候选分子
- 生成 F 取代系列（`--f-range 0-3`）
- 生成自定义计算链（经路由白名单校验，并自动纠正
  `def2-TZVP → def2TZVP`、`M06-2X → M062X`）

模型：`deepseek-v4-flash`（快/便宜，默认）、`deepseek-v4-pro`（最强）。
可在 GUI 或 `--model v4flash|v4pro` 切换。

```bash
python gflow.py ai --prompt "生成 3 个含氟碳酸酯候选" --n-molecules 3 --json
python gflow.py ai-workflow --prompt "B3LYP/6-311+G(d) opt+freq，然后 B2PLYP/def2TZVP 气相 SP，最后 M062X/6-31+G(d) PCM(DMSO) SP" --model v4pro --json
```

---

## 配置

首次运行会自动生成 `config.json`；`config.example.json` 展示所有默认值。

常用环境变量：

| 变量 | 作用 |
|---|---|
| `GFLOW_HOME` | 数据目录（默认项目目录；打包 exe 默认 `%LOCALAPPDATA%\GaussFlow`） |
| `GFLOW_CONFIG` | 指定 `config.json` 路径 |
| `G16_HOME` | Gaussian 16 安装目录（可选；`detect` 可自动找） |
| `DEEPSEEK_API_KEY` | 可选的 DeepSeek Key |

主要目录：

```
runs/<study_id>/              Gaussian 输入/输出/chk/result.json
generated_gjf/<gen_id>/       AI / F 枚举生成的 gjf + index.json
exports/system_<batch_id>/    分析 Excel/CSV/报告/XYZ
logs/results/results_*.jsonl  实时结果日志
logs/worker.log               后台 worker 日志
gaussflow.db                  SQLite 状态库
```

运行数据目录、数据库、本地配置和 API Key 都不会提交到仓库。

---

## 架构

```
GUI / CLI / Agent
        │  受控命令（--json）
        ▼
gaussflow
  gui.py          桌面 GUI（customtkinter）
  cli.py          受控命令接口
  config.py       配置与路径
  locate.py       自动检测 Gaussian / GaussView / Shermo
  recipes.py      通用计算链规划
  pipeline.py     study 创建、阶段串联、重试、汇总
  structure.py    RDKit 3D 结构、离子/配位位点
  gjf.py          .gjf 生成（含 counterpoise BSSE）
  runner.py       g16.exe 执行 + 墙钟看门狗
  daemon.py       后台 worker（claim → run → parse → validate → advance）
  parse.py        cclib + regex（能量/频率/热化学/HOMO-LUMO/BSSE）
  validate.py     科学验收规则
  diagnose.py     失败分类 + 白名单重试建议
  analysis.py     结合能 / 氧化还原 / 前线轨道 / 热化学
  export.py       CSV / Excel / XYZ / Markdown
  resultlog.py    实时 JSONL 结果日志 + 崩溃恢复
  fluorinate.py   F 位点确定性枚举
  llm.py          可选的 DeepSeek 集成
  store.py        SQLite（WAL）
```

---

## 科学口径

- 能量单位 Hartree，同时给出 kcal/mol 和 eV。
- `Eb = E(complex) − E(monomer) − E(Li⁺)`
- Gaussian counterpoise 约定：`corrected = raw + BSSE`（`BSSE > 0`），所以 `Eb(BSSE) = Eb(raw) + BSSE`。
- `G(level) = E(level SP) + [G(opt+freq) − E(opt+freq)]`。
- 氧化还原电位是**相对值**：`E = −ΔG_sol(eV) + reference_offset`；`config.redox.reference_offset_v` 需按文献参考电极设置。
- Gaussian 16 **不支持** PCM 与 `counterpoise` 同时使用，因此 BSSE 只在气相层级计算；溶剂相结合能为未校正值。
- 虚频：极小点 `0`，过渡态 `1`。
- Gaussian 写法：`def2TZVP`（不是 `def2-TZVP`）、`M062X`（不是 `M06-2X`）、`empiricaldispersion=gd3bj`。

### 验收规则

study 只有在关键检查全部通过（或仅有 warning）时才可用：

- 输出存在并以 `Normal termination of Gaussian 16` 结束
- 无 `Error termination`、`Convergence failure`、`QPErr`
- 原子数、电荷、自旋多重度与输入一致
- SCF 收敛
- 几何优化收敛（`Optimization completed` / `Stationary point found`）
- 虚频数在容差内
- 最终能量存在且属于最后阶段
- BSSE 任务解析到 `Counterpoise corrected energy` / `BSSE energy`

---

## 结果

`python gflow.py analysis --batch bt_xxx --json` 返回：

- `binding[]` — 每个 Li⁺ 位点、每个层级的 Eb、BSSE、Eb(BSSE)
- `redox[]` — ΔG_gas、ΔG_sol、相对电位
- `orbitals[]` — HOMO、LUMO、Gap
- `thermo[]` — ZPE、H(298)、G(298)、熵、可选 Shermo G

`export-system` 在 `exports/system_<batch_id>/` 下生成：

- `analysis.xlsx` — 总览 / 结合能 / 氧化还原 / 前线轨道 / 热化学
- `binding_energies.csv`、`redox_potentials.csv`、`frontier_orbitals.csv`
- `analysis.json`、`analysis_report.md`
- `batch/` — 通用 CSV/Excel/XYZ/报告产物

---

## 崩溃恢复

每个完成的 job 和 study 都会以 `flush + fsync` 追加到 `logs/results/results_YYYYMMDD.jsonl`。

- GUI：**结果 → 从日志恢复结果**
- 启动时：若日志中有数据库中缺失的已完成 study，会询问是否恢复
- 启动时：若上次退出留下运行中任务，会询问 **继续未完成任务**
- 命令行：

```bash
python gflow.py logs --json
python gflow.py recover --json
python gflow.py recover --batch bt_xxx --json
```

`recover` 只重建数据库中不存在的研究，重复执行是安全的。

---

## 命令速查

| 命令 | 作用 |
|---|---|
| `doctor` | 环境自检（可选 Gaussian 试算） |
| `detect` | 自动检测 Gaussian / GaussView / Shermo |
| `prepare` | 生成结构 + 路由预览（不计算） |
| `submit` / `submit-batch` | 异步提交 / 后台 worker |
| `status` | study/job 状态 |
| `result` | 结构化结果摘要 |
| `jobs` | 某个 study 的 job 明细 |
| `diagnose` | 失败分类 + 重试建议 |
| `export` | CSV / Excel / XYZ / JSON / Markdown |
| `recipes` | 列出计算方案 |
| `recipe` | 创建计算体系 |
| `analysis` | 结合能 / 氧化还原 / 前线轨道 / 热化学 |
| `export-system` | 导出体系分析 |
| `ai` | AI 生成分子并导出 gjf（需要 API Key） |
| `ai-workflow` | AI 生成并保存计算链（需要 API Key） |
| `fluorinate` | F 位点枚举 |
| `gjf-batch` | 批量导入 gjf |
| `logs` / `recover` | 查看 / 恢复实时日志 |
| `stop` | 终止 Gaussian 进程并取消排队 |
| `purge` | 清理 cancelled 记录 |
| `daemon` | 后台 worker 管理 |
| `preview` | 用 GaussView 打开结构/结果 |

所有命令都支持 `--json`，stdout 只输出一个 JSON 对象。

---

## 测试

```bash
python -m unittest discover -s tests -v
```

覆盖：路由校验、GJF 生成（含 counterpoise）、3D 结构/Li⁺ 位点、cclib 解析、
HOMO/LUMO、BSSE 算术、验收规则、失败诊断/重试、SQLite、Shermo、方案规划、
F 位点枚举、结果日志恢复、AI 返回校验。

---

## 安全与隐私

- **仓库中不包含任何 API Key。**
- `config.json`、`.env`、`logs/`、`runs/`、`exports/`、`generated_gjf/`、`gaussflow.db` 均被 git 忽略。
- API Key 只保存在本机 `config.json`，只发送到你配置的 API 地址。
- Gaussian 输入/输出、分子结构和结果都留在本机，除非你主动分享。
- 分享工作区或提交 issue 前，请检查 `config.json` 和日志。

---

## 贡献

欢迎 issue 和 pull request。提交 PR 前请先跑测试：

```bash
python -m unittest discover -s tests -v
```

新增计算链请保持通用、可配置；不要提交论文专用生产流程、本地路径、API Key 或生成结果。

---

## 许可证

MIT License，见 [LICENSE](LICENSE)。

Gaussian 16 与 GaussView 6 是 Gaussian, Inc. 的商业产品。
GaussFlow 是独立的自动化层，与 Gaussian, Inc. 无隶属关系。
