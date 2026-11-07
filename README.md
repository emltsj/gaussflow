# GaussFlow

**Gaussian 16 automation for computational chemistry — deterministic inputs, scientific validation, structured results.**

GaussFlow turns a molecule (SMILES, a batch of `.gjf` files, or an optional AI-generated structure) into a complete Gaussian 16 workflow: 3D structure generation, controlled input writing, queued execution, output parsing, scientific validation, Li⁺ binding energies with BSSE, redox potentials, frontier orbitals, thermochemistry, and Excel/CSV export.

The repository ships **generic, reusable presets only**. It does not include any paper-specific production chain.

[中文说明](README.zh-CN.md)

---

## Why GaussFlow?

| Problem | GaussFlow approach |
|---|---|
| Hand-writing `.gjf` files is error-prone | Deterministic input builder with a route whitelist and Gaussian-syntax normalisation |
| "Normal termination" does not mean usable results | Validates SCF and geometry convergence, imaginary frequencies, atom/charge/spin consistency, BSSE parsing |
| Long jobs and GUI crashes lose results | Every finished job/study is appended to a fsync'd JSONL result log; recover results after a crash |
| Fluorination sites are ambiguous | Deterministic RDKit H→F enumeration, symmetry de-duplicated, with CH3/CH2/OH site labels |
| Papers use different calculation chains | Generic presets plus optional AI-generated chains validated against the route whitelist |
| Batch comparison is tedious | One batch per system, structured binding/redox/orbital/thermo tables, Excel/CSV export |

---

## Feature overview

- **Three input modes**
  - SMILES (RDKit ETKDGv3 3D embedding + force-field cleanup)
  - One or many `.gjf` / `.com` files, parsed for route, charge, multiplicity
  - Optional AI-assisted molecule generation (DeepSeek; **bring your own API key**)
- **Calculation presets**
  - `opt_freq` — B3LYP/6-31G(d) opt+freq (default, generic starting point)
  - `standard` — opt+freq → gas-phase SP → solvent-phase SP at B3LYP/6-31G(d), Water PCM
  - AI-generated custom chains, validated and saved as `ai_*` recipes
- **Li⁺ binding and BSSE**
  - Places Li⁺ at O/F sites; computes complex, monomer, and Li⁺ reference
  - Gaussian `counterpoise=2` BSSE at gas-phase levels
- **Redox**
  - Cation/anion opt+freq + SP; relative potentials `E = −ΔG_sol(eV) + offset`
- **F substitution site enumeration**
  - Replaces every H systematically; canonical-SMILES symmetry de-duplication
  - Labels such as `C1 (CH3 terminal)`, `C2 (CH2)`, `O3 (OH)`, `C1×2 + O3`
  - Multi-F series (0–3) and one `.gjf` per structure
- **Scientific validation**
  - Normal termination, no Error/QPErr, SCF convergence, optimisation convergence,
    imaginary frequencies, atom/charge/multiplicity consistency, final-stage energy, BSSE parsing
- **Results and export**
  - Li⁺ binding energies (raw and BSSE-corrected), redox potentials, HOMO/LUMO/gap, thermochemistry
  - Excel workbook, CSV files, `analysis.json`, Markdown report, XYZ structures
- **Crash safety**
  - Real-time JSONL log of each job and completed study
  - Restore missing completed studies from the log
  - Resume interrupted queue items after reopening the GUI
- **Two entry points**
  - Desktop GUI: `python GaussFlow.py` (`run_gui.cmd`)
  - CLI/Agent: `python gflow.py <command> --json`

---

## Requirements

- Windows 10/11 64-bit
- **Gaussian 16** (commercial; you must own a licence)
- **GaussView 6** — optional, used only for manual structure inspection
- **Shermo** — optional thermochemistry post-processing
- Python 3.10+ (3.14 tested) for source runs
- Python packages: `numpy`, `scipy`, `rdkit`, `cclib`, `openpyxl`, `customtkinter`, `Pillow`, `psutil`

GaussFlow is **not affiliated with Gaussian, Inc. or GaussView**. You are responsible for complying with their licences.

---

## Installation (from source)

```bash
git clone https://github.com/emltsj/gaussflow.git
cd gaussflow

python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt

python gflow.py doctor --json
python gflow.py detect --json
```

`detect` searches environment variables, PATH, common install locations, the registry, and running processes. If it cannot find Gaussian, set the path explicitly:

```bash
python gflow.py detect --json
# or edit config.json / use the Settings tab
```

Launch the GUI:

```bash
python GaussFlow.py
# or: run_gui.cmd
```

Build a single-file Windows executable (optional):

```bash
build_exe.cmd
# output: dist\GaussFlow.exe
```

---

## Quick start

### 1. Inspect the available recipes

```bash
python gflow.py recipes --json
```

### 2. Prepare a molecule

```bash
python gflow.py prepare --smiles "CCO" --name ethanol --route "B3LYP/6-31G(d)" --workflow opt_freq --json
```

This writes the structure, route preview, warnings, and a study id; it does **not** submit.

### 3. Submit and poll

```bash
python gflow.py submit --study st_xxx --json
python gflow.py status --study st_xxx --json
python gflow.py result --study st_xxx --json
```

### 4. Create a whole calculation system

```bash
python gflow.py recipe --smiles "O" --name water --recipe standard --li --redox --submit --json
python gflow.py analysis --batch bt_xxx --json
python gflow.py export-system --batch bt_xxx --json
```

### 5. Batch-import existing `.gjf` files

```bash
python gflow.py gjf-batch --files "a.gjf;b.gjf" --submit --json
python gflow.py gjf-batch --dir "D:\inputs\gjf" --submit --json
```

### 6. F substitution series

```bash
# list symmetry-unique single-F sites (no structures written)
python gflow.py fluorinate --smiles "CCO" --n-f 1 --describe --json

# write all symmetry-unique F=2 structures and .gjf files
python gflow.py fluorinate --smiles "CCO" --n-f 2 --json
```

---

## AI-assisted features (optional)

AI features are **optional**. The repository contains **no API key**. To use them, provide your own DeepSeek API key:

- GUI: **Settings → DeepSeek** → paste key → **Test connection** → **Save**
- CLI/environment: `DEEPSEEK_API_KEY`
- Or set `llm.api_key` in your local `config.json`

The key is stored only in your local `config.json`, which is **git-ignored**. Never commit it.

AI can:

- Generate one or more candidate molecules from a natural-language description
- Generate an F-substitution series (`--f-range 0-3`)
- Generate a custom calculation chain (validated against the route whitelist, with
  `def2-TZVP → def2TZVP`, `M06-2X → M062X` normalisation)

Model choices: `deepseek-v4-flash` (fast/cheap, default) and `deepseek-v4-pro` (strongest).
Switch in the GUI or with `--model v4flash|v4pro`.

```bash
python gflow.py ai --prompt "generate three fluorinated carbonate candidates" --n-molecules 3 --json
python gflow.py ai-workflow --prompt "B3LYP/6-311+G(d) opt+freq, then B2PLYP/def2TZVP gas SP, then M062X/6-31+G(d) PCM(DMSO) SP" --model v4pro --json
```

---

## Configuration

`config.json` is created automatically at first run. `config.example.json` documents every default.

Key environment variables:

| Variable | Purpose |
|---|---|
| `GFLOW_HOME` | Data home (default: project directory; the packaged exe uses `%LOCALAPPDATA%\GaussFlow`) |
| `GFLOW_CONFIG` | Explicit `config.json` path |
| `G16_HOME` | Gaussian 16 installation directory (optional; `detect` can find it) |
| `DEEPSEEK_API_KEY` | Optional DeepSeek key for AI features |

Important paths:

```
runs/<study_id>/              Gaussian input/output/chk/result.json
generated_gjf/<gen_id>/       AI / F-enumeration .gjf files + index.json
exports/system_<batch_id>/    analysis workbook, CSV, report, XYZ
logs/results/results_*.jsonl  real-time result log
logs/worker.log               background worker log
gaussflow.db                  SQLite state
```

None of the runtime data directories, the database, local config, or API keys are committed to this repository.

---

## Architecture

```
GUI / CLI / Agent
        │  controlled commands (`--json`)
        ▼
gaussflow
  gui.py          desktop GUI (customtkinter)
  cli.py          controlled command interface
  config.py       configuration + path detection
  locate.py       auto-detect Gaussian / GaussView / Shermo
  recipes.py      generic calculation-chain planning
  pipeline.py     study creation, stage chaining, retries, summaries
  structure.py    RDKit 3D embedding, ion/coordination sites
  gjf.py          .gjf writer (including counterpoise BSSE)
  runner.py       g16.exe execution + wall-clock watchdog
  daemon.py       background worker (claim → run → parse → validate → advance)
  parse.py        cclib + regex parsing (energy/frequency/thermo/HOMO-LUMO/BSSE)
  validate.py     scientific acceptance checks
  diagnose.py     failure classification + whitelisted retry advice
  analysis.py     binding energy / redox / orbitals / thermochemistry
  export.py       CSV / Excel / XYZ / Markdown
  resultlog.py    real-time JSONL result log + crash recovery
  fluorinate.py   deterministic F-site enumeration
  llm.py          optional DeepSeek integration
  store.py        SQLite (WAL)
```

---

## Scientific conventions

- Energies in Hartree, with kcal/mol and eV conversions.
- `Eb = E(complex) − E(monomer) − E(Li⁺)`
- Gaussian counterpoise convention: `corrected = raw + BSSE` (`BSSE > 0`), so `Eb(BSSE) = Eb(raw) + BSSE`.
- `G(level) = E(level SP) + [G(opt+freq) − E(opt+freq)]`.
- Redox potentials are **relative**: `E = −ΔG_sol(eV) + reference_offset`; set `config.redox.reference_offset_v` to the literature reference electrode.
- Gaussian 16 does **not** support PCM together with `counterpoise`, so BSSE is calculated only at gas-phase levels; solvent-phase binding energies are uncorrected.
- Imaginary frequencies: `0` for minima, `1` for transition states.
- Gaussian syntax: `def2TZVP` (not `def2-TZVP`), `M062X` (not `M06-2X`), `empiricaldispersion=gd3bj`.

### Validation gates

A study is usable only when all critical checks pass (or only warnings remain):

- Output exists and ends with `Normal termination of Gaussian 16`
- No `Error termination`, `Convergence failure`, `QPErr`
- Atom count, charge, and multiplicity match the input
- SCF converged
- Geometry optimisation converged (`Optimization completed` / `Stationary point found`)
- Imaginary frequency count within tolerance
- Final energy exists and belongs to the last stage
- BSSE task parsed `Counterpoise corrected energy` / `BSSE energy`

---

## Results

`python gflow.py analysis --batch bt_xxx --json` returns:

- `binding[]` — Eb, BSSE, Eb(BSSE), per level and Li⁺ site
- `redox[]` — ΔG_gas, ΔG_sol, relative potential
- `orbitals[]` — HOMO, LUMO, gap
- `thermo[]` — ZPE, H(298), G(298), entropy, optional Shermo G

`export-system` writes `exports/system_<batch_id>/` with:

- `analysis.xlsx` — overview / binding / redox / orbitals / thermochemistry
- `binding_energies.csv`, `redox_potentials.csv`, `frontier_orbitals.csv`
- `analysis.json`, `analysis_report.md`
- `batch/` — generic CSV/Excel/XYZ/report artifacts

---

## Crash safety and recovery

Every completed job and study is appended to `logs/results/results_YYYYMMDD.jsonl` with `flush + fsync`.

- GUI: **Results → Restore from result log**
- Startup: if the log contains completed studies missing from the database, GaussFlow offers to restore them
- Startup: if a previous close left running jobs, GaussFlow offers **Resume unfinished jobs**
- CLI:

```bash
python gflow.py logs --json
python gflow.py recover --json
python gflow.py recover --batch bt_xxx --json
```

`recover` only recreates studies missing from the database; re-running it is safe.

---

## Command reference

| Command | Purpose |
|---|---|
| `doctor` | Environment check (optional Gaussian test run) |
| `detect` | Auto-detect Gaussian / GaussView / Shermo |
| `prepare` | Build structure + route preview (no compute) |
| `submit` / `submit-batch` | Async submission / background worker |
| `status` | Study/job status |
| `result` | Structured result summary |
| `jobs` | Job details for a study |
| `diagnose` | Failure classification + retry advice |
| `export` | CSV / Excel / XYZ / JSON / Markdown |
| `recipes` | List calculation recipes |
| `recipe` | Create a calculation system |
| `analysis` | Binding / redox / orbitals / thermochemistry |
| `export-system` | Export system analysis |
| `ai` | Generate molecules and `.gjf` files (API key required) |
| `ai-workflow` | Generate and save a custom calculation chain (API key required) |
| `fluorinate` | Symmetry-unique F-site enumeration |
| `gjf-batch` | Batch-import `.gjf` files |
| `logs` / `recover` | Inspect / restore the real-time result log |
| `stop` | Kill running Gaussian processes and cancel the queue |
| `purge` | Remove cancelled study records |
| `daemon` | Background worker management |
| `preview` | Open a structure/result in GaussView |

All commands accept `--json` and print exactly one JSON object on stdout.

---

## Tests

```bash
python -m unittest discover -s tests -v
```

The suite covers route validation, `.gjf` writing (including counterpoise), 3D embedding and Li⁺ sites, cclib parsing, HOMO/LUMO, BSSE arithmetic, validation rules, failure diagnosis/retry, SQLite storage, Shermo parsing, recipe planning, F-site enumeration, result-log recovery, and AI response validation.

---

## Security and privacy

- **No API key is included in this repository.**
- `config.json`, `.env`, `logs/`, `runs/`, `exports/`, `generated_gjf/`, and `gaussflow.db` are git-ignored.
- API keys are stored locally in `config.json` and sent only to the configured API endpoint.
- Gaussian input/output files, molecular structures, and results remain on your machine unless you choose to share them.
- Review `config.json` and logs before sharing a workspace or a bug report.

---

## Contributing

Issues and pull requests are welcome. Please run the test suite before opening a PR:

```bash
python -m unittest discover -s tests -v
```

Keep new calculation chains generic and configurable; do not commit paper-specific
production protocols, local paths, API keys, or generated results.

---

## License

MIT License. See [LICENSE](LICENSE).

Gaussian 16 and GaussView 6 are commercial products of Gaussian, Inc.
GaussFlow is an independent automation layer and is not affiliated with Gaussian, Inc.
