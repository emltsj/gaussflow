"""Desktop GUI for GaussFlow (customtkinter).

Tabs
----
计算      molecule input, preset chain, Li⁺/redox options, resources, start/stop
任务/日志  live study/job table, per-job output tail, event log
结果      Li⁺ binding energies (with BSSE), redox potentials, frontier orbitals,
          thermochemistry; export to Excel/CSV and open the folder
设置      auto-detect Gaussian/GaussView/Shermo, resource limits, dependency check

Long calculations run in a background thread (in-process worker) so the window
never blocks; all Gaussian execution still goes through the deterministic core.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path

import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from . import analysis, fluorinate, llm, locate, pipeline, resultlog
from .config import Config
from .store import Store

APP_TITLE = "GaussFlow · Gaussian 16 自动化计算平台"
FONT = "Microsoft YaHei UI"
FONT_MONO = "Consolas"

ROLE_LABELS = {
    "neutral": "中性分子",
    "li_ref": "Li⁺ 参考",
    "redox_cation": "阳离子(氧化)",
    "redox_anion": "阴离子(还原)",
}


def _role_label(role: str) -> str:
    if role in ROLE_LABELS:
        return ROLE_LABELS[role]
    if role.startswith("li_") and "redox" in role:
        return "Li⁺配合物-" + ("氧化" if role.endswith("cation") else "还原")
    if role.startswith("li_"):
        return "Li⁺配合物@" + role[3:]
    return role or ""


def _fmt(v, dec: int = 3) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{dec}f}"
    return str(v)


# ---------------------------------------------------------------------------
class CalcWorker(threading.Thread):
    """Runs pending jobs of the current batch inside the GUI process."""

    def __init__(self, cfg: Config, app: "App"):
        super().__init__(daemon=True)
        self.cfg = cfg
        self.app = app
        self.stop_event = threading.Event()

    def run(self) -> None:
        from .daemon import Worker
        try:
            worker = Worker(self.cfg)
        except Exception:
            self.app.after(0, lambda: self.app._on_run_error(traceback.format_exc()))
            return
        # A forced GUI close leaves jobs in "running" with no live process.
        # Recover them first (parse finished output / mark interrupted), then
        # continue with any pending or auto-retry jobs.
        try:
            worker.recover()
        except Exception:
            worker.store.add_event(
                "GUI worker recovery error: " + traceback.format_exc()[:900],
                level="error")
        while not self.stop_event.is_set():
            try:
                if worker.store.pending_count() == 0:
                    break
                worker.process_one()
            except Exception:
                worker.store.add_event(
                    "GUI worker error: " + traceback.format_exc()[:900], level="error")
                time.sleep(1.0)
        self.app.after(0, self.app._on_run_finished)

    def stop(self) -> None:
        self.stop_event.set()


# ---------------------------------------------------------------------------
class App(ctk.CTk):
    def __init__(self) -> None:
        ctk.set_appearance_mode("light")
        ctk.set_default_color_theme("blue")
        super().__init__()
        self.title(APP_TITLE)
        # fit any screen (small laptops included) and centre the window
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        w, h = min(1240, sw - 60), min(840, sh - 60)
        x, y = max(0, (sw - w) // 2), max(0, (sh - h) // 2)
        self.geometry(f"{w}x{h}+{x}+{y}")
        self.minsize(1000, 640)
        self.cfg = Config.load()
        self.batch_id: str | None = None
        self.runner: CalcWorker | None = None
        self._busy = False
        self._refresh_busy = False
        self._log_seen = 0
        self._gjf_files: list[str] = []
        self._ai_molecules: list[dict] = []
        self._ai_busy = False
        self.batch_ids: list[str] = []
        self.font_scale = 1.2
        self._stopped_by_user = False
        self._build_ui()
        self.after(300, self._startup_detect)
        self.after(900, self._load_latest_batch)
        self.after(1400, self._auto_recover_from_logs)
        self.after(2600, self._check_unfinished_jobs)
        self.after(1200, self._refresh_loop)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _load_latest_batch(self) -> None:
        """Show the most recent batch (e.g. one submitted from the CLI)."""
        try:
            store = Store(self.cfg.db_path)
            studies = store.list_studies(limit=500)
            store.close()
            batch = next((s.get("batch_id") for s in reversed(studies)
                          if s.get("batch_id")), None)
            if batch and not self.batch_id:
                self.batch_id = batch
                self._refresh_batch_choices()
                self._refresh_tasks()
                self._refresh_results(force=True)
                self.run_status.configure(text=f"已载入最近批次：{batch}")
        except Exception:
            pass

    # --- batch selector ----------------------------------------------------
    def _batch_entries(self) -> list[tuple[str, str]]:
        """Return [(display, batch_id)] of recent batches (newest first).

        Batches whose studies are *all* cancelled are hidden from the selector
        (they can be removed permanently via 「清理取消记录」).
        """
        try:
            store = Store(self.cfg.db_path)
            studies = store.list_studies(limit=800)
            store.close()
        except Exception:
            return []
        order: dict[str, str] = {}
        has_active: dict[str, bool] = {}
        for s in studies:  # ascending creation order
            b = s.get("batch_id")
            if not b:
                continue
            if b not in order:
                order[b] = (s.get("name") or "").split("-")[0][:18]
            if s.get("status") != "cancelled":
                has_active[b] = True
        entries = []
        for b in reversed(list(order.keys())):
            if not has_active.get(b) and b != self.batch_id:
                continue
            entries.append((f"{b} | {order[b]}", b))
        return entries[:40]

    def _refresh_batch_choices(self) -> None:
        if not hasattr(self, "batch_combo"):
            return
        entries = self._batch_entries()
        values = [label for label, _ in entries] or ["（无批次）"]
        self.batch_combo.configure(values=values)
        chosen = values[0]
        for label, b in entries:
            if b == self.batch_id:
                chosen = label
                break
        self.batch_combo.set(chosen)

    def _on_batch_combo(self, choice: str) -> None:
        if not choice or choice.startswith("（"):
            return
        batch_id = choice.split(" | ")[0].strip()
        if batch_id == self.batch_id:
            return
        self.batch_id = batch_id
        self._refresh_tasks()
        self._refresh_results(force=True)
        self.run_status.configure(text=f"已载入批次：{batch_id}")

    # --- crash-safe recovery from the real-time result log -----------------
    def _recover_from_logs(self, silent: bool = False) -> None:
        try:
            from .resultlog import recover, summary as log_summary
            store = Store(self.cfg.db_path)
            res = recover(store, self.cfg)
            store.close()
        except Exception as exc:
            messagebox.showerror("恢复失败", f"{type(exc).__name__}: {exc}")
            return
        if res.get("count"):
            entries = self._batch_entries()
            if entries:
                self.batch_id = entries[0][1]
            self._refresh_batch_choices()
            self._refresh_tasks()
            self._refresh_results(force=True)
        if not silent:
            info = log_summary(self.cfg)
            messagebox.showinfo(
                "从日志恢复",
                f"恢复 {res.get('count', 0)} 个研究（已有记录跳过 {res.get('skipped', 0)} 个）。\n"
                f"日志目录：{info.get('log_dir')}\n"
                f"日志内容：{info.get('counts')}\n"
                f"最新记录：{info.get('latest') or '—'}")

    def _auto_recover_from_logs(self) -> None:
        """If completed results exist in the JSONL log but not in the DB, offer
        to restore them (covers a deleted/corrupted DB or a forced close)."""
        try:
            store = Store(self.cfg.db_path)
            from .resultlog import missing_studies
            missing = missing_studies(store, self.cfg)
            store.close()
            if missing.get("count", 0) <= 0:
                return
            if messagebox.askyesno(
                    "发现结果日志",
                    f"实时结果日志里有 {missing['count']} 个已完成研究不在当前数据库中。\n"
                    f"是否从日志恢复这些结果显示？\n\n{missing.get('log_dir')}"):
                self._recover_from_logs(silent=True)
        except Exception:
            pass

    def _check_unfinished_jobs(self) -> None:
        """Offer to resume jobs left running by a forced close of the GUI."""
        if self._busy or self.runner:
            return
        try:
            from .daemon import daemon_status
            if daemon_status(self.cfg).get("running"):
                return
            store = Store(self.cfg.db_path)
            running = len(store.running_jobs())
            pending = max(0, store.pending_count() - running)
            store.close()
            if running <= 0 and pending <= 0:
                return
            if messagebox.askyesno(
                    "继续未完成任务",
                    f"发现上次退出时未完成的任务：\n"
                    f"  运行中（已中断）：{running} 个\n  排队：{pending} 个\n\n"
                    "是否立即回收并继续计算？\n"
                    "（会先解析已完成的 Gaussian 输出，再继续排队任务）"):
                self._resume_jobs()
        except Exception:
            pass

    def _resume_jobs(self) -> None:
        """Start the in-process worker without creating a new batch."""
        if self._busy:
            return
        self._busy = True
        self._stopped_by_user = False
        self.run_btn.configure(state="disabled")
        self.stop_btn.configure(state="normal")
        self.progress.configure(mode="indeterminate")
        self.progress.start()
        self.run_status.configure(text="正在回收上次未完成任务并继续计算…")
        self.tabs.set("任务 / 日志")
        self.runner = CalcWorker(self.cfg, self)
        self.runner.start()

    # ------------------------------------------------------------------ UI
    def _build_ui(self) -> None:
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        header = ctk.CTkFrame(self, corner_radius=0, fg_color="#1F4E78", height=64)
        header.grid(row=0, column=0, sticky="ew")
        header.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(header, text="  GaussFlow", font=(FONT, 22, "bold"),
                     text_color="white").grid(row=0, column=0, padx=(16, 8), pady=12)
        ctk.CTkLabel(header, text="Gaussian 16 计算链 · 自动生成输入 / 提交 / 验收 / 汇总",
                     font=(FONT, 12), text_color="#D6E4F0").grid(
            row=0, column=1, sticky="w")
        self.status_badge = ctk.CTkLabel(header, text="● 检测中…", font=(FONT, 12, "bold"),
                                         text_color="#FFD966")
        self.status_badge.grid(row=0, column=2, padx=16)

        self.tabs = ctk.CTkTabview(self, corner_radius=10)
        self.tabs.grid(row=1, column=0, sticky="nsew", padx=12, pady=(10, 12))
        for name in ("计算", "任务 / 日志", "结果", "设置"):
            self.tabs.add(name)

        self._build_calc_tab(self.tabs.tab("计算"))
        self._build_tasks_tab(self.tabs.tab("任务 / 日志"))
        self._build_results_tab(self.tabs.tab("结果"))
        self._build_settings_tab(self.tabs.tab("设置"))
        self.after(300, self._apply_fonts)

    # ------------------------------------------------------------- 计算页
    def _build_calc_tab(self, parent) -> None:
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        left = ctk.CTkScrollableFrame(parent, corner_radius=10, label_text="分子与计算方案")
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 8), pady=4)
        right = ctk.CTkScrollableFrame(parent, corner_radius=10, label_text="选项与运行")
        right.grid(row=0, column=1, sticky="nsew", padx=(8, 0), pady=4)
        self.left_scroll, self.right_scroll = left, right

        def label(master, text, row, bold=False):
            ctk.CTkLabel(master, text=text, font=(FONT, 12, "bold" if bold else "normal"),
                         anchor="w").grid(row=row, column=0, sticky="w", padx=12,
                                          pady=(10 if bold else 6, 2))

        # --- molecule -------------------------------------------------
        label(left, "① 分子输入", 0, True)
        self.mode_var = tk.StringVar(value="SMILES")
        seg = ctk.CTkSegmentedButton(
            left, values=["SMILES", "GJF 文件（可多选）", "AI 生成（DeepSeek）"],
            variable=self.mode_var, command=self._on_mode_change)
        seg.grid(row=1, column=0, sticky="ew", padx=12, pady=4)

        # SMILES panel
        self.smiles_panel = ctk.CTkFrame(left, fg_color="transparent")
        self.smiles_panel.grid(row=2, column=0, sticky="ew")
        self.smiles_entry = ctk.CTkEntry(self.smiles_panel,
                                        placeholder_text="例如：CCO  或  [Li+]")
        self.smiles_entry.pack(fill="x", padx=12, pady=4)

        # GJF multi-select panel
        self.gjf_panel = ctk.CTkFrame(left, fg_color="transparent")
        self.gjf_panel.grid(row=2, column=0, sticky="ew")
        btn_row = ctk.CTkFrame(self.gjf_panel, fg_color="transparent")
        btn_row.pack(fill="x", padx=12, pady=4)
        ctk.CTkButton(btn_row, text="选择 .gjf 文件（可多选）", width=190,
                      command=self._pick_gjf_files).pack(side="left")
        ctk.CTkButton(btn_row, text="清空", width=60, fg_color="#5A6B7B",
                      command=self._clear_gjf_files).pack(side="left", padx=(6, 0))
        self.gjf_list = ctk.CTkTextbox(self.gjf_panel, height=90, font=(FONT_MONO, 11))
        self.gjf_list.pack(fill="x", padx=12, pady=(0, 4))
        self.gjf_list.insert("end", "（尚未选择；可一次选择多个 gjf，排队计算）")
        self.gjf_list.configure(state="disabled")
        self.gjf_panel.grid_remove()

        # AI generation panel
        self.ai_panel = ctk.CTkFrame(left, fg_color="transparent")
        self.ai_panel.grid(row=2, column=0, sticky="ew")
        ctk.CTkLabel(self.ai_panel,
                     text="用自然语言描述分子（DeepSeek 自动生成结构 → gjf）",
                     font=(FONT, 11), anchor="w").pack(fill="x", padx=12)
        self.ai_prompt = ctk.CTkTextbox(self.ai_panel, height=78, font=(FONT, 12))
        self.ai_prompt.pack(fill="x", padx=12, pady=4)
        self.ai_prompt.insert(
            "end", "例如：生成 LiFSI 电解质的溶剂分子 TEP；或：给出 3 个不同的含氟溶剂候选")
        f_row = ctk.CTkFrame(self.ai_panel, fg_color="transparent")
        f_row.pack(fill="x", padx=12, pady=2)
        ctk.CTkLabel(f_row, text="一次生成分子数", font=(FONT, 11)).pack(side="left")
        self.ai_n_mol = ctk.CTkEntry(f_row, width=52, placeholder_text="1")
        self.ai_n_mol.pack(side="left", padx=(4, 12))
        self.ai_n_mol.insert(0, "1")
        ctk.CTkLabel(f_row, text="F 原子个数（AI 直接决定，可留空）",
                     font=(FONT, 11)).pack(side="left")
        self.ai_f_range = ctk.CTkEntry(f_row, width=80, placeholder_text="0-3")
        self.ai_f_range.pack(side="left", padx=(4, 0))
        model_row = ctk.CTkFrame(self.ai_panel, fg_color="transparent")
        model_row.pack(fill="x", padx=12, pady=2)
        ctk.CTkLabel(model_row, text="AI 模型", font=(FONT, 11)).pack(side="left")
        self.ai_model = ctk.CTkOptionMenu(
            model_row, width=190,
            values=[llm.MODEL_LABELS[m] for m in llm.MODEL_CHOICES],
            command=self._on_ai_model_change)
        self.ai_model.set(llm.model_label(self.cfg.get("llm.model")))
        self.ai_model.pack(side="left", padx=(4, 12))
        ctk.CTkLabel(model_row, text="（V4 Flash 快；V4 Pro 最强）",
                     font=(FONT, 11), text_color="#666").pack(side="left")
        enum_row = ctk.CTkFrame(self.ai_panel, fg_color="transparent")
        enum_row.pack(fill="x", padx=12, pady=2)
        self.ai_enum_f = ctk.CTkCheckBox(
            enum_row, text="按对称不等价位点枚举 F（推荐：确定性生成，AI 只给骨架）",
            font=(FONT, 11), command=self._on_enum_toggle)
        self.ai_enum_f.pack(side="left")
        self.ai_enum_n = ctk.CTkEntry(enum_row, width=48, placeholder_text="1")
        self.ai_enum_n.pack(side="left", padx=(8, 0))
        self.ai_enum_n.insert(0, "1")
        self.ai_enum_n.configure(state="disabled")
        self.ai_gen_btn = ctk.CTkButton(self.ai_panel, text="① 生成分子并导出 gjf",
                                        command=self._ai_generate, height=32)
        self.ai_gen_btn.pack(fill="x", padx=12, pady=4)
        self.ai_status = ctk.CTkLabel(self.ai_panel, text="", font=(FONT, 11),
                                      anchor="w", justify="left", wraplength=520,
                                      text_color="#555")
        self.ai_status.pack(fill="x", padx=12)
        self.ai_list = ctk.CTkTextbox(self.ai_panel, height=96, font=(FONT_MONO, 11))
        self.ai_list.pack(fill="x", padx=12, pady=4)
        self.ai_list.insert(
            "end", "（AI 生成的分子显示在这里，并写入 generated_gjf/ 目录）")
        self.ai_list.configure(state="disabled")
        self.ai_panel.grid_remove()

        self.name_entry = ctk.CTkEntry(left, placeholder_text="体系名称（默认自动生成）")
        self.name_entry.grid(row=4, column=0, sticky="ew", padx=12, pady=4)

        # --- recipe ---------------------------------------------------
        label(left, "② 计算方案（预设 / AI 生成）", 5, True)
        recipes = pipeline_recipe_names(self.cfg)
        self.recipe_var = tk.StringVar(value=self._default_recipe_key())
        self.recipe_menu = ctk.CTkOptionMenu(
            left, values=list(recipes.values()), variable=self.recipe_var,
            command=self._on_recipe_change)
        self.recipe_menu.grid(row=6, column=0, sticky="ew", padx=12, pady=4)
        ctk.CTkButton(left, text="AI 生成计算链…（DeepSeek，按需求自定义）",
                      width=340, fg_color="#5A6B7B",
                      command=self._ai_workflow_dialog).grid(
            row=7, column=0, sticky="w", padx=12, pady=(0, 2))
        self.recipe_desc = ctk.CTkLabel(left, text="", font=(FONT, 11), justify="left",
                                        anchor="w", wraplength=520, text_color="#555")
        self.recipe_desc.grid(row=8, column=0, sticky="ew", padx=12, pady=(0, 6))
        self._recipes = recipes

        # --- stage overview ------------------------------------------
        label(left, "③ 理论级别（当前方案）", 9, True)
        self.stages_box = ctk.CTkTextbox(left, height=150, font=(FONT_MONO, 11))
        self.stages_box.grid(row=10, column=0, sticky="ew", padx=12, pady=4)
        self.stages_box.configure(state="disabled")

        # --- right: options ------------------------------------------
        label(right, "Li⁺ 配位与 BSSE", 0, True)
        self.li_switch = ctk.CTkSwitch(right, text="计算 Li⁺ 配位配合物（含 counterpoise BSSE）")
        self.li_switch.grid(row=1, column=0, sticky="w", padx=12, pady=4)
        self.li_switch.select()
        li_row = ctk.CTkFrame(right, fg_color="transparent")
        li_row.grid(row=2, column=0, sticky="w", padx=24, pady=2)
        self.li_o = ctk.CTkCheckBox(li_row, text="O 位点")
        self.li_o.pack(side="left", padx=(0, 16))
        self.li_f = ctk.CTkCheckBox(li_row, text="F 位点")
        self.li_f.pack(side="left")
        self.li_o.select()

        label(right, "氧化还原（阳离子 / 阴离子）", 3, True)
        self.redox_switch = ctk.CTkSwitch(right, text="计算氧化态与还原态（opt+freq + SP）")
        self.redox_switch.grid(row=4, column=0, sticky="w", padx=12, pady=4)
        self.redox_switch.select()
        redox_row = ctk.CTkFrame(right, fg_color="transparent")
        redox_row.grid(row=5, column=0, sticky="w", padx=24, pady=2)
        self.redox_neutral = ctk.CTkCheckBox(redox_row, text="中性分子")
        self.redox_neutral.pack(side="left", padx=(0, 16))
        self.redox_li = ctk.CTkCheckBox(redox_row, text="Li⁺ 配合物")
        self.redox_li.pack(side="left")
        self.redox_neutral.select()
        self.redox_li.select()

        label(right, "预优化（防 SCF/几何爆炸）", 6, True)
        self.preopt_switch = ctk.CTkSwitch(right, text="先做 HF/3-21G 预优化")
        self.preopt_switch.grid(row=7, column=0, sticky="w", padx=12, pady=4)
        self.preopt_switch.select()

        label(right, "溶剂 / 电荷 / 资源", 8, True)
        form = ctk.CTkFrame(right, fg_color="transparent")
        form.grid(row=9, column=0, sticky="ew", padx=12, pady=4)
        form.grid_columnconfigure(1, weight=1)
        self.solvent_entry = ctk.CTkEntry(form)
        self.solvent_entry.insert(0, self.cfg.get("chemistry.solvent") or "DMSO")
        self._row(form, "隐式溶剂 (PCM)", self.solvent_entry, 0)
        self.charge_entry = ctk.CTkEntry(form, placeholder_text="自动")
        self._row(form, "总电荷", self.charge_entry, 1)
        self.mult_entry = ctk.CTkEntry(form, placeholder_text="自动")
        self._row(form, "自旋多重度", self.mult_entry, 2)
        self.nproc_entry = ctk.CTkEntry(form, width=80)
        self.nproc_entry.insert(0, str(self.cfg.resources["nproc"]))
        self._row(form, "并行核数", self.nproc_entry, 3)
        self.mem_entry = ctk.CTkEntry(form, width=120)
        self.mem_entry.insert(0, str(self.cfg.resources["mem"]))
        self._row(form, "内存 (%mem)", self.mem_entry, 4)

        label(right, "运行", 10, True)
        self.progress = ctk.CTkProgressBar(right)
        self.progress.grid(row=11, column=0, sticky="ew", padx=12, pady=(4, 2))
        self.progress.set(0)
        self.run_status = ctk.CTkLabel(right, text="就绪", anchor="w",
                                       font=(FONT, 11), text_color="#555")
        self.run_status.grid(row=12, column=0, sticky="ew", padx=12)
        btn_row = ctk.CTkFrame(right, fg_color="transparent")
        btn_row.grid(row=13, column=0, sticky="ew", padx=12, pady=10)
        self.run_btn = ctk.CTkButton(btn_row, text="开始计算", height=40,
                                     font=(FONT, 13, "bold"), command=self._start_calc)
        self.run_btn.pack(side="left", fill="x", expand=True)
        self.stop_btn = ctk.CTkButton(btn_row, text="停止", width=100, height=40,
                                      fg_color="#B22222", hover_color="#8B1A1A",
                                      command=self._stop_calc, state="disabled")
        self.stop_btn.pack(side="left", padx=(8, 0))
        self.open_dir_btn = ctk.CTkButton(right, text="打开任务目录", height=32,
                                          fg_color="#5A6B7B", command=self._open_runs_dir)
        self.open_dir_btn.grid(row=14, column=0, sticky="ew", padx=12, pady=(0, 10))

        self._on_recipe_change(self.recipe_var.get())
        self.mode_var.set("SMILES")
        self._on_mode_change()
        self.after(120, self._reset_scroll)

    def _reset_scroll(self) -> None:
        for frame in (getattr(self, "left_scroll", None), getattr(self, "right_scroll", None)):
            try:
                frame._parent_canvas.yview_moveto(0.0)
            except Exception:
                pass

    def _row(self, form, text, widget, row):
        ctk.CTkLabel(form, text=text, font=(FONT, 11), anchor="w").grid(
            row=row, column=0, sticky="w", pady=3)
        widget.grid(row=row, column=1, sticky="w", pady=3)

    def _on_mode_change(self, value: str | None = None) -> None:
        mode = self.mode_var.get()
        for panel in (getattr(self, "smiles_panel", None),
                      getattr(self, "gjf_panel", None),
                      getattr(self, "ai_panel", None)):
            if panel is not None:
                panel.grid_remove()
        if mode.startswith("GJF"):
            self.gjf_panel.grid()
        elif mode.startswith("AI"):
            self.ai_panel.grid()
        else:
            self.smiles_panel.grid()

    # --- .gjf multi-select -------------------------------------------------
    def _pick_gjf_files(self) -> None:
        paths = filedialog.askopenfilenames(
            title="选择 Gaussian 输入文件（可多选）",
            filetypes=[("Gaussian 输入", "*.gjf *.com"), ("全部", "*.*")])
        for p in paths or ():
            if p not in self._gjf_files:
                self._gjf_files.append(p)
        self._refresh_gjf_list()

    def _clear_gjf_files(self) -> None:
        self._gjf_files = []
        self._refresh_gjf_list()

    def _refresh_gjf_list(self) -> None:
        self.gjf_list.configure(state="normal")
        self.gjf_list.delete("1.0", "end")
        if not self._gjf_files:
            self.gjf_list.insert("end", "（尚未选择；可一次选择多个 gjf，排队计算）")
        else:
            self.gjf_list.insert("end", f"已选择 {len(self._gjf_files)} 个文件：\n")
            for i, p in enumerate(self._gjf_files):
                self.gjf_list.insert("end", f"{i + 1}. {p}\n")
        self.gjf_list.configure(state="disabled")

    # --- AI (DeepSeek) molecule generation ---------------------------------
    def _on_ai_model_change(self, choice: str) -> None:
        """Persist the quick model switch (V4 Flash / V4 Pro)."""
        model = llm.model_from_label(choice, self.cfg)
        self.cfg.set("llm.model", model)
        try:
            self.cfg.save()
        except Exception:
            pass
        if getattr(self, "llm_model", None) is not None:
            self.llm_model.set(llm.model_label(model))

    def _on_enum_toggle(self) -> None:
        try:
            self.ai_enum_n.configure(
                state="normal" if self.ai_enum_f.get() else "disabled")
        except Exception:
            pass

    def _ai_generate(self) -> None:
        if self._ai_busy:
            return
        prompt = self.ai_prompt.get("1.0", "end").strip()
        if not prompt or prompt.startswith("例如："):
            messagebox.showwarning("缺少描述", "请先用自然语言描述要生成的分子。")
            return
        try:
            n_mol = max(1, int(self.ai_n_mol.get().strip() or 1))
        except ValueError:
            messagebox.showerror("参数错误", "生成分子数必须是整数。")
            return
        f_range = self.ai_f_range.get().strip() or None
        enum = bool(self.ai_enum_f.get())
        try:
            enum_n = max(1, int(self.ai_enum_n.get().strip() or 1))
        except ValueError:
            messagebox.showerror("参数错误", "枚举 F 个数必须是整数。")
            return
        # use the model selected right above the prompt, and remember it
        self.cfg.set("llm.model", llm.model_from_label(self.ai_model.get(), self.cfg))
        try:
            self.cfg.save()
        except Exception:
            pass
        self._ai_busy = True
        recipe_key = self._recipe_key()
        name_prefix = self.name_entry.get().strip() or "AI"
        self.ai_gen_btn.configure(state="disabled", text="正在调用 DeepSeek…")
        self.ai_status.configure(text="正在调用 DeepSeek API 生成分子结构…")

        def work():
            try:
                res = llm.generate_molecules(
                    self.cfg, prompt, f_range=f_range, n_molecules=n_mol)
                if enum:
                    res = self._enumerate_f_variants(res, enum_n, f_range)
                preview = pipeline.export_gjf_previews(
                    self.cfg, res["molecules"],
                    recipe_name=recipe_key, name_prefix=name_prefix)
                res["preview"] = preview
                self.after(0, lambda: self._on_ai_generated(res))
            except Exception as exc:
                msg = f"{type(exc).__name__}: {exc}"
                self.after(0, lambda: self._on_ai_error(msg))

        threading.Thread(target=work, daemon=True).start()

    def _enumerate_f_variants(self, res: dict, enum_n: int,
                              f_range: str | None) -> dict:
        """Deterministically build one molecule per symmetry-unique F site."""
        maxv = int(self.cfg.get("llm.max_fluorination_variants") or 20)
        counts = llm.parse_f_range(f_range) or [enum_n]
        molecules: list[dict] = []
        site_notes = []
        for base in res["molecules"]:
            built, notes = fluorinate.molecules_for_f_counts(
                base, counts, max_variants=maxv)
            molecules.extend(built)
            site_notes.extend(notes)
        if not molecules:
            raise RuntimeError("按位点枚举 F 没有得到任何结构：" + "；".join(site_notes))
        res = dict(res)
        res["molecules"] = molecules
        res["enumeration"] = {"n_f_values": counts, "max_variants": maxv,
                              "notes": site_notes}
        res["explanation"] = ((res.get("explanation") or "")
                              + f"\n[位点枚举] 对每个骨架按对称不等价位点替换 "
                                f"F×{counts}，共 {len(molecules)} 个结构。").strip()
        return res

    def _on_ai_generated(self, res: dict) -> None:
        self._ai_busy = False
        self.ai_gen_btn.configure(state="normal", text="① 生成分子并导出 gjf")
        self._ai_molecules = res["molecules"]
        lines = []
        files = (res.get("preview") or {}).get("files") or []
        for i, m in enumerate(self._ai_molecules):
            gjf = files[i].get("gjf", "") if i < len(files) else ""
            site = f"  [{m['site_label']}]" if m.get("site_label") else ""
            lines.append(
                f"{i + 1}. {m['name']}  {m['formula']}  F={m['n_F']}  "
                f"q/m={m['charge']}/{m['multiplicity']}{site}\n"
                f"   SMILES: {m['smiles']}\n   gjf: {gjf}")
        self.ai_list.configure(state="normal")
        self.ai_list.delete("1.0", "end")
        self.ai_list.insert("end", "\n".join(lines))
        self.ai_list.configure(state="disabled")
        note = res.get("explanation") or ""
        missing = res.get("missing_f_counts") or []
        status = (f"生成 {len(self._ai_molecules)} 个分子；gjf 已写入 "
                  f"{(res.get('preview') or {}).get('gjf_dir', '')}")
        if res.get("enumeration"):
            status += f"；按位点枚举 F 数 = {res['enumeration'].get('n_f_values')}"
            site_kinds = sorted({s for m in self._ai_molecules
                                 for s in (m.get("sites") or [])})
            if site_kinds:
                status += f"；覆盖位点 {', '.join(site_kinds)}"
        if missing:
            status += f"；缺少 F 个数 {missing}"
        if note:
            status += f"\n{note[:260]}"
        self.ai_status.configure(text=status)

    def _on_ai_error(self, msg: str) -> None:
        self._ai_busy = False
        self.ai_gen_btn.configure(state="normal", text="① 生成分子并导出 gjf")
        self.ai_status.configure(text="生成失败：" + msg[:400])
        messagebox.showerror("AI 生成失败", msg[-1500:])

    # --- AI generated calculation chain -----------------------------------
    def _ai_workflow_dialog(self) -> None:
        win = ctk.CTkToplevel(self)
        win.title("AI 生成计算链（DeepSeek）")
        win.geometry("760x520")
        win.transient(self)
        win.grab_set()
        ctk.CTkLabel(win, text="用自然语言描述需要的计算链（方法/基组/溶剂/Li⁺/氧化还原）",
                     font=(FONT, 13, "bold")).pack(anchor="w", padx=12, pady=(10, 2))
        txt = ctk.CTkTextbox(win, height=150, font=(FONT, 12))
        txt.pack(fill="x", padx=12, pady=4)
        txt.insert("end",
                   "例如：先在 B3LYP/6-311+G(d) 加 gd3bj 优化并算频率，再用 "
                   "B2PLYP/def2TZVP 加 gd3bj 算气相单点，最后在 DMSO 中用 "
                   "M06-2X/6-31+G(d) 算溶剂化单点并输出 HOMO/LUMO；"
                   "不考虑 Li+ 和氧化还原。")
        wf_model_row = ctk.CTkFrame(win, fg_color="transparent")
        wf_model_row.pack(fill="x", padx=12, pady=(0, 2))
        ctk.CTkLabel(wf_model_row, text="模型", font=(FONT, 11)).pack(side="left")
        wf_model = ctk.CTkOptionMenu(
            wf_model_row, width=190,
            values=[llm.MODEL_LABELS[m] for m in llm.MODEL_CHOICES])
        wf_model.set(llm.model_label(self.cfg.get("llm.model")))
        wf_model.pack(side="left", padx=(4, 8))
        ctk.CTkLabel(wf_model_row, text="计算链设计建议用 V4 Pro",
                     font=(FONT, 11), text_color="#666").pack(side="left")
        btns = ctk.CTkFrame(win, fg_color="transparent")
        btns.pack(fill="x", padx=12, pady=4)
        gen_btn = ctk.CTkButton(btns, text="生成并保存计算链", width=170,
                                command=lambda: generate())
        gen_btn.pack(side="left")
        ctk.CTkButton(btns, text="取消", width=80, fg_color="#5A6B7B",
                      command=win.destroy).pack(side="left", padx=8)
        status = ctk.CTkLabel(win, text="", font=(FONT, 11), anchor="w",
                              justify="left", wraplength=700, text_color="#555")
        status.pack(fill="x", padx=12)
        out = ctk.CTkTextbox(win, height=150, font=(FONT_MONO, 10))
        out.pack(fill="both", expand=True, padx=12, pady=(4, 10))
        out.configure(state="disabled")

        def generate():
            prompt = txt.get("1.0", "end").strip()
            if not prompt or prompt.startswith("例如："):
                status.configure(text="请先描述计算链。")
                return
            model = llm.model_from_label(wf_model.get(), self.cfg)
            self.cfg.set("llm.model", model)
            try:
                self.cfg.save()
            except Exception:
                pass
            if getattr(self, "ai_model", None) is not None:
                self.ai_model.set(llm.model_label(model))
            gen_btn.configure(state="disabled", text="生成中…")
            status.configure(text="正在调用 DeepSeek（thinking 模式）设计计算链…")

            def work():
                try:
                    res = llm.generate_workflow(self.cfg, prompt)
                    reg = pipeline.register_ai_recipe(self.cfg, res["recipe"])
                    self.after(0, lambda: done(res, reg))
                except Exception as exc:
                    msg = f"{type(exc).__name__}: {exc}"
                    self.after(0, lambda: fail(msg))

            threading.Thread(target=work, daemon=True).start()

        def fail(msg: str) -> None:
            gen_btn.configure(state="normal", text="生成并保存计算链")
            status.configure(text="生成失败：" + msg[:600])

        def done(res: dict, reg: dict) -> None:
            gen_btn.configure(state="normal", text="生成并保存计算链")
            if not reg.get("ok"):
                status.configure(
                    text="生成的计算链未通过校验：" + "；".join(reg.get("problems") or [])[:600])
                return
            self._recipes = pipeline_recipe_names(self.cfg)
            self.recipe_menu.configure(values=list(self._recipes.values()))
            self.recipe_var.set(self._recipes.get(reg["key"], reg["key"]))
            self._on_recipe_change()
            chain = reg["recipe"]
            lines = [f"方案: {chain['name']}   (key={reg['key']})",
                     f"说明: {chain.get('description', '')}",
                     f"溶剂: {chain.get('solvent')}", "", "阶段:"]
            for st in chain["core_stages"]:
                lines.append(f"  {st['name']:<10} {st['kind']:<8} "
                             f"{st.get('phase', 'gas'):<8} {st['route']}")
            lines += ["", "Li⁺ 配合物: " + ("是" if chain["li"]["on"] else "否"),
                      "氧化还原: " + ("是" if chain["redox"]["on"] else "否"),
                      "", "设计说明: " + (res.get("explanation") or "")]
            out.configure(state="normal")
            out.delete("1.0", "end")
            out.insert("end", "\n".join(lines))
            out.configure(state="disabled")
            status.configure(text="✓ 已保存为计算方案并自动选中；可直接点「开始计算」。")
            win.after(2500, win.destroy)

    # ------------------------------------------------------------- 任务页
    def _build_tasks_tab(self, parent) -> None:
        parent.grid_rowconfigure(0, weight=3)
        parent.grid_rowconfigure(1, weight=2)
        parent.grid_columnconfigure(0, weight=1)
        card = ctk.CTkFrame(parent, corner_radius=10)
        card.grid(row=0, column=0, sticky="nsew", pady=(0, 8))
        card.grid_rowconfigure(1, weight=1)
        card.grid_columnconfigure(0, weight=1)
        bar = ctk.CTkFrame(card, fg_color="transparent")
        bar.grid(row=0, column=0, sticky="ew", padx=10, pady=(8, 4))
        ctk.CTkLabel(bar, text="研究列表", font=(FONT, 14, "bold")).pack(side="left")
        ctk.CTkButton(bar, text="刷新", width=70, command=self._refresh_tasks).pack(side="right")
        ctk.CTkButton(bar, text="GaussView 查看", width=120,
                      command=self._gview_selected).pack(side="right", padx=6)
        ctk.CTkButton(bar, text="A+", width=44,
                      command=lambda: self._zoom(+0.1)).pack(side="right", padx=(2, 0))
        ctk.CTkButton(bar, text="A−", width=44,
                      command=lambda: self._zoom(-0.1)).pack(side="right", padx=(2, 0))
        self.zoom_label = ctk.CTkLabel(bar, text="字体 120%", font=(FONT, 12))
        self.zoom_label.pack(side="right", padx=6)
        ctk.CTkButton(bar, text="清理取消记录", width=110, fg_color="#5A6B7B",
                      command=self._purge_cancelled).pack(side="right", padx=6)
        self.hide_cancelled = ctk.CTkCheckBox(
            bar, text="隐藏已取消", font=(FONT, 12), command=self._refresh_tasks)
        self.hide_cancelled.pack(side="right", padx=6)
        self.hide_cancelled.select()

        cols = ("name", "role", "status", "stage", "validation", "energy", "homo", "nimag")
        heads = ("研究", "角色", "状态", "当前阶段", "验收", "E / Hartree",
                 "HOMO / eV", "虚频")
        frame = tk.Frame(card, bg="#E8EBEE")
        frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=(0, 10))
        self.task_tree = self._make_tree(frame, cols, heads,
                                         widths=(230, 130, 90, 100, 80, 150, 110, 70))

        log_card = ctk.CTkFrame(parent, corner_radius=10)
        log_card.grid(row=1, column=0, sticky="nsew")
        log_card.grid_rowconfigure(1, weight=1)
        log_card.grid_columnconfigure(0, weight=1)
        log_head = ctk.CTkFrame(log_card, fg_color="transparent")
        log_head.grid(row=0, column=0, sticky="ew", padx=10, pady=(8, 2))
        ctk.CTkLabel(log_head, text="运行日志（事件 / 当前任务输出尾部）",
                     font=(FONT, 13, "bold")).pack(side="left")
        ctk.CTkLabel(log_head, text="可上下 / 左右滚动，A+/A− 调整字号",
                     font=(FONT, 11), text_color="#666").pack(side="left", padx=10)
        log_frame = tk.Frame(log_card, bg="white")
        log_frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=(0, 10))
        log_frame.grid_rowconfigure(0, weight=1)
        log_frame.grid_columnconfigure(0, weight=1)
        self.log_text = tk.Text(log_frame, font=(FONT_MONO, 12), wrap="none",
                                bg="white", fg="#222", relief="flat",
                                padx=8, pady=6)
        log_vs = ttk.Scrollbar(log_frame, orient="vertical", command=self.log_text.yview)
        log_hs = ttk.Scrollbar(log_frame, orient="horizontal", command=self.log_text.xview)
        self.log_text.configure(yscrollcommand=log_vs.set, xscrollcommand=log_hs.set)
        self.log_text.grid(row=0, column=0, sticky="nsew")
        log_vs.grid(row=0, column=1, sticky="ns")
        log_hs.grid(row=1, column=0, sticky="ew")
        self.log_text.insert("end", "（尚无日志）")
        self.log_text.configure(state="disabled")

    # --- font zoom ---------------------------------------------------------
    def _apply_fonts(self) -> None:
        size = max(9, int(round(11 * self.font_scale)))
        row_h = max(24, int(round(26 * self.font_scale)))
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure("Treeview", font=(FONT, size), rowheight=row_h,
                        background="white", fieldbackground="white")
        style.configure("Treeview.Heading", font=(FONT, size, "bold"))
        if getattr(self, "log_text", None) is not None:
            self.log_text.configure(font=(FONT_MONO, size))
        if getattr(self, "zoom_label", None) is not None:
            self.zoom_label.configure(text=f"字体 {int(round(self.font_scale * 100))}%")

    def _zoom(self, delta: float) -> None:
        self.font_scale = min(2.0, max(0.8, round(self.font_scale + delta, 2)))
        self._apply_fonts()

    # ------------------------------------------------------------- 结果页
    def _build_results_tab(self, parent) -> None:
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(0, weight=0)
        top = ctk.CTkFrame(parent, corner_radius=10)
        top.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        self.result_title = ctk.CTkLabel(top, text="暂无计算体系", font=(FONT, 14, "bold"))
        self.result_title.pack(side="left", padx=12, pady=10)
        ctk.CTkButton(top, text="导出 Excel / CSV", width=150,
                      command=self._export_results).pack(side="right", padx=8, pady=8)
        ctk.CTkButton(top, text="刷新分析", width=100,
                      command=lambda: self._refresh_results(force=True)).pack(
            side="right", padx=4, pady=8)
        ctk.CTkButton(top, text="打开导出目录", width=120,
                      command=self._open_exports_dir).pack(side="right", padx=4, pady=8)
        ctk.CTkButton(top, text="从日志恢复结果", width=120, fg_color="#5A6B7B",
                      command=self._recover_from_logs).pack(side="right", padx=4, pady=8)
        self.batch_combo = ctk.CTkComboBox(
            top, values=["（无批次）"], width=330, command=self._on_batch_combo)
        self.batch_combo.pack(side="right", padx=6, pady=8)
        self.batch_combo.set("（无批次）")
        ctk.CTkLabel(top, text="批次：", font=(FONT, 12)).pack(side="right")

        self.result_tabs = ctk.CTkTabview(parent, corner_radius=8)
        self.result_tabs.grid(row=1, column=0, sticky="nsew")
        parent.grid_rowconfigure(1, weight=1)
        for name in ("Li⁺ 结合能", "氧化还原", "前线轨道", "热化学"):
            self.result_tabs.add(name)
        b = self.result_tabs.tab("Li⁺ 结合能")
        b.grid_rowconfigure(0, weight=1); b.grid_columnconfigure(0, weight=1)
        self.binding_tree = self._make_tree(
            b, ("complex", "site", "level", "eb", "bsse", "eb_bsse", "eb_ev", "vertical"),
            ("配合物", "位点", "理论级别", "Eb(未校正) kcal/mol", "BSSE kcal/mol",
             "Eb(BSSE) kcal/mol", "Eb(BSSE) eV", "垂直CP kcal/mol"),
            widths=(200, 60, 240, 140, 110, 140, 120, 130))
        r = self.result_tabs.tab("氧化还原")
        r.grid_rowconfigure(0, weight=1); r.grid_columnconfigure(0, weight=1)
        self.redox_tree = self._make_tree(
            r, ("species", "state", "dg_gas", "dg_sol", "dg_sol_kcal", "potential"),
            ("物种", "状态", "ΔG_gas / eV", "ΔG_sol / eV", "ΔG_sol kcal/mol", "E / V (相对)"),
            widths=(240, 70, 120, 120, 140, 120))
        f = self.result_tabs.tab("前线轨道")
        f.grid_rowconfigure(0, weight=1); f.grid_columnconfigure(0, weight=1)
        self.orb_tree = self._make_tree(
            f, ("study", "role", "homo", "lumo", "gap"),
            ("研究", "角色", "HOMO / eV", "LUMO / eV", "Gap / eV"),
            widths=(240, 140, 110, 110, 100))
        t = self.result_tabs.tab("热化学")
        t.grid_rowconfigure(0, weight=1); t.grid_columnconfigure(0, weight=1)
        self.thermo_tree = self._make_tree(
            t, ("study", "role", "zpe", "h", "g", "s", "shermo"),
            ("研究", "角色", "ZPE / Ha", "H(298) / Ha", "G(298) / Ha",
             "S / cal·mol⁻¹K⁻¹", "Shermo G / Ha"),
            widths=(240, 140, 120, 140, 140, 140, 130))

    # ------------------------------------------------------------- 设置页
    def _build_settings_tab(self, parent) -> None:
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(0, weight=1)
        scroll = ctk.CTkScrollableFrame(parent, corner_radius=10)
        scroll.grid(row=0, column=0, sticky="nsew")
        self.settings_scroll = scroll

        ctk.CTkLabel(scroll, text="软件路径（可自动检测，也可手动指定）",
                     font=(FONT, 14, "bold")).grid(row=0, column=0, columnspan=3,
                                                    sticky="w", padx=12, pady=(12, 4))
        self.path_entries: dict[str, ctk.CTkEntry] = {}
        for i, (key, label) in enumerate((
                ("g16_exe", "Gaussian 16 (g16.exe)"),
                ("gview_exe", "GaussView 6 (gview.exe)"),
                ("shermo_exe", "Shermo (shermo.exe，可选)"))):
            ctk.CTkLabel(scroll, text=label, font=(FONT, 11)).grid(
                row=1 + i, column=0, sticky="w", padx=12, pady=4)
            entry = ctk.CTkEntry(scroll, width=560)
            entry.insert(0, str(self.cfg.get(f"__path__.{key}") or
                               self.cfg.data.get(key) or ""))
            entry.grid(row=1 + i, column=1, sticky="ew", padx=8, pady=4)
            self.path_entries[key] = entry
            ctk.CTkButton(scroll, text="浏览", width=60,
                          command=lambda k=key: self._browse_exe(k)).grid(
                row=1 + i, column=2, padx=(0, 12))
        scroll.grid_columnconfigure(1, weight=1)
        ctk.CTkButton(scroll, text="自动检测软件", height=36,
                      command=self._auto_detect).grid(row=4, column=1, sticky="w",
                                                      padx=8, pady=6)
        ctk.CTkButton(scroll, text="保存设置", height=36,
                      command=self._save_settings).grid(row=4, column=2, sticky="e",
                                                        padx=(0, 12), pady=6)

        ctk.CTkLabel(scroll, text="计算资源（硬上限由 config.json 控制）",
                     font=(FONT, 14, "bold")).grid(row=5, column=0, columnspan=3,
                                                    sticky="w", padx=12, pady=(14, 4))
        self.set_nproc = ctk.CTkEntry(scroll, width=100)
        self.set_nproc.insert(0, str(self.cfg.resources["nproc"]))
        ctk.CTkLabel(scroll, text="并行核数 %nproc", font=(FONT, 11)).grid(
            row=6, column=0, sticky="w", padx=12)
        self.set_nproc.grid(row=6, column=1, sticky="w", padx=8)
        self.set_mem = ctk.CTkEntry(scroll, width=100)
        self.set_mem.insert(0, str(self.cfg.resources["mem"]))
        ctk.CTkLabel(scroll, text="内存 %mem", font=(FONT, 11)).grid(
            row=7, column=0, sticky="w", padx=12)
        self.set_mem.grid(row=7, column=1, sticky="w", padx=8)

        ctk.CTkLabel(scroll, text="Python 依赖自检",
                     font=(FONT, 14, "bold")).grid(row=8, column=0, columnspan=3,
                                                    sticky="w", padx=12, pady=(14, 4))
        self.deps_label = ctk.CTkLabel(scroll, text="检查中…", font=(FONT, 11),
                                       justify="left", anchor="w", wraplength=760)
        self.deps_label.grid(row=9, column=0, columnspan=3, sticky="w", padx=12)
        ctk.CTkButton(scroll, text="重新检查", width=110,
                      command=self._check_deps).grid(row=10, column=0, sticky="w",
                                                      padx=12, pady=6)
        self.install_deps_btn = ctk.CTkButton(
            scroll, text="自动安装缺失依赖", width=170, command=self._install_deps)
        self.install_deps_btn.grid(row=10, column=1, sticky="w", padx=8, pady=6)
        self.deps_log = ctk.CTkTextbox(scroll, height=130, font=(FONT_MONO, 10))
        self.deps_log.grid(row=11, column=0, columnspan=3, sticky="ew",
                           padx=12, pady=(0, 12))
        scroll.grid_columnconfigure(1, weight=1)

        # --- DeepSeek LLM (molecule generation) ---
        ctk.CTkLabel(scroll, text="DeepSeek 大模型（AI 生成分子 / 计算链）",
                     font=(FONT, 14, "bold")).grid(
            row=12, column=0, columnspan=3, sticky="w", padx=12, pady=(14, 4))
        self.llm_key = ctk.CTkEntry(scroll, width=560, show="*")
        self.llm_key.insert(0, str(self.cfg.get("llm.api_key") or ""))
        ctk.CTkLabel(scroll, text="API Key", font=(FONT, 11)).grid(
            row=13, column=0, sticky="w", padx=12)
        self.llm_key.grid(row=13, column=1, sticky="ew", padx=8, pady=3)
        raw_models = list(self.cfg.get("llm.model_options") or llm.MODEL_CHOICES)
        model_ids: list[str] = []
        for item in raw_models:
            model = llm.resolve_model(item, self.cfg)
            if model not in model_ids:
                model_ids.append(model)
        cur_model = llm.resolve_model(self.cfg.get("llm.model"), self.cfg)
        if cur_model not in model_ids:
            model_ids.insert(0, cur_model)
        model_labels = [llm.model_label(m) for m in model_ids]
        ctk.CTkLabel(scroll, text="模型", font=(FONT, 11)).grid(
            row=14, column=0, sticky="w", padx=12)
        self.llm_model = ctk.CTkOptionMenu(scroll, values=model_labels, width=240)
        self.llm_model.set(llm.model_label(cur_model))
        self.llm_model.grid(row=14, column=1, sticky="w", padx=8, pady=3)
        ctk.CTkLabel(scroll, text="thinking（推理模式，更慢但更准）",
                     font=(FONT, 11)).grid(row=14, column=2, sticky="w", padx=8)
        self.llm_thinking = ctk.CTkSwitch(scroll, text="开启")
        if str(self.cfg.get("llm.thinking") or "disabled").lower() == "enabled":
            self.llm_thinking.select()
        self.llm_thinking.grid(row=14, column=3, sticky="w", padx=8, pady=3)
        self.llm_base = ctk.CTkEntry(scroll, width=420)
        self.llm_base.insert(0, str(self.cfg.get("llm.base_url") or "https://api.deepseek.com"))
        ctk.CTkLabel(scroll, text="API 地址", font=(FONT, 11)).grid(
            row=15, column=0, sticky="w", padx=12)
        self.llm_base.grid(row=15, column=1, columnspan=2, sticky="w", padx=8, pady=3)
        ctk.CTkButton(scroll, text="测试 DeepSeek 连接", width=170,
                      command=self._test_llm).grid(row=16, column=0, sticky="w",
                                                   padx=12, pady=6)
        self.llm_status = ctk.CTkLabel(scroll, text="", font=(FONT, 11), anchor="w",
                                       justify="left", wraplength=720, text_color="#555")
        self.llm_status.grid(row=16, column=1, columnspan=3, sticky="w", padx=8)
        ctk.CTkLabel(scroll, font=(FONT, 11), text_color="#666", anchor="w",
                     text=("大模型固定为 DeepSeek（OpenAI 兼容接口）；"
                           "V4 Flash = 快 / 便宜（默认），V4 Pro = 最强 / 较慢。"
                           "AI 生成分子页与「AI 生成计算链」窗口里都可以随时切换。"
                           "API Key 只保存在本机 config.json，申请：https://platform.deepseek.com 。"),
                     ).grid(row=17, column=0, columnspan=4, sticky="w", padx=12, pady=(0, 6))

        ctk.CTkLabel(
            scroll, font=(FONT, 11), text_color="#666", justify="left", anchor="w",
            text=("说明：Gaussian 16 / GaussView 6 为商业软件，需自行安装并持有许可；"
                  "本程序只负责调用、解析与汇总。Shermo 为可选的热化学后处理工具（免费）。"
                  "所有任务目录、输入输出与日志均保留在数据目录中用于审计。"),
        ).grid(row=18, column=0, columnspan=4, sticky="w", padx=12, pady=(4, 14))

    # ------------------------------------------------------------- helpers
    def _make_tree(self, parent, cols, heads, widths=None):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        size = max(9, int(round(11 * self.font_scale)))
        style.configure("Treeview", font=(FONT, size),
                        rowheight=max(24, int(round(26 * self.font_scale))),
                        background="white", fieldbackground="white")
        style.configure("Treeview.Heading", font=(FONT, size, "bold"))
        frame = tk.Frame(parent, bg="white")
        frame.pack(fill="both", expand=True)
        tree = ttk.Treeview(frame, columns=cols, show="headings")
        for i, (c, h) in enumerate(zip(cols, heads)):
            tree.heading(c, text=h)
            tree.column(c, width=(widths[i] if widths else 110), anchor="w")
        vs = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        hs = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vs.set, xscrollcommand=hs.set)
        tree.grid(row=0, column=0, sticky="nsew")
        vs.grid(row=0, column=1, sticky="ns")
        hs.grid(row=1, column=0, sticky="ew")
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        return tree

    def _fill_tree(self, tree, rows: list[tuple]) -> None:
        tree.delete(*tree.get_children())
        for row in rows:
            tree.insert("", "end", values=row)

    # --------------------------------------------------------- auto detect
    def _startup_detect(self) -> None:
        self._check_deps()

        def work():
            try:
                found = locate.detect_all({
                    "g16_exe": self.cfg.data.get("g16_exe"),
                    "gview_exe": self.cfg.data.get("gview_exe"),
                    "shermo_exe": self.cfg.data.get("shermo_exe"),
                })
            except Exception:
                found = None
            self.after(0, lambda: self._apply_detection(found, save=True))

        threading.Thread(target=work, daemon=True).start()

    def _auto_detect(self) -> None:
        self._auto_detect_clicked()

    def _auto_detect_clicked(self) -> None:
        self.status_badge.configure(text="● 检测中…", text_color="#FFD966")

        def work():
            try:
                found = locate.detect_all()
                self.after(0, lambda: self._apply_detection(found, save=True))
            except Exception:
                self.after(0, lambda: self._apply_detection(None, save=False))

        threading.Thread(target=work, daemon=True).start()

    def _apply_detection(self, found, save: bool = False) -> None:
        if not found:
            self.status_badge.configure(text="● 未检测到 Gaussian", text_color="#FF6B6B")
            return
        g = found["gaussian"]["best"]
        v = found["gaussview"]["best"]
        s = found["shermo"]["best"]
        for key, item in (("g16_exe", g), ("gview_exe", v), ("shermo_exe", s)):
            if item and key in self.path_entries:
                self.path_entries[key].delete(0, "end")
                self.path_entries[key].insert(0, item["exe"])
        if g:
            gv = f" · GaussView {v['version']}" if v and v.get("version") else ""
            self.status_badge.configure(
                text=f"● Gaussian {g.get('version') or ''} 已就绪{gv}", text_color="#7CE38B")
        else:
            self.status_badge.configure(text="● 未检测到 Gaussian", text_color="#FF6B6B")
        if save:
            self._save_settings(verbose=False)
        self._check_deps()

    def _browse_exe(self, key: str) -> None:
        path = filedialog.askopenfilename(title="选择可执行文件",
                                          filetypes=[("可执行文件", "*.exe"), ("全部", "*.*")])
        if path:
            self.path_entries[key].delete(0, "end")
            self.path_entries[key].insert(0, path)

    def _save_settings(self, verbose: bool = True) -> None:
        try:
            for key, entry in self.path_entries.items():
                self.cfg.set(key, entry.get().strip() or None)
            nproc = max(1, int(self.set_nproc.get().strip() or "1"))
            self.cfg.set("resources.nproc", min(nproc, int(self.cfg.resources.get("max_nproc", 18))))
            self.cfg.set("resources.mem", self.set_mem.get().strip() or "4GB")
            if getattr(self, "llm_key", None) is not None:
                self.cfg.set("llm.api_key", self.llm_key.get().strip())
                self.cfg.set("llm.model",
                             llm.model_from_label(self.llm_model.get(), self.cfg))
                self.cfg.set("llm.thinking",
                             "enabled" if self.llm_thinking.get() else "disabled")
                self.cfg.set("llm.base_url",
                             self.llm_base.get().strip() or "https://api.deepseek.com")
                if getattr(self, "ai_model", None) is not None:
                    self.ai_model.set(llm.model_label(self.cfg.get("llm.model")))
            self.cfg.save()
            if verbose:
                messagebox.showinfo("已保存", "设置已保存。\n新任务将使用这些路径与资源。")
            self._check_deps()
        except Exception as exc:
            messagebox.showerror("保存失败", str(exc))

    def _test_llm(self) -> None:
        if getattr(self, "llm_key", None) is not None:
            self.cfg.set("llm.api_key", self.llm_key.get().strip())
            self.cfg.set("llm.model",
                         llm.model_from_label(self.llm_model.get(), self.cfg))
            self.cfg.set("llm.thinking",
                         "enabled" if self.llm_thinking.get() else "disabled")
            self.cfg.set("llm.base_url",
                         self.llm_base.get().strip() or "https://api.deepseek.com")
            if getattr(self, "ai_model", None) is not None:
                self.ai_model.set(llm.model_label(self.cfg.get("llm.model")))
        self.llm_status.configure(text="正在测试 DeepSeek 连接…", text_color="#555")

        def work():
            try:
                from .llm import chat
                res = chat(self.cfg, [
                    {"role": "user",
                     "content": '只输出 JSON：{"ok": true, "msg": "pong"}'}])
                self.after(0, lambda: self.llm_status.configure(
                    text=f"✓ 连接成功（模型 {res.get('model')}）", text_color="#2E7D32"))
            except Exception as exc:
                msg = f"✗ {type(exc).__name__}: {exc}"
                self.after(0, lambda: self.llm_status.configure(
                    text=msg[:400], text_color="#B22222"))

        threading.Thread(target=work, daemon=True).start()

    # -------------------------------------------------------- dependencies
    def _pip_available(self) -> bool:
        return not getattr(sys, "frozen", False)

    def _check_deps(self) -> None:
        mods = {
            "numpy": "numpy", "scipy": "scipy", "rdkit": "rdkit",
            "cclib": "cclib", "openpyxl": "openpyxl", "customtkinter": "customtkinter",
            "PIL": "Pillow",
        }
        lines, missing = [], []
        for mod, label in mods.items():
            try:
                m = __import__(mod)
                ver = getattr(m, "__version__", "") or getattr(m, "VERSION", "")
                lines.append(f"✓ {label} {ver}")
            except Exception:
                lines.append(f"✗ {label}（缺失）")
                missing.append(label)
        mode = "打包 exe（依赖已内置）" if not self._pip_available() else "源码运行"
        text = f"运行模式：{mode}\n" + "    ".join(lines)
        if missing and not self._pip_available():
            text += "\n（exe 内依赖已打包，若仍缺失请重新构建）"
        self.deps_label.configure(text=text)
        self.install_deps_btn.configure(
            state="normal" if (missing and self._pip_available()) else "disabled")

    def _install_deps(self) -> None:
        if not self._pip_available():
            return
        self.deps_log.delete("1.0", "end")
        self.deps_log.insert("end", "开始安装依赖…\n")

        def work():
            proc = subprocess.run(
                [sys.executable, "-m", "pip", "install", "rdkit", "cclib", "openpyxl",
                 "numpy", "scipy", "customtkinter", "pillow"],
                capture_output=True, text=True, encoding="utf-8", errors="replace")
            out = (proc.stdout or "") + "\n" + (proc.stderr or "")
            self.after(0, lambda: (self.deps_log.insert("end", out[-4000:]),
                                   self.deps_log.see("end"),
                                   self._check_deps()))

        threading.Thread(target=work, daemon=True).start()

    # ---------------------------------------------------------- input read
    def _default_recipe_key(self) -> str:
        recipes = pipeline_recipe_names(self.cfg)
        names = list(recipes.values())
        return names[0] if names else ""

    def _recipe_key(self) -> str:
        label = self.recipe_var.get()
        for k, v in self._recipes.items():
            if v == label:
                return k
        return next(iter(self._recipes))

    def _on_recipe_change(self, _value=None) -> None:
        key = self._recipe_key()
        try:
            from .recipes import get_recipe
            recipe = get_recipe(self.cfg, key)
            self.recipe_desc.configure(text=recipe.get("description", ""))
            self._update_stage_overview()
        except Exception:
            pass

    def _update_stage_overview(self) -> None:
        key = self._recipe_key() if getattr(self, "_recipes", None) else "opt_freq"
        try:
            from .recipes import get_recipe
            recipe = get_recipe(self.cfg, key)
        except Exception:
            return
        lines = []
        if getattr(self, "preopt_switch", None) is not None and self.preopt_switch.get():
            lines.append("# 预优化   HF/3-21G  opt")
        for st in recipe.get("core_stages", []):
            route = st.get("route", "")
            kind = st.get("kind")
            if kind == "optfreq":
                lines.append(f"# {st['name']:<8s} {route}  opt freq")
            elif kind == "sp":
                solvent = None
                if getattr(self, "solvent_entry", None) is not None:
                    solvent = self.solvent_entry.get()
                solv = f"  scrf=(pcm,solvent={st.get('solvent') or solvent})" \
                    if st.get("phase") == "solvent" else ""
                lines.append(f"# {st['name']:<8s} {route}{solv}")
        li = recipe.get("li") or {}
        li_on = (getattr(self, "li_switch", None) is not None
                 and self.li_switch.get())
        if li.get("on") and li_on:
            lines.append("")
            lines.append(f"# Li⁺ 配合物：额外 opt+freq → SP×2，并对 "
                         f"{'、'.join(li.get('bsse_levels') or [])} 做 counterpoise=2 BSSE")
        rd = recipe.get("redox") or {}
        redox_on = (getattr(self, "redox_switch", None) is not None
                    and self.redox_switch.get())
        if rd.get("on") and redox_on:
            lines.append("# 氧化还原：cation / anion 各自 opt+freq → SP×2")
        self.stages_box.configure(state="normal")
        self.stages_box.delete("1.0", "end")
        self.stages_box.insert("end", "\n".join(lines))
        self.stages_box.configure(state="disabled")

    # -------------------------------------------------------------- run
    def _start_calc(self) -> None:
        if self._busy:
            return
        self._stopped_by_user = False
        # avoid two workers competing for the same queue
        try:
            from .daemon import daemon_status, stop_daemon
            if daemon_status(self.cfg).get("running"):
                if not messagebox.askyesno(
                        "后台 worker 正在运行",
                        "检测到后台计算 worker 正在运行（可能是命令行提交的任务）。\n"
                        "为避免两个 worker 抢任务，请先停止它。是否现在停止并继续？"):
                    return
                stop_daemon(self.cfg, timeout=8)
        except Exception:
            pass
        mode = self.mode_var.get()
        smiles = self.smiles_entry.get().strip()
        name = self.name_entry.get().strip()
        options = {
            "li": bool(self.li_switch.get()),
            "li_elements": (["O"] if self.li_o.get() else []) +
                           (["F"] if self.li_f.get() else []),
            "redox": bool(self.redox_switch.get()),
            "redox_states": ["cation", "anion"],
            "redox_targets": (["neutral"] if self.redox_neutral.get() else []) +
                             (["li"] if self.redox_li.get() else []),
            "preopt": bool(self.preopt_switch.get()),
            "solvent": self.solvent_entry.get().strip() or None,
        }
        charge = self.charge_entry.get().strip()
        mult = self.mult_entry.get().strip()

        if mode.startswith("AI") and not self._ai_molecules:
            messagebox.showwarning(
                "缺少分子", "请先点「① 生成分子并导出 gjf」用 DeepSeek 生成结构。")
            return
        if mode.startswith("AI") and len(self._ai_molecules) > 1 and \
                (options["li"] or options["redox"]):
            if not messagebox.askyesno(
                    "任务量提醒",
                    f"已生成 {len(self._ai_molecules)} 个分子，且启用了 Li⁺ 配位 / 氧化还原。\n"
                    "每个分子会展开为 16–22 个研究（数百个 Gaussian 任务），耗时较长。\n"
                    "如果只想比较 F 数系列，建议先关闭 Li⁺ 与氧化还原。\n\n仍要继续吗？"):
                return
        if mode.startswith("GJF") and not self._gjf_files:
            messagebox.showwarning("缺少输入", "请先选择至少一个 .gjf 文件（可多选）。")
            return
        if mode == "SMILES" and not smiles:
            messagebox.showwarning("缺少输入", "请输入分子 SMILES（例如 CCO）。")
            return
        try:
            nproc = int(self.nproc_entry.get().strip() or self.cfg.resources["nproc"])
            mem = self.mem_entry.get().strip() or self.cfg.resources["mem"]
        except ValueError:
            messagebox.showerror("参数错误", "并行核数必须是整数。")
            return
        self.cfg.set("resources.nproc", nproc)
        self.cfg.set("resources.mem", mem)
        self.cfg.save()

        recipe_key = self._recipe_key()
        self._busy = True
        self.run_btn.configure(state="disabled")
        self.stop_btn.configure(state="normal")
        self.progress.configure(mode="indeterminate")
        self.progress.start()
        self.run_status.configure(text="正在生成结构/输入并建立任务…")
        self.tabs.set("任务 / 日志")

        def work():
            try:
                store = Store(self.cfg.db_path)
                if mode.startswith("AI"):
                    res = pipeline.create_ai_systems(
                        store, self.cfg, self._ai_molecules,
                        recipe_name=recipe_key, options=options,
                        name_prefix=name or "AI", submit=True)
                    store.close()
                    self.after(0, lambda: self._on_ai_ran(res))
                elif mode.startswith("GJF"):
                    res = pipeline.create_batch_from_gjfs(
                        store, self.cfg, list(self._gjf_files),
                        name_prefix=name or None, submit=True)
                    store.close()
                    self.after(0, lambda: self._on_batch_created(res))
                else:
                    res = pipeline.create_recipe_system(
                        store, self.cfg, name=name or _auto_name(smiles, ""),
                        smiles=smiles, recipe_name=recipe_key, options=options,
                        charge=int(charge) if charge else None,
                        multiplicity=int(mult) if mult else None)
                    pipeline.submit_batch(store, self.cfg, res["batch_id"])
                    store.close()
                    self.after(0, lambda: self._on_batch_created(res))
            except Exception:
                msg = traceback.format_exc()
                self.after(0, lambda: self._on_run_error(msg))

        threading.Thread(target=work, daemon=True).start()

    def _on_ai_ran(self, res: dict) -> None:
        batches = res.get("batches") or []
        if not batches:
            self._on_run_error("没有创建任何批次（AI 分子列表为空）")
            return
        self.batch_ids = [b["batch_id"] for b in batches]
        self.batch_id = self.batch_ids[-1]
        self.run_status.configure(
            text=(f"已为 {res.get('n_molecules')} 个 AI 分子创建 "
                  f"{len(batches)} 个体系并提交排队，开始计算…"))
        self.runner = CalcWorker(self.cfg, self)
        self.runner.start()
        self._refresh_batch_choices()
        self._refresh_tasks()
        self._refresh_results(force=True)

    def _on_batch_created(self, res: dict) -> None:
        self.batch_ids = [res["batch_id"]]
        self.batch_id = res["batch_id"]
        self.run_status.configure(
            text=f"已创建 {res.get('n_studies', '?')} 个研究（{res.get('batch_id')}），开始计算…")
        self.runner = CalcWorker(self.cfg, self)
        self.runner.start()
        self._refresh_batch_choices()
        self._refresh_tasks()
        self._refresh_results(force=True)

    def _on_run_finished(self) -> None:
        self._busy = False
        self.progress.stop()
        self.progress.configure(mode="determinate")
        self.run_btn.configure(state="normal")
        self.stop_btn.configure(state="disabled")
        if self._stopped_by_user:
            self._stopped_by_user = False
            self.progress.set(0)
            self.run_status.configure(text="已停止：当前 Gaussian 进程已终止，排队任务已取消。")
            self._refresh_batch_choices()
            self._refresh_tasks()
            self._refresh_results(force=True)
            return
        self.progress.set(1.0)
        self.run_status.configure(text="计算结束，正在汇总结果…")
        self._refresh_batch_choices()
        self._refresh_tasks()
        self._refresh_results(force=True)
        self.run_status.configure(text="计算结束。可在“结果”页导出 Excel / CSV。")

    def _on_run_error(self, msg: str) -> None:
        self._busy = False
        self.progress.stop()
        self.progress.configure(mode="determinate")
        self.progress.set(0)
        self.run_btn.configure(state="normal")
        self.stop_btn.configure(state="disabled")
        self.run_status.configure(text="运行失败（见日志）")
        self.log_text.configure(state="normal")
        self.log_text.insert("end", "\n" + msg + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")
        messagebox.showerror("运行失败", msg[-1500:])

    def _stop_calc(self) -> None:
        """Force-stop: kill the running g16.exe tree and cancel the queue."""
        if not self._busy and not self.runner:
            return
        if not messagebox.askyesno(
                "停止计算",
                "将立即终止当前正在运行的 Gaussian 进程（g16.exe 及其 link 子进程），\n"
                "并取消本批次所有排队任务。已经完成的任务结果保留。\n\n确定停止吗？"):
            return
        self._stopped_by_user = True
        self.stop_btn.configure(state="disabled")
        self.run_status.configure(text="正在终止 Gaussian 进程并取消排队任务…")
        # 1) stop the in-process worker loop (it exits after the killed job returns)
        if self.runner:
            self.runner.stop()
        # 2) stop a background daemon if one is running (it may own the jobs)
        daemon_running = False
        try:
            from .daemon import daemon_status, stop_daemon
            if daemon_status(self.cfg).get("running"):
                daemon_running = True
                stop_daemon(self.cfg, timeout=8)
        except Exception:
            pass
        # 3) kill running g16 process trees + cancel pending jobs in the DB
        try:
            store = Store(self.cfg.db_path)
            batches = list(self.batch_ids) or ([self.batch_id] if self.batch_id else [])
            total_run = total_pending = 0
            if batches and not daemon_running:
                for bid in batches:
                    res = pipeline.cancel_jobs(store, self.cfg, batch_ids=[bid])
                    total_run += len(res.get("running_cancelled") or [])
                    total_pending += len(res.get("pending_cancelled") or [])
            else:
                res = pipeline.cancel_jobs(store, self.cfg, all_jobs=True)
                total_run += len(res.get("running_cancelled") or [])
                total_pending += len(res.get("pending_cancelled") or [])
            store.close()
            self.run_status.configure(
                text=(f"已停止：强制终止 {total_run} 个运行中任务，"
                      f"取消 {total_pending} 个排队任务。"))
            self.log_text.configure(state="normal")
            self.log_text.insert(
                "end", f"\n[用户停止] 运行中={total_run} 排队取消={total_pending}\n")
            self.log_text.see("end")
            self.log_text.configure(state="disabled")
        except Exception:
            msg = traceback.format_exc()
            self.run_status.configure(text="停止时出错，详见日志。")
            self.log_text.configure(state="normal")
            self.log_text.insert("end", "\n" + msg + "\n")
            self.log_text.see("end")
            self.log_text.configure(state="disabled")
        # 4) reset the UI immediately (worker thread will exit shortly)
        self._busy = False
        self.progress.stop()
        self.progress.configure(mode="determinate")
        self.progress.set(0)
        self.run_btn.configure(state="normal")
        self._refresh_tasks()
        self._refresh_results(force=True)

    # -------------------------------------------------------- refresh loop
    def _refresh_loop(self) -> None:
        try:
            self._refresh_tasks()
            self._refresh_log()
            if self.batch_id:
                self._refresh_results(force=False)
        except Exception:
            pass
        self.after(2000, self._refresh_loop)

    def _refresh_tasks(self) -> None:
        if not self.batch_id:
            self._fill_tree(self.task_tree, [])
            return
        try:
            store = Store(self.cfg.db_path)
            studies = store.list_studies(batch_id=self.batch_id, limit=500)
            hide_cancelled = (getattr(self, "hide_cancelled", None) is not None
                              and self.hide_cancelled.get())
            rows = []
            self._study_names = {}
            for s in studies:
                if hide_cancelled and s["status"] == "cancelled":
                    continue
                jobs = store.study_jobs(s["id"])
                current = next((j for j in jobs if j["status"] in ("running", "pending")),
                               jobs[-1] if jobs else None)
                last_res = store.get_result(current["id"]) if current else None
                parsed = (last_res or {}).get("parsed") or {}
                val = ((last_res or {}).get("validation") or {}).get("overall")
                meta = s.get("meta") or {}
                rows.append((
                    s["name"], _role_label(meta.get("role") or ""),
                    s["status"], current["stage"] if current else "—",
                    val or "—",
                    _fmt(parsed.get("final_energy_hartree"), 6),
                    _fmt((parsed.get("orbitals") or {}).get("homo_ev"), 3),
                    _fmt((parsed.get("freq") or {}).get("n_imag")),
                ))
                self._study_names[s["id"]] = s["name"]
            self._fill_tree(self.task_tree, rows)
            store.close()
        except Exception:
            pass

    def _purge_cancelled(self) -> None:
        """Delete cancelled studies from the task list (run files kept)."""
        if not messagebox.askyesno(
                "清理取消记录",
                "将从任务列表删除所有「已取消」的研究（数据库记录）。\n"
                "对应的 runs/ 目录文件默认保留（如需彻底删除请手动删除目录）。\n\n"
                "确定清理吗？"):
            return
        try:
            store = Store(self.cfg.db_path)
            batch_ids = [self.batch_id] if self.batch_id else None
            res = pipeline.purge_cancelled(store, self.cfg, batch_ids=batch_ids)
            if self.batch_id and not store.list_studies(batch_id=self.batch_id):
                self.batch_id = None
                self.result_title.configure(text="暂无计算体系")
            store.close()
            self._refresh_batch_choices()
            self._refresh_tasks()
            self._refresh_results(force=True)
            self.run_status.configure(
                text=f"已清理 {res.get('count', 0)} 条取消记录（文件保留在 runs/）。")
        except Exception as exc:
            messagebox.showerror("清理失败", f"{type(exc).__name__}: {exc}")

    def _refresh_log(self) -> None:
        if not self.batch_id:
            return
        try:
            store = Store(self.cfg.db_path)
            studies = store.list_studies(batch_id=self.batch_id, limit=500)
            ids = {s["id"] for s in studies}
            events = [e for e in store.list_events(limit=120) if e.get("study_id") in ids]
            events.reverse()
            text = "\n".join(
                f"[{e['ts'][11:19]}] {e.get('level','info').upper():7s} {e['message']}"
                for e in events[-80:])
            # running job output tail
            running = [j for s in studies for j in store.study_jobs(s["id"])
                       if j["status"] == "running"]
            if running:
                job = running[0]
                out = Path(job.get("output_path") or "")
                if out.exists():
                    try:
                        tail = out.read_text(encoding="utf-8", errors="replace").splitlines()[-8:]
                        text += "\n--- 当前任务输出尾部 ---\n" + "\n".join(tail)
                    except OSError:
                        pass
            self.log_text.configure(state="normal")
            self.log_text.delete("1.0", "end")
            self.log_text.insert("end", text or "（暂无日志）")
            self.log_text.see("end")
            self.log_text.configure(state="disabled")
            store.close()
        except Exception:
            pass

    def _refresh_results(self, force: bool = False) -> None:
        if not self.batch_id or self._refresh_busy:
            return
        if not force and time.time() - getattr(self, "_last_result", 0) < 8:
            return
        self._refresh_busy = True

        def work():
            try:
                store = Store(self.cfg.db_path)
                data = analysis.system_summary(store, self.cfg, self.batch_id)
                store.close()
            except Exception:
                data = None
            self.after(0, lambda: self._apply_results(data))

        threading.Thread(target=work, daemon=True).start()
        self._last_result = time.time()

    def _apply_results(self, data) -> None:
        self._refresh_busy = False
        if not data or not data.get("ok"):
            return
        self._analysis = data
        self.result_title.configure(
            text=f"体系：{data.get('system_name')} · {data.get('n_studies')} 个研究 · "
                 f"方案 {data.get('recipe')} · 溶剂 {data.get('solvent')}")
        rows = []
        for row in data["binding"]:
            rows.append((row["complex"], row["site"], row["level_label"],
                         _fmt(row["E_binding_kcal"], 2),
                         _fmt(row["BSSE_kcal"] if row["BSSE_kcal"] is not None else None, 2),
                         _fmt(row["E_binding_BSSE_kcal"], 2),
                         _fmt(row["E_binding_BSSE_ev"], 3),
                         _fmt(row["vertical_CP_corrected_kcal"], 2)))
        self._fill_tree(self.binding_tree, rows)
        rows = []
        for row in data["redox"]:
            rows.append((row["species"], row["state"],
                         _fmt(row["dG_gas_ev"], 3), _fmt(row["dG_sol_ev"], 3),
                         _fmt(row["dG_sol_kcal"], 1), _fmt(row["E_potential_v"], 3)))
        self._fill_tree(self.redox_tree, rows)
        rows = []
        for row in data["orbitals"]:
            rows.append((row["study"], _role_label(row["role"]),
                         _fmt(row["HOMO_eV"], 3), _fmt(row["LUMO_eV"], 3),
                         _fmt(row["gap_eV"], 3)))
        self._fill_tree(self.orb_tree, rows)
        rows = []
        for row in data["thermo"]:
            rows.append((row["study"], _role_label(row["role"]),
                         _fmt(row["ZPE_hartree"], 6), _fmt(row["H_298_hartree"], 6),
                         _fmt(row["G_298_hartree"], 6), _fmt(row["S_cal_mol_K"], 2),
                         _fmt(row["Shermo_G_hartree"], 6)))
        self._fill_tree(self.thermo_tree, rows)

    def _export_results(self) -> None:
        if not self.batch_id:
            messagebox.showwarning("尚未计算", "请先运行一个计算体系。")
            return
        self.run_status.configure(text="正在导出 Excel / CSV…")

        def work():
            try:
                store = Store(self.cfg.db_path)
                res = analysis.export_system(store, self.cfg, self.batch_id)
                store.close()
            except Exception:
                msg = traceback.format_exc()
                self.after(0, lambda: messagebox.showerror("导出失败", msg[-1500:]))
                return

            def done():
                self.run_status.configure(text="导出完成。")
                if messagebox.askyesno("导出完成",
                                       f"已导出到：\n{res['out_dir']}\n\n是否打开目录？"):
                    self._open_path(res["out_dir"])

            self.after(0, done)

        threading.Thread(target=work, daemon=True).start()

    # ------------------------------------------------------------ misc
    def _open_path(self, path: str) -> None:
        try:
            os.startfile(path)  # type: ignore[attr-defined]
        except Exception as exc:
            messagebox.showerror("打开失败", str(exc))

    def _gview_selected(self) -> None:
        """Open the selected study in GaussView for human review."""
        gview = (self.path_entries.get("gview_exe").get().strip()
                 if getattr(self, "path_entries", None) else "") or \
            str(self.cfg.data.get("gview_exe") or "")
        if not gview or not Path(gview).exists():
            messagebox.showwarning("未配置", "请先在“设置”里指定 gview.exe。")
            return
        if not self.batch_id:
            messagebox.showwarning("无任务", "当前没有可查看的批次。")
            return
        sel = self.task_tree.selection()
        if not sel:
            messagebox.showinfo("提示", "请先在任务列表中选中一个研究。")
            return
        name = self.task_tree.item(sel[0], "values")[0]
        try:
            store = Store(self.cfg.db_path)
            studies = store.list_studies(batch_id=self.batch_id, limit=500)
            study = next((s for s in studies if s["name"] == name), None)
            target = None
            if study:
                initial = self.cfg.work_root / study["id"] / "initial.gjf"
                if initial.exists():
                    target = initial
            store.close()
            if target is None:
                messagebox.showwarning("未找到结构", "该研究还没有可打开的结构文件。")
                return
            subprocess.Popen([gview, str(target)])
        except Exception as exc:
            messagebox.showerror("打开 GaussView 失败", str(exc))

    def _open_runs_dir(self) -> None:
        path = self.cfg.work_root
        if self.batch_id:
            path = self.cfg.work_root
        self._open_path(str(path))

    def _open_exports_dir(self) -> None:
        self._open_path(str(self.cfg.exports_dir))

    def _on_close(self) -> None:
        if self._busy and not messagebox.askyesno(
                "退出", "计算仍在进行，退出会中断当前任务。确定退出？"):
            return
        if self.runner:
            self.runner.stop()
        self.destroy()


# ---------------------------------------------------------------------------
def pipeline_recipe_names(cfg) -> dict[str, str]:
    from .recipes import list_recipes
    try:
        return list_recipes(cfg)
    except Exception:
        return {"opt_freq": "优化 + 频率"}


def _auto_name(smiles: str, file_path: str) -> str:
    import re
    if smiles:
        base = re.sub(r"[^0-9A-Za-z]+", "_", smiles).strip("_")
        return (base[:28] or "molecule")
    if file_path:
        return Path(file_path).stem
    return "molecule"


# ---------------------------------------------------------------------------
def selftest(outfile: str) -> int:
    from .apphome import bootstrap_home
    home = bootstrap_home()
    info: dict = {"ok": True, "frozen": bool(getattr(sys, "frozen", False)),
                  "python": sys.version, "home": str(home),
                  "executable": sys.executable}
    mods = {}
    for mod in ("numpy", "scipy", "rdkit", "cclib", "openpyxl", "customtkinter", "PIL"):
        try:
            m = __import__(mod)
            mods[mod] = getattr(m, "__version__", "ok")
        except Exception as exc:
            mods[mod] = f"MISSING: {exc}"
            info["ok"] = False
    info["modules"] = mods
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        info["rd_embed_atoms"] = mol.GetNumAtoms()
    except Exception as exc:
        info["rd_embed_error"] = str(exc)
        info["ok"] = False
    try:
        found = locate.detect_all()
        info["detected"] = {
            "gaussian": (found["gaussian"]["best"] or {}).get("exe"),
            "gaussview": (found["gaussview"]["best"] or {}).get("exe"),
            "shermo": (found["shermo"]["best"] or {}).get("exe"),
        }
    except Exception as exc:
        info["detect_error"] = str(exc)
    try:
        cfg = Config.load()
        store = Store(cfg.db_path)
        store.close()
        info["config"] = str(cfg.path)
        info["work_root"] = str(cfg.work_root)
    except Exception as exc:
        info["config_error"] = str(exc)
        info["ok"] = False
    try:
        Path(outfile).write_text(json.dumps(info, ensure_ascii=False, indent=2),
                                 encoding="utf-8")
    except OSError:
        return 1
    return 0 if info.get("ok") else 2


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--selftest" in argv:
        idx = argv.index("--selftest")
        out = argv[idx + 1] if idx + 1 < len(argv) else "gaussflow_selftest.json"
        return selftest(out)
    from .apphome import bootstrap_home
    bootstrap_home()
    app = App()
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
