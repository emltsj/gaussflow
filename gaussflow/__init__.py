"""gaussflow — Agent-driven Gaussian 16 / GaussView calculation pipeline.

Deterministic core: structure preparation -> GJF generation -> job queue ->
Gaussian execution -> parsing -> scientific validation -> aggregation/export.
The LLM/Agent is only the orchestrator on top of these controlled tools.

Key design rule: the Agent never edits raw shell commands for Gaussian.
It calls prepare/submit/status/result/diagnose/export, and everything
numeric is parsed and validated deterministically.
"""

__version__ = "0.2.0"
