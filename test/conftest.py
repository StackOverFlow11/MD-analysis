"""Shared pytest fixtures and helpers for MD Analysis tests."""

from __future__ import annotations

from md_analysis.water.WaterAnalysis._common import _parse_abc_from_md_inp as parse_abc_from_md_inp

__all__ = ["parse_abc_from_md_inp"]
