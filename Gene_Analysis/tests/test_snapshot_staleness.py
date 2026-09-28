"""Tests for evidence snapshot age and staleness logic."""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

from Gene_Analysis.evidence import is_snapshot_stale, snapshot_age_days


def test_snapshot_age_days_and_stale_detection() -> None:
    now = datetime.now(UTC)
    fresh = {"metadata": {"refresh_status": "complete", "created_at": (now - timedelta(days=3)).isoformat()}}
    old = {"metadata": {"refresh_status": "complete", "created_at": (now - timedelta(days=45)).isoformat()}}

    fresh_age = snapshot_age_days(fresh)
    old_age = snapshot_age_days(old)
    assert fresh_age is not None and 2.0 < fresh_age < 5.0
    assert old_age is not None and old_age > 40.0

    assert not is_snapshot_stale(fresh, stale_after_days=30)
    assert is_snapshot_stale(old, stale_after_days=30)


def test_snapshot_stale_when_missing_timestamp() -> None:
    assert is_snapshot_stale({"metadata": {}}, stale_after_days=30)
    assert is_snapshot_stale(None, stale_after_days=30)



def test_partial_failed_and_legacy_snapshots_are_not_current():
    now = datetime.now(UTC).isoformat()
    for status in (None, "partial", "failed"):
        assert is_snapshot_stale({"metadata": {"created_at": now, "refresh_status": status}}, 30)
