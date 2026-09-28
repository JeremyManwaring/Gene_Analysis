"""Supplementary API updates cannot advertise failed or modified snapshots as fresh."""

from __future__ import annotations

import io
import json
from datetime import UTC, datetime

import pytest

from Gene_Analysis.evidence import http, snapshot


def test_graphql_200_error_and_empty_response_are_recorded_failures(monkeypatch):
    class Response(io.BytesIO):
        headers = {}
    for raw in (b'{"errors":[{"message":"unknown field"}]}', b'', b'{}'):
        monkeypatch.setattr(http.urllib.request, "urlopen", lambda *a, **kw: Response(raw))
        client = http.HttpClient(retries=0)
        result = client.post_json("https://example.test/query", {"query": "query { field }"})
        assert result["_error"]
        assert client.requests[0]["status"] == "failed"
        assert client.requests[0]["attempted_at"]
        assert "fetched_at" not in client.requests[0]


def test_success_response_has_hash_and_source_date(monkeypatch):
    class Response(io.BytesIO):
        headers = {"Last-Modified": "Thu, 24 Sep 2026 09:06:57 GMT", "ETag": "version-1"}
    monkeypatch.setattr(http.urllib.request, "urlopen", lambda *a, **kw: Response(b'{"data":[1]}'))
    client = http.HttpClient(retries=0)
    assert client.get_json("https://example.test/") == {"data": [1]}
    audit = client.requests[0]
    assert audit["status"] == "success"
    assert audit["source_last_modified"].startswith("Thu")
    assert len(audit["response_sha256"]) == 64
    assert audit["fetched_at"]


def test_partial_refresh_preserves_previous_pointer_and_records_failure(monkeypatch, tmp_path):
    latest = tmp_path / "latest.json"
    latest.write_text('{"previous":"complete"}')
    monkeypatch.setattr(snapshot, "LATEST_PATH", latest)
    monkeypatch.setattr(snapshot, "CACHE_DIR", tmp_path)
    class FailingClient:
        def __init__(self, **kwargs):
            self.requests = []
        def get_json(self, url, **kwargs):
            self.requests.append({"url": url, "status": "failed", "error": "timeout"})
            return {"_error": "timeout"}
    monkeypatch.setattr(snapshot, "HttpClient", FailingClient)
    result = snapshot.update_evidence_snapshot(markers=[])
    assert result["metadata"]["refresh_status"] == "failed"
    assert result["metadata"]["requests_failed"] == 1
    assert json.loads(latest.read_text()) == {"previous": "complete"}
    loaded = snapshot.load_snapshot(result["metadata"]["path"])
    assert result == loaded
    assert snapshot.is_snapshot_stale(result, 30)
    path = result["metadata"]["path"]
    payload = json.loads(open(path).read())
    payload["metadata"]["refresh_status"] = "complete"
    with open(path, "w") as handle:
        json.dump(payload, handle)
    with pytest.raises(ValueError, match="checksum"):
        snapshot.load_snapshot(path)


def test_pinned_snapshot_missing_and_overwrite_fail(tmp_path):
    with pytest.raises(FileNotFoundError):
        snapshot.load_snapshot(str(tmp_path / "missing.json"))
    target = tmp_path / "immutable.json"
    snapshot._json_dump(target, {"a": 1}, immutable=True)
    with pytest.raises(FileExistsError):
        snapshot._json_dump(target, {"a": 2}, immutable=True)
    assert json.loads(target.read_text()) == {"a": 1}


def test_future_snapshot_timestamp_never_counts_as_fresh():
    metadata = {"refresh_status": "complete", "created_at": "2999-01-01T00:00:00+00:00"}
    assert snapshot.is_snapshot_stale({"metadata": metadata}, 30)
