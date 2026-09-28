"""Bounded JSON HTTP retrieval with auditable source status and response hashes."""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from typing import Any, Dict, Optional


class HttpClient:
    """Small retrying client; HTTP 200 GraphQL errors are failures, not empty success."""

    def __init__(self, timeout_seconds: int = 20, retries: int = 2, backoff_seconds: float = 0.7):
        self.timeout_seconds = timeout_seconds
        self.retries = retries
        self.backoff_seconds = backoff_seconds
        self.requests: list[dict[str, Any]] = []

    def get_json(self, url: str, params: Optional[Dict[str, Any]] = None) -> Any:
        if params:
            query = urllib.parse.urlencode(params, doseq=True)
            url = f"{url}?{query}"
        return self._request_json(url=url, method="GET", payload=None)

    def post_json(self, url: str, payload: Dict[str, Any]) -> Any:
        return self._request_json(url=url, method="POST", payload=payload)

    def _request_json(self, url: str, method: str, payload: Optional[Dict[str, Any]]) -> Any:
        body = None
        headers = {"Accept": "application/json", "User-Agent": "Gene-Analysis/1.0"}
        if payload is not None:
            body = json.dumps(payload, sort_keys=True).encode("utf-8")
            headers["Content-Type"] = "application/json"
        audit: dict[str, Any] = {"url": url, "method": method, "status": "failed"}
        if body:
            audit["request_body_sha256"] = hashlib.sha256(body).hexdigest()
        last_error: Optional[Exception] = None
        for attempt in range(self.retries + 1):
            audit["attempts"] = attempt + 1
            audit["attempted_at"] = datetime.now(UTC).isoformat()
            req = urllib.request.Request(url=url, data=body, headers=headers, method=method)
            try:
                with urllib.request.urlopen(req, timeout=self.timeout_seconds) as response:
                    raw = response.read()
                    if not raw.strip():
                        raise ValueError("Empty JSON response")
                    data = json.loads(raw.decode("utf-8"))
                    if data == {}:
                        raise ValueError("Empty JSON object")
                    if not isinstance(data, (dict, list)):
                        raise ValueError("Unexpected JSON response type")
                    if isinstance(data, dict) and (data.get("errors") or data.get("error")):
                        raise ValueError(f"API returned an error: {str(data.get('errors') or data.get('error'))[:500]}")
                    audit.update(status="success", fetched_at=datetime.now(UTC).isoformat(),
                                 response_sha256=hashlib.sha256(raw).hexdigest(),
                                 source_last_modified=response.headers.get("Last-Modified"),
                                 source_etag=response.headers.get("ETag"))
                    self.requests.append(audit)
                    return data
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError,
                    ValueError, UnicodeError) as exc:
                last_error = exc
                if attempt < self.retries:
                    time.sleep(self.backoff_seconds * (attempt + 1))
        audit["error"] = str(last_error) if last_error else "Unknown request error"
        self.requests.append(audit)
        return {"_error": audit["error"], "_url": url, "_method": method}
