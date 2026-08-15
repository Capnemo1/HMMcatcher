"""Visit counter backed by CounterAPI (https://counterapi.dev), v2.

Requires a free CounterAPI account and API key (v1 was free/keyless but is
now deprecated). The key must be provided via Streamlit secrets
(COUNTERAPI_KEY) -- never hardcode it in source, since this repo is public.

Every function here is best-effort: any failure (missing secret, network
error, unexpected response shape) returns None instead of raising, so a
broken or unconfigured counter never breaks the rest of the app.
"""
import requests

BASE_URL = "https://api.counterapi.dev/v2"
DEFAULT_WORKSPACE = "hmmcatcher"
DEFAULT_COUNTER = "visits"
REQUEST_TIMEOUT_SECONDS = 3


def _extract_count(payload):
    data = payload.get("data", payload) if isinstance(payload, dict) else None
    if not isinstance(data, dict):
        return None
    for key in ("up_count", "count", "value"):
        if key in data:
            try:
                return int(data[key])
            except (TypeError, ValueError):
                continue
    return None


def record_visit(api_key, workspace=DEFAULT_WORKSPACE, counter=DEFAULT_COUNTER):
    """Increment the visit counter once and return the new total, or None on any failure."""
    if not api_key:
        return None
    try:
        response = requests.get(
            f"{BASE_URL}/{workspace}/{counter}/up",
            headers={"Authorization": f"Bearer {api_key}"},
            timeout=REQUEST_TIMEOUT_SECONDS,
        )
        response.raise_for_status()
        return _extract_count(response.json())
    except Exception:
        return None
