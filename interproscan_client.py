"""Minimal client for the EBI InterProScan5 REST API.

https://www.ebi.ac.uk/Tools/services/rest/iprscan5
No local database required: sequences are submitted as jobs, polled until
finished, and results are fetched as TSV.
"""
import time

import requests

BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
POLL_INTERVAL_SECONDS = 10
MAX_POLL_ATTEMPTS = 60  # ~10 minutes per job
REQUEST_RETRIES = 3
REQUEST_RETRY_BACKOFF_SECONDS = 5


class InterProScanError(RuntimeError):
    pass


def _request_with_retry(method, url, **kwargs):
    # EBI's REST endpoints intermittently return transient 5xx errors under load.
    last_exc = None
    for attempt in range(REQUEST_RETRIES):
        try:
            response = method(url, **kwargs)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as exc:
            last_exc = exc
            if attempt < REQUEST_RETRIES - 1:
                time.sleep(REQUEST_RETRY_BACKOFF_SECONDS)
    raise InterProScanError(f"Request to {url} failed after {REQUEST_RETRIES} attempts: {last_exc}")


def submit_job(sequence, email, title=None):
    payload = {"email": email, "sequence": sequence, "stype": "p"}
    if title:
        payload["title"] = title
    response = _request_with_retry(requests.post, f"{BASE_URL}/run", data=payload)
    return response.text.strip()


def poll_job(job_id):
    for _ in range(MAX_POLL_ATTEMPTS):
        try:
            response = requests.get(f"{BASE_URL}/status/{job_id}")
            response.raise_for_status()
        except requests.exceptions.RequestException:
            # EBI's endpoint occasionally returns transient 5xx errors mid-poll; retry.
            time.sleep(POLL_INTERVAL_SECONDS)
            continue
        status = response.text.strip()
        if status == "FINISHED":
            return status
        if status in ("ERROR", "FAILURE", "NOT_FOUND"):
            raise InterProScanError(f"InterProScan job {job_id} ended with status {status}")
        time.sleep(POLL_INTERVAL_SECONDS)
    raise InterProScanError(f"InterProScan job {job_id} timed out waiting for completion")


def fetch_result_tsv(job_id):
    response = _request_with_retry(requests.get, f"{BASE_URL}/result/{job_id}/tsv")
    return response.text


def run_interproscan_batch(records, email, output_file, progress_callback=None):
    """Submit one InterProScan job per SeqRecord, poll each, and merge TSV output.

    progress_callback(done, total, record_id) is called after each sequence finishes.
    """
    all_lines = []
    total = len(records)
    for i, record in enumerate(records, start=1):
        sequence = f">{record.id}\n{str(record.seq)}"
        job_id = submit_job(sequence, email, title=record.id)
        poll_job(job_id)
        tsv_text = fetch_result_tsv(job_id)
        all_lines.extend(line for line in tsv_text.splitlines() if line.strip())
        if progress_callback:
            progress_callback(i, total, record.id)

    with open(output_file, "w") as fh:
        fh.write("\n".join(all_lines) + ("\n" if all_lines else ""))
    return output_file
