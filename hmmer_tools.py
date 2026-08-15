"""Subprocess wrappers around ClustalO / HMMER used by the Streamlit app."""
import subprocess
from pathlib import Path

from Bio import SeqIO


class ToolError(RuntimeError):
    """Raised when an external CLI tool exits with a non-zero status."""


def _run(command, description):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise ToolError(f"{description} failed:\n{result.stderr or result.stdout}")
    return result


def run_clustalo(input_fasta, output_sto):
    _run(
        ["clustalo", "-i", str(input_fasta), "-o", str(output_sto), "--outfmt=st", "--force"],
        "Alignment (clustalo)",
    )
    return output_sto


def run_hmmbuild(alignment_sto, output_hmm):
    _run(["hmmbuild", str(output_hmm), str(alignment_sto)], "HMM profile creation (hmmbuild)")
    return output_hmm


def run_hmmsearch(hmm_profile, target_fasta, output_tblout, evalue):
    _run(
        ["hmmsearch", "-E", str(evalue), "--tblout", str(output_tblout), str(hmm_profile), str(target_fasta)],
        "HMM search (hmmsearch)",
    )
    return output_tblout


def parse_hmmsearch_hits(tblout_file):
    """Return a list of dicts: {accession, evalue, score} from a hmmsearch --tblout file."""
    hits = []
    with open(tblout_file, "r") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            columns = line.split()
            hits.append({
                "accession": columns[0],
                "evalue": float(columns[4]),
                "score": float(columns[5]),
            })
    return hits


def extract_sequences(hits, database_fasta, output_fasta):
    wanted_ids = {hit["accession"] for hit in hits}
    extracted = [record for record in SeqIO.parse(database_fasta, "fasta") if record.id in wanted_ids]
    SeqIO.write(extracted, output_fasta, "fasta")
    return output_fasta


def generate_output_filenames(workdir, base_name):
    workdir = Path(workdir)
    return {
        "alignment": workdir / f"{base_name}_aligned.sto",
        "hmm_profile": workdir / f"{base_name}_profile.hmm",
        "search_results": workdir / f"{base_name}_search_results.tsv",
        "extracted_sequences": workdir / f"{base_name}_extracted_sequences.fasta",
        "interproscan_results": workdir / f"{base_name}_interproscan_results.tsv",
    }
