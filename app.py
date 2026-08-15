import shutil
import tempfile
from pathlib import Path

import streamlit as st
from Bio import SeqIO

import hmmer_tools as hmmer
import interproscan_client as interpro
import visit_counter

# --- Initial Configuration ---
st.set_page_config(
    page_title="HMMcatcher | Protein Family Mining Tool",
    page_icon="🧬",
    layout="centered",
)

# --- Header ---
st.markdown(
    """
    <div style="text-align:center;">
        <h1 style="color:#2E86C1;">🧬 HMMcatcher</h1>
        <h3>Protein Family Discovery & HMM Profiling</h3>
        <p style="color:gray; font-size: 1.1em; margin-top: 10px;">
            This tool facilitates the search for <b>protein families, transcription factors, or any protein of interest</b>
            by <b>building your own HMM profile</b> from known orthologous sequences, or by
            <b>using an existing HMM profile</b> (e.g. downloaded from PANTHER), to search sensitively
            within <b>unannotated proteomes</b>.
        </p>
        <p style="color:#5D6D7E; font-style:italic;">
            HMMcatcher integrates Clustal Omega, HMMER, and EBI InterProScan in a streamlined
            Streamlit interface, making complex bioinformatics analyses accessible to non-specialists.
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)
st.divider()

# --- Citation Info (Sidebar) ---
st.sidebar.markdown("## 📚 How to Cite")
st.sidebar.markdown(
    """
    **Cita Sugerida (Suggested Citation):**

    Arroyo-Álvarez, E. (2025). HMMcatcher (v0.1-beta) [Computer software]. Zenodo.

    <a href="https://doi.org/10.5281/zenodo.17266955" target="_blank" style="color: #2E86C1; text-decoration: none;">
        https://doi.org/10.5281/zenodo.17266955
    </a>
    """,
    unsafe_allow_html=True,
)
st.sidebar.markdown("---")

# --- Visit Counter (best-effort; silently absent if not configured) ---
if "visit_count" not in st.session_state:
    try:
        api_key = st.secrets.get("COUNTERAPI_KEY")
        workspace = st.secrets.get("COUNTERAPI_WORKSPACE", visit_counter.DEFAULT_WORKSPACE)
        counter_name = st.secrets.get("COUNTERAPI_COUNTER", visit_counter.DEFAULT_COUNTER)
    except Exception:
        api_key, workspace, counter_name = None, visit_counter.DEFAULT_WORKSPACE, visit_counter.DEFAULT_COUNTER
    st.session_state.visit_count = visit_counter.record_visit(api_key, workspace, counter_name)

if st.session_state.visit_count is not None:
    st.sidebar.caption(f"👁️ {st.session_state.visit_count:,} visits")

EVALUE_OPTIONS = ["1e-3", "1e-4", "1e-5", "1e-10", "1e-20", "1e-30"]
FASTA_EXTENSIONS = ["fasta", "fa", "faa", "fna", "txt"]


def get_workdir():
    if "workdir" not in st.session_state:
        st.session_state.workdir = tempfile.mkdtemp(prefix="hmmcatcher_")
    return Path(st.session_state.workdir)


def reset_session():
    workdir = st.session_state.get("workdir")
    if workdir and Path(workdir).exists():
        shutil.rmtree(workdir, ignore_errors=True)
    for key in list(st.session_state.keys()):
        del st.session_state[key]


def save_upload(uploaded_file, destination):
    with open(destination, "wb") as fh:
        fh.write(uploaded_file.getbuffer())
    return destination


def run_step(step_name, fn, *args, **kwargs):
    """Run fn inside a labeled st.status block, matching the app's existing step UI."""
    with st.status(f"⏳ {step_name} in progress...", expanded=True) as status:
        try:
            result = fn(*args, **kwargs)
            status.update(label=f"✅ {step_name} completed.", state="complete")
            return result
        except (hmmer.ToolError, interpro.InterProScanError) as e:
            status.update(label=f"⚠️ Error during {step_name}", state="error")
            st.code(str(e), language="text")
            raise


# --- Step 1: search strategy ---
st.markdown("### 🧭 Step 1. Choose search strategy")
search_mode = st.radio(
    "Search strategy",
    ["Build HMM from reference proteins", "Use existing HMM profile"],
    label_visibility="collapsed",
    horizontal=True,
)

# --- Step 2: uploads ---
st.markdown("### 📂 Step 2. Upload your files")
reference_file = None
hmm_file = None
col1, col2 = st.columns(2)
if search_mode == "Build HMM from reference proteins":
    with col1:
        reference_file = st.file_uploader(
            "Reference proteins (FASTA)",
            type=FASTA_EXTENSIONS,
            help="Known orthologous protein sequences to build the HMM profile.",
        )
else:
    with col1:
        hmm_file = st.file_uploader(
            "HMM profile (.hmm or .txt)",
            type=["hmm", "txt"],
            help="A pre-built HMM profile, e.g. downloaded from PANTHER (PANTHER exports these as .txt, "
            "but the content is a standard HMM file).",
        )
with col2:
    proteome_file = st.file_uploader(
        "Target proteome (FASTA)",
        type=FASTA_EXTENSIONS,
        help="Protein database where candidates will be searched.",
    )

st.divider()

# --- Step 3: search parameters ---
st.markdown("### ⚙️ Step 3. Search parameters")
evalue_label = st.selectbox("E-value threshold", EVALUE_OPTIONS, index=2)
evalue = float(evalue_label)

st.divider()

# --- Step 4: InterProScan validation (optional) ---
st.markdown("### 🧠 Step 4. Domain validation (InterProScan) — optional")
run_interpro = st.checkbox(
    "Validate extracted candidates with InterProScan",
    value=False,
    help="Not required to run the analysis. The HMM search results and extracted FASTA are always "
    "produced regardless of this option, e.g. for manual genome-wide curation.",
)
interpro_email = None
max_interpro_sequences = 20
if run_interpro:
    interpro_email = st.text_input(
        "Email (required by EBI InterProScan)",
        help="EBI requires an email address to track submitted jobs.",
    )
    max_interpro_sequences = st.number_input("Max candidate sequences to validate", min_value=1, value=20, step=1)
    st.caption(
        "InterProScan jobs run one sequence at a time and typically take 1-5 minutes each, "
        "so validating many candidates can take a while. If this step fails or times out, "
        "the HMM search results and extracted sequences remain available for download."
    )

st.divider()

# --- Main Execution ---
if st.button("🚀 Run Analysis", type="primary"):
    inputs_ok = proteome_file is not None and (
        (search_mode == "Build HMM from reference proteins" and reference_file is not None)
        or (search_mode == "Use existing HMM profile" and hmm_file is not None)
    )
    if not inputs_ok:
        st.warning("⚠️ Please upload all required files for the selected search strategy.")
        st.stop()
    if run_interpro and not interpro_email:
        st.warning("⚠️ Please provide an email address to run InterProScan, or uncheck that option to skip it.")
        st.stop()

    workdir = get_workdir()
    proteome_path = save_upload(proteome_file, workdir / "target_proteome.fasta")

    # --- Essential pipeline: HMM search + extraction. Always run; results are kept even if the
    # optional InterProScan step below fails. ---
    try:
        if search_mode == "Build HMM from reference proteins":
            base_name = Path(reference_file.name).stem
            reference_path = save_upload(reference_file, workdir / "reference_proteins.fasta")
            paths = hmmer.generate_output_filenames(workdir, base_name)

            run_step("Protein alignment (Clustal Omega)", hmmer.run_clustalo, reference_path, paths["alignment"])
            run_step("HMM profile construction (hmmbuild)", hmmer.run_hmmbuild, paths["alignment"], paths["hmm_profile"])
            alignment_path = paths["alignment"]
        else:
            base_name = Path(hmm_file.name).stem
            paths = hmmer.generate_output_filenames(workdir, base_name)
            paths["hmm_profile"] = save_upload(hmm_file, workdir / f"{base_name}.hmm")
            alignment_path = None

        run_step(
            f"HMM database search (E-value < {evalue_label})",
            hmmer.run_hmmsearch,
            paths["hmm_profile"],
            proteome_path,
            paths["search_results"],
            evalue,
        )
        hits = hmmer.parse_hmmsearch_hits(paths["search_results"])
        st.success(f"✅ {len(hits)} candidate(s) found at E-value <= {evalue_label}.")

        run_step(
            "Homologous sequence extraction",
            hmmer.extract_sequences,
            hits,
            proteome_path,
            paths["extracted_sequences"],
        )
    except hmmer.ToolError:
        st.stop()

    st.session_state.results = {
        "hits": hits,
        "paths": {
            "alignment": alignment_path,
            "hmm_profile": paths["hmm_profile"],
            "search_results": paths["search_results"],
            "extracted_sequences": paths["extracted_sequences"],
            "interproscan_results": None,
        },
    }

    # --- Optional InterProScan validation. A failure here does not remove the results above. ---
    if run_interpro and hits:
        candidates = list(SeqIO.parse(paths["extracted_sequences"], "fasta"))
        capped_candidates = candidates[: int(max_interpro_sequences)]
        if len(candidates) > len(capped_candidates):
            st.warning(
                f"⚠️ {len(candidates)} candidates found; validating only the first "
                f"{len(capped_candidates)} with InterProScan (raise the limit above to validate more)."
            )

        with st.status(f"🧠 Validating {len(capped_candidates)} sequence(s) with InterProScan...", expanded=True) as status:
            progress_bar = st.progress(0.0)

            def update_progress(done, total, record_id):
                progress_bar.progress(done / total, text=f"InterProScan: {record_id} ({done}/{total})")

            try:
                interpro.run_interproscan_batch(
                    capped_candidates, interpro_email, paths["interproscan_results"], update_progress
                )
                status.update(label="✅ InterProScan validation completed.", state="complete")
                st.session_state.results["paths"]["interproscan_results"] = paths["interproscan_results"]
            except interpro.InterProScanError as e:
                status.update(label="⚠️ InterProScan validation failed (other results are still available below)", state="error")
                st.code(str(e), language="text")
    elif run_interpro and not hits:
        st.info("ℹ️ No candidates passed the E-value threshold; skipping InterProScan.")

# --- Results (persists across reruns, e.g. clicking a download button) ---
if st.session_state.get("results"):
    results = st.session_state.results
    st.divider()
    st.markdown("### 📁 Results ready for download")

    if results["hits"]:
        st.dataframe(results["hits"], use_container_width=True)
    else:
        st.info("No candidate proteins passed the E-value threshold.")

    paths = results["paths"]
    download_items = [
        ("📊 Alignment (Stockholm)", paths["alignment"]),
        ("📈 HMM Profile", paths["hmm_profile"]),
        ("📋 HMMsearch Results (TSV)", paths["search_results"]),
        ("🧩 Extracted Sequences (FASTA)", paths["extracted_sequences"]),
        ("🧠 InterProScan Results (TSV)", paths["interproscan_results"]),
    ]
    available = [(label, path) for label, path in download_items if path and Path(path).exists()]
    if available:
        download_cols = st.columns(len(available))
        for col, (label, path) in zip(download_cols, available):
            with open(path, "rb") as fh:
                col.download_button(label, fh.read(), file_name=Path(path).name, key=f"download_{label}")

    st.divider()
    if st.button("🔄 Start new analysis"):
        reset_session()
        st.rerun()

# --- Footer ---
st.markdown(
    """
    <hr>
    <div style="text-align:center; color:gray;">
        <small>
            Developed by <b>Erick Arroyo</b> · Center for Scientific Research of Yucatán (CICY) <br>
            Contact: <a href="mailto:erick.arroyo@cicy.mx">erick.arroyo@cicy.mx</a> <br>
            <i>Beta version – Powered by Streamlit, HMMER, Clustal Omega & InterProScan.</i>
        </small>
    </div>
    """,
    unsafe_allow_html=True,
)
