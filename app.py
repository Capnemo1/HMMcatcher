import os
import subprocess
import shutil
from pathlib import Path
from Bio import SeqIO
import streamlit as st

# --- Initial Configuration ---
st.set_page_config(
    page_title="HMMcatcher | Protein Alignment & Profiling Tool",
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
            by utilizing <b>known orthologous sequences</b> provided in <b>FASTA format</b>. 
            The algorithm subsequently constructs a <b>Hidden Markov Model (HMM)</b> profile 
            to enable sensitive searching within <b>unannotated proteomes</b>.
        </p>
        <p style="color:#5D6D7E; font-style:italic;">
            HMMcatcher integrates Clustal Omega and HMMER in a streamlined Streamlit interface,
            making complex bioinformatics analyses accessible to non-specialists.
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
    unsafe_allow_html=True
)
st.sidebar.markdown("---")


# --- Global Variables ---
PFAM_DB_PATH = "Pfam-A.hmm"
TEMP_DIR = "temp_processing"

# --- Run System Commands ---
def run_command(command, step_name):
    """Executes a system command, handling errors and capturing output."""
    with st.status(f"⏳ {step_name} in progress...", expanded=True) as status:
        try:
            result = subprocess.run(
                command, shell=True, check=True, capture_output=True, text=True
            )
            status.update(label=f"✅ {step_name} completed.", state="complete")
            if result.stderr:
                st.code(result.stderr, language="text")
            return True
        except FileNotFoundError:
            status.update(label=f"❌ Tool not found for '{step_name}'", state="error")
            st.error(f"Make sure '{command.split()[0]}' is installed and in PATH.")
            return False
        except subprocess.CalledProcessError as e:
            status.update(label=f"⚠️ Error during {step_name}", state="error")
            st.code(e.stderr, language="text")
            return False
        except Exception as e:
            status.update(label=f"⚠️ Unexpected error in {step_name}", state="error")
            st.error(f"An unexpected error occurred: {e}")
            return False

# --- Generate Dynamic Filenames ---
def generate_output_filenames(input_file):
    base_name = Path(input_file).stem
    os.makedirs(TEMP_DIR, exist_ok=True)
    return {
        "alignment": os.path.join(TEMP_DIR, f"{base_name}_aligned.sto"),
        "hmm_profile": os.path.join(TEMP_DIR, f"{base_name}_profile.hmm"),
        "search_results": os.path.join(TEMP_DIR, f"{base_name}_search_results.tsv"),
        "extracted_sequences": os.path.join(TEMP_DIR, f"{base_name}_extracted_sequences.fasta"),
        "pfam_results": os.path.join(TEMP_DIR, f"{base_name}_pfam_scan.tsv"),
    }

# --- Extract Sequences from hmmsearch results ---
def extract_sequences_from_results(hmm_results, database, output_file):
    protein_names = set()
    try:
        with open(hmm_results, "r") as results_file:
            for line in results_file:
                if not line.startswith("#") and line.strip():
                    columns = line.split()
                    protein_names.add(columns[0])

        extracted_sequences = [
            record for record in SeqIO.parse(database, "fasta") if record.id in protein_names
        ]
        SeqIO.write(extracted_sequences, output_file, "fasta")
        return True
    except Exception as e:
        st.error(f"Error extracting sequences: {e}")
        return False

# --- File Upload Section ---
st.markdown("### 📂 Step 1. Upload your files")
col1, col2 = st.columns(2)
FASTA_EXTENSIONS = ["fasta", "fa", "faa", "fna", "txt"]

with col1:
    protein_file = st.file_uploader(
        "Input sequences (FASTA)",
        type=FASTA_EXTENSIONS,
        help="Protein sequences to build the HMM profile.",
    )
with col2:
    database_file = st.file_uploader(
        "Search database (FASTA)",
        type=FASTA_EXTENSIONS,
        help="Protein database where homologs will be searched.",
    )

st.divider()

# --- Main Execution ---
if st.button("🚀 Run Analysis", type="primary"):
    if not protein_file or not database_file:
        st.warning("⚠️ Please upload both FASTA files to continue.")
        st.stop()

    protein_filename = os.path.join(TEMP_DIR, "temp_protein_input.fasta")
    database_filename = os.path.join(TEMP_DIR, "temp_protein_database.fasta")
    output_files = generate_output_filenames(protein_filename)
    pfam_results_generated = False

    try:
        with open(protein_filename, "wb") as f:
            f.write(protein_file.read())
        with open(database_filename, "wb") as f:
            f.write(database_file.read())
    except Exception as e:
        st.error(f"Error saving uploaded files: {e}")
        st.stop()

    st.info("🧠 Starting bioinformatics workflow...")
    st.divider()

    # Step 1: Clustal Omega
    alignment_command = f"clustalo -i {protein_filename} -o {output_files['alignment']} --outfmt=st --force"
    if not run_command(alignment_command, "Protein alignment (Clustal Omega)"):
        st.stop()

    # Step 2: hmmbuild
    hmmbuild_command = f"hmmbuild {output_files['hmm_profile']} {output_files['alignment']}"
    if not run_command(hmmbuild_command, "HMM profile construction"):
        st.stop()

    # Step 3: hmmsearch
    hmmsearch_command = (
        f"hmmsearch --tblout {output_files['search_results']} -E 1e-4 "
        f"{output_files['hmm_profile']} {database_filename}"
    )
    if not run_command(hmmsearch_command, "HMM database search (E-value < 1e-4)"):
        st.stop()

    # Step 4: Extract sequences
    if not extract_sequences_from_results(
        output_files["search_results"], database_filename, output_files["extracted_sequences"]
    ):
        st.stop()
    st.success("✅ Homologous sequence extraction completed.")

    # Step 5: Optional Pfam validation
    if os.path.exists(PFAM_DB_PATH):
        pfam_scan_command = (
            f"hmmscan --domtblout {output_files['pfam_results']} "
            f"{PFAM_DB_PATH} {output_files['extracted_sequences']}"
        )
        if run_command(pfam_scan_command, "Domain validation with Pfam"):
            pfam_results_generated = True
    else:
        st.warning("⚠️ Pfam-A.hmm not found. Domain validation step skipped.")

    st.divider()
    st.markdown("### 📁 Results ready for download")

    # --- Persistent Download Section ---
    if "download_ready" not in st.session_state:
        st.session_state.download_ready = True

    if st.session_state.download_ready:
        download_files = [
            ("📊 Alignment (Stockholm)", output_files["alignment"], True),
            ("📈 HMM Profile", output_files["hmm_profile"], True),
            ("📋 HMMsearch Results (TSV)", output_files["search_results"], True),
            ("🧩 Extracted Sequences (FASTA)", output_files["extracted_sequences"], True),
            ("🧠 Pfam Results (TSV)", output_files["pfam_results"], pfam_results_generated),
        ]
        filtered_download_files = [item for item in download_files if item[2]]
        download_cols = st.columns(len(filtered_download_files))

        for i, (label, path, _) in enumerate(filtered_download_files):
            try:
                with open(path, "rb") as f:
                    download_cols[i].download_button(
                        label, f.read(), file_name=Path(path).name, key=f"download_{i}"
                    )
            except FileNotFoundError:
                download_cols[i].error("File not found.")

    st.divider()
    try:
        shutil.rmtree(TEMP_DIR)
        st.info(f"🧹 Temporary files removed from '{TEMP_DIR}'.")
    except Exception:
        pass

# --- Footer ---
st.markdown(
    """
    <hr>
    <div style="text-align:center; color:gray;">
        <small>
            Developed by <b>Erick Arroyo</b> · Center for Scientific Research of Yucatán (CICY) <br>
            Contact: <a href="mailto:erick.arroyo@cicy.mx">erick.arroyo@cicy.mx</a> <br>
            <i>Beta version – Powered by Streamlit, HMMER & Clustal Omega.</i>
        </small>
    </div>
    """,
    unsafe_allow_html=True,
)
