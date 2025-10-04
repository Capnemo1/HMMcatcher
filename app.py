import os
import subprocess
import shutil
from pathlib import Path
from Bio import SeqIO
import streamlit as st

# --- Configuración Inicial ---
st.set_page_config(
    page_title="HMMcatcher | Protein Alignment & Profiling Tool",
    page_icon="🧬",
    layout="centered",
)

# --- Encabezado principal ---
st.markdown(
    """
    <div style="text-align:center;">
        <h1 style="color:#2E86C1;">🧬 HMMcatcher</h1>
        <h3>Protein Alignment & HMM Profiling Tool</h3>
        <p style="color:gray;">Analiza secuencias, construye perfiles HMM y valida dominios opcionalmente con Pfam.</p>
    </div>
    """,
    unsafe_allow_html=True,
)
st.divider()

# --- Variables globales ---
PFAM_DB_PATH = "Pfam-A.hmm"
TEMP_DIR = "temp_processing"

# --- Función para ejecutar comandos del sistema ---
def run_command(command, step_name):
    """Ejecuta un comando del sistema, maneja errores y captura la salida."""
    with st.status(f"⏳ {step_name} en progreso...", expanded=True) as status:
        try:
            result = subprocess.run(
                command, shell=True, check=True, capture_output=True, text=True
            )
            status.update(label=f"✅ {step_name} completado.", state="complete")
            if result.stderr:
                st.code(result.stderr, language="text")
            return True
        except FileNotFoundError:
            status.update(label=f"❌ Herramienta no encontrada para '{step_name}'", state="error")
            st.error(f"Asegúrate de que '{command.split()[0]}' esté instalado.")
            return False
        except subprocess.CalledProcessError as e:
            status.update(label=f"⚠️ Error en {step_name}", state="error")
            st.code(e.stderr, language="text")
            return False
        except Exception as e:
            status.update(label=f"⚠️ Error inesperado en {step_name}", state="error")
            st.error(f"Ocurrió un error inesperado: {e}")
            return False

# --- Función para generar nombres de archivos ---
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

# --- Función para extraer secuencias ---
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
        st.error(f"Error al extraer secuencias: {e}")
        return False

# --- Sección de carga de archivos ---
st.markdown("### 📂 1. Sube tus archivos")
col1, col2 = st.columns(2)
FASTA_EXTENSIONS = ["fasta", "fa", "faa", "fna", "txt"]

with col1:
    protein_file = st.file_uploader(
        "Secuencias iniciales (FASTA)",
        type=FASTA_EXTENSIONS,
        help="Secuencias proteicas para construir el perfil HMM.",
    )
with col2:
    database_file = st.file_uploader(
        "Base de datos de búsqueda (FASTA)",
        type=FASTA_EXTENSIONS,
        help="Base de datos de proteínas donde buscarás los homólogos.",
    )

st.divider()

# --- Ejecución principal ---
if st.button("🚀 Iniciar Flujo de Análisis", type="primary"):
    if not protein_file or not database_file:
        st.warning("⚠️ Por favor, sube ambos archivos FASTA para continuar.")
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
        st.error(f"Error al guardar archivos: {e}")
        st.stop()

    st.info("🧠 Iniciando procesamiento bioinformático...")
    st.divider()

    # Paso 1: Clustal Omega
    alignment_command = f"clustalo -i {protein_filename} -o {output_files['alignment']} --outfmt=st --force"
    if not run_command(alignment_command, "Alineación con Clustal Omega"):
        st.stop()

    # Paso 2: hmmbuild
    hmmbuild_command = f"hmmbuild {output_files['hmm_profile']} {output_files['alignment']}"
    if not run_command(hmmbuild_command, "Creación del perfil HMM"):
        st.stop()

    # Paso 3: hmmsearch
    hmmsearch_command = (
        f"hmmsearch --tblout {output_files['search_results']} -E 1e-4 "
        f"{output_files['hmm_profile']} {database_filename}"
    )
    if not run_command(hmmsearch_command, "Búsqueda HMM en la base de datos (E-value < 1e-4)"):
        st.stop()

    # Paso 4: extracción
    if not extract_sequences_from_results(
        output_files["search_results"], database_filename, output_files["extracted_sequences"]
    ):
        st.stop()
    st.success("✅ Extracción de secuencias homólogas completada.")

    # Paso 5: Validación opcional con Pfam
    if os.path.exists(PFAM_DB_PATH):
        pfam_scan_command = (
            f"hmmscan --domtblout {output_files['pfam_results']} "
            f"{PFAM_DB_PATH} {output_files['extracted_sequences']}"
        )
        if run_command(pfam_scan_command, "Validación de dominios con Pfam"):
            pfam_results_generated = True
    else:
        st.warning("⚠️ No se encontró la base de datos Pfam-A.hmm. La validación de dominios se omitió.")

    st.divider()
    st.markdown("### 📁 Resultados listos para descargar")

    # Mostrar los botones de descarga con íconos
    download_files = [
        ("📊 Alineación (Stockholm)", output_files["alignment"], True),
        ("📈 Perfil HMM", output_files["hmm_profile"], True),
        ("📋 Resultados HMMsearch (TSV)", output_files["search_results"], True),
        ("🧩 Secuencias Extraídas (FASTA)", output_files["extracted_sequences"], True),
        ("🧠 Resultados Pfam (TSV)", output_files["pfam_results"], pfam_results_generated),
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
            download_cols[i].error("Archivo no encontrado.")

    st.divider()
    try:
        shutil.rmtree(TEMP_DIR)
        st.info(f"🧹 Archivos temporales eliminados de '{TEMP_DIR}'.")
    except Exception:
        pass

# --- Pie de página ---
st.markdown(
    """
    <hr>
    <div style="text-align:center; color:gray;">
        <small>Developed with ❤️ by Erick Arroyo · Version Beta · Powered by Streamlit & HMMER</small>
    </div>
    """,
    unsafe_allow_html=True,
)
