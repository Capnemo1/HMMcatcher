import os
import subprocess # Importado para la ejecución robusta de comandos
import shutil 
from pathlib import Path
from Bio import SeqIO
import streamlit as st

# --- Configuración Inicial ---
st.set_page_config(
    page_title="Protein Alignment & HMM Profiling Tool",
    page_icon="🧬",
    layout="centered",
)

# Título
st.title("Protein Alignment & HMM Profiling Tool")
st.subheader("Analiza secuencias, construye perfiles HMM y valida dominios.")

# Ruta para la base de datos Pfam (DEBE existir en el entorno de despliegue)
PFAM_DB_PATH = "Pfam-A.hmm"
TEMP_DIR = "temp_processing" # Directorio temporal para todos los archivos

# Función para ejecutar comandos de forma segura (con subprocess)
def run_command(command, step_name):
    """Ejecuta un comando del sistema, maneja errores y captura la salida."""
    st.info(f"Ejecutando: {step_name}...")
    try:
        result = subprocess.run(
            command, 
            shell=True, 
            check=True,  # Lanza CalledProcessError si el comando falla (código de salida != 0)
            capture_output=True, # Captura stdout y stderr
            text=True
        )
        st.success(f"{step_name} completado.")
        # Opcional: mostrar advertencias o información de salida si hay
        if result.stderr:
            st.code(result.stderr, language='text') 
        return True
    except FileNotFoundError:
        st.error(f"Error: La herramienta para '{step_name}' no se encuentra.")
        st.error(f"Asegúrate de que el ejecutable ('{command.split()[0]}') esté instalado y en el PATH.")
        return False
    except subprocess.CalledProcessError as e:
        # Muestra el error exacto reportado por la herramienta (clustalo, hmmbuild, etc.)
        st.error(f"Error en {step_name}. La herramienta reportó un fallo.")
        st.code(e.stderr, language='text')
        return False
    except Exception as e:
        st.error(f"Ocurrió un error inesperado durante {step_name}: {e}")
        return False


# Función para generar nombres dinámicos (usando la carpeta temporal)
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

# Función para extraer secuencias de resultados (Uso de Biopython)
def extract_sequences_from_results(hmm_results, database, output_file):
    protein_names = set()
    try:
        # Leer IDs de las proteínas de los resultados de hmmsearch (columna 1)
        with open(hmm_results, "r") as results_file:
            for line in results_file:
                if not line.startswith("#") and line.strip():
                    columns = line.split()
                    protein_names.add(columns[0])

        # Extraer las secuencias del archivo de la base de datos
        extracted_sequences = []
        for record in SeqIO.parse(database, "fasta"): 
            if record.id in protein_names:
                extracted_sequences.append(record)

        # Escribir las secuencias extraídas
        SeqIO.write(extracted_sequences, output_file, "fasta")
        return True
    except Exception as e:
        st.error(f"Error al extraer secuencias: {e}")
        return False

# --- UI y Flujo Principal ---

st.markdown("### 1. Sube tus archivos")
col1, col2 = st.columns(2)

# Extensiones FASTA ampliadas (¡Excelente mejora!)
FASTA_EXTENSIONS = ["fasta", "fa", "faa", "fna", "txt"]

with col1:
    protein_file = st.file_uploader(
        "Secuencias Iniciales (FASTA)", 
        type=FASTA_EXTENSIONS, 
        help="Secuencias proteicas para construir el HMM."
    )
with col2:
    database_file = st.file_uploader(
        "Base de Datos de Búsqueda (FASTA)", 
        type=FASTA_EXTENSIONS, 
        help="Base de datos de proteínas donde buscarás los homólogos."
    )

if st.button("🧬 Iniciar Flujo de Análisis", type="primary"):
    if not protein_file or not database_file:
        st.warning("Por favor, sube ambos archivos FASTA para continuar.")
    elif not os.path.exists(PFAM_DB_PATH):
           st.error(f"Error de base de datos: Falta el archivo '{PFAM_DB_PATH}'. Este debe estar instalado en el servidor.")
    else:
        # Inicializar archivos de entrada en la carpeta temporal
        protein_filename = os.path.join(TEMP_DIR, "temp_protein_input.fasta")
        database_filename = os.path.join(TEMP_DIR, "temp_protein_database.fasta")
        
        # Generar nombres de archivos de salida
        output_files = generate_output_filenames(protein_filename)
        
        try:
            # 1. Guardar archivos subidos
            with open(protein_filename, "wb") as f:
                f.write(protein_file.read())
            with open(database_filename, "wb") as f:
                f.write(database_file.read())
        except Exception as e:
            st.error(f"Error al guardar los archivos: {e}")
            st.stop()


        with st.spinner('Procesando el flujo de bioinformática...'):
            # Paso 1: Alinear proteínas (clustalo)
            alignment_command = f"clustalo -i {protein_filename} -o {output_files['alignment']} --outfmt=st --force"
            if not run_command(alignment_command, "Alineación con Clustal Omega"):
                st.stop()

            # Paso 2: Crear perfil HMM (hmmbuild)
            hmmbuild_command = f"hmmbuild {output_files['hmm_profile']} {output_files['alignment']}"
            if not run_command(hmmbuild_command, "Creación del perfil HMM"):
                st.stop()

            # Paso 3: Buscar en la base de datos (hmmsearch)
            # *** CAMBIO APLICADO AQUÍ: Filtro de E-value ***
            hmmsearch_command = f"hmmsearch --tblout {output_files['search_results']} -E 1e-4 {output_files['hmm_profile']} {database_filename}"
            if not run_command(hmmsearch_command, "Búsqueda HMM en la base de datos (E-value < 1e-4)"):
                st.stop()

            # Paso 4: Extraer secuencias (Biopython)
            if not extract_sequences_from_results(
                output_files["search_results"], database_filename, output_files["extracted_sequences"]
            ):
                st.stop()
            st.success("Extracción de secuencias homólogas completada.")
            
            # Paso 5: Validar dominios con Pfam (hmmscan)
            pfam_scan_command = f"hmmscan --domtblout {output_files['pfam_results']} {PFAM_DB_PATH} {output_files['extracted_sequences']}"
            if not run_command(pfam_scan_command, "Validación de dominios con Pfam"):
                st.stop()
        
        # --- Resultados y Descarga ---
        st.markdown("### 🎉 Resultados listos para descargar")
        
        download_cols = st.columns(5)
        
        download_files = [
            ("Alineación (Stockholm)", output_files["alignment"]),
            ("Perfil HMM", output_files["hmm_profile"]),
            ("Resultados de Búsqueda (TSV)", output_files["search_results"]),
            ("Secuencias Extraídas (FASTA)", output_files["extracted_sequences"]),
            ("Resultados Pfam (TSV)", output_files["pfam_results"]),
        ]
        
        for i, (label, path) in enumerate(download_files):
            try:
                with open(path, "rb") as f:
                    download_cols[i].download_button(
                        label, 
                        f.read(), 
                        file_name=Path(path).name,
                        key=f"download_{i}"
                    )
            except FileNotFoundError:
                download_cols[i].error("Archivo no encontrado.")
        
        # 6. Limpiar archivos temporales
        st.markdown("---")
        
        try:
            # Elimina la carpeta temporal y su contenido después de la descarga
            shutil.rmtree(TEMP_DIR)
            st.info(f"Archivos temporales en '{TEMP_DIR}' eliminados.")
        except Exception:
            # En caso de error de permisos o si la carpeta no existe
            pass 
