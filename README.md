# **MotSASi 2.0**  

Este repositorio contiene todos los archivos de código necesarios para ejecutar el pipeline **MotSASi 2.0**.  

## **Descripción**  

**MotSASi 2.0** es una estrategia diseñada para mejorar la predicción de **Short Linear Motifs (SLiMs)** funcionales mediante la integración de información sobre variantes genómicas y datos estructurales. Utiliza tanto estructuras cristalográficas depositadas en **PDB** como modelos de **AlphaFold2**.  

Además, permite la generación de **matrices de Sustitución de Un Solo Aminoácido (SAS)**, proporcionando predicciones para todas las variantes missense que ocurren dentro de SLiMs.  

Si esta herramienta te resulta útil, te pedimos que cites nuestro trabajo.  

---

## **Requisitos e Instrucciones de Uso**  

Para asegurar compatibilidad y reproducibilidad, recomendamos ejecutar los scripts dentro de un entorno **Anaconda**.  

### **1. Instalación de Anaconda y Creación del Entorno**  

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2024.02-0-Linux-x86_64.sh  
bash Anaconda3-2024.02-0-Linux-x86_64.sh  
conda create --name motsasi python=3.10  
conda activate motsasi  
```

Instalar **Biopython**:  

```bash
conda install -c conda-forge biopython  
```

### **2. Clonar el Repositorio**  

```bash
git clone https://github.com/fginob1/MotSASi2.0.git  
cd MotSASi2.0  
git lfs pull  
```

### **3. Crear Directorios Necesarios**  

```bash
mkdir tmp Motifs PDB AF_human_proteome repaired_pdbs secondary_structures scratch  
```

---

## **Instalación de Dependencias**  

### **1. FoldX**  

Se requiere una licencia académica para ejecutar **FoldX** localmente.  

- Solicitar la licencia en: [https://foldxsuite.crg.eu/academic-license-info](https://foldxsuite.crg.eu/academic-license-info).  
- Colocar el **ejecutable** de FoldX y la carpeta **molecules** dentro del directorio principal `MotSASi2.0`.  

### **2. FreeSASA**  

Instalar con `pip`:  

```bash
pip install freesasa  
```

### **3. SCRATCH**  

```bash
cd scratch  
wget http://download.igb.uci.edu/SCRATCH-1D_1.3.tar.gz  
tar -zxf SCRATCH-1D_1.3.tar.gz  
cd SCRATCH-1D_1.3  
perl install.pl  
cd ../..  
```

---

## **Descarga de Bases de Datos**  

### **1. Base de datos PDB (estructuras cristalográficas)**  

Si ya tienes la base de datos PDB descargada, solo indica la ruta correcta.  

```bash
cd PDB  
mkdir zipped unzipped  
rsync -rlpt -v -z --delete --port=33444
rsync.rcsb.org::ftp_data/structures/all/pdb/ ./zipped  
cd ..  
```

### **2. Base de datos AlphaFold para el proteoma humano**  

```bash
cd AF_human_proteome  
wget -c https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar  
tar -xvf UP000005640_9606_HUMAN_v4.tar  
rm UP000005640_9606_HUMAN_v4.tar  
cd ..  
```

### **3. Base de datos UniProt (SwissProt & Trembl)**  

```bash
cd UniProt  
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz  
zcat uniprot_sprot.fasta.gz | bgzip > uniprot_sprot.fasta.bgz  

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz  
zcat uniprot_trembl.fasta.gz | bgzip > uniprot_trembl.fasta.bgz  
cd ..  
```

---

## **Ejecutando MotSASi 2.0**  

### **MotSASi_2_cluster.py: Integración de Información Estructural y de Variantes para el Análisis de SLiMs**  

El script **MotSASi_2_cluster.py** realiza las siguientes tareas:  

- **Identificación de motivos:** busca el motivo de interés dentro del proteoma humano.  
- **Recolección de variantes:** obtiene variantes benignas y patogénicas de **ClinVar** y **gnomAD**.  
- **Cálculo de la Matriz de Sustitución:**  
  - Si hay estructuras PDB disponibles, se genera con **FoldX**.  
  - En caso contrario, se usa una predicción de **AlphaFold2**.  
- **Generación de la Matriz ClinSig y de Frecuencia Alélica:** basada en instancias de motivos de **ELM** (control positivo).  
- **Cálculo de conservación** para el control positivo.  
- **Cálculo de probabilidad de estructura secundaria** usando **SCRATCH**.  
- **Cálculo de Solvent Accessible Surface Area (SASA)**.  
- **Obtención de términos GO** desde **UniProt**.  

### **Ejemplo de Ejecución**  

```bash
python3 MotSASi_2_cluster.py [MOTIF] [DOT-SEPARATED MOTIF] [ELM MOTIF NAME]  
```

#### **Ejemplo concreto:**  

```bash
python3 MotSASi_2_cluster.py [RK]P[^P][^P]L.[LIVMF] RK.P.^P.^P.L.x.LIVMF DOC_MAPK_JIP1_4  
```

### **Incorporando Modelos de AlphaFold2 en el Pipeline**  

Los modelos AlphaFold deben ubicarse en una carpeta **Seed** dentro del directorio correspondiente al motivo en la carpeta **Motifs**.  

Ejemplo de ejecución:

```bash
python3 alphafold_parameters.py DOC_MAPK_JIP1_4  
```

Este paso es necesario para calcular métricas de **Energía de Interacción** y **Confianza**, y seleccionar los modelos óptimos.  

