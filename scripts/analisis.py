# analisis.py
# Anabel Yu Flores Moral
# Herramientas y Algoritmos en Bioinformática
# 2025-11-05

# Importaciones necesarias para el análisis funcional:
import gseapy as gp  # Librería para análisis de enriquecimiento génico con Enrichr.
import pandas as pd  # Para manipulación y exportación de resultados en DataFrames.
import numpy as np   # Para cálculos numéricos en métricas de significancia (e.g., log p-valores).
import matplotlib.pyplot as plt  # Para generar visualizaciones de términos enriquecidos.
import os  # Para crear directorios y manejar paths de salida.

def analisis_go_ora(genes: list[str], salida: str = "resultados_gseapy", umbral_fdr: float = 0.05):
    """
    Realiza un análisis de sobrerrepresentación génica (ORA) para identificar procesos
    biológicos, funciones moleculares y componentes celulares enriquecidos en una lista
    de genes humanos, utilizando Enrichr de la librería GSEApy.

    Método:
        ORA compara la lista de genes contra conjuntos de anotaciones de Gene Ontology (GO)
        para detectar términos sobrerrepresentados estadísticamente (p-valor ajustado < umbral).
        Adecuado para genes mitocondriales como COX4I2 (subunidad citocromo c oxidasa),
        ND1 y ATP6 (complejo NADH deshidrogenasa y ATP sintasa), evaluando vías como
        fosforilación oxidativa.

    Bases de datos utilizadas:
        - GO_Biological_Process_2021: Procesos biológicos (e.g., transporte de electrones).
        - GO_Molecular_Function_2021: Funciones moleculares (e.g., actividad ATPasa).
        - GO_Cellular_Component_2021: Componentes celulares (e.g., membrana mitocondrial interna).
        Fuente: Enrichr (integración de GO Consortium, actualizado 2021).

    Parámetros:
        genes (list[str]): Lista de símbolos génicos en formato HUGO (e.g., ["COX4I2", "ND1", "ATP6"]).
        salida (str): Directorio de salida para CSV de resultados y PNG del gráfico (default: "resultados_gseapy").
        umbral_fdr (float): Umbral para p-valor ajustado por FDR (False Discovery Rate); default 0.05.

    Devuelve:
        pd.DataFrame: Tabla con términos GO enriquecidos, incluyendo Term, Adjusted P-value, Genes y Odds Ratio.
        Si no hay términos significativos, devuelve el DataFrame vacío.
    """
    # Crear directorio de salida si no existe, para organizar archivos generados.
    os.makedirs(salida, exist_ok=True)
    print(f"Ejecutando ORA en {len(genes)} genes usando Enrichr...")

    # Ejecutar análisis ORA: Enriquecimiento de genes contra bases GO via Enrichr API.
    enr = gp.enrichr(
        gene_list=genes,  # Lista de símbolos HUGO (COX4I2, ND1, ATP6).
        gene_sets=["GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"],
        # Subontologías GO para procesos, funciones y componentes celulares.
        organism="Human",  # Especifica Homo sapiens para anotaciones humanas/mitocondriales.
        cutoff=umbral_fdr,  # Filtra resultados por p-valor ajustado (FDR < 0.05 por defecto).
    )

    # Mostrar resultados
    df = enr.results
    if df.empty:
        print("No se encontraron términos GO significativos.")
        return df

    # Guardar resultados
    archivo_csv = os.path.join(salida, "resultados_enrichr.csv")
    df.to_csv(archivo_csv, index=False)
    print(f"Resultados guardados en: {archivo_csv}")

# --------- REPRESENTACIÓN GRÁFICA ---------

    # Top 15 términos por p-valor ajustado.
    df_plot = df.sort_values("Adjusted P-value").head(15).iloc[::-1].copy()
    df_plot['Gene_Count'] = df_plot['Genes'].apply(lambda x: len(x.split(';')) if pd.notna(x) else 0)

    plt.figure(figsize=(8, 10))
    y_pos = np.arange(len(df_plot))

    # Dot plot: x = -log10(p-valor), y = términos, size = conteo de genes.
    scatter = plt.scatter(
        -np.log10(df_plot["Adjusted P-value"]),
        y_pos,
        s=df_plot['Gene_Count'] * 50,
        c=df_plot["Adjusted P-value"],
        cmap='Blues_r',
        alpha=0.7
    )

    plt.yticks(y_pos, df_plot["Term"])  # Términos en y.
    plt.xlabel("-log10(Adjusted P-value)", fontsize=12)  # Significancia en x.
    plt.title("Términos GO más enriquecidos (Dot Plot)", fontsize=14)
    plt.grid(axis='x', alpha=0.3)
    plt.gca().invert_yaxis()  # Invertir y para que top sea el más significativo.

    # Añadir colorbar para p-valor.
    plt.colorbar(scatter, label="Adjusted P-value")

    plt.tight_layout()
    plt.savefig(os.path.join(salida, "grafico_enriquecimiento_go_dotplot.png"), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Gráfico guardado en: {os.path.join(salida, 'grafico_enriquecimiento_go.png')}")
    return df

# --------- EJECUCIÓN ---------

genes = ["COX4I2", "ND1", "ATP6"]
resultados = analisis_go_ora(genes)
