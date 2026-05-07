# 1. Planteamiento del Problema y Contexto de los Datos

[← Inicio](index.md) | [Siguiente: Exploración y preprocesamiento →](02_preprocesamiento.md)

---

## El Cáncer de Páncreas: Un Problema de Alta Relevancia Clínica

El **adenocarcinoma ductal pancreático (PDAC)**, clasificado en TCGA como **PAAD** (*Pancreatic Adenocarcinoma*), representa una de las neoplasias malignas con peor pronóstico a nivel mundial. A pesar de los avances en oncología molecular, la tasa de supervivencia a 5 años sigue siendo de apenas el **12%**, principalmente debido a su diagnóstico tardío y a la falta de marcadores moleculares fiables para estratificar a los pacientes.

Desde una perspectiva molecular, el cáncer de páncreas es extraordinariamente heterogéneo: distintos pacientes con el mismo diagnóstico clínico pueden presentar perfiles de expresión génica completamente diferentes, lo que se traduce en respuestas distintas a los tratamientos disponibles. Identificar **subgrupos moleculares** dentro de esta enfermedad es, por lo tanto, un paso fundamental para avanzar hacia una medicina de precisión.

---

## Objetivo del Análisis

El objetivo central de este proyecto es aplicar técnicas de **aprendizaje no supervisado** para descubrir estructura latente en los datos de expresión génica del cáncer de páncreas. Específicamente, se busca:

1. **Reducir la dimensionalidad** de los datos mediante PCA para explorar la estructura global del dataset.
2. **Identificar subgrupos de pacientes** con perfiles moleculares similares mediante algoritmos de clustering (K-Means y Clustering Jerárquico).
3. **Caracterizar clínicamente** cada subgrupo encontrado, relacionándolo con variables como estadio de la enfermedad, supervivencia, género y edad.
4. **Identificar genes diferencialmente expresados** entre los clusters, que potencialmente puedan servir como biomarcadores o dianas terapéuticas.

La pregunta central que guía este análisis es:

> **¿Existen subtipos moleculares bien definidos en el cáncer de páncreas que se asocien con diferencias clínicas o de supervivencia significativas?**

---

## Fuente de los Datos

Los datos utilizados en este proyecto provienen de **The Cancer Genome Atlas (TCGA)**, accedidos a través de la plataforma **UCSC Xena** ([xenabrowser.net](https://xenabrowser.net)), específicamente de la cohorte:

> 📂 **GDC TCGA Pancreatic Cancer (PAAD)**  
> [https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Pancreatic%20Cancer%20(PAAD)](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Pancreatic%20Cancer%20(PAAD))

TCGA es un consorcio que ha caracterizado genómicamente más de 20,000 tumores primarios de 33 tipos de cáncer distintos. Los datos de la plataforma GDC (*Genomic Data Commons*) han sido re-analizados de forma uniforme utilizando el ensamblaje de referencia del genoma humano **hg38**, garantizando consistencia y reproducibilidad.

Se trabajó con **tres archivos principales**:

| Dataset | Descripción |
|---|---|
| `STAR - TPM` | Datos de expresión génica (RNA-seq), cuantificados en TPM (*Transcripts Per Million*) mediante el alineador STAR |
| `GDC Phenotype` | Variables clínicas del paciente: estadio, género, edad al diagnóstico, histología, entre otras |
| `Survival data` | Datos de supervivencia: tiempo de seguimiento y estado vital (vivo/muerto) para análisis de pronóstico |

> **Nota sobre los datos de expresión génica:** Los valores de TPM representan la abundancia relativa de transcritos, normalizando tanto por la longitud del gen como por la profundidad de secuenciación, lo que facilita la comparación entre muestras.

---

## Justificación del Enfoque

El análisis de expresión génica a escala genómica genera datos de **alta dimensionalidad** (decenas de miles de genes por muestra), lo que plantea desafíos computacionales y estadísticos considerables. Las técnicas de aprendizaje no supervisado permiten abordar este problema sin depender de etiquetas previas, revelando estructura que puede ser clínicamente relevante y que de otro modo permanecería oculta.

Estudios previos utilizando los mismos datos de TCGA-PAAD han identificado entre 2 y 4 subtipos moleculares con diferencias significativas en supervivencia y respuesta a quimioterapia, lo que valida la pertinencia del enfoque utilizado en este proyecto.

---

## Metodología General

El flujo del análisis sigue los siguientes pasos:

1. Exploración del conjunto de datos
2. Tratamiento de datos
3. Reducción de dimensionalidad (PCA)
4. Construcción y comparación de modelos de clustering
5. Interpretación de clusters / análisis de genes
6. Interpretación de clusters / análisis clínico
7. Comparación con un artículo científico
8. Conclusión

---

## Nota Preliminar

> Cabe recalcar que no somos genetistas y mucho menos médicos, de modo que no podemos explicar cada uno de los miles de genes con los que trabajamos — tendríamos un glosario demasiado grande. Trabajaremos sin ignorar el contexto de los genes, pero con lógica y un buen enfoque para tratar datos. Haremos el mejor esfuerzo posible.

---

[← Inicio](index.md) | [Siguiente: Exploración y preprocesamiento →](02_preprocesamiento.md)
