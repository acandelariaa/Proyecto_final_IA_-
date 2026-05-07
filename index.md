# Análisis de Clustering en Cáncer de Páncreas

## Identificación de subtipos moleculares mediante aprendizaje no supervisado

> Proyecto Parcial 3 — Inteligencia Artificial I  
> Datos: TCGA-PAAD · The Cancer Genome Atlas

---

## ¿De qué trata este proyecto?

El **adenocarcinoma ductal pancreático (PDAC)** tiene una tasa de supervivencia a 5 años de apenas el 12%. Una de las razones principales es la enorme heterogeneidad molecular entre pacientes: dos personas con el mismo diagnóstico pueden tener perfiles genéticos completamente distintos y responder de forma diferente al tratamiento.

Este proyecto aplica técnicas de **aprendizaje no supervisado** sobre datos de expresión génica del TCGA-PAAD para descubrir si existen subgrupos moleculares con diferencias clínicas o de supervivencia reales.

---

## Estructura del proyecto

| Sección | Descripción |
|---|---|
| [1. Planteamiento y datos](01_planteamiento.md) | Contexto clínico, fuente de datos y metodología general |
| [2. Exploración y preprocesamiento](02_preprocesamiento.md) | Limpieza, valores faltantes, outliers y transformación de variables |
| [3. Reducción de dimensionalidad (PCA)](03_pca.md) | Aplicación de PCA a los 60,660 genes del dataset TPM |
| [4. K-Means Clustering](04_kmeans.md) | Determinación de K óptimo, clustering y evaluación |
| [5. Clustering Jerárquico](05_jerarquico.md) | Dendrograma, corte y comparación con K-Means |
| [6. Análisis de genes y supervivencia](06_analisis_genes.md) | Genes diferenciales, curvas Kaplan-Meier y scorecard comparativo |
| [7. Conclusiones](07_conclusion.md) | Hallazgos, limitaciones y reflexión final |
| [EXTRA: Filtrado MAD](08_extra_mad.md) | ¿Mejora el clustering si filtramos genes por varianza? |

---

## Resultado principal

El **Clustering Jerárquico con 4 clusters** fue el método más efectivo, logrando separar grupos con hasta **1,278 días de diferencia en mediana de supervivencia** (~3.5 años) y explicando el 27% de la variabilidad en tiempo de sobrevida. K-Means con 2 clusters capturó apenas el 8%.

El **procedimiento extra** demostró que filtrar genes por MAD (Median Absolute Deviation) antes de aplicar los modelos mejora el K-Means en ~23 puntos porcentuales.

---

## Datos utilizados

Todos los datos provienen de **TCGA-PAAD** vía [UCSC Xena](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Pancreatic%20Cancer%20(PAAD)):

| Archivo | Contenido |
|---|---|
| `TCGA-PAAD.star_tpm.tsv.gz` | Expresión génica (RNA-seq, 60,660 genes) |
| `TCGA-PAAD.clinical.tsv.gz` | Variables clínicas (estadio, género, edad, etc.) |
| `TCGA-PAAD.survival.tsv.gz` | Tiempo de seguimiento y estado vital |

---

## Stack tecnológico

`Python` · `pandas` · `scikit-learn` · `scipy` · `matplotlib` · `seaborn` · `lifelines`
