# 7. Conclusión y Reflexión Final

[← Análisis de Genes](06_analisis_genes.md) | [Siguiente: EXTRA — Filtrado MAD →](08_extra_mad.md)

---

## ¿Qué encontramos?

Aplicamos K-Means y Clustering Jerárquico sobre datos de expresión génica de pacientes con cáncer de páncreas (TCGA-PAAD), buscando subgrupos moleculares que pudieran relacionarse con diferencias en supervivencia.

El resultado más concreto: el **Clustering Jerárquico ganó**. Con 4 clusters logró separar grupos con hasta **1,278 días de diferencia en la mediana de supervivencia** — aproximadamente 3.5 años — explicando el 27% de la variabilidad en OS.time. K-Means, con solo 2 grupos, capturó apenas el 8%. En una enfermedad donde la supervivencia global a 5 años es del 12%, esa diferencia no es menor.

Dicho eso, las correlaciones entre los genes encontrados y el tiempo de supervivencia fueron bajas en ambos métodos (~0.15 y ~0.12 en promedio). Eso no necesariamente significa que los genes estén mal elegidos — en biología del cáncer, el pronóstico rara vez depende de un solo gen. Lo que sí indica es que estos genes deben interpretarse como parte de un **programa transcripcional más amplio**, no como biomarcadores individuales listos para usar.

---

## Lo que aprendimos de los datos

El cáncer de páncreas es genuinamente heterogéneo, y los datos lo reflejan. El scatter de PCA mostró una nube compacta sin separaciones obvias, lo que desde el principio anticipó que el clustering no iba a ser limpio. Aun así, al darle más granularidad al modelo (4 clusters en lugar de 2), la señal de supervivencia apareció con más claridad.

Eso sugiere que **forzar el problema a solo 2 subtipos pierde información real**.

---

## Limitaciones Honestas

- Trabajamos con ~183 muestras y más de 60,000 genes. Es un escenario de **alta dimensionalidad con poca muestra**, lo que hace frágil cualquier conclusión.
- La reducción a 10 componentes principales fue una decisión pragmática, no óptima. Es posible que información biológicamente importante haya quedado descartada.
- No validamos los clusters contra subtipos moleculares ya descritos en la literatura (como los subtipos *Classical* y *Basal-like* del PDAC). Sin eso, no podemos afirmar con certeza que los grupos que encontramos corresponden a entidades biológicas reales.
- No somos médicos ni genetistas, y eso importa: hay decisiones de preprocesamiento y selección de genes que alguien con dominio clínico hubiera tomado de forma más informada.

---

## ¿Qué se podría mejorar?

Tres cosas concretas que harían este análisis más sólido:

**1. Filtrar genes antes de clustering usando varianza (MAD)**  
En lugar de usar los 60,000 genes crudos, seleccionar los más variables reduce el ruido de entrada y mejora la calidad de los clusters. El [análisis extra](08_extra_mad.md) demuestra exactamente esto.

**2. Análisis de vías biológicas (Pathway Analysis)**  
Herramientas como GSEA o Enrichr convertirían la lista de genes en conocimiento interpretable: ¿están relacionados con metabolismo? ¿con respuesta inmune? ¿con señalización KRAS? Una lista de identificadores Ensembl no comunica nada clínico por sí sola.

**3. Validación externa**  
Reproducir el análisis en una cohorte independiente (por ejemplo, ICGC) para saber si los subtipos encontrados son reproducibles o específicos de este dataset.

---

## Reflexión Final

Este proyecto fue, en el fondo, un ejercicio de **humildad computacional**. Los datos de expresión génica del cáncer de páncreas son ruidosos, complejos y difíciles de interpretar sin contexto clínico profundo. Aun así, con herramientas de aprendizaje no supervisado fue posible detectar estructura latente que tiene correlación real con el pronóstico de los pacientes.

No encontramos biomarcadores listos para la clínica, pero sí confirmamos algo importante: **la biología del PDAC no cabe en dos grupos**. La granularidad importa, y explorarla con más rigor — más muestras, más integración de datos, más validación — es exactamente hacia donde apunta la oncología de precisión.

---

| | K-Means | Jerárquico |
|---|---|---|
| Clusters | 2 | 4 |
| % Varianza OS | 8.26% | **27.08%** |
| Diff. mediana OS.time | 214 días | **1,278 días** |
| Score | 0/3 | **3/3** |
----
# Comparación con Literatura Científica

---

## Artículo de Referencia

> Ji, H., Zhang, Q., Wang, X.-X., Li, J., Wang, X., Pan, W., Zhang, Z., Ma, B., & Zhang, H.-M. (2022).
> **Identification of stromal microenvironment characteristics and key molecular mining in pancreatic cancer.**
> *Discover Oncology*, 13, 83. https://doi.org/10.1007/s12672-022-00532-y

---

## ¿Por qué este artículo?

Trabaja exactamente con los mismos datos (TCGA-PAAD), aplica clustering jerárquico sobre expresión genetica y busca subgrupos con diferencias en supervivencia. La diferencia clave es el enfoque: ellos parten de conocimiento biológico previo sobre la **matriz extracelular (ECM)**; nosotros partimos de los datos crudos sin suposiciones temáticas.

---

## ¿Qué hicieron ellos?

Seleccionaron 14 vías biológicas relacionadas con la ECM, aplicaron **GSVA** para puntuar cada vía por paciente, y realizaron clustering jerárquico obteniendo **4 grupos** con diferencias significativas en supervivencia (p = 0.0021). A partir de ahí construyeron el modelo **PECMS** con 8 genes usando regresión Lasso, logrando un AUC de 0.90–0.93, validado en una cohorte independiente (CPTAC-3) y en 20 pacientes clínicos propios.

---

## Similitudes

El hallazgo más directo: **ambos estudios convergen en 4 clusters como número óptimo** para TCGA-PAAD. En el paper, el mejor grupo no alcanzó mediana de supervivencia y el peor tuvo 418 días. En nuestro caso, el clustering jerárquico con 4 grupos produjo una diferencia de hasta **1,278 días entre medianas**, explicando el 27% de la variabilidad en OS.time — frente al 8% de K-Means con 2 grupos.

Ambos estudios confirman lo mismo: **la biología del PDAC no se simplifica bien a 2 subtipos**, y la granularidad importa.

- **Interpretación biológica**: identificaron a **KLHL32** como gen protector y a **SLC2A1 (GLUT1)** como marcador de peor pronóstico, con validación en tejido (inmunohistoquímica). Nosotros detectamos genes diferenciales pero no llegamos a interpretar su función.
- **Validación externa**: confirmaron sus resultados en dos cohortes adicionales. aqui se trabajó solo con TCGA-PAAD.

---

## En resumen

El artículo funciona como una validación indirecta de que la señal que encontramos es real. Al mismo tiempo, muestra el camino natural hacia un análisis más robusto: incorporar conocimiento biológico en la selección de genes, análisis de vías (GSEA, Enrichr), y validación en datos independientes.



---

[← Análisis de Genes](06_analisis_genes.md) | [Siguiente: EXTRA — Filtrado MAD →](08_extra_mad.md)




