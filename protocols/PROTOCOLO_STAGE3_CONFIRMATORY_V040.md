# PROTOCOLO CONGELADO — PAPER A, STAGE 3
## Focused confirmatory synthetic experiment v0.40

**Fecha de congelación:** 2026-08-08  
**Versión:** v0.40  
**Estado:** congelado después de Stage 2 y antes del primer cálculo confirmatorio  
**Amazon:** permanece completamente congelado; Stage 3 es exclusivamente sintético.

---

## 1. Motivación

Stage 1 mostró una discrepancia reproducible entre un null construido directamente sobre el índice anual y un null construido en la resolución mensual nativa. Stage 2 mostró que:

1. la discrepancia no desaparece al eliminar el block-conditioning de 50 años;
2. aumenta en el régimen de fuerte concentración estacional;
3. la fase anual la modula, pero no la explica por sí sola;
4. `count` y `maxrun` muestran el patrón con buena fidelidad espectral del surrogate anual;
5. `nruns` produce discrepancias muy grandes, pero con peor fidelidad IAAFT anual, por lo que no será outcome confirmatorio principal.

Stage 3 deja de explorar mecanismos nuevos. Su finalidad es **confirmar o refutar** el efecto principal con mucha mayor precisión Monte Carlo y muchas más trayectorias generadoras independientes.

---

## 2. Claim confirmatorio principal

En series finitas generadas por procesos mensuales de memoria corta conocida, un null IAAFT construido directamente sobre el índice anual `count` produce más rechazos que un null construido en resolución mensual nativa y posteriormente sometido al mismo operador threshold + agregación.

Este claim se refiere a una **discrepancia inferencial de muestra finita**. No afirma que el proceso mensual subyacente tenga índice extremal asintótico menor que uno ni que la agregación “cree memoria física”.

---

## 3. Distinción conceptual

Un rechazo del null anual puede ser correcto para la pregunta estadística definida sobre el índice anual. El evento de interés mecanístico es:

\[
M_i = 1\{\text{annual reject}\}\,1\{\text{native does not reject}\}.
\]

En los generadores de Stage 3 sabemos por construcción que no existe memoria fraccional ni memoria no-Markoviana. Por ello este evento se denomina **mechanistic misattribution event**, no “false positive” estadístico en sentido genérico.

---

## 4. Diseño congelado

### Familias de dependencia

1. Gaussian AR(1)
2. periodic-innovation Gaussian AR(1)
3. Student-t innovation AR(1), \(\nu=5\)

### Persistencia

\[
\phi\in\{0.2,0.5,0.8\}.
\]

### Concentración estacional

\[
\lambda\in\{0,0.5,1,2,4\}.
\]

La suma esperada de meses-evento permanece fija:

\[
\sum_{m=1}^{12}p_m=3.
\]

### Longitud

Se simulan 252 años mensuales. El outcome anual de fase 0 utiliza los primeros 251 años, preservando exactamente el contrato de Stage 2.

### Réplicas

\[
N=300
\]

trayectorias generadoras independientes por celda.

Número total:

\[
3\times3\times5\times300=13\,500
\]

trayectorias.

### Surrogates

\[
B=249
\]

por null.

Con

\[
p=(1+r)/(B+1)
\]

y la regla congelada \(p<0.05\), el tamaño Monte Carlo alcanzable bajo exchangeability es:

\[
12/250=0.048.
\]

Esta elección aproxima mejor el 5% que Stage 1/2 sin requerir el costo de \(B=499\).

---

## 5. Outcome y fase

### Primario

`count`: número de meses-evento en una ventana anual no solapada.

### Secundario confirmatorio

`maxrun`: máxima racha de meses-evento dentro del año.

### Excluido de los headline claims

`nruns`: se excluye del Stage 3 confirmatorio porque Stage 2 mostró una fidelidad espectral IAAFT anual materialmente peor.

### Fase

Solo fase 0 es confirmatoria.

Las restantes fases ya cumplieron su función mecanística en Stage 2 y no se vuelven a explorar en Stage 3.

---

## 6. Contratos del native null

### Primario

\[
\text{block}=252\ \text{años},
\]

es decir, grupos por mes calendario sobre todo el registro, sin block-conditioning multidecadal adicional.

### Robustez prespecificada

\[
\text{block}=50\ \text{años}.
\]

No se probarán otros bloques en Stage 3.

---

## 7. Hipótesis H1 — primaria

Para `count`, block = 252, se define en cada trayectoria:

- \(A_i=1\): annual null rechaza;
- \(N_i=1\): native null rechaza.

Se cuentan los pares discordantes:

\[
n_{A\bar N}=\#\{A_i=1,N_i=0\},
\]

\[
n_{\bar A N}=\#\{A_i=0,N_i=1\}.
\]

### Test primario

Test binomial exacto unilateral equivalente al McNemar exacto condicional:

\[
H_0:P(A\bar N)\le P(\bar A N),
\]

\[
H_1:P(A\bar N)>P(\bar A N).
\]

Nivel:

\[
\alpha=0.05.
\]

### Criterio de éxito H1

H1 pasa si:

1. el test exacto unilateral da \(p<0.05\); y
2. la diferencia pareada de tasas es positiva.

H1 es el único test requerido para declarar éxito confirmatorio primario, siempre que los gates técnicos pasen.

---

## 8. Hipótesis H2 — concentración estacional

H2 se evalúa **jerárquicamente solo después de H1**.

Se compara el evento mechanistic misattribution entre:

\[
\lambda=4
\]

y

\[
\lambda=0.
\]

mediante Fisher exacto unilateral, agregando el diseño balanceado sobre familias y \(\phi\).

### Criterio de éxito H2

\[
p<0.05
\]

y diferencia de riesgo \(>0\).

Este contraste respalda la afirmación direccional de que el problema es mayor bajo fuerte concentración estacional. No se exige monotonicidad perfecta en todos los niveles intermedios.

---

## 9. Trend test secundario

Se reportará un Cochran-Armitage-style score test unilateral usando los scores congelados:

\[
\lambda=(0,0.5,1,2,4).
\]

Es evidencia secundaria. No reemplaza H2 ni se utilizará para rescatar H2 si el contraste \(\lambda=4\) vs \(0\) falla.

---

## 10. Robustez block=50

El mismo test pareado se reportará para block = 50.

No forma parte del gate primario porque Stage 2 ya estableció que el block-conditioning no es la explicación dominante. Su función es demostrar que la dirección no depende del contrato menos condicionado.

---

## 11. Gates técnicos prespecificados

Antes de interpretar H1/H2 deben satisfacerse simultáneamente:

### G1 — conservación exacta

\[
\texttt{native\_count\_exact\_fraction}=1
\]

y

\[
\texttt{native\_group\_count\_exact\_fraction}=1
\]

en todas las realizaciones.

### G2 — validez del índice extremal

Fracción válida mínima para el outcome primario:

\[
\ge0.99.
\]

### G3 — fidelidad espectral para `count`

Annual IAAFT:

- mediana \(\le0.005\);
- percentil 99 \(\le0.02\).

Native IAAFT:

- mediana \(\le0.01\);
- percentil 99 \(\le0.03\).

### G4 — native null no anticonservador

Para `count`, block = 252, la tasa pooled de rechazo nativo debe ser:

\[
\le0.065.
\]

El valor nominal Monte Carlo de referencia es 0.048. Una tasa inferior se interpreta como conservadora, no como invalidación del gate.

---

## 12. Jerarquía de decisión

### Stage 3 primary success

Requiere:

1. G1–G4 pasan;
2. H1 pasa.

### Stage 3 full directional success

Requiere además:

3. H2 pasa.

`maxrun`, block=50 y el trend test son secundarios y no rescatan un H1 fallido.

---

## 13. Regla frente a un resultado negativo

Si H1 falla con los gates técnicos válidos:

> se abandona el claim confirmatorio general de que el index-resolution null produce sistemáticamente mayor mechanistic misattribution en este benchmark.

No se modificarán:

- \(B\);
- número de réplicas;
- thresholds;
- \(\lambda\);
- \(\phi\);
- familias;
- fase;
- regla \(p<0.05\);
- estimador de índice extremal;
- block primario.

Si H1 pasa pero H2 falla:

> se retiene la discrepancia annual/native, pero se abandona el claim confirmatorio de que aumenta específicamente con fuerte concentración estacional.

---

## 14. Qué no se hará en Stage 3

- no se toca Amazon;
- no se exploran nuevas fases;
- no se cambian bloques;
- no se añaden familias;
- no se prueban nuevos thresholds;
- no se ejecuta todavía la curva de potencia con \(\kappa>0\);
- no se usa `nruns` como headline outcome.

---

## 15. Próxima etapa si Stage 3 pasa

Solo después de un Stage 3 exitoso se autoriza el estudio de potencia:

> **Stage 4 — power under prespecified injected alternatives**.

La potencia se reportará por mecanismo y no como función universal de \(\Delta\theta\).

---

## 16. Implementación y reproducibilidad

El script confirmatorio:

`paperA_stage3_confirmatory_v040.py`

usa:

- semillas deterministas por familia, \(\phi\), \(\lambda\), réplica y surrogate;
- escritura incremental por chunks;
- opción `--resume` para continuar un run interrumpido sin cambiar las realizaciones;
- metadata automática;
- tests confirmatorios automáticos;
- gates técnicos automáticos;
- archivo final `CONFIRMATORY_DECISION.json`.

La malla confirmatoria no puede modificarse desde la línea de comandos.
