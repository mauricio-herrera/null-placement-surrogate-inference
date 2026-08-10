# PROTOCOLO CONGELADO — PAPER A, STAGE 5
## Mitigation and independent closure audit v0.60

**Fecha de congelación:** 2026-08-09  
**Versión:** v0.60  
**Estado:** congelado después del Stage 4 y antes de ejecutar Stage 5.  
**Regla de cierre:** Stage 5 es la última etapa computacional planificada para Paper A. No se autoriza Stage 6 salvo fallo técnico demostrado.

---

## 1. Objetivo

Stage 5 tiene tres módulos y ninguna finalidad de rescate estadístico:

- **5A — estimator-offset diagnostic:** determinar cuánto del offset finito de \(\Delta\theta\) puede asociarse al techo \(\hat\theta\le 1\), cuánto persiste con el estimador no truncado y cuánto puede quedar asociado a sesgo finito/conditioning del surrogate.
- **5B — percentile-threshold mitigation:** comprobar si igualar la probabilidad mensual de evento reduce la *mechanistic misattribution* observada bajo thresholds absolutos fuertemente estacionales, y cuantificar la degradación cuando los percentiles se estiman desde períodos de referencia finitos.
- **5C — independent implementation:** verificar de manera independiente los números titulares del Stage 3 y reproducir la dirección inferencial sobre una submuestra congelada con un segundo pipeline.

Amazon permanece completamente congelado.

---

# MÓDULO 5A — ESTIMATOR-OFFSET DIAGNOSTIC

## 2. Estimador raw y truncado

Para cada serie anual se calculará el estimador de intervalos de Ferro–Segers antes del truncamiento superior:

\[
\hat\theta_{\rm raw}=\frac{\text{num}}{\text{den}},
\]

y la versión usada hasta ahora:

\[
\hat\theta=\min\{1,\max(0,\hat\theta_{\rm raw})\}.
\]

Se registrará por trayectoria y contrato:

- \(\hat\theta_{\rm obs,raw}\);
- \(\hat\theta_{\rm obs}\);
- mediana raw y truncada de cada null;
- fracción de surrogates con \(\hat\theta_{\rm raw}\ge1\);
- offset raw y truncado;
- indicador de trayectoria completamente libre de ceiling: observación raw < 1 y todos los surrogates raw < 1.

## 3. Estadístico no truncado auxiliar

Se calculará además la longitud media de rachas consecutivas de excedencias por encima del cuantil anual \(q=0.90\). Este estadístico será diagnóstico, no un nuevo outcome confirmatorio.

## 4. Interpretación permitida

Se podrá afirmar que el ceiling contribuye al offset solo si el desplazamiento disminuye materialmente al usar \(\hat\theta_{\rm raw}\) y/o si existe una diferencia marcada en la masa raw por encima de uno.

No se atribuirá automáticamente el offset residual a binarización. Puede contener sesgo de muestra finita y conditioning del surrogate.

---

# MÓDULO 5B — PERCENTILE-THRESHOLD MITIGATION

## 5. Generadores

Se mantienen las nueve combinaciones congeladas:

\[
\{\text{Gaussian},\text{Periodic innovation},\text{Student-t innovation}\}
\times
\{\phi=0.2,0.5,0.8\}.
\]

Cada trayectoria tiene 252 años mensuales y produce 251 índices anuales no solapados de fase 0.

## 6. Réplicas y surrogates

Por celda familia × \(\phi\):

\[
N=150.
\]

Total de trayectorias generadoras:

\[
9\times150=1350.
\]

Por null:

\[
B=149.
\]

Con la regla estricta \(p<0.05\), el tamaño Monte Carlo alcanzable es

\[
7/150=0.046\overline{6}.
\]

## 7. Cinco contratos de threshold

Todos se aplican a la **misma trayectoria de evaluación**.

### T0 — absolute_lambda4
Threshold físico sintético con concentración estacional \(\lambda=4\) y \(\sum_m p_m=3\).

### T1 — oracle_percentile
Probabilidad mensual uniforme:

\[
p_m=3/12=0.25,\qquad m=1,\ldots,12.
\]

### T2 — estimated_percentile_15y
Percentil 25 mensual estimado desde una serie de calibración independiente de 15 años.

### T3 — estimated_percentile_30y (primario aplicado)
Percentil 25 mensual estimado desde una serie de calibración independiente de 30 años.

### T4 — estimated_percentile_60y
Percentil 25 mensual estimado desde una serie de calibración independiente de 60 años.

Las series de calibración son independientes de la trayectoria de evaluación.

## 8. Gate de consistencia P0

En la misma trayectoria se construirán también los thresholds definidos por \(\lambda=0\). Por el lema de rangos/cópula:

\[
\lambda=0\quad\Longleftrightarrow\quad p_m=0.25\ \forall m.
\]

El vector indicador producido por T1 debe coincidir **exactamente** con el vector indicador producido por la implementación \(\lambda=0\).

Gate técnico:

\[
\text{indicator identity fraction}=1.0.
\]

## 9. Nulls

Para cada contrato se comparan:

1. **index-resolution annual IAAFT** sobre el conteo anual detrendido;
2. **native-resolution constrained monthly IAAFT**, con block = 252 años (grupos por mes calendario en todo el registro), seguido por exactamente el mismo threshold y agregación.

Los mismos surrogates mensuales nativos se reutilizan entre los cinco contratos de threshold de una misma trayectoria. Esto preserva el emparejamiento y reduce costo computacional sin alterar la hipótesis.

Outcome primario:

\[
M=1\{p_{\rm annual}<0.05,\ p_{\rm native}\ge0.05\},
\]

la *mechanistic misattribution*.

## 10. Hipótesis congeladas

### H_P0 — consistency
T1 reproduce exactamente la rama \(\lambda=0\) a nivel del indicador.

### H_P1 — oracle seasonal mitigation (primaria)

\[
P(M=1\mid T1)<P(M=1\mid T0).
\]

Se usará comparación pareada exacta unilateral sobre discordancias T0-only versus T1-only.

PASS si:

- diferencia de tasa T0 − T1 > 0; y
- \(p_{\rm paired}<0.05\).

### H_P2 — 30-year estimated mitigation (secundaria confirmatoria)

\[
P(M=1\mid T3)<P(M=1\mid T0).
\]

Mismo criterio pareado unilateral.

### H_P3 — calibration-length behavior (secundaria descriptiva)
Se reportarán, para 15/30/60 años:

- RMSE del vector de probabilidades reales respecto de 0.25;
- rango \(\max p_m-\min p_m\);
- \(N_{\rm eff,p}\);
- mechanistic misattribution.

No se exige monotonicidad perfecta como gate confirmatorio.

## 11. Claims permitidos

Si H_P1 pasa:

> Equalizing event probability across the seasonal cycle substantially reduces the mechanistic misattribution associated with a strongly concentrated absolute threshold.

Si H_P2 pasa:

> The mitigation remains present when monthly percentiles are estimated from an independent 30-year reference period.

No se podrá afirmar:

- que percentiles eliminan universalmente el problema;
- que el residuo es causado exclusivamente por binarización;
- que todos los índices percentílicos climatológicos están protegidos;
- que T1 debe tener exactamente 5% de rechazo anual.

---

# MÓDULO 5C — IMPLEMENTACIÓN INDEPENDIENTE

## 12. Submuestra congelada

Se seleccionan **10 replicates por cada una de las 45 celdas** de Stage 3 mediante ranking SHA256 del identificador `(family, phi, lambda, replicate)` con namespace fijo `paperA_stage5c_selection_v060`.

Total:

\[
450\ \text{trayectorias}.
\]

La selección se guarda antes de ejecutar el pipeline B en `STAGE5C_SELECTION_V060.csv` y no usa ningún resultado estadístico.

## 13. Pipeline B

El script independiente:

- no importa funciones de Stage 3;
- regenera las trayectorias desde la especificación y semillas congeladas de Stage 3;
- implementa por separado la agregación mensual→anual;
- implementa por separado detrending cúbico;
- implementa Ferro–Segers raw/truncado;
- implementa IAAFT anual y constrained monthly IAAFT con estructura de código distinta;
- usa semillas de surrogate propias de Stage 5C;
- usa \(B=99\) para el rerun independiente de inferencia.

## 14. Verificaciones exactas y estadísticas

### C1 — deterministic theta identity
Sobre las 450 trayectorias seleccionadas, el \(\hat\theta_{obs}\) regenerado debe coincidir con `theta_obs` de Stage 3 dentro de tolerancia:

\[
\max |\Delta\hat\theta|\le10^{-10}.
\]

### C2 — exact recomputation from frozen CSV
Un agregador independiente debe reproducir desde `replicate_results.csv`, para `count`, block=252:

- annual rejections = 1697;
- native rejections = 629;
- annual-only = 1084;
- native-only = 16.

### C3 — independent surrogate direction
En la submuestra de 450 trayectorias con \(B=99\):

- native rejection rate \(\le0.07\);
- annual rejection rate > native rejection rate;
- annual-only > native-only;
- paired rate difference \(\ge0.03\).

C3 verifica robustez de dirección, no identidad exacta de p-values con Stage 3.

---

# 15. Diagnósticos técnicos comunes

Se exigirán:

- ausencia de NaN masivos;
- exactitud de conteos mensuales bajo native surrogates = 1.0;
- fracción válida de \(\theta\) ≥ 0.99;
- fidelidad espectral comparable a stages previos;
- hashes del protocolo, scripts y selección registrados antes del run completo.

---

# 16. Regla final de cierre

Después de Stage 5:

- **no** Stage 6;
- **no** nuevos valores de \(\lambda\);
- **no** nuevos \(\kappa\);
- **no** segundo índice climático antes de redactar;
- **no** cambios en Amazon;
- **no** escalera retrospectiva de nulls.

El siguiente artefacto será la matriz final `claim–evidence–limitation`, seguida del manuscrito.
