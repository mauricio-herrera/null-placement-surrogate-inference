# PROTOCOLO CONGELADO — STAGE 2 MECHANISM AUDIT v0.30

**Fecha:** 2026-08-08  
**Estado:** congelado después del Stage-1 pilot v0.21 y antes de observar resultados Stage 2.  
**Objeto:** identificar qué parte de la discrepancia entre el null anual y el null mensual puede atribuirse a (i) block-conditioning, (ii) fase de agregación, y (iii) forma del observable derivado.

## 1. Lo que NO se modifica

- No se toca Amazon.
- No se modifica el threshold de 3.3 mm/day de Amazon.
- No se recalcula el gate Amazon con otros bloques.
- No se añaden familias, phi o lambda en respuesta a resultados del pilot.
- No se modifica q = 0.90.
- Ferro–Segers sigue siendo el estimador primario.

## 2. Malla congelada

Familias:

- `gaussian`
- `periodic`
- `tinnov`

Dependencia:

\[
\phi\in\{0.2,0.5,0.8\}.
\]

Concentración:

\[
\lambda\in\{0,0.5,1,2,4\}.
\]

Stage-2 mechanism pilot:

- 12 réplicas por celda;
- B = 30;
- 45 celdas;
- 540 trayectorias generadoras.

La simulación usa 252 años mensuales. Cada una de las 12 fases produce exactamente 251 ventanas no solapadas de 12 meses sin circular wrapping.

## 3. Auditoría de block-conditioning

Se comparan cuatro contratos del null mensual:

\[
25,\ 50,\ 100,\ 252\ \text{años}.
\]

`252` corresponde a un único grupo por mes calendario durante todo el registro; por tanto elimina el condicionamiento de frecuencia a bloques multidecadales manteniendo la marginal por mes calendario y el conteo total de eventos por mes.

El contrato de 50 años sigue siendo el primario para comparación con el diseño que motivó el estudio.

### Explicación competidora H_block

Parte de la reducción observada al pasar del null anual al mensual podría originarse en la conservación exacta de frecuencias dentro de bloques largos, y no en thresholding + aggregation por sí solos.

Se reportará cómo cambian:

- `native_rejection_rate`;
- `mechanistic_misattribution_rate`;
- `median_null_gap`;
- `delta_native`;

a medida que el bloque pasa 25 -> 50 -> 100 -> 252 años.

No se elegirá retrospectivamente un bloque “preferido”.

## 4. Auditoría de fase

Se evalúan las 12 fases posibles del operador de agregación anual:

\[
s=0,1,\ldots,11.
\]

Cada fase agrupa 12 meses consecutivos sin superposición. El análisis principal de fase se reporta bajo el contrato mensual congelado de 50 años.

Se cuantificará:

- rango de `delta_native` entre fases;
- desviación estándar entre fases;
- rango de `null_gap`;
- tasa de misattribution por fase.

Una dependencia sistemática respecto de la fase se interpretará como evidencia de que la ubicación del corte temporal respecto del ciclo estacional forma parte del mecanismo.

No se afirmará que el corte anual cree dependencia en un proceso mensual independiente.

## 5. Ventana móvil de 12 meses

Se calcula adicionalmente el mismo observable para todas las ventanas móviles solapadas de 12 meses.

**Este cálculo es descriptivo únicamente.** Las ventanas solapadas inducen dependencia por construcción, por lo que su índice extremal no se tratará como un contraste inferencial equivalente al de las ventanas anuales no solapadas.

Su función es mostrar cómo cambia la organización del observable cuando se elimina una única frontera anual fija.

## 6. Observables

Primario:

\[
C_y=\text{número de meses-evento dentro de cada ventana de 12 meses}.
\]

Secundario mecanístico 1:

\[
L_y=\text{máxima longitud de una racha dentro de la ventana}.
\]

Secundario mecanístico 2:

\[
R_y=\text{número de rachas distintas dentro de la ventana}.
\]

Los headline claims siguen dependiendo de `count`. `maxrun` y `nruns` se utilizan para discriminar persistencia de fragmentación.

## 7. Null anual y null mensual

Para cada fase y observable:

- se calcula el índice extremal observado después del detrending cúbico;
- se genera un IAAFT directamente sobre el observable anual;
- se genera el constrained monthly IAAFT sobre el proceso mensual para cada contrato de bloque;
- el surrogate mensual se transforma después mediante exactamente el mismo threshold y operador de agregación.

La semilla del surrogate mensual no depende de fase ni observable: la misma realización nativa se proyecta a todas las transformaciones. Esto mantiene una comparación pareada limpia.

## 8. Validaciones obligatorias

Para cada contrato de bloque se exige:

- `native_count_exact_fraction = 1`;
- `native_group_count_exact_fraction = 1`;
- errores espectrales finitos y pequeños;
- ausencia de NaN sistemáticos;
- calibración threshold de Gaussian/Periodic a precisión numérica.

Con B=30 y regla `p < 0.05`, el menor p Monte Carlo es 1/31 = 0.03226; el tamaño discreto alcanzable bajo exchangeability es por tanto 1/31, no exactamente 0.05.

## 9. Outputs obligatorios

- `phase_observed_and_annual.csv`
- `native_block_results.csv`
- `rolling_diagnostics.csv`
- `phase_summary_block50.csv`
- `block_summary_phase0.csv`
- `mechanism_overview_count_phase0.csv`
- `phase_overview_count_block50.csv`
- `phase_dispersion_count_block50.csv`
- `rolling_summary.csv`
- `RUN_METADATA.json`

## 10. Regla de interpretación

Stage 2 no confirma aún el paper. Debe decidir cuál de estas explicaciones permanece plausible:

1. threshold + aggregation;
2. block-conditioning;
3. fase de la frontera anual;
4. discretización/forma del observable;
5. interacción de los anteriores.

Solo después de esta auditoría se congelará el diseño confirmatorio y la curva de potencia.

## 11. Regla de detención

Si la discrepancia anual-vs-nativa desaparece al eliminar block-conditioning (252 años), la explicación por bloques pasa al centro y el claim de coarse graining deberá reducirse.

Si la discrepancia persiste en el contrato 252 años y además varía sistemáticamente con la concentración y/o fase, se justifica una fase confirmatoria específicamente diseñada para separar esos mecanismos.

No se realizarán nuevos cálculos Amazon en ninguno de los dos casos.
