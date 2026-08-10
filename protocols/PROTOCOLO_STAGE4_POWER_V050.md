# PROTOCOLO CONGELADO — PAPER A, STAGE 4
## Power under prespecified injected alternatives v0.50

**Fecha de congelación:** 2026-08-09  
**Versión:** v0.50  
**Estado:** congelado después del PASS confirmatorio de Stage 3 y antes del cálculo Stage 4 del usuario.  
**Amazon:** permanece completamente congelado.

---

## 1. Pregunta de Stage 4

Stage 3 ya estableció que, bajo generadores mensuales de memoria corta conocida, el null construido directamente sobre el índice anual rechaza sustancialmente más que el null construido a resolución mensual nativa y pasado por el mismo operador threshold + agregación.

Stage 4 responde una pregunta diferente:

> **¿Qué potencia tiene el gate mensual nativo, con ocho unidades inferenciales y B = 500, para detectar formas y magnitudes prespecificadas de organización temporal adicional?**

El objetivo es interpretar la ausencia de significación del caso Amazon sin recalcular Amazon.

La potencia **no** se considerará una función universal de \(\Delta\theta\). Se reportará separadamente por mecanismo inyectado.

---

## 2. Diseño tipo “ocho modelos”

Cada réplica de potencia contiene exactamente:

\[
8
\]

unidades inferenciales independientes.

El universo sintético tiene las nueve combinaciones:

\[
\{\text{Gaussian},\text{Periodic},\text{Student-t}\}
\times
\{0.2,0.5,0.8\}.
\]

En cada cohorte se seleccionan **8 de las 9 combinaciones sin reemplazo**, mediante una semilla fijada por el número de réplica. La misma composición de cohorte se reutiliza para todos los mecanismos y niveles de intensidad de esa réplica.

Esto evita escoger arbitrariamente una única cohorte de ocho celdas y produce heterogeneidad entre unidades manteniendo exactamente ocho unidades inferenciales.

Stage 4 **no** reproduce las 24 ramas modelo–escenario de CMIP6. Deliberadamente trabaja en el nivel inferencial de ocho unidades, evitando tratar ramas que comparten historia como réplicas independientes.

---

## 3. Generador base

Se mantienen las tres familias auditadas en Stage 3:

1. Gaussian AR(1);
2. periodic-innovation Gaussian AR(1);
3. Student-t innovation AR(1), \(\nu=5\).

Persistencia:

\[
\phi\in\{0.2,0.5,0.8\}.
\]

Concentración estacional primaria:

\[
\lambda=4.
\]

Se conserva:

\[
\sum_{m=1}^{12}p_m=3.
\]

Para \(\lambda=4\), el número efectivo de meses contribuyentes está en el régimen fuertemente concentrado confirmado en Stage 3.

Longitud:

\[
252\times12
\]

meses simulados, produciendo 251 índices anuales de fase 0.

---

## 4. Gate mensual nativo

El único null utilizado en Stage 4 es el null nativo mensual auditado:

- block = 50 años;
- preservación por mes calendario × bloque;
- IAAFT constrained;
- reconstrucción del índice anual después del surrogate;
- detrending cúbico anual;
- índice extremal de Ferro–Segers;
- \(q=0.90\);
- \(B=500\).

Para una unidad \(i\):

\[
\Delta_i
=
\operatorname{median}(\theta_{i,\mathrm{null}})
-
\theta_{i,\mathrm{obs}}.
\]

La unidad tiene dirección positiva si:

\[
\Delta_i>0.
\]

---

## 5. Gate agregado de ocho unidades

El gate reproduce la estructura lógica usada para la interpretación de Amazon:

### Condición 1 — dirección consistente

\[
\#\{i:\Delta_i>0\}\ge6/8.
\]

### Condición 2 — contraste agregado

Se define:

\[
T_{\mathrm{obs}}
=
\frac{1}{8}\sum_{i=1}^8
\left[
\operatorname{median}(\theta_{i,\mathrm{null}})
-
\theta_{i,\mathrm{obs}}
\right].
\]

Para cada índice de surrogate \(b=1,\ldots,B\):

\[
T_b
=
\frac{1}{8}\sum_{i=1}^8
\left[
\operatorname{median}(\theta_{i,\mathrm{null}})
-
\theta_{i,b}^{\mathrm{null}}
\right].
\]

El valor Monte Carlo es:

\[
p_{\mathrm{agg}}
=
\frac{1+\#\{b:T_b\ge T_{\mathrm{obs}}\}}{B+1}.
\]

El gate completo pasa si y solo si:

\[
\boxed{
\#\{\Delta_i>0\}\ge6
\quad\text{y}\quad
p_{\mathrm{agg}}<0.05.
}
\]

Este contraste agregado es una definición prespecificada para Stage 4. Se usa para cuantificar potencia de un diseño de ocho unidades; no se presenta como reconstrucción retrospectiva exacta de cada detalle computacional del análisis Amazon.

---

## 6. Número de réplicas

Por punto del diseño:

\[
N_{\mathrm{cohort}}=150.
\]

Con ocho unidades por cohorte y \(B=500\), el costo está deliberadamente en el mismo orden que Stage 3.

Se reportarán intervalos de Wilson del 95 % para la potencia del gate.

---

## 7. Mecanismo A — persistent regime alignment

### Propósito

Representar organización por un estado latente persistente compartido entre meses, sin cambiar las marginales mensuales.

Se genera un score anual latente:

\[
R_y=0.85R_{y-1}+\sqrt{1-0.85^2}\,\varepsilon_y.
\]

Dentro de cada mes calendario, los valores originales se reordenan por rangos usando:

\[
S_{y,m}
=(1-\kappa_R)Z^{\mathrm{base}}_{y,m}
+
\kappa_R R_y.
\]

Los valores mensuales originales se asignan según el orden de \(S_{y,m}\).

Esto preserva exactamente:

- el multiset de valores de cada mes calendario;
- el número de cruces del threshold por mes calendario;
- la marginal estacional.

Intensidades congeladas:

\[
\kappa_R\in\{0,0.2,0.4,0.6,0.8\}.
\]

\(\kappa_R=0\) es el generador base sin inyección.

Este mecanismo puede modificar fuertemente la dependencia de segundo orden. Por ello no se espera necesariamente que el gate, que preserva el espectro, tenga la misma potencia que frente al mecanismo B. Esa diferencia es parte del resultado científico.

---

## 8. Mecanismo B — history feedback

### Propósito

Representar una organización history-dependent/endógena que actúa sobre la cópula temporal manteniendo exactamente las marginales mensuales.

A partir del indicador actual se construye:

\[
H_t
=
\sum_{\ell=1}^{24}w_\ell I_{t-\ell},
\]

con

\[
w_\ell
\propto
\exp(-\ell/6).
\]

Dentro de cada mes calendario, el score de rango se modifica como:

\[
S_{y,m}
=
Z^{\mathrm{base}}_{y,m}
-
\kappa_H\widetilde H_{y,m}.
\]

La operación se itera cuatro veces, remapeando siempre sobre el multiset mensual original.

Intensidades congeladas:

\[
\kappa_H\in\{0,0.25,0.50,0.75,1.00\}.
\]

También preserva exactamente:

- la marginal de cada mes calendario;
- el número de eventos por mes calendario.

---

## 9. Punto nulo compartido

Se calcula una única condición baseline:

\[
\kappa=0,
\]

compartida por las dos curvas de potencia.

Por tanto, Stage 4 contiene nueve puntos computacionales:

- 1 baseline;
- 4 niveles de persistent regime alignment;
- 4 niveles de history feedback.

---

## 10. Outcomes primarios y secundarios

### Primario

\[
\text{Power}_{8}
=
P(\text{gate completo pasa}).
\]

### Secundarios

Se reportarán:

1. potencia del contraste agregado sin el requisito 6/8;
2. frecuencia de cumplir solamente el requisito 6/8;
3. tasa de rechazo de una unidad individual;
4. mediana de \(\Delta\theta\) inducida;
5. número mediano de unidades positivas;
6. AC(1) del conteo anual;
7. correlación del indicador mensual a lag 1 y lag 12.

---

## 11. Referencia Amazon congelada

Se registra solamente como referencia descriptiva:

\[
\Delta\theta_{\mathrm{Amazon}}=0.085009.
\]

**Prohibido:** escoger o modificar la malla de \(\kappa\) para acercarse a 0.085.

Después de ejecutar Stage 4, para cada mecanismo se reportará únicamente:

> el punto de la malla prespecificada cuyo \(\Delta\theta\) mediano esté más cerca de 0.085 y la potencia observada en ese punto.

No se autoriza una interpolación universal \(\text{Power}=f(\Delta\theta)\).

---

## 12. Gates técnicos

Antes de interpretar potencia deben cumplirse:

### G1 — exactitud de la inyección

En todas las unidades:

- multiset mensual original = multiset mensual inyectado;
- conteo de eventos por mes calendario idéntico.

### G2 — exactitud del surrogate

\[
\min\texttt{native\_count\_exact\_fraction}=1,
\]

\[
\min\texttt{native\_group\_count\_exact\_fraction}=1.
\]

### G3 — validez

\[
\min\texttt{native\_null\_valid\_fraction}\ge0.99.
\]

### G4 — fidelidad espectral

Mediana del error espectral nativo:

\[
\le0.01,
\]

y percentil 99:

\[
\le0.03.
\]

### G5 — calibración single-unit del baseline

Con \(B=500\), la tasa alcanzable bajo \(p<0.05\) es aproximadamente:

\[
25/501=0.04990.
\]

Se exige:

\[
\text{baseline single-unit rejection}\le0.065.
\]

### G6 — calibración del gate completo

Se exige:

\[
\text{baseline full-gate false-positive rate}\le0.10.
\]

Este límite se fija antes del cálculo considerando solo 150 cohortes y la composición conjunta de dos condiciones del gate.

---

## 13. Lectura de potencia

La potencia se reporta separadamente para cada mecanismo.

Se identificarán descriptivamente los primeros puntos de la malla que alcancen:

\[
50\%
\]

y

\[
80\%
\]

de potencia, si existen.

No alcanzar 80 % dentro de la malla no autoriza extender retrospectivamente \(\kappa\).

---

## 14. Interpretación permitida respecto de Amazon

Después de Stage 4 se podrá escribir, dependiendo de los resultados, una frase del tipo:

> Under a persistent-regime alternative producing a deficit comparable to the frozen Amazon value, the eight-unit native-resolution gate had X% power, whereas under a history-feedback alternative of comparable induced deficit it had Y% power.

Esto permitiría distinguir:

- poca potencia frente a efectos de esa forma/magnitud;
- potencia razonable para excluir efectos mayores de esa clase.

No permitirá afirmar ausencia universal de memoria no lineal.

---

## 15. Claims prohibidos

Stage 4 no permitirá afirmar que:

- la potencia sea una función universal de \(\Delta\theta\);
- los dos mecanismos representen exhaustivamente la dinámica amazónica;
- un resultado negativo en Amazon pruebe ausencia de memoria física;
- las ocho unidades sintéticas sean equivalentes a ocho modelos CMIP6 independientes en sentido físico;
- todos los mecanismos con igual \(\Delta\theta\) tengan igual detectabilidad;
- el gate mensual sea óptimo.

---

## 16. Regla de detención

Stage 4 queda completo cuando existan:

1. 150 cohortes por los nueve puntos del diseño;
2. \(B=500\) por unidad;
3. gates técnicos completos;
4. curvas de potencia por mecanismo;
5. comparación descriptiva con \(\Delta\theta_{Amazon}=0.085009\);
6. umbrales de 50 % y 80 % si se alcanzan.

No se añadirán mecanismos ni niveles de \(\kappa\) después de observar resultados.

---

## 17. Próximo paso después de Stage 4

Si los gates técnicos pasan, no habrá Stage 5 exploratorio automático.

El siguiente paso será construir la **matriz final claim–evidence del Paper A**, integrando:

- Gate 0 bibliográfico;
- Stage 3 confirmatorio;
- Stage 4 de potencia;
- resultado Amazon congelado.

Solo entonces se decidirá si hace falta un segundo índice climático real.
