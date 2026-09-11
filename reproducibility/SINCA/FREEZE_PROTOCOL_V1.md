# SINCA PM2.5 AOAS — Frozen primary cohort V1

**Status: FROZEN**

Freeze ID: `SINCA_PM25_AOAS_FROZEN_COHORT_V1`  
Primary cohort: **31 stations**  
Analysis window: **2018-01-01 to 2025-12-31**  
Station-list SHA-256: `59b019cd344b1cdf9c0f6fb90b381ad12d416d72461d94d31d6ce0039ff010e6`

## Principle

The primary station cohort is fixed before any SINCA station-level native/index surrogate or null-transport inference. No subsequent station may be added or removed on the basis of surrogate p-values, null discrepancies, effect directions, or visual inspection of inferential results. Any later change must be versioned as a new cohort freeze and justified independently of the inferential outcomes.

## Frozen eligibility rule

A station entered V1 only if all conditions were met:

1. unique short-code to station-catalog mapping;
2. >=95% PM2.5 daily coverage over 2018–2025;
3. >=95% of the 417 complete Monday–Sunday weeks have at least 6 observed days;
4. >=100 observed days with PM2.5 > 50;
5. <=0.1% observations outside the diagnostic range [0,500];
6. no gross recent scale-shift flag.

The range rule is a QC diagnostic, not a claim about physical impossibility. This freeze does not certify instrument homogeneity or source-unit provenance.

## Frozen summary

- Regions represented: **8**
- Total observed exceedance days >50: **10,738**
- N_eff,p range: **3.398–5.673**
- Median N_eff,p: **4.297**
- Minimum fixed-window daily coverage: **97.365%**
- Minimum >=6/7-week fraction: **96.403%**

## Frozen station roster

- `RIX/901` — TE — Las Encinas Temuco — IX
- `RIX/902` — PLCII — Padre Las Casas II — IX
- `RIX/905` — TEII — Nielol Temuco — IX
- `RM/D12` — LF — La Florida — RM
- `RM/D14` — SA — Parque O'Higgins — RM
- `RM/D15` — PU — Pudahuel — RM
- `RM/D18` — CN — Cerro Navia — RM
- `RM/D30` — QU — Quilicura — RM
- `RVI/609` — RGAI — Rancagua I — VI
- `RVI/611` — RNG — Rengo — VI
- `RVI/615` — RGAII — Rancagua II — VI
- `RVII/709` — CRCI — Curicó — VII
- `RVII/711` — TALIII — Universidad de Talca — VII
- `RVIII/802` — THNI — Consultorio - San Vicente — VIII
- `RVIII/804` — HPNIII — JUNJI — VIII
- `RVIII/806` — THNV — Inpesca — VIII
- `RVIII/810` — CHLI — INIA. Chillán — VIII
- `RVIII/832` — CNH — Balneario Curanilahue — VIII
- `RVIII/837` — THNVI — Nueva Libertad — VIII
- `RVIII/841` — HQI — Hualqui — VIII
- `RVIII/854` — CHYI — Punteras — VIII
- `RVIII/873` — CHLII — Puren — VIII
- `RVIII/874` — LAII — Los Ángeles Oriente — VIII
- `RVIII/875` — LAI — 21 de mayo — VIII
- `RX/A01` — OS — Osorno — X
- `RX/A07` — PMII — Mirasol — X
- `RX/A08` — PMI — Alerce — X
- `RXI/B03` — COI — Coyhaique I — XI
- `RXI/B04` — COII — Coyhaique II — XI
- `RXI/B05` — VI — Vialidad — XI
- `RXIV/E04` — LU — La Union — XIV

## Deliberately not frozen at this stage

- handling of weeks with one missing day;
- any imputation;
- exact native/index surrogate contracts for the applied analysis;
- primary station-level statistic and across-station aggregation/multiplicity rule.

Those items should be prespecified in the next protocol before running the substantive SINCA null-transport analysis.
