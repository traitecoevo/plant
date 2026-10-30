# July 2026 ATLS results (pre-port)

These files were computed on 2026-07-22 to 24 on `feature/tf24-atls-thermal-damage` (TF24t merged revision, before #591 moved the leaf into phylloptim). They were gitignored `.rds` caches until the port; they are kept as the record of what the July code gave, not as a baseline for the ported code.

They predate phylloptim's leaf-temperature-in-the-optimiser (#84), the Q10 respiration trait (#86), the leaf-to-air deficit (#93) and the birth-date coordinate (#590), so the ported code is expected to give different numbers.

| File | What it holds |
| --- | --- |
| `scm-viability-k0.00125.csv` | Low-density R0 (birth rate 1) across constant air temperatures 23 to 31 °C: TF24 with the energy balance vs default TF24t, at `k_mat = 0.00125 /day` |
| `scm-viability-strat-k0.00125.csv` | The same for four TF24t strategies (Generalist, Tolerant, Repairer, Acclimator), 23 to 33 °C |
| `running-mean.csv` | The running-mean leaf-temperature demo trajectory: air, operating and running-mean leaf temperature through time |
| `*-call.R` | The call that produced each file, as stored in its cache |
| `../gate_atls_vs_pm_A.csv`, `_B.csv` | The verification gate's two tables, computed on the superseded N-switch version |

The calls use helpers from the July branch's `overstorey_staging/atls_demo_helpers.R`.
