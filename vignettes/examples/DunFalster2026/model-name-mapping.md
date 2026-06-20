# Canopy-model naming: paper <-> plant code

| Paper name | Meaning | `plant` `shading_model` code |
|---|---|---|
| Deep Crown | smooth (Yokozawa) shade; assimilation integrated over crown depth | deep-crown |
| PPA (smoothed) | layered, stepped light profile, boundaries smoothed (runnable; all PPA results) | ppa, ppa_layer_smoothing > 0 (default 0.3) |
| PPA (hard step) | literal field discretisation; discontinuous, does not run | ppa, ppa_layer_smoothing = 0 |
| Flat Top | box / thin-layer shade (field's usual sense; = Stage-1 MAESPA "Flat Top") | flat-top-box (hard step), flat-top-soft-box (smoothed) |
| Mean-Light | smooth shade; light averaged over crown depth, one assimilation evaluation | mean-light  (plant-only, new) |
| Crown-Centre | smooth shade; one assimilation evaluation at the crown centre | crown-centre  (plant-only, new) |

- The code label `crown-centre` matches the paper name **Crown-Centre** (smooth shade, crown-centre assimilation) — it is **not** the field's "Flat Top", which is the box / thin-layer shade (plant's `flat-top-box` / `flat-top-soft-box`).
- **Flat Top** = the box family: "Flat Top (hard step)" (does not run) vs "Flat Top (smoothed)" (runs but biased).
- **PPA** is one shading model (`ppa`) at two settings of `ppa_layer_smoothing`, not two `shading_model` strings: "PPA (smoothed)" (`> 0`, default 0.3 — runnable; all reported PPA results) vs "PPA (hard step)" (`= 0` — the literal field discretisation; discontinuous, does not run in the adaptive solver).
- **Mean-Light** and **Crown-Centre** appear only in the dynamic (Stage 2) results.
- Crown-Centre and Flat Top use the *same* crown-centre assimilation and differ only in the shade cast — so Crown-Centre giving usable fitness while Flat Top fails isolates the shade cast, not the assimilation rule.
