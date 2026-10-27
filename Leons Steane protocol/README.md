# Phase 5: Vollständige iQPE auf zwei Steane-Logikqubits

Lokal erstellt, nicht auf GitHub.

## Protokoll (`iqpe_full.py`)

- **C (iQPE-Ancilla):** in jeder Iteration frisch |+_L⟩, am Ende X-Auslese.
- **S (System):** |1_L⟩, Eigenzustand von U = P(2πφ). Wird einmal präpariert und bleibt über alle 3 Iterationen erhalten. Fehler auf S wirken also in allen folgenden Iterationen nach.
- **Kontrollierte Potenzen** wie in Dustins `cu`: CP(2a) = P(a)_C · P(a)_S · CX · P(−a)_S · CX, mit a = j·2^m·π/8.
  - Für j = 2, 4, 6 ist das in Clifford+T exakt: CZ und CS. CS enthält echte T-Gadgets, darunter eines auf dem verschränkten Paar.
  - Für ungerade j bräuchte die letzte Iteration √T. Das ist hier nicht umgesetzt.
- **Bausteine und Rauschen** wie in `phase4_qiskit_hpc`:
  - SD6, alle Gates verrauscht, kein Pauli-Frame;
  - Flag-QEC nach jedem logischen Gate, FT-Magic-State;
  - transversales CX zwischen den Blöcken, danach QEC auf beiden.
- **Simulation:** gebündelte Zustandsvektor-Trajektorien (numpy).
  - C und S sind nur zwischen den beiden CX eines kontrollierten Gates verschränkt. Dafür werden sie zu 16 Qubits zusammengelegt, sonst läuft jeder Block mit 9 Qubits.
  - Die Faktorisierung wird pro Shot per SVD geprüft. Bislang war sie immer exakt, mit Rest ≈ 2·10⁻⁷, und sie ist bitgenau identisch zur vollen 16-Qubit-Simulation (`run_full`).
  - Die Teleportation des T-Gadgets wird exakt ohne 21-Qubit-Zustand gerechnet.

## Tests

1. **`--selftest`** (rauschfrei): kontrollierte Gates, vollständige iQPE für j = 0, 2, 4, 6 exakt, S bleibt in |1_L⟩.
2. **`--ftcheck`:** exakter Einzelfehler-Test für eine Iteration mit CZ, mit CS + S† und mit CS³ + S.
   - 6970, 10966 bzw. 12829 Einzelfehler.
   - Kein Fehler ändert das Ergebnis von C oder kippt S. **Bestanden.**
3. **Gegenprobe Ramsey (ein Block)** gegen `phase4_qiskit_hpc`: Alle Bausteine und die bedingte Rückkopplung stimmen innerhalb von etwa 1σ überein.
   - Offen: Ramsey j = 3 weicht beim Qiskit-Lauf um 3.6σ ab (Qiskit hatte dort nur 8000 Shots). Die übrigen Winkel passen.

## Sweep (`sweep.py`)

- Winkel j = 2, 4, 6.
- p ∈ {2e-3, 1e-3, 5e-3, 5e-4, 2e-4, 1e-4}, in Chunks von 512 Shots.
- Dazu das nackte 2-Qubit-iQPE (gleiches Rauschen, P(a) als ein verrauschtes Gate) mit 2·10⁶ Shots pro Punkt.
- Ergebnisdateien: eine pro (Lauf, j, p, Seed) in `results/`. Bestehende Dateien werden übersprungen.

```bash
python sweep.py --procs 7                        # lokal
sbatch run_hpc.sh                                # HPC: 32 Tasks x 16 Kerne, jeder Task jeder 32. Job
python sweep.py --summary                        # Tabelle + results/sweep_full.png
```

**Aufwand:** etwa 4800 Steane-Jobs × 512 Shots × ~0.15 s ≈ 100 Kernstunden.

## Ergebnisse (HPC + lokal, 26.10.2026)

Datensatz: `results/sweep_full.json`, Abbildung `results/compare.png` (zum Vergleich mit Ramsey aus `phase4_qiskit_hpc`).

- **Umfang:** 4824 von 4962 geplanten Jobs. Es fehlen 138 Steane-Chunks bei kleinem p, vermutlich Tasks am Zeitlimit. Deshalb sind es bei 10⁻⁴ und 2·10⁻⁴ etwas weniger Shots als geplant.
- **Darstellung:** mittleres Δφ über j = 2, 4, 6.

| p | Steane, vollständiges iQPE | Shots | nacktes 2-Qubit-iQPE | Steane / nackt | S am Ende gekippt |
|---|---|---|---|---|---|
| 1e-4 | 9.13e-5 ± 0.55e-5 | 1.14 M | 5.88e-4 | **0.16** | 1.1e-4 |
| 2e-4 | 3.28e-4 ± 0.13e-4 | 762 k | 1.18e-3 | **0.28** | 3.8e-4 |
| 5e-4 | 2.01e-3 ± 0.05e-3 | 385 k | 2.93e-3 | **0.69** | 2.3e-3 |
| 1e-3 | 7.82e-3 ± 0.17e-3 | 98 k | 5.82e-3 | 1.34 | 8.8e-3 |
| 2e-3 | 2.91e-2 ± 0.05e-2 | 49 k | 1.15e-2 | 2.52 | 3.5e-2 |
| 5e-3 | 1.26e-1 ± 0.01e-1 | 25 k | 2.77e-2 | 4.55 | 1.7e-1 |

- **Break-even des vollständigen iQPE bei p ≈ 7.4·10⁻⁴** (log-log interpoliert). Ramsey liegt bei ≈ 1.0·10⁻³.
- **Darunter gewinnt die Codierung deutlich:** Faktor 3.6 bei 2·10⁻⁴, Faktor 6.4 bei 10⁻⁴.
- **Skalierung:** Die Steigung des logischen Fehlers ist 1.84–1.98, also quadratisch, wie es bei Fehlertoleranz zu erwarten ist. Das nackte iQPE wächst linear.
- **System-Block:** Die logischen Kipp-Ereignisse von S skalieren ebenfalls etwa quadratisch, sie sind ein wesentlicher Teil des Fehlerbudgets.
- **j = 4** (nur ein CZ, keine T-Gates) ist am besten, **j = 6** (CS³ mit 4 T-Gadgets) am schlechtesten.
- **Faktorisierungs-Kontrolle:** Der zweite Singulärwert liegt relativ zum ersten bei ~10⁻¹⁵, die Faktorisierung ist also exakt. Die in den JSON-Dateien gespeicherte Größe `non_product`/`max_split_residual` hat bis 26.10. versehentlich die float32-Normdrift gemessen (√(1−σ₁²) ≈ 10⁻³ bei einer Normabweichung von 10⁻⁶). Sie ist inzwischen korrigiert und hat die Ergebnisse nicht beeinflusst.
