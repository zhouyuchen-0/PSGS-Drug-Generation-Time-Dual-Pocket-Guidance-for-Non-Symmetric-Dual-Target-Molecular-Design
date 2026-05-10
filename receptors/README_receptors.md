# Receptor and docking-box files

This folder contains prepared receptor files and Vina configuration records used for the PSGS-Drug downstream docking protocol.

- `3fappro.pdbqt`: prepared 3FAP receptor PDBQT file.
- `7pqvpro_clean_fix2.pdbqt`: prepared 7PQV receptor PDBQT file.
- `3fap.pkl` and `7pqv.pkl`: pocket representation files used by the pocket-prior module.
- `3fapconfig.txt` and `7pqvconfig.txt`: Vina box configuration records.

Docking boxes:

| Target | Center (x, y, z) | Size (x, y, z) |
|---|---|---|
| 3FAP | -11.655, 24.756, 33.139 | 37.5, 37.5, 45.0 |
| 7PQV | -5.223, 67.587, 35.180 | 15.144, 16.921, 10.290 |

Docking scores are used as computational prioritization signals and not as experimental evidence of binding affinity.
