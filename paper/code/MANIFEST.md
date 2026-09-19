# Code manifest

Synchronised from the project repository at commit `ccc719c` on 2026-09-19.

Each entry is the SHA-256 prefix of the file as shipped. Regenerate and re-verify with:

```bash
python paper/sync_code.py --check
```

| File | SHA-256 (first 16) | Bytes |
|---|---|---:|
| `code/configs/default.yaml` | `5d03c8b204c34a69` | 1700 |
| `code/configs/experiments.yaml` | `9d61a25470c5bdc5` | 1205 |
| `code/data/external/aca.all.cfg` | `4c1d48465343cb0a` | 1211 |
| `code/data/external/alma.all.cfg` | `1396bc73e796621e` | 10059 |
| `code/data/external/SOURCES.md` | `5692346fba0036dd` | 4015 |
| `code/data/external/vla.a.cfg` | `e958cbdf312f056f` | 1545 |
| `code/experiments/06_layout_shootout.py` | `6df6e2a95035faa0` | 11421 |
| `code/experiments/07_formulation_comparison.py` | `ae1fd7d390211eac` | 8607 |
| `code/experiments/08_optimizer_benchmark.py` | `758a8c79dc2b2120` | 19450 |
| `code/experiments/09_science_cases.py` | `bdc94b559264c9d5` | 17146 |
| `code/experiments/10_multiepoch.py` | `2cdd755021f7cdcb` | 14170 |
| `code/experiments/11_real_alma.py` | `beeaff0623433a10` | 17288 |
| `code/experiments/12_qubo_solver.py` | `23bddcc64eb286ce` | 17844 |
| `code/experiments/13_science_case.py` | `5ed874608c79f313` | 14100 |
| `code/experiments/common.py` | `6074ef5db54ecc6e` | 6643 |
| `code/pytest.ini` | `3b3ffd97382ae1ec` | 61 |
| `code/requirements.txt` | `a36ae851cbdc6027` | 665 |
| `code/src/thz_opt/__init__.py` | `6af959cffcfd100f` | 987 |
| `code/src/thz_opt/arrays/__init__.py` | `e3b0c44298fc1c14` | 0 |
| `code/src/thz_opt/arrays/fibonacci.py` | `766b1655ec8eff4a` | 4988 |
| `code/src/thz_opt/arrays/geometries.py` | `f29d30261f1e6b2e` | 11929 |
| `code/src/thz_opt/arrays/golden_spiral.py` | `7bcd83a5a5a2b756` | 2406 |
| `code/src/thz_opt/arrays/random_array.py` | `7d0ce673dbcc957f` | 2027 |
| `code/src/thz_opt/arrays/real_arrays.py` | `11071d5bc91541ea` | 7329 |
| `code/src/thz_opt/arrays/validation.py` | `57236f8e8a4bd8df` | 1528 |
| `code/src/thz_opt/constraints/__init__.py` | `e3b0c44298fc1c14` | 0 |
| `code/src/thz_opt/constraints/atmosphere_data.py` | `b0997c86b211b4b0` | 6800 |
| `code/src/thz_opt/constraints/cable.py` | `506316e1fae8e5dd` | 3845 |
| `code/src/thz_opt/constraints/coherence.py` | `d55a2375957c31c2` | 9308 |
| `code/src/thz_opt/constraints/phase.py` | `5425e803066341d0` | 5293 |
| `code/src/thz_opt/constraints/pwv.py` | `cd9b4254f9c9dbee` | 2404 |
| `code/src/thz_opt/constraints/separation.py` | `992653d9252b5a0b` | 2397 |
| `code/src/thz_opt/interferometry/__init__.py` | `e3b0c44298fc1c14` | 0 |
| `code/src/thz_opt/interferometry/baselines.py` | `ab0062b5e24699be` | 1815 |
| `code/src/thz_opt/interferometry/earth_rotation.py` | `1d8914658e118c50` | 6002 |
| `code/src/thz_opt/interferometry/gridding.py` | `863330944eefc5a1` | 6444 |
| `code/src/thz_opt/interferometry/multifrequency.py` | `a44f762310953d83` | 5716 |
| `code/src/thz_opt/interferometry/psf.py` | `6da76f65d5aa87ae` | 3112 |
| `code/src/thz_opt/interferometry/uv.py` | `28f1b2d8ca3836e7` | 5577 |
| `code/src/thz_opt/metrics/__init__.py` | `e3b0c44298fc1c14` | 0 |
| `code/src/thz_opt/metrics/baseline_distribution.py` | `2b10f1909395bed6` | 3969 |
| `code/src/thz_opt/metrics/density_matching.py` | `802c447afb53600f` | 7347 |
| `code/src/thz_opt/metrics/psf_metrics.py` | `b71185bb84fe7d69` | 3524 |
| `code/src/thz_opt/metrics/redundancy.py` | `2427146ffed8afc8` | 1625 |
| `code/src/thz_opt/metrics/uv_coverage.py` | `2f537c0fba66cc72` | 5568 |
| `code/src/thz_opt/optimize/__init__.py` | `a0b1d6e2b7b63be9` | 2126 |
| `code/src/thz_opt/optimize/exact.py` | `eb9b69bb26d45471` | 6876 |
| `code/src/thz_opt/optimize/heuristics.py` | `3952be466cc5fdb9` | 13790 |
| `code/src/thz_opt/optimize/objectives.py` | `19fcea8a0d73828e` | 4121 |
| `code/src/thz_opt/optimize/pareto.py` | `9bec9752d1e5b3e0` | 3584 |
| `code/src/thz_opt/optimize/science_cases.py` | `adbdc88dd97221d4` | 7261 |
| `code/src/thz_opt/optimize/state.py` | `75cbf0bfe64b6337` | 7249 |
| `code/src/thz_opt/qubo/__init__.py` | `e3b0c44298fc1c14` | 0 |
| `code/src/thz_opt/qubo/baseline_qubo.py` | `4595b8cef7fba6ed` | 25300 |
| `code/src/thz_opt/qubo/coefficients.py` | `7fe72f9ba4d5b08c` | 8604 |
| `code/src/thz_opt/qubo/exhaustive.py` | `cbf3169072fff3e4` | 2496 |
| `code/src/thz_opt/qubo/multiepoch.py` | `fd854ce1dc88ee0c` | 10036 |
| `code/src/thz_opt/qubo/objective.py` | `32dc41314b196a0a` | 12472 |
| `code/src/thz_opt/qubo/robust.py` | `fd4e85ad24ad4919` | 3403 |
| `code/src/thz_opt/qubo/validation.py` | `0342605d7e1fc721` | 5235 |
| `code/src/thz_opt/science/__init__.py` | `1da8feb5cff81a92` | 58 |
| `code/src/thz_opt/science/disk_gap.py` | `f86b60b0c14b3937` | 15902 |
| `code/tests/test_arrays.py` | `4140cb46ac58e48b` | 3659 |
| `code/tests/test_baseline_qubo.py` | `991c8853c68316d0` | 9226 |
| `code/tests/test_baselines.py` | `ed6032641f46f51b` | 2116 |
| `code/tests/test_coherence.py` | `1ca3ca70be3bedf3` | 6361 |
| `code/tests/test_extensions.py` | `4b91a71bc2cbb436` | 9752 |
| `code/tests/test_geometries.py` | `7513a8ab9a7738b1` | 6444 |
| `code/tests/test_ladakh_design.py` | `7a987f3ef8731fd2` | 11144 |
| `code/tests/test_metrics.py` | `18035b79405c7a68` | 6679 |
| `code/tests/test_optimize.py` | `019c131fe3eb422a` | 2528 |
| `code/tests/test_qubo.py` | `a947af5ae1b34591` | 6969 |
| `code/tests/test_qubo_solver.py` | `ae07699c14045141` | 5291 |
| `code/tests/test_real_arrays.py` | `007aef074c486fb6` | 3845 |
| `code/tests/test_robust_cable_sparse.py` | `4a8378fd0f020f4d` | 6094 |
| `code/tests/test_science.py` | `66f96805fa227621` | 4980 |
| `code/tests/test_uv.py` | `ca502923bd599d7a` | 7825 |
| `data/exp06_shootout.csv` | `7b7a2b3cf1815451` | 3579 |
| `data/exp07_formulations.csv` | `139af6d1b1a656aa` | 1278 |
| `data/exp08_certified_instance.csv` | `1e3894a8b3e343e3` | 496 |
| `data/exp08_optimizer_benchmark.csv` | `7ea0f83867c442bf` | 2116 |
| `data/exp08_summary.json` | `bd89c89bd617ac47` | 11823 |
| `data/exp09_multifrequency.csv` | `540a5a351afd9ac5` | 353 |
| `data/exp09_science_cases.csv` | `3382b7eb044bc546` | 493 |
| `data/exp09_soft_vs_hard.csv` | `760966e4fd178d35` | 506 |
| `data/exp09_summary.json` | `3d063fb633f25d29` | 4634 |
| `data/exp10_multiepoch.csv` | `aca128e93630ac47` | 471 |
| `data/exp10_summary.json` | `e67293c4f4ac329e` | 2944 |
| `data/exp10_window_sweep.csv` | `2e160e493d8531fd` | 280 |
| `data/exp11_real_alma.csv` | `8bf092e2b7b07bcf` | 1243 |
| `data/exp11_summary.json` | `cb8e11ad44bf5edc` | 10207 |
| `data/exp12_gap_decomposition.csv` | `713c6a315e6abf1f` | 23258 |
| `data/exp12_solver_benchmark.csv` | `4c75c30e79357277` | 1678 |
| `data/exp12_summary.json` | `603602269fa95ca8` | 7934 |
| `data/exp13_science_case.csv` | `00bf6682fe4e578c` | 980 |
| `data/exp13_summary.json` | `f9b2f8560c24b785` | 4443 |
| `figures/fig16_formulation_comparison.png` | `a96072f98354724f` | 120249 |
| `figures/fig17_optimizer_comparison.png` | `43ae55addb5d1c09` | 102214 |
| `figures/fig18_optimized_layout.png` | `027c2465288440cc` | 94933 |
| `figures/fig19_science_case_layouts.png` | `dae506d90279c9a9` | 70033 |
| `figures/fig21_multifrequency_gain.png` | `b0da7cd4e8bfc186` | 82519 |
| `figures/fig22_multiepoch_tradeoff.png` | `7ae287e4888cdfb2` | 111684 |
| `figures/fig23_real_alma_pads.png` | `b01dbbb8211caa05` | 67423 |
| `figures/fig24_real_alma_selection.png` | `f3d5f80cffc00acc` | 65613 |
| `figures/fig25_qubo_gap_decomposition.png` | `f8bf750d74b7f13e` | 130000 |
| `figures/fig26_qubo_solver_convergence.png` | `9432aed491bbfe92` | 93752 |
| `figures/fig27_science_weight_derivation.png` | `7041e2d9ed199705` | 146103 |
| `figures/fig28_science_vs_cellcount.png` | `dc3683ea6c406de8` | 70796 |

108 files.
