# maf2con demo (coverage-aware, chrY example)

A tiny, self-contained chrY-style MAF that exercises `realignpro maf2con`'s
coverage-aware behavior — no external binaries and no large data required.

## Files
- `chrY_demo.maf` — 2-block MAF; reference `hg38.chrY`, targets `hs1..hs4`.
  - Block 1 (`chrY:0-8`): all assemblies present, fully conserved.
  - Block 2 (`chrY:20-28`): `hs4` **missing** (lower depth) and one **polymorphic**
    column at `chrY:24` (where `hs3` carries `T` instead of `A`).
- `run_maf2con_demo.sh` — runs `maf2con` in several modes and checks the output.

## Run
```bash
cd examples/maf2con
./run_maf2con_demo.sh
```

## What it shows
| Command | Result | Why |
|---|---|---|
| *default* | `chrY 0 8` | `N_exp` is auto-derived per chromosome (= 5 here); with `--min-call-rate 0.99` the floor is `ceil(0.99·5) = 5`, so block 2 (depth 4) is below the floor and is dropped. |
| `--min-call-rate 0` | `chrY 0 8`, `chrY 20 24`, `chrY 25 28` | No relative floor and **no whole-block gate**: block 2 is assessed on its 4 present assemblies; the polymorphic column at `chrY:24` splits the constrained interval. |
| `--fixed-only --emit-depth` | same intervals + `N_cov`, major-allele frequency | Only 100 %-fixed (monomorphic) columns are kept; the extra columns report coverage and frequency. |

Genome IDs are the text before the first `.` in each `s` line's source
(`hg38.chrY` → `hg38`), matching the rest of ReAlignPro.

Outputs are written to `examples/maf2con/work/` (git-ignored).
