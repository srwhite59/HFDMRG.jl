# Frozen-Be Stationarity Diagnostic

Date: 2026-07-20

Owner: `hfdmrg-manager`

Branch/base: `experiment/be-stationarity-20260720` / `06023eace6b08c4087f7b760c8206c4cadb6f881`

Status: **calculation complete; stop for paper-manager review**

## Outcome

Tightening or exhausting the local HFDMRG microiteration budget does not relieve the
sweep-140 residual floor. Ten fresh restarted sweeps at `scf_cutoff=1e-12`, `1e-14`, and `0`
all remain near
`6.33e-7`, even though the routes use 2,940, 3,155, and 4,900 local Fock builds. Their
occupied spaces move by at most `1.03e-9` from the common input at the requested sample
sweeps.

One unconditional full-space Fock update reduces the residual by `1011.17x`, from
`6.3275583e-7` to `6.2576393e-10`, while changing the physical energy by only
`+3.55e-14 Ha`. The retained monotone full-space UHF procedure accepts no update: its first
damped trial changes energy by `+1.78e-14 Ha`, so the historical `1e-10 Ha` equality test
stops immediately and returns the input. Thus energy flatness hides a removable
nonstationarity from that full-space polisher. This does not establish the cause of the
HFDMRG sweep floor, which persists even when every local window exhausts its iteration
budget.

This is diagnostic evidence, not a production algorithm decision. No HFDMRG interface,
source, test, accepted cache account, March/history code, or stopping policy changed.
It does not show that repeated raw Fock updates are stable or energy-monotone and does not
validate a stationarity stopping rule.

## Convention And Frozen Input

For each spin,

```text
D_sigma = C_sigma C_sigma'
R_sigma = F_sigma C_sigma - C_sigma(C_sigma' F_sigma C_sigma)
R_OV = hypot(norm(R_alpha), norm(R_beta))
E = 1/2 Tr[D_alpha(F_alpha+H)] + 1/2 Tr[D_beta(F_beta+H)].
```

The Fock matrices and energy use the accepted independent six-index physical contraction.
The residual uses the raw Fock. For the occupied Fock eigenspace, only the diagonalization
copy is averaged as `(F+F')/2`; if `U` is its lowest occupied eigenspace,

```text
theta_sigma = atan(opnorm(U-C(C'U)), sigma_min(C'U)).
```

Occupied-space movement is the symmetric two-sided projector distance from the common input.
`S2 = Sz(Sz+1)+Nbeta-norm(Calpha' Cbeta)^2`, and Z-spin is the alpha expectation
minus the beta expectation of the saved `Zop`.

The 468-by-2 spin orbitals came from the forced sweep-140 checkpoint. Every route received a
fresh copy with identical raw-byte hashes:

```text
packet SHA-256:     fba29651c02cde3cdde90850d90ba30a1abd120ebfa8ac8b64cf2e5a92fea629
checkpoint SHA-256: 4b17397a83e97d90667cac1210e6b1bd95d7e263118764e4e1e073484e4f3112
alpha orbital SHA:  75a2a5571d9390edb46dcf3e607d2407b59f7e61a0b153f4ebc4b8d816e50e41
beta orbital SHA:   32f762d51c6df9b15296c842308c8574efd0a793556aee17eb7d037391a03a63
```

The independently recomputed starting values are:

| Quantity | Sweep-140 value |
|---|---:|
| Physical energy | `-14.573351154440026 Ha` |
| Checkpoint reported minus physical energy | `+9.41e-14 Ha` |
| Raw occupied-virtual residual | `6.327558312251281e-7` |
| Fock angle, alpha / beta | `1.257440084e-8 / 1.257444645e-8 rad` |
| Occupied Fock gap, alpha / beta | `0.3199864739 / 0.3199864737 Ha` |
| `S2` | `0.13165093167387432` |
| Z-spin | `-1.1194310437547075` |

## Fresh-Restart Comparison

Each row is independently diagnosed after the observer return. `dP` is the larger spin
movement. Fock calls and wall times are `this sweep / cumulative`; diagnostics are excluded
from solver timing. `cutoff=0` disables the outer energy stop. These runs are not literal
continuations: each rebuilds cached environments and damping history from the same orbitals.

| Local cutoff | Sweep | Physical energy (Ha) | `R_OV` | `theta_max` (rad) | `dP` | `S2` | Z-spin | Fock calls | Wall (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1e-12` | 1 | `-14.573351154439997` | `6.3274494e-7` | `1.2573557e-8` | `1.6920e-10` | `0.131650931552007` | `-1.119431043192306` | `294 / 294` | `0.1949 / 0.1949` |
| `1e-12` | 2 | `-14.573351154439916` | `6.3290631e-7` | `1.2572917e-8` | `3.8793e-10` | `0.131650931330194` | `-1.119431042209692` | `294 / 588` | `0.1210 / 0.3159` |
| `1e-12` | 5 | `-14.573351154440008` | `6.3283997e-7` | `1.2571514e-8` | `6.1829e-10` | `0.131650931182394` | `-1.119431041490022` | `294 / 1470` | `0.1212 / 0.7586` |
| `1e-12` | 10 | `-14.573351154439989` | `6.3274302e-7` | `1.2570646e-8` | `9.2710e-10` | `0.131650930926687` | `-1.119431040306306` | `294 / 2940` | `0.1215 / 1.4798` |
| `1e-14` | 1 | `-14.573351154439859` | `6.3294821e-7` | `1.2573617e-8` | `2.2796e-10` | `0.131650931519991` | `-1.119431043055687` | `330 / 330` | `0.1269 / 0.1269` |
| `1e-14` | 2 | `-14.573351154439866` | `6.3283323e-7` | `1.2572886e-8` | `3.3192e-10` | `0.131650931401092` | `-1.119431042510459` | `336 / 666` | `0.1632 / 0.2900` |
| `1e-14` | 5 | `-14.573351154440035` | `6.3276027e-7` | `1.2571721e-8` | `7.4700e-10` | `0.131650931077871` | `-1.119431041042491` | `311 / 1628` | `0.1234 / 0.6995` |
| `1e-14` | 10 | `-14.573351154440005` | `6.3280476e-7` | `1.2570628e-8` | `9.9714e-10` | `0.131650930855698` | `-1.119431040004272` | `299 / 3155` | `0.1554 / 1.4096` |
| `0` | 1 | `-14.573351154440061` | `6.3275427e-7` | `1.2573638e-8` | `1.1575e-10` | `0.131650931615407` | `-1.119431043461416` | `490 / 490` | `0.1270 / 0.1270` |
| `0` | 2 | `-14.573351154439926` | `6.3284345e-7` | `1.2572933e-8` | `3.1444e-10` | `0.131650931513735` | `-1.119431042990684` | `490 / 980` | `0.1367 / 0.2638` |
| `0` | 5 | `-14.573351154439937` | `6.3285442e-7` | `1.2571578e-8` | `5.5694e-10` | `0.131650931229397` | `-1.119431041692109` | `490 / 2450` | `0.1730 / 0.7468` |
| `0` | 10 | `-14.573351154440008` | `6.3288830e-7` | `1.2570840e-8` | `1.0323e-9` | `0.131650930805450` | `-1.119431039792378` | `490 / 4900` | `0.1288 / 1.4993` |

Across all twelve samples, the residual range is `6.3274302e-7` to
`6.3294821e-7`, the Fock-angle range is `1.2570628e-8` to `1.2573638e-8 rad`,
and the largest reported-versus-physical energy mismatch is `2.49e-13 Ha`. Every sweep uses
98 incremental cached absorptions and 98 local windows, with zero fallback or projection
windows.

| Local cutoff | Initialization (s) | Ten sweeps (s) | Solve call (s) | Total Fock calls | Cumulative sweep allocation |
|---:|---:|---:|---:|---:|---:|
| `1e-12` | `0.0435` | `1.4798` | `1.5234` | `2940` | `1.788 GiB` |
| `1e-14` | `0.0445` | `1.4096` | `1.4541` | `3155` | `1.796 GiB` |
| `0` | `0.0405` | `1.4993` | `1.5398` | `4900` | `1.886 GiB` |

The fixed-budget route proves that all four local updates were taken: each of its 98 windows
uses five Fock builds. The process maximum resident size was `1.237 GiB`; cumulative allocation
is not resident storage.

## Full-Space Comparisons From The Same Input

| Route | Physical energy (Ha) | `Delta E` (Ha) | `R_OV` | `theta_max` | Movement alpha / beta | Iterations and result | Algorithm time |
|---|---:|---:|---:|---:|---:|---|---:|
| One full Fock update | `-14.573351154439990` | `+3.55e-14` | `6.2576393e-10` | `1.0375199e-9` | `1.7783167e-8 / 1.7783232e-8` | one applied update; one physical Fock contraction | `0.1338 s` |
| Monotone UHF polish | `-14.573351154440026` | `0` | `6.3275583e-7` | `1.2574446e-8` | `2.43e-15 / 2.54e-15` | one outer iteration, one trial, zero accepted; equality stop | `0.2506 s` |

The full update gives `S2=0.1316509316575183` and
`Zspin=-1.119431043423488`. The polish returns the input values. Post-route physical
diagnostics required a separate contraction and are excluded from the algorithm times.

The polish reproduces the retained `hf_damped!` control procedure—initial alpha/beta damping
`0.9/0.89`, SVD subspace mixing, strict energy-decrease acceptance, ten backtracks, persistent
damping, and the `1e-10 Ha` equality stop—using the accepted physical oracle and averaged
symmetric Fock. It is not byte-for-byte execution of the historical direct builder. The frozen
reference file SHA-256 is
`937c45c91fbe4c949a96aee673df14bcf4c9884efc30d374a4580fd51bc65147`.

## Validation And Evidence

The run used Julia 1.12.6, one Julia thread, and eight-thread ILP64 OpenBLAS on
`rh310l.ps.uci.edu`. A discarded warmup preceded the matched measurements. The driver enforced
all per-sweep route/count gates. A separate audit recomputed physical energies and residuals
from every one of the 15 serialized determinant states and rechecked
orthogonality/idempotency below `1e-11`, table schemas, representative route counts, totals,
and input hashes.

Machine-local evidence is under
`/Users/srw/dmrgtmp/hfdmrg_be_stationarity_20260720/`:

- driver: `run_be_stationarity.jl`, SHA-256
  `3a6b3ddc57d64c61075421487113f2ae7bc1ad1dcb497086ba5068aa17fe4439`;
- audit: `audit_be_stationarity.jl` and `audit.log` (`serialized_state_audit_ok=true`);
- reduced tables: `restart_samples.tsv`, `restart_summary.tsv`, `fullspace.tsv`, and
  `polish_history.tsv`;
- raw states: `restart_states.jls` and `fullspace_states.jls`;
- run log/PID: `run.log` and `run.pid` (completed PID 63735).

The log was monitored directly with
`sed -n '1,200p' /Users/srw/dmrgtmp/hfdmrg_be_stationarity_20260720/run.log`.
The next action is paper-manager review of this evidence before any convergence design or March
port.

-- hfdmrg-manager@macmini
