# Packet 10A7E public physical-input audit

Packet 10A7E is a source-free audit based on HFDMRG commit
`68364562d78c80da31ae351417bddc88138f155f`. It performed no solve and changes
no production source, test, dependency, control, or preserved Packet 10A6D
documentation. The frozen scientific source is PPP commit
`0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c`; it is an external oracle and is
not a runtime dependency of the proposed public reproduction.

## Exact physical producer chain

The accepted H-chain fixture is created by `sliced_hydrogen_chain` from the
GaussletBases revision pinned by frozen PPP,
`c680bbbcce5e7ff9ef6405c60fec0d3c2fabb2f0`. Its resolved inputs are:

- open H chains at separation `R=3.6` bohr;
- 18 longitudinal sites per atom (spacing `0.2` bohr), atom-centered phase
  `0.0`, two-bohr left and right padding, and `:atom_centered` alignment;
- H/STO-6G transverse data with `:periodic_template` reuse;
- one-body tolerance `1e-8` and interaction tolerance `1e-5`;
- basis fingerprint
  `73cc87798a9c78e7cc0f41ec56fc078bf8ccb8500fa6acb9655da8941bc788c7`.

GaussletBases constructs the open-boundary one-body bands by summing the
attraction from every nucleus. Therefore H1 is length- and boundary-dependent:
the H10 bands cannot be reused for H20. The interaction compressor is also
trained on the target chain in the accepted fixture. Its retained ranks differ:
H10 has 183 sites, tail/residual ranks 62/11 and maximum scalar rank 82; H20
has 363 sites, ranks 64/12 and maximum scalar rank 85. A fixture must carry two
separately certified operators rather than changing only `sites`.

The compact solve needs, per length, the exact `BandedOneBody.bands`, the
`UnitCellInteraction` values `sites`, `phase_at_first`, `tail_transfer`,
`residual_transfer`, `phase_source`, `phase_sink`, and `local_triangle`, and
the explicit producer-derived atom intervals. Those intervals must be stored,
not regenerated from a presumed integral-cell layout: the accepted open chain
has padded terminal fragments and Packet 10A7B correctly rejects an omitted
partition for this ambiguous geometry. The dimer RHF and atomic-Neel UHF starts
then require no coefficient checkpoint.

The physical producer audit additionally needs the uncompressed sliced H1 and
Vee values and the nuclear-repulsion constant. The accepted constants are
`5.358245149911816 Ha` for H10 and `14.431886984131566 Ha` for H20. These are
not inputs to `solve_hfdmrg`; they enter only a `ProducerHamiltonian` audit (or
are added when comparing an electronic represented energy with a total-energy
reference).

## Existing provenance and reference boundary

Frozen PPP's `validation/packet6a4/run_anchor.jl` owns the physical constructor
and `src/SlicedHydrogenModels.jl` owns the GaussletBases-to-banded/MPO adapter.
The later block-native runner imports that closure. Current HFDMRG large-chain
runners still import it by path, so they are audit records, not clean-checkout
public reproductions.

The Packet 6A4 H10/H20 records retained at frozen PPP contain total-energy
anchors (electronic plus nuclear repulsion): RHF/UHF are respectively
`-4.2015409365128`/`-4.751294669000901 Ha` for H10 and
`-8.458420930773865`/`-9.559260520716624 Ha` for H20. Their embedded generation
commit is `507cd8493d87e8bd4f1b1f02c618e321860aae2f`; the frozen `0bcd4f1` tree is
the immutable source of the retained records. They are provenance anchors, not
new recipient measurements. The accepted recipient H20 UHF electronic
represented energy is independently recorded as approximately
`-23.991147504847 Ha`; adding `14.431886984131566 Ha` explains its agreement
with the PPP total-energy basin.

`represented_energy` is the electronic energy of the compressed
`UnitCellInteraction` used by the compact solver. `producer_energy` is a
separate streamed observation of a caller-supplied uncompressed Hamiltonian;
its `constant` can include nuclear repulsion. It neither changes the compact
state nor proves producer-Hamiltonian stationarity, a residual bound, or a
global HF minimum.

## Distribution alternatives

### Pinned serialized input (recommended bounded gate)

The smallest self-contained route is one versioned H10/H20 data bundle plus a
short public runner. Serialize a NamedTuple of primitive arrays and scalars,
not private PPP/GaussletBases structs. Pin Julia 1.12.6 and little-endian
Float64, record every shape and SHA-256 in a text manifest, and include:

1. the two length-specific compact H1/MPO payloads and exact atom intervals;
2. chain parameters, units, tolerances, basis fingerprint, ranks, source
   revisions, construction-report facts, and nuclear constant;
3. the full ordered producer Vee (and physical H1 bands) if public
   represented-versus-producer comparison is claimed; Packet 10A7F found that
   independently evaluated opposite triangles are not bitwise identical, so
   triangle packing would not preserve the pinned producer input; and
4. reference results as metadata, never as construction inputs.

Julia `Serialization` keeps the loader dependency-free but is intentionally
version-pinned and opaque. A portable explicit binary-array schema is the
better later archival form; it costs a small loader but avoids binding the
fixture to Julia's serialization format. Packet 10A7F selected that portable
form and retained the full ordered producer Vee to preserve exact source
parity; the bounded H10/H20 sizes make materialization by the runner
reasonable.

### Public construction

GaussletBases `c680bbbc...` has an HTTPS GitHub source URL and an MIT license
with Steven R. White's notice. It can construct the physical H1/Vee from the
listed parameters. It cannot by itself construct HFDMRG's accepted compact
MPO: that adapter and compression closure lives in PPP, whose frozen commit has
no `LICENSE`, `COPYING`, or `NOTICE`. Copying that broad closure would also
violate the bounded-fixture objective. Consequently a public-construction
route is not currently a clean HFDMRG workflow, even if GaussletBases is
publicly fetchable. It would require a separately authorized and packaged
compact-operator constructor, versioned as scientific production code.

HFDMRG and GaussletBases use matching MIT terms and their notices can accompany
a generated fixture. Before publication, the owner should explicitly approve
distribution of the generated H10/H20 compact coefficients, producer arrays,
and retained PPP reference records under the HFDMRG distribution terms. That
small data authorization is the only rights decision identified here; no
external repository needs relicensing.

## Public runner contract

The clean-checkout runner should load a selected `h10` or `h20` record, verify
revision/schema/shape/hash metadata before constructors, then call:

```julia
solve_hfdmrg(one_body, interaction;
    spin = :rhf, initialization = :dimer,
    state_maxdim = 16, state_cutoff = 1e-12,
    energy_tolerance = 1e-11, maximum_half_sweeps = 20,
    starting_orientation = :physical, exact = false,
    atom_intervals = fixture.atom_intervals)
```

For UHF, use `spin=:uhf, initialization=:atomic_neel` with the same controls;
`spin_swapped=false` is the physical default. The runner must report
`converged`, `reason`, `represented_energy`, `half_sweeps`, accepted/rejected
update counts, ranks, discarded weights, replay energy/projector changes,
`requested_energy_tolerance`, `effective_energy_tolerance`, and
`stopping_tolerance_reason`. UHF reports separate alpha/beta ranks, pair ranks,
discarded weights, and projector changes. A producer audit should print the
decomposed `ProducerEnergy` fields and label its `total` separately.

## Pending Packet 10A6D integration

Do not rewrite the six preserved 10A6D files during fixture extraction. After
the fixture passes review, update its existing reproducibility guide in one
docs-only reconciliation: replace the current external-fixture limitation with
the public fixture path and hashes; add exact H10/H20 commands; retain the PPP
oracle-versus-recipient distinction and all H100/H1000 timing attribution; and
link the new command from the already-pending README/manual indices and release
checklist. The guide must continue to label clean GitHub checkout/CI execution
as future evidence until it actually runs.

## Recommendation and next gate

Proceed with one bounded fixture packet, not a producer-framework port:

1. obtain explicit owner authorization for the generated data and retained
   reference metadata;
2. generate H10 and H20 independently from the exact pinned constructor;
3. store only primitive payloads and manifest provenance, and prove dense
   H1/Vee and compact interaction entry parity against the producer before any
   solve;
4. run one clean temporary-checkout RHF and UHF reproduction per length under
   Julia 1.12.6, checking requested/effective tolerance reporting and the
   represented/total-energy offset; and
5. only after that gate, reconcile and commit the preserved 10A6D guide.

No current distributable runner/input satisfies this gate. The missing bounded
fixture is the remaining public paper-reproduction item; it does not invalidate
the package examples, accepted compact solver evidence, or historical
backends.
