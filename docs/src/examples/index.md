# Examples

The checked-in examples are small deterministic workflows that validate their
own results. Run them from the repository root in this order:

1. **Basic density-density RHF** -- `examples/quickstart.jl`

   Establishes the input shapes, return tuple, and orbital orthogonality.

   ```sh
   julia --project=. examples/quickstart.jl
   ```

2. **Basic density-density UHF** -- `examples/basic_uhf.jl`

   Uses independent alpha/beta occupations and a spin-asymmetric seed, then
   reconstructs the full physical densities, Focks, electronic energy,
   occupied--virtual residuals, gaps, and Aufbau projectors.

   ```sh
   julia --project=. examples/basic_uhf.jl
   ```

3. **Sweep observation and restart** -- `examples/observer_checkpoint.jl`

   Records complete-sweep history, writes an atomic local checkpoint, requests
   a clean stop, verifies readback, and performs a one-sweep restart. Julia
   `Serialization` is used only for the self-contained example; applications
   should choose an archival format appropriate to their project.

   ```sh
   julia --project=. examples/observer_checkpoint.jl
   ```

4. **Target-space residual UHF** -- `examples/target_residual.jl`

   Builds a signed residual in a two-orbital target space and exercises the
   target-residual backend without a global four-index tensor.

   ```sh
   julia --project=. examples/target_residual.jl
   ```

5. **Cached sliced RHF** -- `examples/cached_sliced.jl`

   Constructs a symmetric fixed-slice tensor, aligns the physical partition,
   and explicitly selects the cached production backend.

   ```sh
   julia --project=. examples/cached_sliced.jl
   ```

## Scope

These examples demonstrate package mechanics rather than benchmark claims.
They use small self-contained fixtures and omit private research controls,
heavyweight consumer artifacts, and production datasets.
