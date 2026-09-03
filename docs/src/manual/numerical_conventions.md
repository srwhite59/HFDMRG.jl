# Numerical conventions

HFDMRG uses real occupied-orbital matrices with orthonormal columns and real
symmetric one-body operators. It performs no unit conversion. The returned
energy has the units of the supplied operators and excludes nuclear repulsion
or any consumer-defined scalar offset.

For RHF, `psiup0` contains one occupied spatial orbital per column and denotes
one spin density. Thus ``n`` columns represent ``2n`` electrons. UHF supplies
``C_\alpha`` and ``C_\beta`` separately, and their column counts remain fixed.

## Density and trace convention

Let

```math
D_\sigma=C_\sigma C_\sigma^T.
```

During density damping, the matrix used inside a local SCF step may instead be
a non-idempotent mixture of successive determinant densities. Every backend
must support this case. With Fock matrices that already include the one-body
terms, HFDMRG uses

```math
E_{\mathrm{UHF}}=
\tfrac12\operatorname{Tr}[D_\alpha(F_\alpha+H_\alpha)]+
\tfrac12\operatorname{Tr}[D_\beta(F_\beta+H_\beta)],
```

and

```math
E_{\mathrm{RHF}}=\operatorname{Tr}[D(F+H)].
```

## Density-density interaction

The matrix backend means

```math
(ij\mid kl)=\delta_{ij}\delta_{kl}V_{ik}.
```

For UHF, define ``q=\operatorname{diag}(D_\alpha+D_\beta)``. Then

```math
F_\sigma=H_\sigma+\operatorname{Diag}(Vq)-V\mathbin{\odot}D_\sigma,
```

```math
E_{\mathrm{UHF}}=
\operatorname{Tr}(D_\alpha H_\alpha)+\operatorname{Tr}(D_\beta H_\beta)
+\tfrac12q^TVq
-\tfrac12\sum_{ij}V_{ij}\left(D_{\alpha,ij}^2+D_{\beta,ij}^2\right).
```

For RHF, with ``d=\operatorname{diag}(D)``:

```math
F=H+2\operatorname{Diag}(Vd)-V\mathbin{\odot}D,
```

```math
E_{\mathrm{RHF}}=2\operatorname{Tr}(DH)+2d^TVd-
\sum_{ij}V_{ij}D_{ij}^2.
```

Sliced and target-residual backends obey the same trace convention. Their
storage and contraction conventions are detailed in [Interaction backends](@ref).

## Compact represented and producer energies

The compact typed solver reports `represented_energy`, the objective evaluated
with its working unit-cell representation. The optional
`producer_energy(result, producer)` operation evaluates the same final compact
determinant with a caller-supplied uncompressed density-density Hamiltonian and
explicit constant. Its total follows the UHF decomposition above (with equal
spin projectors in RHF).

These quantities have deliberately different names. Agreement between them is
a same-state energy comparison, not a bound on the occupied--virtual residual
of the producer Hamiltonian. The producer audit supplies no residual and does
not alter solver convergence status.

## Residual and stationarity

For one spin, a gauge-invariant occupied--virtual residual can be written

```math
R=(I-CC^T)FC.
```

Equivalently, the density commutator ``FD-DF`` vanishes at stationarity. A
small sweep-energy change does not bound either quantity in a soft response
direction. Occupied rotations within `C` change its column gauge but not its
projector or determinant; stationarity concerns occupied-to-unoccupied
rotations.
