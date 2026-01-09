"""
Initialize interaction backend state for a new block.

Arguments:
- side::Symbol: :left or :right. :left means the block grows from low indices
  (ra is 1:l). :right means the block grows from high indices (ra is r:N).
- ra::UnitRange{Int}: physical site range represented by the block.
- raV::UnitRange{Int}: physical site range outside the block that couples via V.
  For :left blocks, raV == ra[end]+1:N. For :right blocks, raV == 1:ra[1]-1.
- phi::AbstractMatrix: block basis on ra with orthonormal columns, size
  (length(ra), m).
- Vee: backend-specific interaction object (density-density uses an N x N matrix).

Returns backend state for the block. The density-density backend caches
pp (ra x m x m), Vpp (raV x m x m), and Vijkl (m x m x m x m).
"""
vee_init_block(side, ra, raV, phi, Vee) =
    error("vee_init_block not implemented for backend")

"""
Update backend state when absorbing a center chunk into an existing block.

Arguments:
- side::Symbol: :left or :right, indicating which block is being grown.
- vee_old: backend state for the old block.
- cra::UnitRange{Int}: physical center range being absorbed (absolute indices).
- Phi_old::AbstractMatrix: part of the isometry for the old block basis; for
  :left it maps the old-left block into the new basis, for :right it maps the
  old-right block into the new basis. Size (m_old, m_new).
- Phi_C::AbstractMatrix: part of the isometry for the absorbed center chunk,
  size (length(cra), m_new).
- phi_new::AbstractMatrix: new block basis on the expanded range, size
  (length(new ra), m_new).
- raV_new::UnitRange{Int}: new outside range after absorption.
- Vee: backend-specific interaction object.

Returns new backend state. The core assumes Phi_old/Phi_C come from an SVD of
occupied orbitals on (old block + center).
"""
vee_absorb_block(side, vee_old, cra, Phi_old, Phi_C, phi_new, raV_new, Vee) =
    error("vee_absorb_block not implemented for backend")

"""
Build a lightweight window interaction object for the current L-C-R partition.

Arguments:
- Lvee/Rvee: backend states for left/right blocks.
- Cra::UnitRange{Int}: physical center range between L and R blocks
  (absolute indices).
- Vee: backend-specific interaction object.

Returns win, a lightweight object containing interaction data in the current
superblock basis [L_basis; C_physical; R_basis]. For density-density backends,
win holds VLL (ml,ml,ml,ml), VLC (ml,ml,mc), VLR (ml,ml,mr,mr),
VCC (mc,mc), VCR (mc,mr,mr), and VRR (mr,mr,mr,mr).
"""
vee_window(Lvee, Rvee, Cra, Vee) =
    error("vee_window not implemented for backend")

"""
Add interaction mean-field contributions to unrestricted (UHF) Fock matrices.

Arguments:
- Fup/Fdn: Fock matrices for spin up/down in the current superblock basis.
- rhoup/rhodn: density matrices from occupied spin orbitals in that basis
  (psi_up * psi_up', psi_dn * psi_dn').
- win: window interaction object from vee_window.

The direct term uses the total density (rhoup + rhodn); exchange uses like-spin
densities. The core computes energy as
0.5*tr(rhoup*(Fup+H1B)) + 0.5*tr(rhodn*(Fdn+H1B)), so conventions must match.
"""
vee_add_fock!(Fup, Fdn, rhoup, rhodn, win) =
    error("vee_add_fock! not implemented for backend")

"""
Add interaction mean-field contributions to a restricted (RHF) Fock matrix.

Arguments:
- F::AbstractMatrix: Fock matrix in the current superblock basis (modified
  in place).
- rho::AbstractMatrix: density matrix from occupied orbitals in that basis
  (psi * psi').
- win: window interaction object from vee_window.

For density-density interactions this corresponds to F += (2J - K). The core
computes energy as tr(rho*(F+H1B)), so conventions must match.
"""
vee_add_fock_r!(F, rho, win) =
    error("vee_add_fock_r! not implemented for backend")
