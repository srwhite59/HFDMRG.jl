# Backend API: interaction-specific code plugs into the generic sweep engine.

# Initialize backend state for a new block with basis `phi` on orbital/site range `ra`.
vee_init_block(side, ra, raV, phi, Vee) =
    error("vee_init_block not implemented for backend")

# Update backend state when absorbing center range `cra` into an existing block.
vee_absorb_block(side, vee_old, cra, Phi_old, Phi_C, phi_new, raV_new, Vee) =
    error("vee_absorb_block not implemented for backend")

# Build a lightweight window object for the current L-C-R partition.
vee_window(Lvee, Rvee, Cra, Vee) =
    error("vee_window not implemented for backend")

# Add interaction mean-field to unrestricted Fock matrices in-place.
vee_add_fock!(Fup, Fdn, rhoup, rhodn, win) =
    error("vee_add_fock! not implemented for backend")

# Add interaction mean-field to restricted Fock matrix in-place.
vee_add_fock_r!(F, rho, win) =
    error("vee_add_fock_r! not implemented for backend")
