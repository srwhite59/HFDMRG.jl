"""
Sliced-basis backend skeleton.

This backend will eventually manage interaction data in a slice-adapted basis.
For now, it only validates obvious inputs and throws a clear error.
"""
struct SlicedBasisBackend{TL, TV}
    layout::TL
    Vee::TV
end

struct SlicedBasisWindow end

function _check_side(side)
    (side === :left || side === :right) || error("side must be :left or :right")
end

function _check_range(ra, name)
    isa(ra, UnitRange) || error("$name must be UnitRange")
    eltype(ra) <: Integer || error("$name must be integer range")
    isempty(ra) && error("$name must be non-empty")
end

function _sliced_not_implemented()
    error("SlicedBasisBackend not implemented yet")
end

function vee_init_block(side, ra, raV, phi, backend::SlicedBasisBackend)
    _check_side(side)
    _check_range(ra, "ra")
    _check_range(raV, "raV")
    _sliced_not_implemented()
end

function vee_absorb_block(side, vee_old, cra, Phi_old, Phi_C, phi_new, raV_new,
    backend::SlicedBasisBackend)
    _check_side(side)
    _check_range(cra, "cra")
    _check_range(raV_new, "raV_new")
    _sliced_not_implemented()
end

function vee_window(Lvee, Rvee, Cra, backend::SlicedBasisBackend)
    _check_range(Cra, "Cra")
    _sliced_not_implemented()
end

function vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::SlicedBasisWindow)
    _sliced_not_implemented()
end

function vee_add_fock_r!(F, rho, win::SlicedBasisWindow)
    _sliced_not_implemented()
end
