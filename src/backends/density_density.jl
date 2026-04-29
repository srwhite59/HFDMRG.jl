# Density-density backend (V is an N x N matrix).

struct DensityDensityBackend{T<:AbstractMatrix}
    V::T
end

struct DDBlockState{T}
    ra::UnitRange{Int}
    raV::UnitRange{Int}
    pp::Array{T,3}
    Vpp::Array{T,3}
    Vijkl::Array{T,4}
end

struct DDWindow{T}
    VLL::Array{T,4}
    VLC::Array{T,3}
    VLR::Array{T,4}
    VCC::Array{T,2}
    VCR::Array{T,3}
    VRR::Array{T,4}
end

function vee_init_block(side, ra, raV, phi, backend::DensityDensityBackend)
    V = backend.V
    lb = length(ra)
    m = size(phi, 2)
    pp = [phi[k, i] * phi[k, j] for k = 1:lb, i = 1:m, j = 1:m]
    Vpp = contract(V[raV, ra], pp)
    Vijkl = contract(pp, contract(V[ra, ra], pp), true, false)
    DDBlockState(ra, raV, pp, Vpp, Vijkl)
end

function vee_absorb_block(side, vee_old::DDBlockState, cra, Phi_old, Phi_C, phi_new, raV_new,
    backend::DensityDensityBackend)
    V = backend.V
    n = size(phi_new, 1)
    m = size(phi_new, 2)
    pp = [phi_new[j, k] * phi_new[j, l] for j = 1:n, k = 1:m, l = 1:m]
    Vijkl = trans4(vee_old.Vijkl, Phi_old)
    lc = length(cra)
    if side == :right
        ppcc = pp[1:lc, :, :]
        VppC = trans2(vee_old.Vpp[cra, :, :], Phi_old)
    elseif side == :left
        ppcc = pp[cra, :, :]
        interabs, lr, cr, nonnull = intersectrange(vee_old.raV, cra)
        VppC = trans2(vee_old.Vpp[lr, :, :], Phi_old)
    else
        error("unknown block side: $side")
    end
    Vijkl += contract(ppcc, VppC, true, false)
    Vijkl += contract(VppC, ppcc, true, false)
    Vijkl += contract(ppcc, contract(V[cra, cra], ppcc), true, false)
    VC = contract(V[raV_new, cra], ppcc)
    if side == :right
        newVpp = trans2(vee_old.Vpp[raV_new, :, :], Phi_old)
    else
        interabs, lr, l2r, nonnull = intersectrange(vee_old.raV, raV_new)
        newVpp = trans2(vee_old.Vpp[lr, :, :], Phi_old)
    end
    Vpp = newVpp + VC
    ra_new = side == :right ? (cra[1]:vee_old.ra[end]) : (vee_old.ra[1]:cra[end])
    DDBlockState(ra_new, raV_new, pp, Vpp, Vijkl)
end

function vee_window(Lvee::DDBlockState, Rvee::DDBlockState, Cra, backend::DensityDensityBackend)
    V = backend.V
    Lra = Lvee.ra
    Rra = Rvee.ra
    lc = length(Cra)
    VLL = Lvee.Vijkl
    VCC = V[Cra, Cra]
    VRR = Rvee.Vijkl
    VLC = permutedims(Lvee.Vpp[1:lc, :, :], (2, 3, 1))
    VCR = Rvee.Vpp[Cra, :, :]
    ml, mr = size(VLL, 1), size(VRR, 1)
    lenR = length(Rra)
    idxL = Rra .- Lra[end]
    VLR = reshape(reshape(Lvee.Vpp[idxL, :, :], lenR, ml^2)' *
                  reshape(Rvee.pp, lenR, mr^2), ml, ml, mr, mr)
    DDWindow(VLL, VLC, VLR, VCC, VCR, VRR)
end

function vee_add_fock_r!(F, rho, win::DDWindow)
    ml, mc, mr = size(win.VLL, 1), size(win.VCC, 1), size(win.VRR, 1)
    sc, sr = ml, ml + mc
    if size(rho) != size(F)
        error("inconsistent sizes rho and F")
    end
    for i = 1:ml, j = 1:ml, l = 1:ml, k = 1:ml
        F[i, j] += rho[k, l] * (2 * win.VLL[i, j, k, l] - win.VLL[i, l, k, j])
    end
    for j = 1:mc, i = 1:mc
        ii, jj = sc + i, sc + j
        v = win.VCC[i, j]
        F[ii, ii] += rho[jj, jj] * 2 * v
        F[ii, jj] -= rho[ii, jj] * v
    end
    for j = 1:mr, l = 1:mr, k = 1:mr, i = 1:mr
        ii, jj, kk, ll = i + sr, j + sr, k + sr, l + sr
        F[ii, jj] += rho[kk, ll] * (2 * win.VRR[i, j, k, l] - win.VRR[i, l, k, j])
    end
    for j = 1:ml, l = 1:mr, k = 1:mr, i = 1:ml
        ii, jj, kk, ll = i, j, k + sr, l + sr
        v = win.VLR[i, j, k, l]
        F[ii, jj] += rho[kk, ll] * 2 * v
        F[kk, ll] += rho[ii, jj] * 2 * v
        F[ii, kk] -= rho[jj, ll] * v
        F[kk, ii] -= rho[ll, jj] * v
    end
    for i = 1:ml, j = 1:ml, k = 1:mc
        kk = k + sc
        v = win.VLC[i, j, k]
        F[i, j] += rho[kk, kk] * 2 * v
        F[kk, kk] += rho[i, j] * 2 * v
        F[i, kk] -= rho[j, kk] * v
        F[kk, i] -= rho[kk, j] * v
    end
    for i = 1:mr, j = 1:mr, k = 1:mc
        ii, jj, kk = i + sr, j + sr, k + sc
        v = win.VCR[k, i, j]
        F[ii, jj] += rho[kk, kk] * 2 * v
        F[kk, kk] += rho[ii, jj] * 2 * v
        F[ii, kk] -= rho[jj, kk] * v
        F[kk, ii] -= rho[kk, jj] * v
    end
end

function vee_add_fock!(Fup, Fdn, rhoup, rhodn, win::DDWindow)
    ml, mc, mr = size(win.VLL, 1), size(win.VCC, 1), size(win.VRR, 1)
    sc, sr = ml, ml + mc
    if size(rhoup) != size(Fup) || size(rhodn) != size(Fup)
        error("inconsistent sizes rho and F")
    end
    for i = 1:ml, j = 1:ml, l = 1:ml, k = 1:ml
        direct = (rhoup[k, l] + rhodn[k, l]) * win.VLL[i, j, k, l]
        Fup[i, j] += direct - rhoup[k, l] * win.VLL[i, l, k, j]
        Fdn[i, j] += direct - rhodn[k, l] * win.VLL[i, l, k, j]
    end
    for j = 1:mc, i = 1:mc
        ii, jj = sc + i, sc + j
        v = win.VCC[i, j]
        direct = (rhoup[jj, jj] + rhodn[jj, jj]) * v
        Fup[ii, ii] += direct
        Fdn[ii, ii] += direct
        Fup[ii, jj] -= rhoup[ii, jj] * v
        Fdn[ii, jj] -= rhodn[ii, jj] * v
    end
    for j = 1:mr, l = 1:mr, k = 1:mr, i = 1:mr
        ii, jj, kk, ll = i + sr, j + sr, k + sr, l + sr
        direct = (rhoup[kk, ll] + rhodn[kk, ll]) * win.VRR[i, j, k, l]
        Fup[ii, jj] += direct - rhoup[kk, ll] * win.VRR[i, l, k, j]
        Fdn[ii, jj] += direct - rhodn[kk, ll] * win.VRR[i, l, k, j]
    end
    for j = 1:ml, l = 1:mr, k = 1:mr, i = 1:ml
        ii, jj, kk, ll = i, j, k + sr, l + sr
        v = win.VLR[i, j, k, l]
        rkl = rhoup[kk, ll] + rhodn[kk, ll]
        rij = rhoup[ii, jj] + rhodn[ii, jj]
        Fup[ii, jj] += rkl * v
        Fdn[ii, jj] += rkl * v
        Fup[kk, ll] += rij * v
        Fdn[kk, ll] += rij * v
        Fup[ii, kk] -= rhoup[jj, ll] * v
        Fdn[ii, kk] -= rhodn[jj, ll] * v
        Fup[kk, ii] -= rhoup[ll, jj] * v
        Fdn[kk, ii] -= rhodn[ll, jj] * v
    end
    for i = 1:ml, j = 1:ml, k = 1:mc
        kk = k + sc
        v = win.VLC[i, j, k]
        rkk = rhoup[kk, kk] + rhodn[kk, kk]
        rij = rhoup[i, j] + rhodn[i, j]
        Fup[i, j] += rkk * v
        Fdn[i, j] += rkk * v
        Fup[kk, kk] += rij * v
        Fdn[kk, kk] += rij * v
        Fup[i, kk] -= rhoup[j, kk] * v
        Fdn[i, kk] -= rhodn[j, kk] * v
        Fup[kk, i] -= rhoup[kk, j] * v
        Fdn[kk, i] -= rhodn[kk, j] * v
    end
    for i = 1:mr, j = 1:mr, k = 1:mc
        ii, jj, kk = i + sr, j + sr, k + sc
        v = win.VCR[k, i, j]
        rkk = rhoup[kk, kk] + rhodn[kk, kk]
        rij = rhoup[ii, jj] + rhodn[ii, jj]
        Fup[ii, jj] += rkk * v
        Fdn[ii, jj] += rkk * v
        Fup[kk, kk] += rij * v
        Fdn[kk, kk] += rij * v
        Fup[ii, kk] -= rhoup[jj, kk] * v
        Fdn[ii, kk] -= rhodn[jj, kk] * v
        Fup[kk, ii] -= rhoup[kk, jj] * v
        Fdn[kk, ii] -= rhodn[kk, jj] * v
    end
end
