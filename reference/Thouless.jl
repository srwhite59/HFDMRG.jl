module Thouless

using LinearAlgebra
using Optim
using FastEigs, Util

# Lightweight trace helper for 2D matrices.
trAB(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64}) = dot(vec(A), vec(B))

# J/K builder interface for DiagYlm-style Coulomb/exchange builds.
# build_J!(J, Ptot) and build_K!(K, Pspin) should overwrite outputs in-place.
struct JKOperators
    build_J!::Function
    build_K!::Function
end

struct JKWorkspace
    Ptot::Matrix{Float64}
    J::Matrix{Float64}
    Kup::Matrix{Float64}
    Kdn::Matrix{Float64}
end

JKWorkspace(nn::Int) = JKWorkspace(zeros(nn, nn), zeros(nn, nn), zeros(nn, nn), zeros(nn, nn))

# V_{ijkl}^{nn'} -> cdag(i,n) cdag(j,n') c(k,n') c(l,n)
# hamup/dn: nj,npnt,nj,npnt
# V: nj,nj,nj,nj,npnt,npnt
# rho: nj,npnt,nj,npnt
# F: (nj * npnt), (nj * npnt)

# getF and methods
#--------------------------------------------
function build_fock!(rho_up::AbstractArray{Float64,4}, rho_dn::AbstractArray{Float64,4},
		V::AbstractArray{Float64,6}, F_up::Array{Float64,2}, F_dn::Array{Float64,2})

	nj = size(rho_up, 1)
	npnt = size(rho_up, 2)

	tmp_Fup = zeros(nj, npnt, nj, npnt)
	tmp_Fdn = zeros(nj, npnt, nj, npnt)
	@inbounds for n = 1:npnt, np = 1:npnt, j = 1:nj, l = 1:nj, k = 1:nj, i = 1:nj
		# Terms diagonal in slices
		tmp_Fup[i,n,j,n] += V[i,l,k,j,n,np] * (rho_up[l,np,k,np] + rho_dn[l,np,k,np])
		tmp_Fdn[i,n,j,n] += V[i,l,k,j,n,np] * (rho_dn[l,np,k,np] + rho_up[l,np,k,np])

		# Remaining Terms
		tmp_Fup[i,n,j,np] -= V[i,k,j,l,n,np] * rho_up[k,np,l,n]
		tmp_Fdn[i,n,j,np] -= V[i,k,j,l,n,np] * rho_dn[k,np,l,n]
	end

	F_up[:,:] += reshape(tmp_Fup, (nj * npnt), (nj * npnt))
	F_dn[:,:] += reshape(tmp_Fdn, (nj * npnt), (nj * npnt))
	symmetrize!(F_up)
    symmetrize!(F_dn)
end

function build_fock!(rhomat_up::AbstractMatrix{Float64}, rhomat_dn::AbstractMatrix{Float64},
		V::AbstractArray{Float64,6}, F_up::Array{Float64,2}, F_dn::Array{Float64,2})

	nj = size(V,1)
	npnt = size(V,6)

	rho_up = reshape(rhomat_up, nj, npnt, nj, npnt)
	rho_dn = reshape(rhomat_dn, nj, npnt, nj, npnt)
	build_fock!(rho_up, rho_dn, V, F_up, F_dn)
end

function build_fock!(rhomat_up::AbstractMatrix{Float64}, rhomat_dn::AbstractMatrix{Float64},
		v2::AbstractMatrix{Float64}, F_up::Array{Float64,2}, F_dn::Array{Float64,2})

	nn = size(rhomat_up, 1)
	rhou = diag(rhomat_up)
	rhod = diag(rhomat_dn)
	occ = rhou .+ rhod
	pot = v2 * occ
	@inbounds for i = 1:nn
		F_up[i, i] += pot[i]
		F_dn[i, i] += pot[i]
	end
	@inbounds for k = 1:nn, i = 1:nn
		F_up[i, k] -= rhomat_up[i, k] * v2[i, k]
		F_dn[i, k] -= rhomat_dn[i, k] * v2[i, k]
	end
	symmetrize!(F_up)
    symmetrize!(F_dn)
end

function build_fock!(rho_up::AbstractArray{Float64,4}, rho_dn::AbstractArray{Float64,4},
		v2::AbstractMatrix{Float64}, F_up::Array{Float64,2}, F_dn::Array{Float64,2})

	nj = size(rho_up, 1)
	npnt = size(rho_up, 2)
	rho_up_mat = reshape(rho_up, nj * npnt, nj * npnt)
	rho_dn_mat = reshape(rho_dn, nj * npnt, nj * npnt)
	build_fock!(rho_up_mat, rho_dn_mat, v2, F_up, F_dn)
end

function build_fock!(rho_up::AbstractArray{Float64,4}, rho_dn::AbstractArray{Float64,4},
		jk::JKOperators, F_up::Array{Float64,2}, F_dn::Array{Float64,2};
		work::Union{Nothing,JKWorkspace}=nothing)

	nj = size(rho_up, 1)
	npnt = size(rho_up, 2)
	rhomat_up = reshape(rho_up, nj * npnt, nj * npnt)
	rhomat_dn = reshape(rho_dn, nj * npnt, nj * npnt)
	build_fock!(rhomat_up, rhomat_dn, jk, F_up, F_dn; work=work)
end

function build_fock!(rhomat_up::AbstractMatrix{Float64}, rhomat_dn::AbstractMatrix{Float64},
		jk::JKOperators, F_up::Array{Float64,2}, F_dn::Array{Float64,2};
		work::Union{Nothing,JKWorkspace}=nothing)

	if work === nothing
		Ptot = rhomat_up + rhomat_dn
		J = similar(rhomat_up)
		Kup = similar(rhomat_up)
		Kdn = similar(rhomat_up)
	else
		Ptot = work.Ptot
		@. Ptot = rhomat_up + rhomat_dn
		J = work.J
		Kup = work.Kup
		Kdn = work.Kdn
	end
	jk.build_J!(J, Ptot)
	jk.build_K!(Kup, rhomat_up)
	jk.build_K!(Kdn, rhomat_dn)

	F_up .+= J
	F_up .-= Kup
	F_dn .+= J
	F_dn .-= Kdn
	symmetrize!(F_up)
    symmetrize!(F_dn)
end
#--------------------------------------------
#--------------------------------------------


# HFenergy and methods
#--------------------------------------------
function hf_energy(rho_up::AbstractArray{Float64,4}, rho_dn::AbstractArray{Float64,4},
		hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		V::AbstractArray{Float64,6})

	nj = size(rho_up, 1)
	npnt = size(rho_up, 2)
	energy = 0
	# Single-patricle terms
	@inbounds for n = 1:npnt, np = 1:npnt, j = 1:nj, i = 1:nj
		energy += (hamup[i,n,j,np] * rho_up[i,n,j,np]) +
			(hamdn[i,n,j,np] * rho_dn[i,n,j,np])
	end

	# Two-particle terms
	@inbounds for n = 1:npnt, np = 1:npnt, j = 1:nj, k = 1:nj, l = 1:nj, i = 1:nj
		energy += 0.5 * V[i,j,k,l,n,np] * (2 * (rho_up[i,n,l,n] * rho_dn[j,np,k,np]) +
			(rho_up[j,np,k,np] * rho_up[i,n,l,n]) - (rho_up[i,n,k,np] * rho_up[j,np,l,n]) +
			(rho_dn[j,np,k,np] * rho_dn[i,n,l,n]) - (rho_dn[i,n,k,np] * rho_dn[j,np,l,n]))
	end
	energy
end

function hf_energy(rhomat_up::AbstractMatrix{Float64}, rhomat_dn::AbstractMatrix{Float64},
		hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		V::AbstractArray{Float64,6})

	nj = size(V,1)
	npnt = size(V,6)

	rho_up = reshape(rhomat_up, nj, npnt, nj, npnt)
	rho_dn = reshape(rhomat_dn, nj, npnt, nj, npnt)
	hf_energy(rho_up, rho_dn, hamup, hamdn, V)
end

function hf_energy(rho_up::AbstractArray{Float64,4}, rho_dn::AbstractArray{Float64,4},
		ham::AbstractArray{Float64,4}, V::AbstractArray{Float64,6})
	hf_energy(rho_up, rho_dn, ham, ham, V)
end

function hf_energy(rhomat_up::AbstractMatrix{Float64}, rhomat_dn::AbstractMatrix{Float64},
		ham::AbstractArray{Float64,4}, V::AbstractArray{Float64,6})

	nj = size(V,1)
	npnt = size(V,6)

	rho_up = reshape(rhomat_up, nj, npnt, nj, npnt)
	rho_dn = reshape(rhomat_dn, nj, npnt, nj, npnt)
	hf_energy(rho_up, rho_dn, ham, ham, V)
end

function hf_energy(rhomat_up::AbstractMatrix{Float64}, rhomat_dn::AbstractMatrix{Float64},
		hammat_up::AbstractMatrix{Float64}, hammat_dn::AbstractMatrix{Float64},
		v2::AbstractMatrix{Float64})

	nn = size(rhomat_up, 1)
	energy = 0.0
	energy += trAB(rhomat_up, hammat_up)
	energy += trAB(rhomat_dn, hammat_dn)

	rhou = diag(rhomat_up)
	rhod = diag(rhomat_dn)
	v2u = v2 * rhou
	v2d = v2 * rhod
	energy += dot(rhou, v2d)
	energy += 0.5 * dot(rhou, v2u)
	energy += 0.5 * dot(rhod, v2d)
	@inbounds for k = 1:nn, i = 1:nn
		energy += 0.5 * v2[i, k] * (-rhomat_up[i, k] * rhomat_up[i, k] -
			rhomat_dn[i, k] * rhomat_dn[i, k])
	end
	energy
end

function hf_energy(rhomat_up::AbstractMatrix{Float64}, rhomat_dn::AbstractMatrix{Float64},
		hammat::AbstractMatrix{Float64}, v2::AbstractMatrix{Float64})
	hf_energy(rhomat_up, rhomat_dn, hammat, hammat, v2)
end

function hf_energy(rho_up::AbstractArray{Float64,4}, rho_dn::AbstractArray{Float64,4},
		hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		v2::AbstractMatrix{Float64})

	nj = size(hamup, 1)
	npnt = size(hamup, 2)
	nn = nj * npnt
	rhomat_up = reshape(rho_up, nn, nn)
	rhomat_dn = reshape(rho_dn, nn, nn)
	hammat_up = reshape(hamup, nn, nn)
	hammat_dn = reshape(hamdn, nn, nn)
	hf_energy(rhomat_up, rhomat_dn, hammat_up, hammat_dn, v2)
end

function hf_energy(rho_up::AbstractArray{Float64,4}, rho_dn::AbstractArray{Float64,4},
		ham::AbstractArray{Float64,4}, v2::AbstractMatrix{Float64})
	hf_energy(rho_up, rho_dn, ham, ham, v2)
end

function hf_energy(rho_up::AbstractArray{Float64,4}, rho_dn::AbstractArray{Float64,4},
		hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		jk::JKOperators; work::Union{Nothing,JKWorkspace}=nothing)

	nj = size(hamup, 1)
	npnt = size(hamup, 2)
	nn = nj * npnt
	rhomat_up = reshape(rho_up, nn, nn)
	rhomat_dn = reshape(rho_dn, nn, nn)
	hammat_up = reshape(hamup, nn, nn)
	hammat_dn = reshape(hamdn, nn, nn)
	hf_energy(rhomat_up, rhomat_dn, hammat_up, hammat_dn, jk; work=work)
end

function hf_energy(rhomat_up::AbstractMatrix{Float64}, rhomat_dn::AbstractMatrix{Float64},
		hammat_up::AbstractMatrix{Float64}, hammat_dn::AbstractMatrix{Float64},
		jk::JKOperators; work::Union{Nothing,JKWorkspace}=nothing)

	if work === nothing
		Ptot = rhomat_up + rhomat_dn
		J = similar(rhomat_up)
		Kup = similar(rhomat_up)
		Kdn = similar(rhomat_up)
	else
		Ptot = work.Ptot
		@. Ptot = rhomat_up + rhomat_dn
		J = work.J
		Kup = work.Kup
		Kdn = work.Kdn
	end
	jk.build_J!(J, Ptot)
	jk.build_K!(Kup, rhomat_up)
	jk.build_K!(Kdn, rhomat_dn)

	Eone = trAB(rhomat_up, hammat_up) + trAB(rhomat_dn, hammat_dn)
	EJ = 0.5 * trAB(Ptot, J)
	EK = -0.5 * (trAB(rhomat_up, Kup) + trAB(rhomat_dn, Kdn))
	return Eone + EJ + EK
end

function hf_energy(rho_up::AbstractArray{Float64,4}, rho_dn::AbstractArray{Float64,4},
		ham::AbstractArray{Float64,4}, jk::JKOperators; work::Union{Nothing,JKWorkspace}=nothing)
	hf_energy(rho_up, rho_dn, ham, ham, jk; work=work)
end

function hf_energy(rhomat_up::AbstractMatrix{Float64}, rhomat_dn::AbstractMatrix{Float64},
		ham::AbstractMatrix{Float64}, jk::JKOperators; work::Union{Nothing,JKWorkspace}=nothing)
	hf_energy(rhomat_up, rhomat_dn, ham, ham, jk; work=work)
end
#--------------------------------------------
#--------------------------------------------


# Standard HF with damping
#--------------------------------------------------
function hf_damped!(psi_up::Array{Float64,2}, psi_dn::Array{Float64,2},
		hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		V; maxiter::Int64=200)
    nj = size(hamup,1)
    npnt = size(hamup,2)
    nn = nj * npnt
    Nup = size(psi_up, 2)
    Ndn = size(psi_dn, 2)

    hammat_up = reshape(hamup, nn, nn)
    hammat_dn = reshape(hamdn, nn, nn)
    work = V isa JKOperators ? JKWorkspace(nn) : nothing

    lambda_up = 0.9
    lambda_dn = 0.89 # break symmetry
    psiupbest = psi_up
    psidnbest = psi_dn
    function getrho(psiup,psidn)
        pup = reshape(psiup,nj*npnt,Nup)
        pdn = reshape(psidn,nj*npnt,Ndn)
        rhoup = reshape(pup*pup',nj,npnt,nj,npnt)
        rhodn = reshape(pdn*pdn',nj,npnt,nj,npnt)
        rhoup,rhodn
    end
    function getenergy(psiup,psidn)
        rhoup,rhodn = getrho(psiup,psidn)
        if V isa JKOperators
            return hf_energy(rhoup, rhodn, hamup, hamdn, V; work=work)
        end
        return hf_energy(rhoup, rhodn, hamup, hamdn, V)
    end
    function getF(psiup,psidn)
        rhoup,rhodn = getrho(psiup,psidn)
        Fup,Fdn = copy(hammat_up),copy(hammat_dn)
        if V isa JKOperators
            build_fock!(rhoup, rhodn, V, Fup, Fdn; work=work)
        else
            build_fock!(rhoup, rhodn, V, Fup, Fdn)
        end
        Fup,Fdn
    end
    function geteigF(psiup,psidn)
        Fup,Fdn = getF(psiup,psidn)
        evals,vup = eig(Fup)
        evals,vdn = eig(Fdn)
        vup[:,1:Nup],vdn[:,1:Ndn]
    end
    minenergy = getenergy(psiupbest,psidnbest)
    @show minenergy
    vecup,vecdn = geteigF(psiupbest,psidnbest)
    for s = 1:maxiter
        done = false
        for lamiter = 1:10
            uppsi = hcat((1-lambda_up) * psiupbest,lambda_up * vecup)
            dnpsi = hcat((1-lambda_dn) * psidnbest,lambda_dn * vecdn)
            u,d,v = svd(uppsi)
            psiup = u[:,1:Nup]
            u,d,v = svd(dnpsi)
            psidn = u[:,1:Ndn]
            energy = getenergy(psiup,psidn)
            flush(stdout)
            if energy < minenergy
                psiupbest,psidnbest = psiup,psidn
                minenergy = energy
                vecup,vecdn = geteigF(psiup,psidn)
                #@show s,lamiter,lambda_up,energy,minenergy
                break
            elseif abs(energy-minenergy) < 1.0e-10
                done = true
                break
            else
                lambda_up *= 0.5
                lambda_dn *= 0.5
            end
        end
        done && break
    end
    psi_up .= psiupbest
    psi_dn .= psidnbest
    if Nup == Ndn
        de = det(psi_up' * psi_dn)
        println("det of psi_up' * psi_dn is ",de)
    end
        #rho_up = reshape(psiupbest * psiupbest', nj, npnt, nj, npnt)
        #rho_dn = reshape(psidnbest * psidnbest', nj, npnt, nj, npnt)
        #energyfinal = hf_energy(rho_up, rho_dn, hamup, hamdn, V)
        #@show energyfinal,minenergy
    minenergy
end

function hf_damped!(psi_up::Array{Float64,2}, psi_dn::Array{Float64,2},
		hammat_up::AbstractMatrix{Float64}, hammat_dn::AbstractMatrix{Float64},
		jk::JKOperators; maxiter::Int64=200)
	nn = size(hammat_up, 1)
	hamup4 = reshape(hammat_up, nn, 1, nn, 1)
	hamdn4 = reshape(hammat_dn, nn, 1, nn, 1)
	hf_damped!(psi_up, psi_dn, hamup4, hamdn4, jk; maxiter=maxiter)
end

function hf_damped!(psi_up::Array{Float64,2}, psi_dn::Array{Float64,2},
		hammat_up::AbstractMatrix{Float64}, hammat_dn::AbstractMatrix{Float64},
		v2::AbstractMatrix{Float64}; maxiter::Int64=200)
	nn = size(hammat_up, 1)
	hamup4 = reshape(hammat_up, nn, 1, nn, 1)
	hamdn4 = reshape(hammat_dn, nn, 1, nn, 1)
	hf_damped!(psi_up, psi_dn, hamup4, hamdn4, v2; maxiter=maxiter)
end
#--------------------------------------------------
#--------------------------------------------------


# HF with DIIS
#--------------------------------------------------
function hf_diis!(psi_up::Array{Float64,2}, psi_dn::Array{Float64,2},
		hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
        V; maxiter = 200, NDIIS = 5) # changes psi_up, psi_dn

	nj = size(hamup,1)
	npnt = size(hamup,2)
	nn = nj * npnt
	Nup = size(psi_up, 2)
	Ndn = size(psi_dn, 2)

    rho_up = psi_up * psi_up'
    rho_dn = psi_dn * psi_dn'
    work = V isa JKOperators ? JKWorkspace(nn) : nothing
    if V isa JKOperators
        energylast = energy = hf_energy(rho_up, rho_dn, hamup, hamdn, V; work=work)
    else
        energylast = energy = hf_energy(rho_up, rho_dn, hamup, hamdn, V)
    end
    @show energy

	# Damping factor
    lambda = 1.0

	# Build Fock operator with current rho
	hammat_up = reshape(hamup, nn, nn)
	hammat_dn = reshape(hamdn, nn, nn)
    Fup = deepcopy(hammat_up)
    Fdn = deepcopy(hammat_dn)
    if V isa JKOperators
        build_fock!(rho_up, rho_dn, V, Fup, Fdn; work=work)
    else
        build_fock!(rho_up, rho_dn, V, Fup, Fdn)
    end

	# Prepare DIIS
    startDIIS = 1
    errVec = zeros(2*nn*nn,NDIIS)
    FupVec = zeros(nn,nn,NDIIS)
    FdnVec = zeros(nn,nn,NDIIS)

	# Start SCF
    avtime = 0.0
    cnt = 0
    for s=1:maxiter
      start = time()
	  # Diagonalize Fock to get new density matrices
      if Nup <= nn/10
          ps = reshape(psi_up[:,1],nn)
          eigup = fasteigs!(Fup,psi_up,4,1e-12)
          vup = psi_up
      else
          eigup = eigfact(Fup)
          vup = eigup[:vectors][:,1:Nup]
      end
      if Ndn <= nn/10
          ps = reshape(psi_dn[:,1],nn)
          eigdn = fasteigs!(Fdn,psi_dn,4,1e-12)
          vdn = psi_dn
      else
          eigdn = eigfact(Fdn)
          vdn = eigdn[:vectors][:,1:Ndn]
      end
# New density matrix with damping
#=
      rho_up = (1-lambda) * rho_up + lambda * vup * vup'
      rho_dn = (1-lambda) * rho_dn + lambda * vdn * vdn'
=#
      rho_up = vup * vup'
      rho_dn = vdn * vdn'

      if V isa JKOperators
          energy = hf_energy(rho_up, rho_dn, hamup, hamdn, V; work=work)
      else
          energy = hf_energy(rho_up, rho_dn, hamup, hamdn, V)
      end
#     (s <= 3 || s == 10 || s%20 == 0) && (@show s,energy)
#     if energy > energylast
#         lambda *= 0.5
#         @show s,lambda
#     end
# Build Fock operator with current rho
	  hammat_up = reshape(hamup, nn, nn)
	  hammat_dn = reshape(hamdn, nn, nn)
      Fup = deepcopy(hammat_up)
      Fdn = deepcopy(hammat_dn)
      if V isa JKOperators
          build_fock!(rho_up, rho_dn, V, Fup, Fdn; work=work)
      else
          build_fock!(rho_up, rho_dn, V, Fup, Fdn)
      end
# Start DIIS
#     @show size(FupVec)
      for i = 1: NDIIS-1
        FupVec[:,:,i] = FupVec[:,:,i+1]
        FdnVec[:,:,i] = FdnVec[:,:,i+1]
        errVec[:,i] = errVec[:,i+1]
      end
      FupVec[:,:,NDIIS] = Fup
      FdnVec[:,:,NDIIS] = Fdn
      tmp1 = Fup*vup
      tmp2 = vup'*Fup
      comm = tmp1*vup' - vup*tmp2
      errVec[1:nn*nn,NDIIS] = reshape(comm,:,1)
      tmp1 = Fdn*vdn
      tmp2 = vdn'*Fdn
      comm = tmp1*vdn' - vdn*tmp2
      errVec[nn*nn+1:2*nn*nn,NDIIS] = reshape(comm,:,1)
      @show s,energy,maximum(abs.(errVec[:,NDIIS]))
#     @show typeof(FupVec), typeof(FdnVec), typeof(errVec)
      avtime += time() - start
      cnt += 1
# Construct the B-matrix in DIIS
      if ( s > NDIIS+startDIIS)
        BMat = ones(NDIIS+1,NDIIS+1)
        BMat[NDIIS+1,NDIIS+1] = 0.0
        rhs = zeros(NDIIS+1,1)
        rhs[NDIIS+1] = 1.0
        for i = 1: NDIIS, j = 1:NDIIS
          BMat[i,j] = dot(errVec[:,i],errVec[:,j])
        end
#       @show BMat
#       @show rhs
        coeff = BMat \ rhs
#       @show typeof(BMat), typeof(rhs), typeof(coeff)
# DIIS extrapolate
        Fup = zeros(nn,nn)
        Fdn = zeros(nn,nn)
        for i = 1: NDIIS
          Fup += coeff[i] * FupVec[:,:,i]
          Fdn += coeff[i] * FdnVec[:,:,i]
        end
      end
# Check convergence
      if (abs(energy-energylast) < 1.0e-9 && s > 4 ) || s == maxiter
          psi_up[:,:] = vup
          psi_dn[:,:] = vdn
          @show (s,energy)
          println("Finished")
          if Nup == Ndn
            de = det(vup' * vdn)
            println("det of vup' * vdn is ",de)
          end
          break
      end
      energylast = energy
      flush(stdout)
    end
    @show avtime
    @show cnt
    @show nn
    @show avtime/cnt
    energy
end

function hf_diis!(psi_up::Array{Float64,2}, psi_dn::Array{Float64,2},
		hammat_up::AbstractMatrix{Float64}, hammat_dn::AbstractMatrix{Float64},
        jk::JKOperators; maxiter = 200, NDIIS = 5)
	nn = size(hammat_up, 1)
	hamup4 = reshape(hammat_up, nn, 1, nn, 1)
	hamdn4 = reshape(hammat_dn, nn, 1, nn, 1)
	hf_diis!(psi_up, psi_dn, hamup4, hamdn4, jk; maxiter=maxiter, NDIIS=NDIIS)
end

function hf_diis!(psi_up::Array{Float64,2}, psi_dn::Array{Float64,2},
		hammat_up::AbstractMatrix{Float64}, hammat_dn::AbstractMatrix{Float64},
        v2::AbstractMatrix{Float64}; maxiter = 200, NDIIS = 5)
	nn = size(hammat_up, 1)
	hamup4 = reshape(hammat_up, nn, 1, nn, 1)
	hamdn4 = reshape(hammat_dn, nn, 1, nn, 1)
	hf_diis!(psi_up, psi_dn, hamup4, hamdn4, v2; maxiter=maxiter, NDIIS=NDIIS)
end
#--------------------------------------------------
#--------------------------------------------------



# HF Thouless
#--------------------------------------------------
F1(x) = cos.(sqrt.(x))
dF1(x) = - 0.5 * sin.(sqrt.(x))./sqrt.(x)
F2(x) = sin.(sqrt.(x))./sqrt.(x)
dF2(x) = 0.5 * ( cos.(sqrt.(x)) - sin.(sqrt.(x))./sqrt.(x) )./x

function mo_from_thouless(Kvo::Array{Float64,2})
    nv, no = size(Kvo)
    nn = no + nv
    L = Kvo' * Kvo
    val, vec = eigen(L)
    Doo = vec * Diagonal(F1(val)) * vec'
    Dvo = Kvo * (vec * Diagonal(F2(val)) * vec')
    return [Doo;Dvo]
end

function thouless_energy(Kup::Array{Float64,2}, Kdn::Array{Float64,2},
		hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		V; work::Union{Nothing,JKWorkspace}=nothing)

    psi_up = mo_from_thouless(Kup)
    psi_dn = mo_from_thouless(Kdn)
    rho_up = psi_up * psi_up'
    rho_dn = psi_dn * psi_dn'

    if V isa JKOperators
        return hf_energy(rho_up, rho_dn, hamup, hamdn, V; work=work)
    end
    return hf_energy(rho_up, rho_dn, hamup, hamdn, V)
end

function getdFmat(lambda,F,dF)
    n = length(lambda)
    cutoff = 1e-6
    dFmat = zeros(n,n)
    Fi = F(lambda)
    dFi = dF(lambda)
    for i = 1:n
    for j = i:n
        if (abs(lambda[i]-lambda[j])<cutoff)
            dFmat[i,j] = (dFi[i] + dFi[j])/2.0
        else
            dFmat[i,j] = (Fi[i]-Fi[j])/(lambda[i]-lambda[j])
        end
        dFmat[j,i] = dFmat[i,j]
    end
    end
    return dFmat
end

function gradpsi_to_gradk(Gpsi::Array{Float64,2}, Kvo::Array{Float64,2})
    nv, no = size(Kvo)
    nn = no + nv
    L = Kvo' * Kvo
    val, vec = eigen(L)
    dF1mat = getdFmat(val,F1,dF1)
    dF2mat = getdFmat(val,F2,dF2)

	# first term = K * (M1+M1')
    Goo = Gpsi[1:no,1:no]
    M1 = vec*((vec'*Goo*vec).*dF1mat)*vec'
    GK = Kvo*(M1+M1')

	# second term = Gvo*F2 + K * (M2+M2')
    Gvo = Gpsi[no+1:nn,1:no]
    GK += Gvo * (vec * Diagonal(F2(val)) * vec')
    M2 = vec*((vec'*(Kvo'*Gvo)*vec).*dF2mat)*vec'
    GK += Kvo*(M2+M2')
    return GK
end

function thouless_energy_grad(Kup::Array{Float64,2}, Kdn::Array{Float64,2},
		hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		V; work::Union{Nothing,JKWorkspace}=nothing)

	# d E / d D_pi =  ( (F+F') * D )_pi
    psiup = mo_from_thouless(Kup)
    psidn = mo_from_thouless(Kdn)
    rhoup = psiup * psiup'
    rhodn = psidn * psidn'

	# build Fock
	N = size(hamup,1) * size(hamup,2)
	hamup_mat = reshape(hamup, N, N)
	hamdn_mat = reshape(hamdn, N, N)
    Fup = deepcopy(hamup_mat)
    Fdn = deepcopy(hamdn_mat)
    if V isa JKOperators
        build_fock!(rhoup, rhodn, V, Fup, Fdn; work=work)
    else
        build_fock!(rhoup, rhodn, V, Fup, Fdn)
    end
    dEdpsiup = (Fup+Fup') * psiup
    dEdpsidn = (Fdn+Fdn') * psidn

	# convert dE/dD to dE/dK
    dEdKup = gradpsi_to_gradk(dEdpsiup,Kup)
    dEdKdn = gradpsi_to_gradk(dEdpsidn,Kdn)
    return (dEdKup,dEdKdn)
end

function pack_thouless(Kup::Array{Float64,2}, Kdn::Array{Float64,2})
    return [reshape(Kup,:,1);reshape(Kdn,:,1)]
end

function unpack_thouless(x,noup,nvup,nodn,nvdn)
    Kup = reshape(x[1:nvup*noup],nvup,noup)
    Kdn = reshape(x[nvup*noup+1:end],nvdn,nodn)
    return (Kup,Kdn)
end

function funVal(x, hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		V, noup, nvup, nodn, nvdn; work::Union{Nothing,JKWorkspace}=nothing)
	Kup, Kdn = unpack_thouless(x,noup,nvup,nodn,nvdn)
    return thouless_energy(Kup, Kdn, hamup, hamdn, V; work=work)
end

function funGrad(x, hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		V, noup, nvup, nodn, nvdn; work::Union{Nothing,JKWorkspace}=nothing)
	Kup, Kdn = unpack_thouless(x,noup,nvup,nodn,nvdn)
    Gup, Gdn = thouless_energy_grad(Kup, Kdn, hamup, hamdn, V; work=work)
    return G = pack_thouless(Gup,Gdn)
end

function thouless_from_mo(psi)
    nn, no = size(psi)
    nv = nn - no
    Doo = psi[1:no,1:no]
    Dvo = psi[no+1:nn,1:no]
    F = svd(Doo)
    theta = acos.(F.S)
    return Dvo * ( F.V * Diagonal(theta./sin.(theta)) * F.U' )
end

function hf_thouless!(psi_up::Array{Float64,2}, psi_dn::Array{Float64,2},
		hamup::AbstractArray{Float64,4}, hamdn::AbstractArray{Float64,4},
		V; maxiter::Int64=200)
	# Get initial K from psi
	Kup = thouless_from_mo(psi_up)
	Kdn = thouless_from_mo(psi_dn)

	# Wrapper of E and dE/dK
	nvup, noup = size(Kup)
	nvdn, nodn = size(Kdn)
	work = nothing
	if V isa JKOperators
		nn = size(hamup, 1) * size(hamup, 2)
		work = JKWorkspace(nn)
	end
	funValLoc(x) = funVal(x, hamup, hamdn, V, noup, nvup, nodn, nvdn; work=work)
	funGradLoc(x) = funGrad(x, hamup, hamdn,  V, noup, nvup, nodn, nvdn; work=work)
	energy = 0.0

	# Optimize with L-BFGS
	res = optimize(funValLoc, funGradLoc,
                   pack_thouless(Kup,Kdn),
                   method = LBFGS(),
                   g_tol = 1e-6,
                   iterations = maxiter,
                   show_trace = true;
                   inplace = false)
    @show summary(res)
    @show Optim.converged(res)
    @show Optim.iterations(res), Optim.f_calls(res), Optim.g_calls(res)
    energy = Optim.minimum(res)
    x = Optim.minimizer(res)
    Kup, Kdn = unpack_thouless(x,noup,nvup,nodn,nvdn)
    psi_up[:,:] = mo_from_thouless(Kup)
    psi_dn[:,:] = mo_from_thouless(Kdn)
    return energy
end

function hf_thouless!(psi_up::Array{Float64,2}, psi_dn::Array{Float64,2},
		hammat_up::AbstractMatrix{Float64}, hammat_dn::AbstractMatrix{Float64},
		jk::JKOperators; maxiter::Int64=200)
	nn = size(hammat_up, 1)
	hamup4 = reshape(hammat_up, nn, 1, nn, 1)
	hamdn4 = reshape(hammat_dn, nn, 1, nn, 1)
	hf_thouless!(psi_up, psi_dn, hamup4, hamdn4, jk; maxiter=maxiter)
end

function hf_thouless!(psi_up::Array{Float64,2}, psi_dn::Array{Float64,2},
		hammat_up::AbstractMatrix{Float64}, hammat_dn::AbstractMatrix{Float64},
		v2::AbstractMatrix{Float64}; maxiter::Int64=200)
	nn = size(hammat_up, 1)
	hamup4 = reshape(hammat_up, nn, 1, nn, 1)
	hamdn4 = reshape(hammat_dn, nn, 1, nn, 1)
	hf_thouless!(psi_up, psi_dn, hamup4, hamdn4, v2; maxiter=maxiter)
end
#--------------------------------------------------
#--------------------------------------------------

# V: nj,nj,nj,nj,Nb,Nb,  nj == Nperblock,  npnt == Nb
function block_vee(V,blocks,blockvecs)
    Nb = length(blocks)
    nj = size(blockvecs[1],2)
    Vb = zeros(nj,nj,nj,nj,Nb,Nb)
    lba = [length(blocks[a]) for a=1:Nb]
    UU = [ reshape([blockvecs[a][i,m]*blockvecs[a][i,n] for i=1:lba[a],m=1:nj,n=1:nj],lba[a],nj^2) 
                        for a=1:Nb]
    @views for a=1:Nb, ap=1:Nb
        v = reshape(UU[a]' * V[blocks[a],blocks[ap]] * UU[ap],nj,nj,nj,nj)
        Vb[:,:,:,:,a,ap] = permutedims(v,[1,3,4,2])
    end
# V_{ijkl}^{nn'} -> cdag(i,n) cdag(j,n') c(k,n') c(l,n)
    Vb
end

# Hb: nj,Nb,nj,Nb
function block_h1(H,blocks,blockvecs)
    Nb = length(blocks)
    nj = size(blockvecs[1],2)
    Hb = zeros(nj,nj,Nb,Nb)
    @views for a=1:Nb, ap=1:Nb
        #term = blockvecs[a]' * H[blocks[a],blocks[ap]] * blockvecs[ap]
        Hb[:,:,a,ap] = blockvecs[a]' * H[blocks[a],blocks[ap]] * blockvecs[ap]
    end
    Hb = permutedims(Hb,[1,3,2,4])
end

# psi: nj,Nb,Nup
function block_orbitals(psi,blocks,blockvecs)
    Nb = length(blocks)
    nj = size(blockvecs[1],2)
    nup = size(psi,2)
    psib = zeros(nj*Nb,nup)
    for a=1:Nb
        psib[(a-1)*nj+1:a*nj,:] = blockvecs[a]' * psi[blocks[a],:]
    end
    psib
end

# psi: nj,Nb,Nup
function unblock_orbitals(psib,blocks,blockvecs)
    Nb = length(blocks)
    nj = size(blockvecs[1],2)
    nup = size(psib,2)
    nn = sum(length.(blocks))
    psi = zeros(nn,nup)
    for a=1:Nb
        psi[blocks[a],:] = blockvecs[a] * psib[(a-1)*nj+1:a*nj,:]
    end
    psi
end

export JKOperators, JKWorkspace, build_fock!, hf_energy, hf_damped!, hf_diis!, hf_thouless!,
       block_vee, block_h1, block_orbitals, unblock_orbitals

end
