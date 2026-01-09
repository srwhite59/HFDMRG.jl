module HF_dmrg

using LinearAlgebra, HFnn, Util   #, BasisGen, AutoVectors, Arpack 

mutable struct LRBlock		# left or right block
    ra::UnitRange{Int64}	# block represents this range, i.e. 1:l for left block
    raH1::UnitRange{Int64}	# range connected to by H1 outside block, i.e. l+1:l+p for some p
    raV::UnitRange{Int64}	# range for all outside sites, connected to by V, i.e. l+1:N
    m::Int64	# number of vectors phi kept
    phi		# ra x m
    H1ij	# m x m
    H1phi	# raH1 x m
    Vijkl	# m x m x m x m
    pp		# ra x m x m
    Vpp		# raV x m x m
end
                                        # Input: ranges for A and B in absolute sites
function intersectrange(Ara,Bra)        # Output: intersection in abs. sites, as part of A range, then B, then isempty?
    first,last = max(Ara[1],Bra[1]),min(Ara[end],Bra[end])
    absra = first:last
    ara = (absra) .- (Ara[1]-1)
    bra = (absra) .- (Bra[1]-1)
    absra,ara,bra,(last >= first)	# last return value is true or false--is there an intersection?
end

function contract(A,B,transpA=false,transpB=false)	# contract over last and first index
    As,Bs = size(A),size(B)
    @views begin
	AA = reshape(A,prod(As[1:end-1]),As[end])
	transpA && (AA = reshape(A,As[1],prod(As[2:end])))
	BB = reshape(B,Bs[1],prod(Bs[2:end]))
	transpB && (BB = reshape(B,prod(Bs[1:end-1]),Bs[end]))
	if !transpA && !transpB
	    res = reshape(AA*BB,As[1:end-1]...,Bs[2:end]...)
	elseif transpA && !transpB
	    res = reshape(AA'*BB,As[2:end]...,Bs[2:end]...)
	elseif !transpA && transpB
	    res = reshape(AA*BB',As[1:end-1]...,Bs[1:end-1]...)
	elseif transpA && transpB
	    res = reshape(AA'*BB',As[2:end]...,Bs[1:end-1]...)
	end
    end
    res
end

function trans2(V,U)	# V has two final indices needing transforming by U, plus some initial indices
    n = length(size(V))
    swaplast = [1:n-2...,n,n-1]
    V2 = permutedims(contract(V,U),swaplast)
    permutedims(contract(V2,U),swaplast)
end

function trans4(V,U)	# V has exactly 4 identical indices needing transforming by U
    V2 = contract(U,contract(V,U),true,false)
    V2 = permutedims(V2,[2,1,4,3])
    V4 = contract(U,contract(V2,U),true,false)
    permutedims(V4,[2,1,4,3])
end

function transrangeleft(A,ra::UnitRange,O)	# transform left index of A in range ra by O
    n = size(A,1)
    oa = O * A[ra,:]
    ra[1] == 1 && return vcat(oa,A[ra[end]+1:n,:])
    ra[end] == n && return vcat(A[1:ra[1]-1,:],oa)
    vcat(A[1:ra[1]-1,:],oa,A[ra[end]+1:n,:])
end
function transrangeleft(A,l::Integer,O)	# automatically set ra
    transrangeleft(A,l:l+size(O,2)-1,O)
end

function getraH1(ra,firstindsH,finalindsH)	# get the range of H1 outside block with range ra
    if ra[1] == 1
	return ra[end]+1:maximum(finalindsH[ra])
    else
	return minimum(firstindsH[ra]):ra[1]-1
    end
end
function getraV(ra,N)	# get the range of V outside block with range ra
    ra[1] != 1 && (return 1:ra[1]-1)
    ra[end]+1:N
end

function getphi(psira)   #old way Matrix(qr(psi[range,:]).Q)
    u,d,v = svd(psira)
    i = findlast(y->y>1.0e-10,d)
    i == nothing && (i = 1)
    mkeep = max(1,i)
    u[:,1:mkeep],mkeep
end

function getblock(range,psi,H,V,firstindsH,finalindsH,dofull)
    N,m = size(psi)
    lb = length(range)
    raH1 = getraH1(range,firstindsH,finalindsH)
    raV = getraV(range,N)
    phi,m = getphi(psi[range,:])    
    if dofull
	H1ij = phi' * H[range,range] * phi
	H1phi = H[raH1,range] * phi

	#@inbounds pp = [phi[k,i] * phi[k,j] for k=1:lb,i=1:m,j=1:m]
	pp = [phi[k,i] * phi[k,j] for k=1:lb,i=1:m,j=1:m]

	Vpp = contract(V[raV,range],pp)
	Vijkl = contract(pp,contract(V[range,range],pp),true,false)
    else
	H1phi = pp = Vpp = Vijkl = 0
    end
    LRBlock(range,raH1,raV,m,phi,H1ij,H1phi,Vijkl,pp,Vpp)
end

# add C block (defined by range cra) to right-hand block
function addblockright(cra,O,rblock,N,H,V,firstindsH)	# O contains the rotation,  lc+m x m, column-orthogonal
    ra = cra[1]:rblock.ra[end]			# range of new block, absolute sites
    mold = rblock.m
    m, n = size(O,2), length(ra)				# size of new block
    raH1 = minimum(firstindsH[ra]):cra[1]-1
    raV = 1:cra[1]-1
    lc = length(cra)
    PhiC, PhiR = O[1:lc,:], O[lc+1:end,:]
    phi = vcat(PhiC,rblock.phi * PhiR)
    Ophi = phi' * phi
    phi = phi * inv(sqrt(Symmetric(Ophi)))
    H1ij = PhiR' * rblock.H1ij * PhiR + PhiC' * H[cra,cra] * PhiC

    interabs,cr,rr,nonnull = intersectrange(cra,rblock.raH1)
    offterm = PhiC[cr,:]' * rblock.H1phi[rr,:] * PhiR
    H1ij += offterm + offterm'

    H1phi = H[raH1,cra] * PhiC
    interabs,cr,rr,nonnull = intersectrange(raH1,rblock.raH1)
    nonnull && (H1phi[cr,:] += rblock.H1phi[rr,:] * PhiR)

    Vijkl = trans4(rblock.Vijkl,PhiR)

    n > size(phi,1) && error("bad size phi")
    m > size(phi,2) && error("bad size phi")

    # This way of calculating pp has been optimized
    phit = phi'
    @inbounds pp = [phit[k,j] * phit[l,j] for j=1:n, k=1:m, l=1:m]

    ppcc = pp[1:lc,:,:]
    VppC = trans2(rblock.Vpp[cra,:,:],PhiR)
    Vijkl += contract(ppcc,VppC,true,false)
    Vijkl += contract(VppC,ppcc,true,false)
    Vijkl += contract(ppcc,contract(V[cra,cra],ppcc),true,false)
    VC = contract(V[raV,cra],ppcc)
    newVpp = trans2(rblock.Vpp[raV,:,:],PhiR)
    Vpp = newVpp + VC
    LRBlock(ra,raH1,raV,m,phi,H1ij,H1phi,Vijkl,pp,Vpp)
end

# add C block (defined by range cra) to left-hand block
function addblockleft(cra,O,lblock,N,H,V,finalindsH)	# O contains the rotation,  m+lc x m, column-orthogonal
    mold = lblock.m
    m = size(O,2)
    ra = lblock.ra[1]:cra[end]			# range of new block, absolute sites
    n = length(ra)				# size of new block
    raH1 = cra[end]+1:maximum(finalindsH[ra])
    raV = cra[end]+1:N
    lc = length(cra)
    PhiC, PhiL = O[mold+1:end,:], O[1:mold,:]
    phi = vcat(lblock.phi * PhiL,PhiC)
    Ophi = phi' * phi
    phi = phi * inv(sqrt(Symmetric(Ophi)))
    H1ij = PhiL' * lblock.H1ij * PhiL + PhiC' * H[cra,cra] * PhiC

    interabs,cr,lr,nonnull = intersectrange(cra,lblock.raH1)
    offterm = PhiC[cr,:]' * lblock.H1phi[lr,:] * PhiL
    H1ij += offterm + offterm'

    H1phi = H[raH1,cra] * PhiC
    interabs,cr,lr,nonnull = intersectrange(raH1,lblock.raH1)
    nonnull && (H1phi[cr,:] += lblock.H1phi[lr,:] * PhiL)

    Vijkl = trans4(lblock.Vijkl,PhiL)

    pp = [phi[j,k] * phi[j,l] for j=1:n, k=1:m, l=1:m]
    #@inbounds pp = [phi[j,k] * phi[j,l] for j=1:n, k=1:m, l=1:m]

    ppcc = pp[cra,:,:]
    interabs,lr,cr,nonnull = intersectrange(lblock.raV,cra)
    VppC = trans2(lblock.Vpp[lr,:,:],PhiL)
    Vijkl += contract(ppcc,VppC,true,false)
    Vijkl += contract(VppC,ppcc,true,false)
    Vijkl += contract(ppcc,contract(V[cra,cra],ppcc),true,false)
    VC = contract(V[raV,cra],ppcc)
    interabs,lr,l2r,nonnull = intersectrange(lblock.raV,raV)
    newVpp = trans2(lblock.Vpp[lr,:,:],PhiL)
    Vpp = newVpp + VC
    LRBlock(ra,raH1,raV,m,phi,H1ij,H1phi,Vijkl,pp,Vpp)
end

function getblocksizes(N,m,blocksize; verbose = false)
    nblocks = div(N-2*m,blocksize)+2
    while nblocks <= 4
	blocksize -= 1
        if blocksize == 0
            @show N,m,blocksize
        end
	nblocks = div(N-2*m,blocksize)+2
    end
    aveblocksize = div(N-2*m,nblocks-2)
    blocksizes = [aveblocksize for i=1:nblocks]
    blocksizes[1] = m
    blocksizes[end] = m
    missing = N-sum(blocksizes)
    #@show missing
    blocksizes[2] += missing
    Cranges = [1:m for b=1:nblocks]
    for b=2:nblocks
	Cranges[b] = Cranges[b-1][end]+1:Cranges[b-1][end]+blocksizes[b]
    end
    verbose && @show nblocks, blocksizes
    nblocks,blocksizes,Cranges
end

function getfinal(H,N)
    function fl(x,r) 
        res = findlast(x,r)
        res == nothing && (res = 1)
        res
    end
    function ff(x,r) 
        res = findfirst(x,r)
        res == nothing && (res = N)
        res
    end
    #@inbounds finalindsH = [findlast(i->abs(H[j,i])>1e-13,1:N) for j=1:N]
    @inbounds finalindsH = [fl(i->abs(H[j,i])>1e-13,1:N) for j=1:N]
    #@inbounds firstindsH = [findfirst(i->abs(H[j,i])>1e-13,1:N) for j=1:N]
    @inbounds firstindsH = [ff(i->abs(H[j,i])>1e-13,1:N) for j=1:N]
    finalindsH,firstindsH
end

function getinitialblocks(nblocks,blocksizes,Cranges,psi,H,V,firstindsH,finalindsH; verbose = false)
    N = size(H,1)
    block = [LRBlock(1:1,1:1,1:1,1,0,0,0,0,0,0) for i=1:nblocks]
    block[1] = getblock(1:blocksizes[1],psi,H,V,firstindsH,finalindsH,true)
    verbose && @show block[1].ra
    block[nblocks] = getblock(N-blocksizes[nblocks]+1:N,psi,H,V,firstindsH,finalindsH,true)

    for b = nblocks-1:-1:3
	cra = Cranges[b]
	lc = length(cra)
	rra = cra[end]+1:N
	crra = cra[1]:N
	phi,m = getphi(psi[crra,:])
	phi1 = block[b+1].phi
	O = vcat(phi[1:lc,:],phi1' * phi[lc+1:end,:])
	block[b] = addblockright(cra,O,block[b+1],N,H,V,firstindsH)
    end
    block
end

#L, C, R: the three blocks.  C is bare center sites
#|   Lra      | Cra |      Rra|  bare sites, contracted to
#         |lra| Cra |rra|

function getH1(LB::LRBlock,RB::LRBlock,H)
    N = size(H,1)
    Cra, Rra = LB.ra[end]+1:RB.ra[1]-1, RB.ra
    lc = length(Cra)
    ml,mr = LB.m, RB.m
    n = ml+lc+mr
    lra, cra, rra = 1:ml, ml+1:ml+lc, ml+lc+1:n

    Hb = zeros(n,n)	# 3x3 in blocks
    Hb[lra,lra] = LB.H1ij		# diagonal blocks
    Hb[cra,cra] = H[Cra,Cra]
    Hb[rra,rra] = RB.H1ij

    interabs,lr,cr,nonnull = intersectrange(LB.raH1,Cra)	# L-C blocks
    if nonnull
	hcr = ml .+ cr
	Hb[lra,hcr] = LB.H1phi[lr,:]'
	Hb[hcr,lra] = Hb[lra,hcr]'
    end

    interabs,lr,rr,nonnull = intersectrange(LB.raH1,Rra)        # L-R blocks
    if nonnull
	Hb[lra,rra] = LB.H1phi[lr,:]' * RB.phi[rr,:]
	Hb[rra,lra] = Hb[lra,rra]'
    end

    interabs,cr,rr,nonnull = intersectrange(Cra,RB.raH1)	# C-R blocks
    if nonnull
	hcr = ml .+ cr
	Hb[hcr,rra] = RB.H1phi[rr,:]
	Hb[rra,hcr] = Hb[hcr,rra]'
    end
    Hb
end

mutable struct VLCR	# contains the two particle interactions, for the L,C,R superblock
    Lra		# ranges
    Cra
    Rra
    VLL::Array{Float64,4}
    VLC::Array{Float64,3}
    VLR::Array{Float64,4}
    VCC::Array{Float64,2}
    VCR::Array{Float64,3}
    VRR::Array{Float64,4}
end

function getVLCR(LB,RB,V)
    Lra = LB.ra; Cra = LB.ra[end]+1:RB.ra[1]-1; Rra = RB.ra
    lc = length(Cra)
    VLL = LB.Vijkl
    VCC = V[Cra,Cra]
    VRR = RB.Vijkl

    lenR = length(Rra)
    VLC = permutedims(LB.Vpp[1:lc,:,:],[2,3,1])
    VCR = RB.Vpp[Cra,:,:]
    ml,mr = LB.m,RB.m
    VLR = reshape(reshape(LB.Vpp[Rra .- Lra[end],:,:],lenR,ml^2)' * 
		  			reshape(RB.pp,lenR,mr^2),ml,ml,mr,mr)
    VLCR(Lra,Cra,Rra,VLL,VLC,VLR,VCC,VCR,VRR)
end

# restricted HF
function dmrggetFr!(rho::Array{Float64,2},Vlcr::VLCR, F::Array{Float64,2})
    ml,mc,mr = size(Vlcr.VLL,1),size(Vlcr.VCC,1),size(Vlcr.VRR,1)
    #@show ml,mc,mr
    sc, sr = ml, ml + mc
    if size(rho) != size(F)
	error("inconsistent sizes rho and F")
    end
    # Same block terms
    for i=1:ml, j=1:ml, l=1:ml, k=1:ml		# LL
	F[i,j] += rho[k,l] * (2*Vlcr.VLL[i,j,k,l] - Vlcr.VLL[i,l,k,j])
    end
    for j=1:mc, i=1:mc				# CC
	ii,jj = sc+i,sc+j
	v = Vlcr.VCC[i,j]
	F[ii,ii] += rho[jj,jj] * 2*v
	F[ii,jj] -= rho[ii,jj] * v
    end
    for j=1:mr, l=1:mr, k=1:mr, i=1:mr		# RR
	ii,jj,kk,ll = i+sr,j+sr,k+sr,l+sr
	F[ii,jj] += rho[kk,ll] * (2*Vlcr.VRR[i,j,k,l] - Vlcr.VRR[i,l,k,j])
    end
    # Different block terms
    for j=1:ml, l=1:mr, k=1:mr, i=1:ml
	ii,jj,kk,ll = i,j,k+sr,l+sr
	v = Vlcr.VLR[i,j,k,l]
	F[ii,jj] += rho[kk,ll] * 2*v
	F[kk,ll] += rho[ii,jj] * 2*v
	F[ii,kk] -= rho[jj,ll] * v
	F[kk,ii] -= rho[ll,jj] * v
    end
    for i=1:ml, j=1:ml, k=1:mc
	kk = k+sc
	v = Vlcr.VLC[i,j,k]
	F[i,j] += rho[kk,kk] * 2*v
	F[kk,kk] += rho[i,j] * 2*v
	F[i,kk] -= rho[j,kk] * v
	F[kk,i] -= rho[kk,j] * v
    end
    for i=1:mr, j=1:mr, k=1:mc
	ii,jj,kk = i+sr,j+sr,k+sc
	v = Vlcr.VCR[k,i,j]
	F[ii,jj] += rho[kk,kk] * 2*v
	F[kk,kk] += rho[ii,jj] * 2*v
	F[ii,kk] -= rho[jj,kk] * v
	F[kk,ii] -= rho[kk,jj] * v
    end
end
# unrestricted HF
function dmrggetF!(rhoup::Array{Float64,2},rhodn::Array{Float64,2},Vlcr::VLCR, 
		   Fup::Array{Float64,2}, Fdn::Array{Float64,2})
    ml,mc,mr = size(Vlcr.VLL,1),size(Vlcr.VCC,1),size(Vlcr.VRR,1)
    sc, sr = ml, ml + mc
    if size(rhoup) != size(Fup) || size(rhodn) != size(Fup)
	error("inconsistent sizes rho and F")
    end
    # Same block terms
    for i=1:ml, j=1:ml, l=1:ml, k=1:ml		# LL
	direct = (rhoup[k,l]+rhodn[k,l]) * Vlcr.VLL[i,j,k,l]
	Fup[i,j] += direct - rhoup[k,l] * Vlcr.VLL[i,l,k,j]
	Fdn[i,j] += direct - rhodn[k,l] * Vlcr.VLL[i,l,k,j]
    end
    for j=1:mc, i=1:mc				# CC
	ii,jj = sc+i,sc+j
	v = Vlcr.VCC[i,j]
	direct = (rhoup[jj,jj]+rhodn[jj,jj]) * v
	Fup[ii,ii] += direct
	Fdn[ii,ii] += direct
	Fup[ii,jj] -= rhoup[ii,jj] * v
	Fdn[ii,jj] -= rhodn[ii,jj] * v
    end
    for j=1:mr, l=1:mr, k=1:mr, i=1:mr		# RR
	ii,jj,kk,ll = i+sr,j+sr,k+sr,l+sr
	direct = (rhoup[kk,ll] + rhodn[kk,ll]) * Vlcr.VRR[i,j,k,l]
	Fup[ii,jj] += direct - rhoup[kk,ll] * Vlcr.VRR[i,l,k,j]
	Fdn[ii,jj] += direct - rhodn[kk,ll] * Vlcr.VRR[i,l,k,j]
    end
    # Different block terms
    for j=1:ml, l=1:mr, k=1:mr, i=1:ml
	ii,jj,kk,ll = i,j,k+sr,l+sr
	v = Vlcr.VLR[i,j,k,l]
	rkl = rhoup[kk,ll] + rhodn[kk,ll]
	rij = rhoup[ii,jj] + rhodn[ii,jj]
	Fup[ii,jj] += rkl * v
	Fdn[ii,jj] += rkl * v
	Fup[kk,ll] += rij * v
	Fdn[kk,ll] += rij * v
	Fup[ii,kk] -= rhoup[jj,ll] * v
	Fdn[ii,kk] -= rhodn[jj,ll] * v
	Fup[kk,ii] -= rhoup[ll,jj] * v
	Fdn[kk,ii] -= rhodn[ll,jj] * v
    end
    for i=1:ml, j=1:ml, k=1:mc
	kk = k+sc
	v = Vlcr.VLC[i,j,k]
	rkk = rhoup[kk,kk] + rhodn[kk,kk]
	rij = rhoup[i,j] + rhodn[i,j]
	Fup[i,j] += rkk * v
	Fdn[i,j] += rkk * v
	Fup[kk,kk] += rij * v
	Fdn[kk,kk] += rij * v
	Fup[i,kk] -= rhoup[j,kk] * v
	Fdn[i,kk] -= rhodn[j,kk] * v
	Fup[kk,i] -= rhoup[kk,j] * v
	Fdn[kk,i] -= rhodn[kk,j] * v
    end
    for i=1:mr, j=1:mr, k=1:mc
	ii,jj,kk = i+sr,j+sr,k+sc
	v = Vlcr.VCR[k,i,j]
	rkk = rhoup[kk,kk] + rhodn[kk,kk]
	rij = rhoup[ii,jj] + rhodn[ii,jj]
	Fup[ii,jj] += rkk*v
	Fdn[ii,jj] += rkk*v
	Fup[kk,kk] += rij*v
	Fdn[kk,kk] += rij*v
	Fup[ii,kk] -= rhoup[jj,kk] * v
	Fdn[ii,kk] -= rhodn[jj,kk] * v
	Fup[kk,ii] -= rhoup[kk,jj] * v
	Fdn[kk,ii] -= rhodn[kk,jj] * v
    end
end

function getpsiexpanded(psiup,Lra,Lphi,Cra,Rra,Rphi)
    ml = size(Lphi,2)
    mr = size(Rphi,2)
    lra = 1:ml
    lc = length(Cra)
    cra = ml+1:ml+lc
    rra = ml+lc+1:ml+lc+mr
    psiL = Lphi * psiup[lra,:]
    psiC = psiup[cra,:]
    psiR = Rphi * psiup[rra,:]
    vcat(psiL,psiC,psiR)
end

function checkenergy(psi,lra,lphi,rra,rphi,H,V)
return		# effectively commented out
    N = size(H,1)
    ml = size(lphi,2)
    mr = size(rphi,2)
    cra = lra[end]+1:rra[1]-1
    lc = length(cra)
    psiwhole = vcat(lphi * psi[1:ml,:],psi[ml+1:ml+lc,:],rphi*psi[ml+lc+1:end,:])
    lblock = getblock(lra,psiwhole,H,V,firstindsH,finalindsH,true)
    rblock = getblock(rra,psiwhole,H,V,firstindsH,finalindsH,true)
    psinew = vcat(lblock.phi' * psiwhole[lra,:],
		  psiwhole[cra,:],rblock.phi'*psiwhole[rra,:])
    Vlcr = getVLCR(lblock,rblock,V)
    rho = psinew * psinew'
    H1B = getH1(lblock,rblock,H)
    F = copy(H1B)
    dmrggetFr!(rho,Vlcr,F)
    energy = tr(rho*(F+H1B))
    @show "check: ",energy
end

function getpsireduced(psiup,Lra,Lphi,Cra,Rra,Rphi)
    psiL = Lphi' * psiup[Lra,:]
    psiC = psiup[Cra,:]
    psiR = Rphi' * psiup[Rra,:]
    vcat(psiL,psiC,psiR)
end

function dohfdmrg(H,V,psiup0,psidn0,restricted=false,nblockcenter=1,blocksize=200,maxiter=1000,cutoff=1e-11; verbose=false)
    Nup,Ndn,N = size(psiup0,2),size(psidn0,2),size(H,1)
    finalindsH,firstindsH = getfinal(H,N)
    m = restricted ? Nup : Nup+Ndn
    psiall = restricted ? psiup0 : hcat(psiup0,psidn0)
    verbose && @show size(psiall)
    nblocks,blocksizes,Cranges = getblocksizes(N,m,blocksize;verbose=verbose)
    block = getinitialblocks(nblocks,blocksizes,Cranges,psiall,H,V,firstindsH,finalindsH)

    Vlcr = getVLCR(block[1],block[2+nblockcenter],V)
    psiup = getpsireduced(psiup0,Vlcr.Lra,block[1].phi,Vlcr.Cra,Vlcr.Rra,block[2+nblockcenter].phi)
    psidn = restricted ? psiup : getpsireduced(psidn0,Vlcr.Lra,block[1].phi,
				      Vlcr.Cra,Vlcr.Rra,block[2+nblockcenter].phi)
    if restricted
	ir = nblockcenter+2
	checkenergy(psiup,block[1].ra,block[1].phi,block[ir].ra,block[ir].phi,H,V)
    else
    end

    # step tuple is leftblock, direction, transformation-direction
    # transformation-direction is which direction is the next step from here
    left2right = [(b,1,1) for b=1:nblocks-1-nblockcenter]
    left2right[end] = (nblocks-1-nblockcenter,1,-1)
    right2left = [(b,-1,-1) for b=nblocks-2-nblockcenter:-1:2]
    fullsweep = vcat(left2right,right2left)
    psiallup = 0
    psialldn = 0
    energyiter = 10000.0
    lambda = [1.0 for i=1:nblocks]
    for iter=1:maxiter
	energy = 0.0
	verbose && println()
	verbose && @show iter
	for (b,dir,transdir) in fullsweep
	    H1B = getH1(block[b],block[b+1+nblockcenter],H)
	    Vlcr = getVLCR(block[b],block[b+1+nblockcenter],V)
	    rbl = block[b+1+nblockcenter]
	    psiallup = getpsiexpanded(psiup,block[b].ra,block[b].phi,block[b].ra[end]+1:rbl.ra[1]-1, 
										    rbl.ra,rbl.phi)
	    !restricted && (psialldn = getpsiexpanded(psidn,block[b].ra,block[b].phi,block[b].ra[end]+1:rbl.ra[1]-1, 
                                               rbl.ra,rbl.phi))
	    energy = energylast = 1e10
	    rhoup = psiup[:,1:Nup] * psiup[:,1:Nup]'
	    !restricted && ( rhodn = psidn[:,1:Ndn] * psidn[:,1:Ndn]')
	    maxiter = 4
	    for s=1:maxiter
		if restricted
		    Fup = copy(H1B)
		    dmrggetFr!(rhoup,Vlcr, Fup)
                    #=
		    for ii=1:size(Fup,1), jj=1:size(Fup,2)
			if isnan(Fup[ii,jj])
			    @show ii,jj,Fup[ii,jj],"nan"
			end
			if isinf(Fup[ii,jj])
			    @show ii,jj,Fup[ii,jj],"inf"
			end
		    end
                    =#
		    evals,evecs = eigsym(Fup)
		    psiup = evecs[:,1:Nup]
		    rhoup = (1-lambda[b]) * rhoup + lambda[b] * psiup * psiup'
		    Fup = copy(H1B)
		    dmrggetFr!(rhoup,Vlcr,Fup)
		    energy = tr(rhoup*(Fup+H1B))
		    #checkenergy(psiup,block[b].ra,block[b].phi,block[b+2].ra,block[b+2].phi)
		else
		    Fup,Fdn = copy(H1B),copy(H1B)
		    dmrggetF!(rhoup,rhodn,Vlcr, Fup,Fdn)
		    for ii=1:size(Fup,1), jj=1:size(Fup,2)
			if isnan(Fup[ii,jj])
			    @show ii,jj,Fup[ii,jj],"nan"
			end
			if isinf(Fup[ii,jj])
			    @show ii,jj,Fup[ii,jj],"inf"
			end
		    end
		    evals,evecs = eigsym(Fup)
		    psiup = evecs[:,1:Nup]
		    rhoup = (1-lambda[b]) * rhoup + lambda[b] * psiup * psiup'
		    evals,evecs = eigsym(Fdn)
		    psidn = evecs[:,1:Ndn]
		    rhodn = (1-lambda[b]) * rhodn + lambda[b] * psidn * psidn'
		    Fup,Fdn = copy(H1B),copy(H1B)
		    dmrggetF!(rhoup,rhodn,Vlcr, Fup,Fdn)
		    energy = 0.5 * tr(rhoup*(Fup+H1B))  + 0.5 * tr(rhodn*(Fdn+H1B))
		    enH1 = tr((rhoup+rhodn)*H1B)
		    if rand() < 0.01
			verbose && @show enH1,energy-enH1,energy
		    end
		    #checkenergy(psiup,block[b].ra,block[b].phi,block[b+2].ra,block[b+2].phi)
		end
		if energy > energylast
		    lambda[b] *= 0.5
		    #@show s,lambda[b]
		end
		if abs(energylast-energy) < cutoff || s == maxiter
    #		@show b,s,energy
		    break
		end
		energylast = energy
	    end
	    #@show b,energy
	    flush(stdout)
	    psi = 0
	    if restricted
		psi = psiup
	    else
		psi = hcat(psiup,psidn)
	    end

	    if nblockcenter < 1
		Cra = Vlcr.Cra          # range in terms of bare sites
		oldlc = length(Cra)     # number of states in C block

		if transdir == 1	# From b:b+1:b+2 blocks to b+1:b+2:b+3
		    #O = Matrix(qr(psi[1:m+oldlc,:]).Q)	# contraction of b and b+1
		    oldm = block[b].m
		    O,m = getphi(psi[1:oldm+oldlc,:])	# contraction of b and b+1
		    block[b+1] = addblockleft(Cra,O,block[b],N,H,V,finalindsH)
		    lc = length(Cranges[b+2])
		    psibb1 = transrangeleft(psi,1,O')                   # contract b and b+1
		    psiR = transrangeleft(psibb1,m+1,block[b+2].phi)	# trans b+2 to bare sites
		    psi = transrangeleft(psiR,m+lc+1,block[b+3].phi')   # contract bare into b+3

		else			# From b:b+1:b+2 blocks to b-1:b:b+1
		    oldm = block[b].m
		    O,m = getphi(psi[oldm+1:end,:])	# contraction of b and b+1
		    #O = Matrix(qr(psi[m+1:end,:]).Q)		# contracts b+1 and b+2 to m states, b+1
		    block[b+1] = addblockright(Cra,O,block[b+2],N,H,V,firstindsH)
		    psiL = transrangeleft(psi,m+1,O')	                # contract b+1 and b+2
		    psiLL = transrangeleft(psiL,1:oldm,block[b].phi)	# expand b to bare
		    psi = transrangeleft(psiLL,1,block[b-1].phi')       # contract bare into b-1
		end
	    else
		d = nblockcenter
		if transdir == 1	# From b: ... :b+1+d blocks to b+1: ... :b+2+d
		    oldm = block[b].m
		    O,m = getphi(psi[1:oldm+blocksizes[b+1],:])	# contraction of b and b+1
		    #O = Matrix(qr(psi[1:m+blocksizes[b+1],:]).Q)        # contraction of b and b+1
		    block[b+1] = addblockleft(Cranges[b+1],O,block[b],N,H,V,finalindsH)
		    Cra = Cranges[b+2][1]:Cranges[b+1+d][end]
		    lc = length(Cra)
		    psibb1 = transrangeleft(psiup,1,O')                   # contract b and b+1
		    posright = (d == 1 ? m+1 : m+blocksizes[b+2]+1)
		    psiR = transrangeleft(psibb1,posright,block[b+1+d].phi)  # trans b+1+d to bare sites
		    psiup = transrangeleft(psiR,m+lc+1,block[b+2+d].phi')   # contract bare into b+1+d
		    if !restricted
			psibb1 = transrangeleft(psidn,1,O')                   # contract b and b+1
			psiR = transrangeleft(psibb1,posright,block[b+1+d].phi)  # trans b+1+d to bare sites
			psidn = transrangeleft(psiR,m+lc+1,block[b+2+d].phi')   # contract bare into b+1+d
		    end
		else			# From b:...:b+1+d blocks to b-1: ... :b+d
		    startO = size(psi,1) - block[b+1+d].m - blocksizes[b+d] + 1
		    oldm = block[b].m
		    #O = Matrix(qr(psi[startO:end,:]).Q)                # contracts b+d and b+1+d to m states, b+1
		    O,m = getphi(psi[startO:end,:])		   	 #  contracts b+d and b+1+d to m states, b+1
		    block[b+d] = addblockright(Cranges[b+d],O,block[b+1+d],N,H,V,firstindsH)
		    psiL = transrangeleft(psiup,startO,O')               # contract b+d and b+1+d
		    psiLL = transrangeleft(psiL,1:oldm,block[b].phi)      # expand b to bare
		    psiup = transrangeleft(psiLL,1,block[b-1].phi')      # contract bare into b-1
		    if !restricted
			psiL = transrangeleft(psidn,startO,O')               # contract b+d and b+1+d
			psiLL = transrangeleft(psiL,1:oldm,block[b].phi)      # expand b to bare
			psidn = transrangeleft(psiLL,1,block[b-1].phi')      # contract bare into b-1
		    end
		end
	    end
	end
	verbose && @show iter,energy,energyiter
	flush(stdout)
	abs(energyiter-energy) < cutoff && break
	energyiter = energy
    end
    psiallup,psialldn,energyiter
end

export dohfdmrg
end

#=
@show size(psiallup)
@show size(psialldn)
serialize("hfpsi.ser",(psiallup,psialldn))
Hm = Matrix(H)
Fup = deepcopy(Hm)
Fdn = deepcopy(Hm)
rhoup = psiallup * psiallup'
rhodn = psialldn * psialldn'
Vm = Matrix(V)

getF!(rhoup,rhodn,Vm,Fup,Fdn)

serialize("HF.ser",(psiallup,psialldn,Fup,Fdn))
fo = open(basedir*"/hfpsi.dat","w")
psiup=psiallup
psidn=psialldn
println(fo,"$N $Nup $Ndn")
for i=1:N
    print(fo,i,"  ")
    global psiup,psidn
    for j=1:Nup
	print(fo,psiup[i,j],"  ")
    end
    for j=1:Ndn
	print(fo,psidn[i,j],"  ")
    end
    println(fo)
end
close(fo)

=#
