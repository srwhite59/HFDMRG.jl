using HFDMRG
using LinearAlgebra

const width = 18
const atoms = 2
const sites = width * atoms
const decay = 0.24
const amplitude = 0.018
const onsite = 0.08

bands = zeros(2, sites)
for site = 1:sites
    phase = 1 + (site - 1) % width
    bands[1, site] = 0.001 * (phase - 9.5)^2
    site < sites && (bands[2, site] = site % width == 0 ? -0.035 : -0.12)
end
one_body = BandedOneBody(bands)

transfer = [exp(-width * decay)]
source = reshape([amplitude * exp(decay * phase)
    for phase = 1:width], width, 1)
sink = reshape([exp(-decay * (width + phase))
    for phase = 1:width], 1, width)
triangle = [first == second ? onsite :
    amplitude * exp(-decay * (second - first))
    for second = 1:width for first = 1:second]
interaction = UnitCellInteraction(sites, 1, transfer, zeros(0, 0),
    source, sink, triangle)

result = solve_hfdmrg(one_body, interaction; spin = :rhf,
    initialization = :dimer, state_maxdim = 1, state_cutoff = 1e-12,
    energy_tolerance = 1e-7, maximum_half_sweeps = 20)

h = Matrix(Diagonal(bands[1, :]))
for site = 1:sites-1
    h[site, site + 1] = h[site + 1, site] = bands[2, site]
end
v = [i == j ? onsite : amplitude * exp(-decay * abs(i - j))
    for i = 1:sites, j = 1:sites]
producer = ProducerHamiltonian(h, v; constant = 0.0, h1_bandwidth = 1)
audit = producer_energy(result, producer)

@assert isfinite(result.represented_energy)
@assert isfinite(audit.total)
@assert audit == producer_energy(result, producer)
println("represented energy = ", result.represented_energy)
println("producer energy = ", audit.total)
