module Util
using LinearAlgebra
export eigsym

function eigsym(A)
    E = eigen(Symmetric(A))
    E.values, E.vectors
end
end
