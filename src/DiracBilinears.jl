module DiracBilinears

import EzXML
import FFTW
import OMEinsum as ES
import LinearAlgebra as LA
import SpecialFunctions as SF
import Printf as PF
import Dates
# import MPI


function make_σ()
    σ = zeros(ComplexF64, (2, 2, 3))
    σ[:, :, 1] = [
        (0.0) (1.0)
        (1.0) (0.0)
    ]
    σ[:, :, 2] = [
        (0.0) (-im)
        (im) (0.0)
    ]
    σ[:, :, 3] = [
        (1.0) (0.0)
        (0.0) (-1.0)
    ]
    return σ
end
const σ = make_σ()

function make_ϵijk()
    ϵ = zeros(3, 3, 3)
    ϵ[1,2,3] = 1.0
    ϵ[2,3,1] = 1.0
    ϵ[3,1,2] = 1.0
    ϵ[2,1,3] = -1.0
    ϵ[1,3,2] = -1.0
    ϵ[3,2,1] = -1.0
    return ϵ
end
const ϵijk = make_ϵijk()
δ(x, y) = ==(x, y)*1
step(x) = (x > 0.0)*1.0

const au2ev = 27.211396127707

include("read_files.jl")
include("density.jl")
include("for_wan2respack.jl")
include("wannier_matrix.jl")
include("module.jl")

end
