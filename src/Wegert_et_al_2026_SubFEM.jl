module Wegert_et_al_2026_SubFEM

using Gridap, Gridap.MultiField, Gridap.ReferenceFEs,
    Gridap.FESpaces, Gridap.TensorValues, Gridap.Geometry, Gridap.CellData
using SemiAnalyticWECs
using SparseArrays
using GridapGmsh
using Gmsh: gmsh

using Gridap.CellData: get_cell_quadrature

include("Utilities.jl")
include("ReflectionTransmission.jl")
include("IntegralOperators.jl")

end
