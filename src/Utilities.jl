using Gridap.Geometry, Gridap.ReferenceFEs, Gridap.TensorValues
using GridapGmsh
using Gmsh: gmsh

# Meshes
function annulus_3d(;
  Lf=12.5,
  Bf=15.0,
  H=10.0,
  Rout=10.0,
  Rin=7.0,
  h=1.0,
  msh_bottom=1.5,
  msh_top=0.25
)
  gmsh.initialize()

  gmsh.option.setNumber("Mesh.SaveAll", 1)
  gmsh.model.add("my_model")

  # Build
  gmsh.model.occ.addBox(-Lf, -Bf, -H, 2Lf, 2Bf, H, 1)
  gmsh.model.occ.addCircle(0, 0, -h, Rout, 13)
  gmsh.model.occ.addCircle(0, 0, -h, Rin, 14)
  gmsh.model.occ.addCurveLoop([14], 7)
  gmsh.model.occ.addCurveLoop([13], 8)
  gmsh.model.occ.addPlaneSurface([8, 7], 7)

  gmsh.model.occ.extrude((2,7),0,0,h)
  gmsh.model.occ.extrude((2,7),0,0,-H+h)
  gmsh.model.occ.fragment([(3, 1)], [(3, 2), (3, 3)])

  # Mesh sizes and tags
  gmsh.model.occ.synchronize()
  gmsh.model.mesh.setSize([(0, 18), (0, 22), (0, 19), (0, 16), (0, 13), (0, 14)], msh_bottom)
  gmsh.model.mesh.setSize([(0, 17), (0, 21), (0, 20), (0, 15), (0, 12), (0, 11), (0, 10), (0, 9)], msh_top)
  gmsh.model.addPhysicalGroup(2, [14], 29)
  gmsh.model.setPhysicalName(2, 29, "Γ_in")
  gmsh.model.addPhysicalGroup(2, [19], 30)
  gmsh.model.setPhysicalName(2, 30, "Γ_out")
  gmsh.model.addPhysicalGroup(2, [16, 10, 20], 31)
  gmsh.model.setPhysicalName(2, 31, "Γ_fs")
  gmsh.model.addPhysicalGroup(2, [7], 32)
  gmsh.model.setPhysicalName(2, 32, "Γ")
  gmsh.model.addPhysicalGroup(1, [13, 14], 33)
  gmsh.model.setPhysicalName(1, 33, "Γ_end")
  gmsh.model.addPhysicalGroup(3, [4, 5], 34)
  gmsh.model.setPhysicalName(3, 34, "interior_OUT")
  gmsh.model.addPhysicalGroup(2, [8, 9, 11, 12], 37)
  gmsh.model.setPhysicalName(2, 37, "interior_OUT")
  gmsh.model.addPhysicalGroup(3, [2], 35)
  gmsh.model.setPhysicalName(3, 35, "interior_T")
  gmsh.model.addPhysicalGroup(3, [3], 36)
  gmsh.model.setPhysicalName(3, 36, "interior_B")

  # for visualisation
  # gmsh.model.occ.synchronize()
  # gmsh.fltk.run()
  # gmsh.finalize()
  # return

  # Generate mesh
  gmsh.model.mesh.generate(3)
  model = GmshDiscreteModel(gmsh)
  gmsh.finalize()
  return model
end

function rectangle_3d(;
  Lf=12.5,
  Bf=15.0,
  H=10.0,
  L=10.0,
  B=5.0,
  y_offset=0.0,
  h=1.0,
  msh_bottom=1.5,
  msh_top=0.25,
  msh_plate=msh_top
)
  gmsh.initialize()

  gmsh.option.setNumber("Mesh.SaveAll", 1)
  gmsh.model.add("plate_model")

  # Build
  gmsh.model.occ.addBox(-Lf, -Bf, -H, 2Lf, 2Bf, H, 1)
  gmsh.model.occ.addRectangle(-L, -B + y_offset, -h, 2*L, 2*B, 7)
  gmsh.model.occ.extrude((2,7),0,0,h)
  gmsh.model.occ.extrude((2,7),0,0,-H+h)
  gmsh.model.occ.fragment([(3, 1)], [(3, 2), (3, 3)])

  # Mesh sizes and tags
  gmsh.model.occ.synchronize()
  gmsh.model.mesh.setSize([(0,24),(0,20),(0,17),(0,22),(0,25),(0,18),(0,19),(0,28)], msh_bottom)
  gmsh.model.mesh.setSize([(0, 9), (0, 10), (0, 11), (0, 12)], msh_plate)
  gmsh.model.mesh.setSize([(0,23),(0,16),(0,13),(0,21),(0,26),(0,14),(0,15),(0,27)], msh_top)

  gmsh.model.addPhysicalGroup(2, [18], 29)
  gmsh.model.setPhysicalName(2, 29, "Γ_in")
  gmsh.model.addPhysicalGroup(2, [23], 30)
  gmsh.model.setPhysicalName(2, 30, "Γ_out")
  gmsh.model.addPhysicalGroup(2, [20,12], 31)
  gmsh.model.setPhysicalName(2, 31, "Γ_fs")
  gmsh.model.addPhysicalGroup(2, [7], 32)
  gmsh.model.setPhysicalName(2, 32, "Γ")
  gmsh.model.addPhysicalGroup(1, [16, 15, 14, 13], 33)
  gmsh.model.setPhysicalName(1, 33, "Γ_end")
  gmsh.model.addPhysicalGroup(3, [4], 34)
  gmsh.model.setPhysicalName(3, 34, "interior_OUT")
  gmsh.model.addPhysicalGroup(2, [16,13,14,15,8,9,10,11], 37)
  gmsh.model.setPhysicalName(2, 37, "interior_OUT")
  gmsh.model.addPhysicalGroup(3, [2], 35)
  gmsh.model.setPhysicalName(3, 35, "interior_T")
  gmsh.model.addPhysicalGroup(3, [3], 36)
  gmsh.model.setPhysicalName(3, 36, "interior_B")

  # for visualisation
  # gmsh.model.occ.synchronize()
  # gmsh.fltk.run()
  # gmsh.finalize()
  # return

  # Generate mesh
  gmsh.model.mesh.generate(3)
  model = GmshDiscreteModel(gmsh)
  gmsh.finalize()
  return model
end

function mesh_2d(;
  Lf=15.0,
  H=10.0,
  h=0.6,
  p1=(-10,-h),
  p2=(10,-h),
  msh_bottom=0.25,
  msh_surface=msh_bottom,
  msh_plate=msh_bottom
)
  gmsh.initialize()

  gmsh.option.setNumber("Mesh.SaveAll", 1)
  gmsh.model.add("segment_model")

  # Build
  gmsh.model.occ.addRectangle(-Lf, -H, 0.0, 2Lf, H, 1)
  gmsh.model.occ.addPoint(p1[1], p1[2], 0.0, msh_plate, 5)
  gmsh.model.occ.addPoint(p2[1], p2[2], 0.0, msh_plate, 6)
  gmsh.model.occ.addLine(5, 6, 5)
  gmsh.model.occ.addPoint(p1[1], -H, 0.0, msh_bottom, 7)
  gmsh.model.occ.addPoint(p1[1], 0, 0.0, msh_surface, 8)
  gmsh.model.occ.addLine(7, 8, 6)
  gmsh.model.occ.addPoint(p2[1], -H, 0.0, msh_bottom, 9)
  gmsh.model.occ.addPoint(p2[1], 0, 0.0, msh_surface, 10)
  gmsh.model.occ.addLine(9, 10, 7)
  gmsh.model.occ.fragment([(2, 1)], [(1, 5),(1, 6),(1, 7)])

  # Mesh sizes and tags
  gmsh.model.occ.synchronize()
  gmsh.model.mesh.setSize([(0, 2), (0, 10)], msh_bottom)
  gmsh.model.mesh.setSize([(0, 9), (0, 1)], msh_surface)

  gmsh.model.addPhysicalGroup(1, [6], 10)
  gmsh.model.setPhysicalName(1, 10, "Γ_in")
  gmsh.model.addPhysicalGroup(1, [15], 11)
  gmsh.model.setPhysicalName(1, 11, "Γ_out")
  gmsh.model.addPhysicalGroup(1, [7,17,14], 12)
  gmsh.model.setPhysicalName(1, 12, "Γ_fs")
  gmsh.model.addPhysicalGroup(1, [5], 13)
  gmsh.model.setPhysicalName(1, 13, "Γ")
  gmsh.model.addPhysicalGroup(0, [4, 7], 14)
  gmsh.model.setPhysicalName(0, 14, "Γ_end")
  gmsh.model.addPhysicalGroup(2, [1,3], 15)
  gmsh.model.setPhysicalName(2, 15, "interior_OUT")
  gmsh.model.addPhysicalGroup(1, [8,9,13,11], 18)
  gmsh.model.setPhysicalName(1, 18, "interior_OUT")
  gmsh.model.addPhysicalGroup(2, [4], 16)
  gmsh.model.setPhysicalName(2, 16, "interior_T")
  gmsh.model.addPhysicalGroup(2, [2], 17)
  gmsh.model.setPhysicalName(2, 17, "interior_B")

  # for visualisation
  # gmsh.model.occ.synchronize()
  # gmsh.fltk.run()
  # gmsh.finalize()
  # return

  # Generate mesh
  gmsh.model.mesh.generate(2)
  model = GmshDiscreteModel(gmsh)
  gmsh.finalize()
  return model
end

# Adapted from MonolithicFEMVLFS.jl
function build_kl_tensor(
  ambient_dim::Int,
  manifold_dim::Int,
  E::T,
  ν::Float64,
  d::Float64
) where T
  ambient_dim in (2, 3) ||
    error("ambient_dim must be 2 or 3, got $ambient_dim")
  manifold_dim in (1, 2) ||
    error("manifold_dim must be 1 or 2, got $manifold_dim")
  manifold_dim <= ambient_dim ||
    error("manifold_dim ($manifold_dim) must be <= ambient_dim ($ambient_dim)")

  I = d^3 / 12
  δ(x, y) = ==(x, y)
  C_type = SymFourthOrderTensorValue{ambient_dim, T}

  max_index = 0
  for i in 1:ambient_dim, j in 1:ambient_dim, k in 1:ambient_dim, l in 1:ambient_dim
    max_index = max(max_index, data_index(C_type, i, j, k, l))
  end
  Cvals = zeros(T, max_index)

  if manifold_dim == 1
    Cvals[data_index(C_type, 1, 1, 1, 1)] = E * I
  else
    μ = E / (2 * (1 + ν))
    λ = ν * E / (1 - ν^2)
    for i in 1:2, j in 1:2, k in 1:2, l in 1:2
      Cvals[data_index(C_type, i, j, k, l)] =
        I * (
          μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) +
          λ * (δ(i, j) * δ(k, l))
        )
    end
  end

  SymFourthOrderTensorValue(Cvals...)
end

# Semi-analytic approach from Wegert et al. 2026 (DOI: 10.1016/j.oceaneng.2026.124792)
function solve_semi_analytic(output_dir,write_output;ω,D,Ib,H,h,L,bc_case,g,modes=100,nx=1000)
  semianalytic_data = solve_submerged_plate_2d(ω,D,Ib,0,0,0,H,h,L,modes,nx;
    bc_case,return_displacements=write_output,g);

  !write_output && return semianalytic_data.R, semianalytic_data.T
  # Reflection and transmission
  open("$output_dir/energy_bal_semianalytic.txt","w") do io
    println(io,"abs(R) = $(abs(semianalytic_data.R))")
    println(io,"abs(T) = $(abs(semianalytic_data.T))")
    println(io,"1-abs(R)^2-abs(T)^2 = $(1-abs(semianalytic_data.R)^2-abs(semianalytic_data.T)^2)")
  end

  w_sa = (semianalytic_data.w);
  X = semianalytic_data.x'

  # For comparison in Paraview
  model_semianalytic = CartesianDiscreteModel((X[1],X[end]),(length(X)-1,))
  @assert maximum(abs,norm.(model_semianalytic.grid.node_coords - X))/maximum(abs,norm.(X)) < 1e-14
  V_sa = TestFESpace(model_semianalytic,ReferenceFE(lagrangian,Float64,1),vector_type=Vector{ComplexF64})
  wh_sa = FEFunction(V_sa,vec(w_sa))
  Γ_sa = Triangulation(model_semianalytic)
  writevtk(Γ_sa,"$output_dir/semianalytic_solution_w",cellfields=["re(w)"=>real(wh_sa),"im(w)"=>imag(wh_sa)])

  η_sa = semianalytic_data.η_s
  X = semianalytic_data.XF

  model_semianalytic = CartesianDiscreteModel((X[1],X[end]),(length(X)-1,))
  @assert maximum(abs,norm.(model_semianalytic.grid.node_coords - X))/maximum(abs,norm.(X)) < 1e-14
  V_sa = TestFESpace(model_semianalytic,ReferenceFE(lagrangian,Float64,1),vector_type=Vector{ComplexF64})
  ηh_sa = FEFunction(V_sa,vec(η_sa))
  Γ_sa = Triangulation(model_semianalytic)

  writevtk(Γ_sa,"$output_dir/semianalytic_solution_η",cellfields=["re(η)"=>real(ηh_sa),"im(η)"=>imag(ηh_sa)])
  return semianalytic_data.R, semianalytic_data.T
end