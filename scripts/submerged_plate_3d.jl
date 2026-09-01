using Gridap, Gridap.MultiField, Gridap.ReferenceFEs
using SemiAnalyticWECs
using IncompleteLU
using IterativeSolvers
using Wegert_et_al_2026_SubFEM
import Wegert_et_al_2026_SubFEM: annulus_3d, rectangle_3d

function solve_submerged_plate_3d(;
  model = annulus_3d(),         # Mesh
  M::Int = 5,                   # Modes in z to inclulde in integral operator bc
  N::Int = 5,                   # Modes in y to inclulde in integral operator bc
  p::Int = 0,                   # Incoming wave mode in y direction (p=0 → plane wave)
  bc_case = :free,              # Boundary condition at end of plate, currently only :free is supported
  E = 3e+12,                    # Plate Young's modulus
  nu = 0.392,                   # Plate Poisson's ratio
  d = 0.01,                     # Plate thickness
  ρ = 1780,                     # Plate density
  ρw = 1025,                    # Fluid density
  g = 9.81,                     # Gravity
  A = 1.0,                      # Wave amplitude
  T = 4.0,                      # Wave period
  τ = 0.0001,                   # Drop tolerance for iLU(τ)
  write_triangulations = false, # Write triangulations for debug
  output_dir = "./results/3d_$(string(bc_case))/",
  output_suffix = ""
)
  @assert M>0 && N>0 "Number of modes M and N must be greater than 0"
  @assert bc_case ∈ (:free,) "Currently only support bc_case = :free"
  mkpath(output_dir)

  # Plate properties
  D = Wegert_et_al_2026_SubFEM.build_kl_tensor(3,2,E,nu,d)
  Ib = ρ*d;

  # Fluid properties
  ω = 2π/T
  α = ω^2/g

  # Background model
  write_triangulations && writevtk(model,"$output_dir/model")
  coords = get_node_coordinates(model)
  Lf = abs(maximum(getindex.(coords,1)) - minimum(getindex.(coords,1)))/2
  Bf = abs(maximum(getindex.(coords,2)) - minimum(getindex.(coords,2)))/2
  H = abs(maximum(getindex.(coords,3)) - minimum(getindex.(coords,3)))

  # Incident wave
  km = conj(dispersion_free_surface(α, M, H)/im);
  @assert all(km.*tanh.(km.*H).≈α)
  kmn = stack(map(m->@.(sqrt(km^2 - (m*π/(2*Bf))^2)),0:N))
  k0 = first(km);
  k0p = sqrt(k0^2-(p*π/(2*Bf))^2)
  @assert isreal(k0p) "k₀ₚ must be real, you should decrease p"
  @assert !isreal(kmn[1,:]) "k₀ᵢ is real for all i∈{0,...,N}. You should increase N
  so that at least k_0N is imaginary, otherwise some of the propagating modes will
  be missed in the computation."
  β = -im*g*A/ω
  φ_in((x,y,z),) = β*exp(im*k0p*x)*cosh(k0*(z+H))/cosh(k0*H)*cos(p*π*(y-Bf)/(2*Bf))

  # Triangulations
  order = 2
  Ω = Triangulation(model)
  hmin = mean(CellField(lazy_map(vol->(vol)^(1/3),get_cell_measure(Ω)),Ω))
  Ω_Λ = Skeleton(model,tags=["interior_OUT","interior_T","interior_B"])
  Γ = BoundaryTriangulation(model,tags=["Γ","Γ_end"])
  Ω_T = Triangulation(model,tags=["interior_T"])
  Ω_B = Triangulation(model,tags=["interior_B"])
  Γ_i = InterfaceTriangulation(Ω_B,Ω_T) # Note: plus -> below Γ & minus -> above Γ
  Γ_fs = BoundaryTriangulation(model,tags=["Γ_fs"])
  Γ_in = BoundaryTriangulation(model,tags=["Γ_in"])
  Γ_out = BoundaryTriangulation(model,tags=["Γ_out"])
  Λ = Skeleton(Γ)
  dΩ = Measure(Ω,3*order)
  dΩ_Λ = Measure(Ω_Λ,3*order)
  dΓ = Measure(Γ,3*order)
  dΓ⁺ = Measure(Γ_i.plus,3*order)
  dΓ⁻ = Measure(Γ_i.minus,3*order)
  dΓ_fs = Measure(Γ_fs,3*order)
  dΓ_in = Measure(Γ_in,10*order)
  dΓ_out = Measure(Γ_out,10*order)
  dΛ = Measure(Λ,3*order)
  n_Ω_Λ = get_normal_vector(Ω_Λ)
  n_Λ = get_normal_vector(Λ)

  if write_triangulations
    writevtk(Ω_Λ,"$output_dir/Ω_Λ")
    writevtk(Γ,"$output_dir/Γ")
    writevtk(Γ_i,"$output_dir/Γ_i")
    writevtk(Λ,"$output_dir/Λ")
  end

  # FE spaces
  V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order),vector_type=Vector{ComplexF64},conformity=:L2)
  R = TestFESpace(Γ ,ReferenceFE(lagrangian,Float64,order),vector_type=Vector{ComplexF64})

  VR = MultiFieldFESpace([V,R])

  ## Weak form
  μₚ = 1/ρw; # Result of coercivity
  μ₁ = order*(order+1)
  μ₂ = order*(order+1)

  a((φ,w),(v,r)) = ∫(∇(φ)⋅∇(v))dΩ -
      ∫(ω^2/g*φ*v)dΓ_fs -
      ∫(mean(∇(φ))⋅jump(v*n_Ω_Λ))dΩ_Λ -
      ∫(mean(∇(v))⋅jump(φ*n_Ω_Λ))dΩ_Λ +
      ∫(μ₁/hmin*jump(φ)*jump(v))dΩ_Λ +
      ∫(im*ω*w*v.plus)dΓ⁺ - ∫(im*ω*w*v.minus)dΓ⁻ +
      ∫(-μₚ*im*ω*ρw*φ.plus*r)dΓ⁺ + ∫(μₚ*im*ω*ρw*φ.minus*r)dΓ⁻ -
      ∫(μₚ*Ib*ω^2*w*r + 0im*(w*r))dΓ +
      ∫(μₚ*((D ⊙ ∇∇(w)) ⊙ ∇∇(r)) + 0im*(w*r))dΓ -
      ∫(μₚ*mean(D ⊙ ∇∇(w)) ⊙ jump(∇(r) ⊗ n_Λ) + 0im*mean(w*r))dΛ -
      ∫(μₚ*mean(D ⊙ ∇∇(r)) ⊙ jump(∇(w) ⊗ n_Λ) + 0im*mean(w*r))dΛ +
      ∫(μₚ*maximum(D)*μ₂/hmin*jump(∇(w) ⊗ n_Λ) ⊙ jump(∇(r) ⊗ n_Λ) + 0im*mean(w*r))dΛ
  l((v,r)) = ∫(- 2im*k0p*(φ_in*v))dΓ_in

  op = AffineFEOperator(a,l,VR,VR)

  # Add integral operator -K(φ,v) to op
  Wegert_et_al_2026_SubFEM.enforce_integral_op_3d!(model,op,dΓ_in,dΓ_out,km,kmn,H,Bf)

  # Solve
  K = get_matrix(op); l = get_vector(op)
  Pl = ilu(K, τ = τ)
  x = IterativeSolvers.idrs(K,l;Pl,verbose=true,reltol=1e-13)
  φh, wh = FEFunction(VR,x)

  # Write results
  η = im*ω/g*φh
  η_in = im*ω/g*CellField(φ_in,Γ_fs)

  _output_suffix = isempty(output_suffix) ? output_suffix : "_"*output_suffix;
  writevtk(Ω,"$output_dir/phi$_output_suffix",cellfields=["re(φ)"=>real(φh),"im(φ)"=>imag(φh)])
  writevtk(Γ_fs,"$output_dir/fs$_output_suffix",cellfields=["re(η)"=>real(η),"im(η)"=>imag(η),"re(η_in)"=>real(η_in),"im(η_in)"=>imag(η_in)])
  writevtk(Γ,"$output_dir/w$_output_suffix",cellfields=["re(w)"=>real(wh),"im(w)"=>imag(wh)])

  # Reflection and transmission
  Rj, Tj, q = Wegert_et_al_2026_SubFEM.compute_R0j_and_T0j_3d(p,φh,km,kmn,H,Lf,Bf,β,dΓ_in,dΓ_out)

  if p == 0
    E = abs(Rj[1])^2+abs(Tj[1])^2 +
      1/2*sum((kmn[1,j]/kmn[1,1])*(abs(Rj[j])^2+abs(Tj[j])^2) for j ∈ 2:q)
  else
    E = (2kmn[1,1]/kmn[1,p+1])*(abs(Rj[1])^2+abs(Tj[1])^2) +
      sum((kmn[1,j]/kmn[1,p+1])*(abs(Rj[j])^2+abs(Tj[j])^2) for j ∈ 2:q)
  end

  open("$output_dir/energy_bal$(_output_suffix).txt","w") do io
    println(io,"k0j = $(kmn[1,:])")
    println(io,"Rj = $(Rj)")
    println(io,"Tj = $(Tj)")
    println(io,"1 - E = $(abs(1 - E))")
  end

  return φh, wh, η, Rj, Tj
end

solve_submerged_plate_3d(;model=rectangle_3d(),p=0,output_dir = "./results/3d_rect/")
solve_submerged_plate_3d(;model=rectangle_3d(;y_offset=3.0),p=0,output_dir = "./results/3d_rect_offset/")

solve_submerged_plate_3d(;model=annulus_3d(),p=0,output_dir = "./results/3d_annulus_p=0/")
solve_submerged_plate_3d(;model=annulus_3d(),p=1,output_dir = "./results/3d_annulus_p=1/")