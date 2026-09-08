using Gridap, Gridap.MultiField, Gridap.ReferenceFEs, Gridap.FESpaces
using SemiAnalyticWECs
using Wegert_et_al_2026_SubFEM
import Wegert_et_al_2026_SubFEM: mesh_2d

function solve_submerged_plate_2d(;
  model = mesh_2d(),            # Mesh
  M::Int = 5,                   # Modes in z to inclulde in integral operator bc
  bc_case = :free,              # Boundary condition at end of plate, :free, :simply_supported, or :clamped
  E = 3e+12,                    # Plate Young's modulus
  nu = 0.392,                   # Plate Poisson's ratio
  d = 0.01,                     # Plate thickness
  ρ = 1780,                     # Plate density
  ρw = 1025,                    # Fluid density
  g = 9.81,                     # Gravity
  A = 1.0,                      # Wave amplitude
  T = 3.0,                      # Wave period
  semi_analytic = false,        # Run semi-analytic approach for comparison, only works for horizontal case
  semi_analytic_nx = 1000,      # Number of x points in discretisation for semi-analytic method
  write_triangulations = false, # Write triangulations for debug
  write_output = true,          # Whether or not to write results to disk
  output_dir = "./results/2d_$(string(bc_case))/",
  output_suffix = ""
)
  @assert M>0 "Number of modes must be greater than 0"
  @assert bc_case ∈ (:free,:simply_supported,:clamped) "Currently only support bc_case = :free, :simply_supported, or :clamped"
  write_output && mkpath(output_dir)

  # Plate properties
  D = d^3/12*E/(1-nu^2);
  Ib = ρ*d;

  # Fluid properties
  ω = 2π/T
  α = ω^2/g

  # Background model
  write_output && write_triangulations && writevtk(model,"$output_dir/model")
  coords = get_node_coordinates(model)
  Lf = abs(maximum(getindex.(coords,1)) - minimum(getindex.(coords,1)))/2
  H = abs(maximum(getindex.(coords,2)) - minimum(getindex.(coords,2)))

  # Incident wave
  km = conj(dispersion_free_surface(α, M, H)/im);
  @assert all(km.*tanh.(km.*H).≈α)
  k = first(km);
  β = -im*g*A/ω
  φ_in((x,z),) = β*exp(im*k*x)*cosh(k*(z+H))/cosh(k*H)

  # Triangulations
  order = 2
  Ω = Triangulation(model)
  hmin = mean(CellField(lazy_map(vol->(vol)^(1/2),get_cell_measure(Ω)),Ω))
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

  ∂Γ = BoundaryTriangulation(Γ,tags=["Γ_end"])
  d∂Γ = Measure(∂Γ,3*order)
  n_∂Γ = get_normal_vector(∂Γ)

  if write_output && write_triangulations
    writevtk(Ω_Λ,"$output_dir/Ω_Λ")
    writevtk(Γ,"$output_dir/Γ")
    writevtk(Γ_i,"$output_dir/Γ_i")
    writevtk(Λ,"$output_dir/Λ")
    writevtk(∂Γ,"$output_dir/∂Γ",cellfields=["n_∂Γ"=>n_∂Γ])
  end

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(model,reffe,vector_type=Vector{ComplexF64},conformity=:L2)
  if bc_case == :simply_supported || bc_case == :clamped
    R = TestFESpace(Γ ,reffe,vector_type=Vector{ComplexF64},dirichlet_tags=["Γ_end"])
  elseif bc_case == :free
    R = TestFESpace(Γ ,reffe,vector_type=Vector{ComplexF64})
  end

  VR = MultiFieldFESpace([V,R])

  # Weak form
  μₚ = 1/ρw;
  μ₁ = order*(order+1)
  μ₂ = order*(order+1)
  c = bc_case == :clamped ? 1 : 0;

  a((φ,w),(v,r)) = ∫(∇(φ)⋅∇(v))dΩ -
      ∫(ω^2/g*φ*v)dΓ_fs -
      ∫(mean(∇(φ))⋅jump(v*n_Ω_Λ))dΩ_Λ -
      ∫(mean(∇(v))⋅jump(φ*n_Ω_Λ))dΩ_Λ +
      ∫(μ₁/hmin*jump(φ)*jump(v))dΩ_Λ +
      ∫(im*ω*w*v.plus)dΓ⁺ - ∫(im*ω*w*v.minus)dΓ⁻ +
      ∫(-μₚ*im*ω*ρw*φ.plus*r)dΓ⁺ + ∫(μₚ*im*ω*ρw*φ.minus*r)dΓ⁻ -
      ∫(μₚ*Ib*ω^2*w*r + 0im*(w*r))dΓ +
      ∫(μₚ*D*Δ(w)*Δ(r) + 0im*(w*r))dΓ -
      c*∫(μₚ*(D*Δ(w))*(∇(r) ⋅ n_∂Γ) + 0im*(w*r))d∂Γ - # clamped
      c*∫(μₚ*(D*Δ(r))*(∇(w) ⋅ n_∂Γ) + 0im*(w*r))d∂Γ + # clamped
      c*∫(μₚ*D*μ₂/hmin*(∇(w) ⋅ n_∂Γ)*(∇(r) ⋅ n_∂Γ) + 0im*(w*r))d∂Γ - # clamped
      ∫(μₚ*mean(D*Δ(w))*jump(∇(r) ⋅ n_Λ) + 0im*mean(w*r))dΛ -
      ∫(μₚ*mean(D*Δ(r))*jump(∇(w) ⋅ n_Λ) + 0im*mean(w*r))dΛ +
      ∫(μₚ*D*μ₂/hmin*jump(∇(w) ⋅ n_Λ)*jump(∇(r) ⋅ n_Λ) + 0im*mean(w*r))dΛ
  l((v,r)) = ∫(- 2im*k*(φ_in*v))dΓ_in

  op = AffineFEOperator(a,l,VR,VR)

  # Add integral operator -K(φ,v) to op
  Wegert_et_al_2026_SubFEM.enforce_integral_op_2d!(model,op,dΓ_in,dΓ_out,km,H)

  # Solve
  φh, wh = solve(op)

  # Write results
  η = im*ω/g*φh
  η_in = im*ω/g*CellField(φ_in,Γ_fs)

  _output_suffix = isempty(output_suffix) ? output_suffix : "_"*output_suffix;
  if write_output
    writevtk(Ω,"$output_dir/phi$_output_suffix",cellfields=["re(φ)"=>real(φh),"im(φ)"=>imag(φh)])
    writevtk(Γ_fs,"$output_dir/fs$_output_suffix",cellfields=["re(η)"=>real(η),"im(η)"=>imag(η),"re(η_in)"=>real(η_in),"im(η_in)"=>imag(η_in)])
    writevtk(Γ,"$output_dir/w$_output_suffix",cellfields=["re(w)"=>real(wh),"im(w)"=>imag(wh)])
  end

  # Reflection and transmission
  R,T = Wegert_et_al_2026_SubFEM.compute_R0_and_T0_2d(φh,k,H,Lf,β,dΓ_in,dΓ_out)

  if write_output
    open("$output_dir/energy_bal$(_output_suffix).txt","w") do io
      println(io,"km = $km")
      println(io,"abs(R) = $(abs(R))")
      println(io,"abs(T) = $(abs(T))")
      println(io,"1-abs(R)^2-abs(T)^2 = $(1-abs(R)^2-abs(T)^2)")
    end
  end

  # Approach from Wegert et al. 2026 (DOI: 10.1016/j.oceaneng.2026.124792)
  # This only works for a horizontal plate
  if semi_analytic
    p0 = first(first(get_cell_coordinates(Γ)))
    h = abs(p0[2])
    _X = getindex.(unique(reduce(vcat,stack(get_cell_coordinates(Γ)))),1)
    get_cell_coordinates(Γ)
    L = (maximum(_X) - minimum(_X))/2
    R_semi,T_semi = Wegert_et_al_2026_SubFEM.solve_semi_analytic(
      output_dir,write_output;ω,D,Ib,H,h,L,bc_case,g,nx=semi_analytic_nx)
    return φh, wh, η, R, T, R_semi, T_semi
  end

  return φh, wh, η, R, T
end

# Generate results for different mesh sizes and compare to semi-analytic code
solve_submerged_plate_2d(;model=mesh_2d(;msh_bottom=1.5),
  output_dir="results/2d_mesh_compare",output_suffix="msh=1.5",semi_analytic=true)
solve_submerged_plate_2d(;model=mesh_2d(;msh_bottom=1.0),
  output_dir="results/2d_mesh_compare",output_suffix="msh=1.0")
solve_submerged_plate_2d(;model=mesh_2d(;msh_bottom=0.5),
  output_dir="results/2d_mesh_compare",output_suffix="msh=0.5")
solve_submerged_plate_2d(;model=mesh_2d(;msh_bottom=0.25),
  output_dir="results/2d_mesh_compare",output_suffix="msh=0.25")
solve_submerged_plate_2d(;model=mesh_2d(;msh_bottom=0.1),
  output_dir="results/2d_mesh_compare",output_suffix="msh=0.1")

# Generate plots of |R| and |T| as a function of the period
#  Note that this takes awhile because the semi-analytic method is run
#  at high resolution
using CairoMakie, DelimitedFiles
# Ts = collect(1:0.01:5)
# model = mesh_2d(;msh_bottom=0.1)
# R = zeros(ComplexF64,length(Ts))
# T = zeros(ComplexF64,length(Ts))
# R_semi = zeros(ComplexF64,length(Ts))
# T_semi = zeros(ComplexF64,length(Ts))
# mkpath("./results/R_and_T_data")
# for i ∈ eachindex(Ts)
#   println("Running case $(i)/$(length(Ts))")
#   _, _, _, Ri, Ti, Ri_semi, Ti_semi = solve_submerged_plate_2d(;
#     model,T=Ts[i],write_output=false,semi_analytic=true,semi_analytic_nx=1000)
#   R[i] = Ri; T[i] = Ti
#   R_semi[i] = Ri_semi; T_semi[i] = Ti_semi
# end
# writedlm("./results/R_and_T_data/R_and_T.txt",hcat(Ts,R,T,R_semi,T_semi),',')
# Outputs
data = readdlm("./results/R_and_T_data/R_and_T.txt",',',ComplexF64)
Ts = data[:,1]; R = data[:,2]; T = data[:,3]; R_semi = data[:,4]; T_semi = data[:,5];
fig = with_theme(theme_latexfonts()) do
    fig = Figure(fontsize = 24,size=(900,500))
    ax = Axis(fig[1, 1], xlabel = L"T_s", ylabel=L"|R|,~|T|",
        yminorticksvisible = true,yminorgridvisible = true,
        xminorticksvisible = true,xminorgridvisible = true)
    lines!(ax,abs.(Ts),abs.(R),label=L"|R_{\mathrm{FEM}}|",linestyle=:solid, linewidth = 4)
    lines!(ax,abs.(Ts),abs.(R_semi),label=L"|R|",linestyle=:dash, linewidth = 4)
    lines!(ax,abs.(Ts),abs.(T),label=L"|T_{\mathrm{FEM}}|",linestyle=:solid, linewidth = 4)
    lines!(ax,abs.(Ts),abs.(T_semi),label=L"|T|",linestyle=:dash, linewidth = 4)
    Legend(fig[1,2],ax,orientation=:vertical, )
    ax2 = Axis(fig[2, 1], yscale=log10, xlabel = L"T_s", ylabel=L"|\xi|",
        )
    lines!(ax2,abs.(Ts),abs.(1 .-abs2.(R)-abs2.(T)),label=L"\xi_{\mathrm{FEM}}",linestyle=:solid, linewidth = 4)
    lines!(ax2,abs.(Ts),abs.(1 .-abs2.(R_semi)-abs2.(T_semi)),label=L"\xi",linestyle=:solid, linewidth = 4)
    ylims!(ax2,5e-17,5e-8)
    Legend(fig[2,2],ax2,orientation=:vertical, )
    fig
end;
save("results/R_and_T_data/R_and_T.png",fig, px_per_unit = 4)

# Slanted plate
model = mesh_2d(;p1=(-10,-0.5),p2=(10,-7),msh_bottom=0.25)
solve_submerged_plate_2d(;model,output_dir="results/slanted_simply_supp",bc_case=:simply_supported)
solve_submerged_plate_2d(;model,output_dir="results/slanted_free",bc_case=:free)
