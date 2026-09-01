# Compute R0 and T0 for a 2d wave problem
function compute_R0_and_T0_2d(φh,k,H,Lf,β,dΓ_in,dΓ_out)
  ψ0((x,z),) = cosh(k*(H+z))/cosh(k*H)
  sq(x) = x^2
  T = exp(-im*k*Lf)*∑(∫(φh*ψ0)dΓ_out)/∑(∫(sq ∘ ψ0)dΓ_out)
  R = exp(-im*k*Lf)*(∑(∫(φh*ψ0)dΓ_in)/∑(∫(sq ∘ ψ0)dΓ_in) - β*exp(-im*k*Lf))
  # Scale by incoming wave coeff so that abs2(R) + abs2(T) ≈ 1
  return [R/β,T/β]
end

# Compute Rij and Tij for a 3d wave problem
function compute_Rij_and_Tij_3d(i,j,p,φh,km,kmn,H,Lf,Bf,β,dΓ_in,dΓ_out)
  Γ_in = get_triangulation(get_cell_quadrature(dΓ_in))
  Γ_out = get_triangulation(get_cell_quadrature(dΓ_out))
  ψi((x,y,z),) = cosh(km[i+1]*(H+z))/cosh(km[i+1]*H)
  χj((x,y,z),) = cos(j*π*(y-Bf)/(2*Bf))
  kij = kmn[i+1,j+1]
  ψi_out = CellField(ψi,Γ_out); χj_out = CellField(χj,Γ_out)
  ψi_in = CellField(ψi,Γ_in); χj_in = CellField(χj,Γ_in)
  sq(x) = x^2
  CiBj_in = ∑(∫((sq ∘ ψi_in)*(sq ∘ χj_in))dΓ_in)
  CiBj_out = ∑(∫((sq ∘ ψi_out)*(sq ∘ χj_out))dΓ_out)
  T = exp(-im*kij*Lf)*∑(∫(φh*ψi_out*χj_out)dΓ_out)/CiBj_out
  if i == 0 && j == p
    R = exp(-im*kij*Lf)*(∑(∫(φh*ψi_in*χj_in)dΓ_in)/CiBj_in - β*exp(-im*kij*Lf))
  else
    R = exp(-im*kij*Lf)* ∑(∫(φh*ψi_in*χj_in)dΓ_in)/CiBj_in
  end
  [R/β,T/β] # Normalise by incident wave amplitude
end

# Compute R0j and T0j for a 3d wave problem
function compute_R0j_and_T0j_3d(p,φh,km,kmn,H,Lf,Bf,β,dΓ_in,dΓ_out)
  q = count(isreal.(kmn[1,:]))
  N = size(kmn,2)
  _compute_R0j_T_0jp(j) = compute_Rij_and_Tij_3d(0,j,p,φh,km,kmn,H,Lf,Bf,β,dΓ_in,dΓ_out)
  RjTj = stack(map(_compute_R0j_T_0jp,0:N-1))
  Rj = RjTj[1,:]
  Tj = RjTj[2,:]
  return Rj, Tj, q
end