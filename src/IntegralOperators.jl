# Add the following integral operator to the matrix in an AffineFEOperator:
#   K(φ,v) = -∑ᵢ{ ikᵢ/Cᵢ²( (ψᵢ,φ)_Γout*(ψᵢ,v)_Γout + (ψᵢ,φ)_Γin*(ψᵢ,v)_Γin )}.
# We construct (ψᵢ,φ)_Γout*(ψᵢ,v)_Γout and (ψᵢ,φ)_Γin*(ψᵢ,v)_Γin using the outer product
# of two sparse vectors, i.e.,
#   (ψᵢ,φ)_Γout*(ψᵢ,v)_Γout → (b ⊗ b)φₕ where b = (ψᵢ,v)_Γout.
# For the 3D case, the same approach is used.
function enforce_integral_op_2d!(::UnstructuredDiscreteModel{2},op,dΓ_in,dΓ_out,km,H)
  K_op = get_matrix(op)
  VR = get_test(op)

  ψn(n) = ((x,z),) -> cosh(km[n]*(z+H))
  sq(x) = x^2
  Cn_in = map(i->sqrt(sum(∫(sq ∘ ψn(i))dΓ_in)),1:length(km))
  Cn_out = map(i->sqrt(sum(∫(sq ∘ ψn(i))dΓ_out)),1:length(km))
  b0 = sparse(assemble_vector(((v,r),)-> ∫((v*ψn(1))/Cn_out[1])dΓ_out,VR))
  K_op .-= im*km[1]*b0*b0'
  map(2:length(km)) do i
    b = sparse(assemble_vector(((v,r),)-> ∫((v*ψn(i))/Cn_out[i])dΓ_out,VR))
    K_op .-= im*km[i]*b*b'
  end
  b0 = sparse(assemble_vector(((v,r),)-> ∫((v*ψn(1))/Cn_in[1])dΓ_in,VR))
  K_op .-= im*km[1]*b0*b0'
  map(2:length(km)) do i
    b = sparse(assemble_vector(((v,r),)-> ∫((v*ψn(i))/Cn_in[i])dΓ_in,VR))
    K_op .-= im*km[i]*b*b'
  end

  return op
end

function enforce_integral_op_3d!(::UnstructuredDiscreteModel{3},op,dΓ_in,dΓ_out,km,kmn,H,Bf)
  K_op = get_matrix(op)
  VR = get_test(op)
  Γ_in = get_triangulation(get_cell_quadrature(dΓ_in))
  Γ_out = get_triangulation(get_cell_quadrature(dΓ_out))
  knm_inds = CartesianIndices(kmn)
  N = size(kmn,2)

  ψn(n) = ((x,y,z),) -> cosh(km[n]*(z+H)) # Eigenfunction in x
  χm(m) = ((x,y,z),) -> cos((m-1)*π*(y-Bf)/(2*Bf)) # Eigenfunction in y
  sq(x) = x^2
  ψn²_in = map(n->CellField(sq ∘ ψn(n),Γ_in), 1:length(km))
  χm²_in = map(m->CellField(sq ∘ χm(m),Γ_in), 1:N+1)
  CnBm_in = map(knm_inds) do I
    sqrt(sum(∫(ψn²_in[I[1]]*χm²_in[I[2]])dΓ_in))
  end
  ψn²_out = map(n->CellField(sq ∘ ψn(n),Γ_out), 1:length(km))
  χm²_out = map(m->CellField(sq ∘ χm(m),Γ_out), 1:N+1)
  CnBm_out = map(knm_inds) do I
    sqrt(sum(∫(ψn²_out[I[1]]*χm²_out[I[2]])dΓ_out))
  end

  b0 = sparse(assemble_vector(((v,r),)-> ∫((χm(1)*v*ψn(1))/CnBm_out[1])dΓ_out,VR))
  K_op .-= im*kmn[1]*b0*b0'
  map(knm_inds[2:end]) do I
    n = I[1]; m = I[2]
    b = sparse(assemble_vector(((v,r),)-> ∫((χm(m)*v*ψn(n))/CnBm_out[I])dΓ_out,VR))
    K_op .-= im*kmn[I]*b*b'
  end
  b0 = sparse(assemble_vector(((v,r),)-> ∫((χm(1)*v*ψn(1))/CnBm_in[1])dΓ_in,VR))
  K_op .-= im*kmn[1]*b0*b0'
  map(knm_inds[2:end]) do I
    n = I[1]; m = I[2]
    b = sparse(assemble_vector(((v,r),)-> ∫((χm(m)*v*ψn(n))/CnBm_in[I])dΓ_in,VR))
    K_op .-= im*kmn[I]*b*b'
  end

  return op
end