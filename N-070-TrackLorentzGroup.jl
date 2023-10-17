### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ d175d7aa-96b2-4823-848b-4db22a523322
begin
	using LinearAlgebra
	using StaticArrays
end

# ╔═╡ a36f3998-a483-4f39-b1d4-712e9f7ece56
md"""
# Tracking Lorentz Transformations using combined SL(2,C)&SO(3,1) representations

Ensuring a coherent framework for amplitude analysis in particle physics, particularly in modeling multibody decays, demands a scrupulous implementation of Lorentz transformations and alignment rotations. With particles possessing intrinsic spin, maintaining the uniformity of quantization axes across diverse decay channels and topologies is paramount. This consistency, pivotal for combining distinct decay trees while preserving the physical integrity of the model, introduces the necessity to adeptly navigate through Lorentz transformations and alignment rotations.

Mathematically embodying this, our work finds its foundation in the Lorentz group, focusing on the systematic application and decomposition of Lorentz transformations in Minkowski spacetime. In terms of mathematical representations, we primarily utilize the representation:

$M = R_z(\psi)R_y(\theta)B_z(\xi)R_z(\phi_\text{rf})R_y(\theta_\text{rf})R_z(\psi_\text{rf})$

where each term signifies a specific rotation or boost in the respective axis, \(Rz(\psi)\), \(Ry(\theta)\), and \(Bz(\xi)\) represent rotations around the z-axis, y-axis, and a boost along the z-direction, parameterized by angles \(\psi\), \(\theta\), and rapidity \(\xi\), respectively. The subscript "rf" refers to parameters in the rest frame, which is crucial for the transformations in the particle's rest frame.

The seem to be no simple way to decode physical parameters of the transformation from the `SL(2,C)` matrix. However, it is rather straightforward to decode the SO(3,1) matrix. The method employed here is to track the both representations. The decoding is done for SO(3,1) part, leaving 2 ambiguities related to possible 2π shift in the azimuthal angle to be checked for SL(2,C) part.

The algorithm for decomposing 4x4 matrix intricately applies the Lorentz transformation to specified vectors, followed by an adroit sequence of inverse transformations and parameter extractions to dissect the attributes of the initiating transformation. By applying the Lorentz transformation to a time-like vector - a vector with components [0, 0, 0, 1] - we ensure that the initial three rotations do not influence the state. Subsequent transformations and reductions allow for the decoding of the angles from the resultant transformation and the extraction of the remaining rotations and boost parameters.

Ensuring the physical relevance of the parameters derived from the 4x4 representations, a meticulous validation against the 2x2 representation is conducted. This involves comparing the 4x4-derived parameters against the 2x2 representation, employing adjustments to azimuthal angles across different configurations.

In culmination, this endeavor promises to navigate through the complexities of multi-body decay processes in particle physics, facilitating researchers in constructing and analyzing physically sound models across varied decay topologies and channels.
"""

# ╔═╡ d9b5a692-10fe-40b3-8ee5-c185ce96159b
md"""
## Matrices of SO(3,1)
"""

# ╔═╡ 4c8b20b3-998f-43b9-89bb-5f97d89267da
begin # Define 4x4 transformation
	Ry_4x4(θ::Float64) = 
	    @SMatrix [
			cos(θ) 0 sin(θ) 0
			0 1 0 0
			-sin(θ) 0 cos(θ) 0
			0 0 0 1
		]
	
	Rz_4x4(φ::Float64) =
		@SMatrix [
			cos(φ) -sin(φ) 0 0
			sin(φ) cos(φ) 0 0
			0 0 1 0
			0 0 0 1
		]
	
	function Bz_4x4(χ::Float64)
	    γ = cosh(χ)
	    βγ = sinh(χ)
	    return @SMatrix [
			1 0 0 0
			0 1 0 0
			0 0 γ βγ
			0 0 βγ γ]
	end
end

# ╔═╡ d170e4a7-4dbd-4b24-acab-2bcb3b4e407e
# Decode rotation parameters from a 4x4 matrix
function decode_rotation(R::AbstractArray)
	φ = atan(R[2,3], R[1,3])
	θ = acos(R[3,3])
	ψ = atan(R[3,2], -R[3,1])
	return φ, θ, ψ
end

# ╔═╡ 2a76d690-c970-4f8d-bf5a-dec96c9f42c0
function build_M_4x4(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf)
	return Rz_4x4(ψ)*Ry_4x4(θ)*Bz_4x4(ξ)*Rz_4x4(φ_rf)*Ry_4x4(θ_rf)*Rz_4x4(ψ_rf)
end

# ╔═╡ e5065104-c5ce-4b51-ade3-8e307ee40f81
md"""
## Matrices of SL(2,C)
"""

# ╔═╡ 08844a92-1c61-45fc-acdf-e9fadfbba3f7
begin
	# Pauli Matrices
	const σx = @SMatrix [0.0 1.0; 1.0 0.0]
	const σy = @SMatrix [0.0 -im; im 0.0]
	const σz = @SMatrix [1.0 0.0; 0.0 -1.0]
	
	Ry(θ) = cos(θ/2)*I - 1im*sin(θ/2)*σy
	Rz(φ) = cos(φ/2)*I - 1im*sin(φ/2)*σz
	Bz(ξ) = cosh(ξ/2)*I + sinh(ξ/2)*σz
end

# ╔═╡ 3e4d4e6b-b6aa-4de2-b3cf-8ea942bc077e
function build_M_2x2(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf)
	return Rz(ψ)*Ry(θ)*Bz(ξ)*Rz(φ_rf)*Ry(θ_rf)*Rz(ψ_rf)
end

# ╔═╡ 31a85762-e5d7-470a-963d-4af8ec202e14
md"""
## Joined representation
"""

# ╔═╡ e33ac4f8-f167-4fe9-a994-58124e3c08fd
begin
	struct LorentzTransformation
	    M_2x2::SMatrix{2,2,Complex{Float64}}
	    M_4x4::SMatrix{4,4,Float64}
	end
	LorentzTransformation(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf) =
		LorentzTransformation(
			build_M_2x2(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf),
			build_M_4x4(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf)
		)
	function Base.:*(A::LorentzTransformation, B::LorentzTransformation)
	    new_M_2x2 = A.M_2x2 * B.M_2x2
	    new_M_4x4 = A.M_4x4 * B.M_4x4
	    return LorentzTransformation(new_M_2x2, new_M_4x4)
	end
end

# ╔═╡ df796c39-bbe1-4714-8594-c98165f8a95c
function decompose_4x4_matrix(M::SMatrix{4,4,Float64})
    # Step 1: Apply M to a time-like vector V0
    m = 1.0  # Assuming mass is 1
    V0 = @SVector [0.0, 0.0, 0.0, m]
    V = M * V0

    # Step 2: Obtain ψ, θ, and ξ from V
    w = V[4]
    p = sqrt(V[1]^2 + V[2]^2 + V[3]^2)
    γ = w / m
    ξ = acosh(γ)  # hyperbolic arccosine
    
    ψ = atan(V[2], V[1])
    θ = acos(V[3] / p)

    # Step 4: Compute M in the rest frame (M_rf)
    M_rf = Bz_4x4(-ξ) * Ry_4x4(-θ) * Rz_4x4(-ψ) * M
    
    # Step 5: Decode the remaining rotation angles (φ_rf, θ_rf, ψ_rf)
    φ_rf, θ_rf, ψ_rf = decode_rotation(@view M_rf[1:3, 1:3])

	pars = ψ, θ, ξ, φ_rf, θ_rf, ψ_rf
	
	@assert build_M_4x4(pars...) ≈ M
	
	return pars
end

# ╔═╡ 8998cfa2-61bb-41a0-a0dd-47091b35fec0
let
	ψ = -1.2
	θ = 2.3
	ξ = 1.4
	φ_rf = 2.5
	θ_rf = 2.6
	ψ_rf = -2.7
	# 
	M = build_M_4x4(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf)
	@assert all((ψ, θ, ξ, φ_rf, θ_rf, ψ_rf) .≈ decompose_4x4_matrix(M))
	(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf) .=> decompose_4x4_matrix(M)
end

# ╔═╡ 370edcd6-c6a2-4813-8818-4681475b05ec
function adjust_2π(M_original_2x2, ψ, θ, ξ, φ_rf, θ_rf, ψ_rf)
	options = [
		(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf),
		(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf+2π)]

	M_2x2 = build_M_2x2(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf)

	if M_original_2x2 ≈ M_2x2
		return (ψ, θ, ξ, φ_rf, θ_rf, ψ_rf)
	elseif M_original_2x2 ≈ -M_2x2
		return (ψ, θ, ξ, φ_rf, θ_rf, ψ_rf + 2π)
	else
		error("No options found")
	end
end

# ╔═╡ acdb8b1a-2069-4f00-b169-d7ccb906ba3b
function decompose(L::LorentzTransformation)
    ψ, θ, ξ, φ_rf, θ_rf, ψ_rf = decompose_4x4_matrix(L.M_4x4)
    return adjust_2π(L.M_2x2, ψ, θ, ξ, φ_rf, θ_rf, ψ_rf)
end

# ╔═╡ 68587656-5eb4-4c68-ac92-f214f216de22
md"""
# Test: decoding matrix composition
"""

# ╔═╡ 84fa9b29-e748-4b16-a860-b9701d8c4a1c
let
	ψ = 2.2
	θ = 1.3
	ξ = 2.4
	φ_rf = 2.6
	θ_rf = 0.6
	ψ_rf = 5.7
	# 
	M = LorentzTransformation(ψ, θ, ξ, φ_rf, θ_rf, ψ_rf)
	X = M*M*M
	d_4x4 = decompose_4x4_matrix(X.M_4x4)
	d_4x4 .=> decompose(X)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
StaticArrays = "~1.6.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "b26880843961472df3a68ffad9c4b192b6d2e634"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "0adf069a2a490c47273727e029371b31d44b72b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.5"

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"
"""

# ╔═╡ Cell order:
# ╟─a36f3998-a483-4f39-b1d4-712e9f7ece56
# ╠═d175d7aa-96b2-4823-848b-4db22a523322
# ╟─d9b5a692-10fe-40b3-8ee5-c185ce96159b
# ╠═4c8b20b3-998f-43b9-89bb-5f97d89267da
# ╠═d170e4a7-4dbd-4b24-acab-2bcb3b4e407e
# ╠═df796c39-bbe1-4714-8594-c98165f8a95c
# ╠═2a76d690-c970-4f8d-bf5a-dec96c9f42c0
# ╠═8998cfa2-61bb-41a0-a0dd-47091b35fec0
# ╟─e5065104-c5ce-4b51-ade3-8e307ee40f81
# ╠═08844a92-1c61-45fc-acdf-e9fadfbba3f7
# ╠═3e4d4e6b-b6aa-4de2-b3cf-8ea942bc077e
# ╟─31a85762-e5d7-470a-963d-4af8ec202e14
# ╠═e33ac4f8-f167-4fe9-a994-58124e3c08fd
# ╠═370edcd6-c6a2-4813-8818-4681475b05ec
# ╠═acdb8b1a-2069-4f00-b169-d7ccb906ba3b
# ╟─68587656-5eb4-4c68-ac92-f214f216de22
# ╠═84fa9b29-e748-4b16-a860-b9701d8c4a1c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
