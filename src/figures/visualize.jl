### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6180733f-4499-47d3-bcc0-a40643c2da43
begin
	# using Images
	using PlutoUI
	# using Plots
	using Printf
	using HypertextLiteral
end

# ╔═╡ 28e2fcef-5d82-4887-a764-1b8bd5399de2
md"""# Julia & Pluto Setup"""

# ╔═╡ 9f55fecb-55d4-4521-a4ac-fc97c409bf00
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 1100px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 20669689-536b-4788-8247-5241668bae40
md"""
# Slider plot
"""

# ╔═╡ f8e27694-2ebc-4973-84ec-703643ac204e
begin
	models = ["DRPM", "sPPM","Gaussian PPMx","Curve PPMx"]
	plots = ["Salso", "Graph", "ARI", "Graph and Boxplot"]
	root_folder = "./"
	ari_format = "png"
end;

# ╔═╡ fee5c53c-c77d-4cdf-b1f9-9872356d4b29
md"""Select the model, the plot, and the week."""

# ╔═╡ dae34474-9b71-49c3-bc9c-79f2d6d72e7a
@bind model_selected Select(models)

# ╔═╡ ac1e3050-7fe2-438c-9036-be90b1842425
@bind plot_selected Select(plots)

# ╔═╡ 194de13e-d3a0-48ae-8b94-b824bbd0817c
begin
	dirs = Dict()
	dirs[models[1]] = "./$(root_folder)/$(models[1])/$(plot_selected)/"
	dirs[models[2]] = "./$(root_folder)/$(models[1])/$(plot_selected)/"
	dirs[models[3]] = "./$(root_folder)/$(models[1])/$(plot_selected)/"
	dirs[models[4]] = "./$(root_folder)/$(models[1])/$(plot_selected)/"
	
	imgs = Dict()
	imgs[models[1]] = readdir(dirs[models[1]]) 
	imgs[models[2]] = readdir(dirs[models[1]]) 
	imgs[models[3]] = readdir(dirs[models[1]]) 
	imgs[models[4]] = readdir(dirs[models[1]]) 
end

# ╔═╡ 0094f0fc-0168-4b27-b1e6-6f11da6211f9
# Creare uno slider e visualizzare l'immagine corrispondente
# @bind week Slider(1:53,default=45,show_value=true)

# wider selection bar, if we will want it
@bind week html"""
<input type=range min=1 max=53 step=1 value=1 style='width: 40%;' oninput='this.nextElementSibling.value=this.value;'>
<output>1</output>
"""

# ╔═╡ 958f1005-26a1-48c1-a0e5-cf1596c54a3c
# when we make public our github this page will run automatically (without the need of running julia) as the images will become Resources, no more LocalResources
if plot_selected == "ARI"
	PlutoUI.LocalResource("$(root_folder)/$(model_selected)/$(plot_selected)/ari.$(ari_format)")
else
	num = @sprintf "%02d" week
	# Images.load(folder ecc)
	PlutoUI.LocalResource("$(root_folder)/$(model_selected)/$(plot_selected)/$(plot_selected)-$num.png")
end

# ╔═╡ 0120f572-84f7-4854-ba92-a8afaa5f8b32
# when moving to images in github repository check this
# https://discourse.julialang.org/t/how-to-render-figures-in-pluto-markdown/67900/7

# ╔═╡ 473d5755-eaa6-4918-bfbd-5690213d0da7
md"""
# Grid 4x4
"""

# ╔═╡ cd2a5f4e-61c8-477e-8020-8b5d6394736d
# begin
# 	drpm_img = load(joinpath(dirs[models[1]],imgs[models[1]][week_grid]))
# 	sppm_img = load(joinpath(dirs[models[2]],imgs[models[2]][week_grid]))
# 	gpmx_img = load(joinpath(dirs[models[3]],imgs[models[3]][week_grid]))
# 	cpmx_img = load(joinpath(dirs[models[4]],imgs[models[4]][week_grid]))

# 	drpm = plot(drpm_img);
# 	sppm = plot(sppm_img);
# 	gpmx = plot(gpmx_img);
# 	cpmx = plot(cpmx_img);
# end;

# ╔═╡ d4d6f55e-1a3d-43e4-8d52-348a99a74704
@bind week_grid html"""
<input type=range min=1 max=53 step=1 value=0 style='width: 40%;' oninput='this.nextElementSibling.value=this.value;'>
<output>1</output>
"""

# ╔═╡ a991b8aa-ac7c-4e6f-aa5a-91fa8679d202
# plot(drpm,sppm,gpmx,cpmx,layout=4,
# 	axis=nothing, showaxis = false,
# 	size=(900,500),
# 	# plot_title="Clusters comparison",
# 	title=["DRPM" "sPPM" "Gaussian PPMx" "Curve PPMx"],
# 	right_margin = 0Plots.mm,
# 	left_margin = 0Plots.mm,
# 	top_margin = 0Plots.mm,
# 	bottom_margin = 0Plots.mm,
# )

# ╔═╡ 9227c396-72f5-4f66-95b7-a42458db357f
# img_test = "https://download.unsplash.com/photo-1429616588302-fec569e203ce"
img_test = PlutoUI.LocalResource("$(root_folder)/DRPM/ARI/ari.$(ari_format)");

# ╔═╡ 0245bb01-526f-4bc6-9974-3c284033bfc1
@htl """
<table>
    <tr>
		<td>$(img_test) $(models[1])</td>
		<td>$(img_test) $(models[2])</td>            
	</tr>
	<tr>
		<td>$(img_test) $(models[3])</td>
		<td>$(img_test) $(models[4])</td>
	</tr>
</table>
"""

# ╔═╡ f95e1f3b-210b-41b8-a378-0c0ecf090c39
md"""
# Gif (Evolving Plots)
"""

# ╔═╡ 963e1a96-3d09-4aaf-9473-2342be1672b6
@bind clocks PlutoUI.Clock(;start_running=false,interval=0.4,fixed=false)

# ╔═╡ b60d286b-277a-4269-9a79-d58b93aff78d
clock_week = Int(mod(clocks,53)+1)

# ╔═╡ f5204fdc-0046-4434-8874-69af85a7dd24
if plot_selected == "ARI"
	PlutoUI.LocalResource("$(root_folder)/$(model_selected)/$(plot_selected)/ari.$(ari_format)")
else
	clock_num = @sprintf "%02d" clock_week
	# Images.load(folder ecc)
	PlutoUI.LocalResource("$(root_folder)/$(model_selected)/$(plot_selected)/$(plot_selected)-$clock_num.png")
end

# ╔═╡ a7565b38-4e09-45bd-b97d-778631aa0c13
# TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
HypertextLiteral = "~0.9.5"
PlutoUI = "~0.7.54"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "6a8aa985a57dbcb2d34133a96d159601c9ff8773"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "793501dcd3fa7ce8d375a2c878dca2296232686e"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "bd7c69c7f7173097e7b5e1be07cee2b8b7447f51"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.54"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─28e2fcef-5d82-4887-a764-1b8bd5399de2
# ╠═6180733f-4499-47d3-bcc0-a40643c2da43
# ╠═9f55fecb-55d4-4521-a4ac-fc97c409bf00
# ╟─20669689-536b-4788-8247-5241668bae40
# ╠═f8e27694-2ebc-4973-84ec-703643ac204e
# ╠═194de13e-d3a0-48ae-8b94-b824bbd0817c
# ╟─fee5c53c-c77d-4cdf-b1f9-9872356d4b29
# ╟─dae34474-9b71-49c3-bc9c-79f2d6d72e7a
# ╟─ac1e3050-7fe2-438c-9036-be90b1842425
# ╟─0094f0fc-0168-4b27-b1e6-6f11da6211f9
# ╠═958f1005-26a1-48c1-a0e5-cf1596c54a3c
# ╠═0120f572-84f7-4854-ba92-a8afaa5f8b32
# ╟─473d5755-eaa6-4918-bfbd-5690213d0da7
# ╠═cd2a5f4e-61c8-477e-8020-8b5d6394736d
# ╠═d4d6f55e-1a3d-43e4-8d52-348a99a74704
# ╠═a991b8aa-ac7c-4e6f-aa5a-91fa8679d202
# ╠═0245bb01-526f-4bc6-9974-3c284033bfc1
# ╠═9227c396-72f5-4f66-95b7-a42458db357f
# ╟─f95e1f3b-210b-41b8-a378-0c0ecf090c39
# ╠═963e1a96-3d09-4aaf-9473-2342be1672b6
# ╟─b60d286b-277a-4269-9a79-d58b93aff78d
# ╠═f5204fdc-0046-4434-8874-69af85a7dd24
# ╠═a7565b38-4e09-45bd-b97d-778631aa0c13
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
