using Documenter, QuadraticModelsXpress

makedocs(
  modules = [QuadraticModelsXpress],
  doctest = true,
  # linkcheck = true,
  strict = true,
  format = Documenter.HTML(assets = ["assets/style.css"], prettyurls = get(ENV, "CI", nothing) == "true"),
  sitename = "QuadraticModelsXpress.jl",
  pages = Any["Home" => "index.md",
              "Tutorial" => "tutorial.md",
              "Reference" => "reference.md"]
)

deploydocs(repo = "github.com/JuliaSmoothOptimizers/QuadraticModelsXpress.jl.git")
