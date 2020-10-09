# QuadraticModelsXpress

A package to use Xpress to optimize linear and quadratic problems in QuadraticModel
format (see QuadraticModels.jl and QPSReader.jl)

# Usage

```julia
using QPSReader, QuadraticModels, QuadraticModelsXpress
qps = readqps("AFIRO.SIF")
qm = QuadraticModel(qps)
stats = xpress(qm)
```
