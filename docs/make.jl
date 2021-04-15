push!(LOAD_PATH,"../src/")
using Documenter, Brunnermeier_Sannikov_Replication

makedocs(modules = [Brunnermeier_Sannikov_Replication], sitename = "Brunnermeier_Sannikov_Replication.jl")

deploydocs(repo = "github.com/TancredePolge/Brunnermeier_Sannikov_Replication.jl.git", devbranch = "main")
