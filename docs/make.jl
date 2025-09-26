using Literate
using Documenter, DocumenterVitepress
using SphericalSpatialTrees

rootdir = dirname(@__DIR__)
Literate.markdown(
    joinpath(rootdir, "scripts", "convert_tiled.jl"), 
    joinpath(@__DIR__, "src");
    documenter = true
)

DocMeta.setdocmeta!(SphericalSpatialTrees, :DocTestSetup, :(using SphericalSpatialTrees); recursive=true)
# Make the docs
makedocs(;
    modules=[SphericalSpatialTrees],
    authors="Fabian Gans <fgans@bgc-jena.mpg.de> and contributors",
    sitename="SphericalSpatialTrees.jl",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/meggart/SphericalSpatialTrees.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Introduction" => "convert_tiled.md",
        # "Methods" => "methods.md",
        # "API" => "api.md",
        # "Source code" => literate_pages
    ],
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo="github.com/meggart/SphericalSpatialTrees.jl",
    devbranch="main",
)