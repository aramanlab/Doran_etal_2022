using DrWatson
@quickactivate "SPIdemo"
using Distributed

@everywhere using DrWatson
@everywhere @quickactivate "SPIdemo"
@everywhere using SPI
@everywhere using StatsBase: sample

using ArgMacros 
using Gotree_jll
using Muon

@structarguments false Args begin
    @argumentrequired String inputfile "-i" "--inputfile"
    @argtest inputfile isfile "Couldn't find the input file."
    @argtest inputfile (f)->(split(f, ".")[end] ∈ ["h5ad"]) "extension not recognized; should be h5ad file"
    @argumentrequired String outputdir "-o" "--outputdir"
    @argumentdefault Int 100 nboot "-b" "--nboot"
    @argtest nboot n->0≤n "nboot must be positive"
end

function julia_main()::Cint
    # parse arguments
    args = Args()    
    @info "Starting SPI inference"
    @info "Setting up workspace"
    # setup output dir
    mkpath(args.outputdir) != nothing || 
        ErrorException("Could not create outputdir: $(args.outputdir)")
    name = first(split(basename(args.inputfile), "."))

    @info "Running SPI" 
    adata = readh5ad(args.inputfile)
    M = adata.X[:, :]
    @time spitree = distributed_calc_spi_tree(M, adata.obs_names.vals)
    
    @info "Writing out SPI Tree"
    open(joinpath(args.outputdir, name * "-tree.nw"), "w") do io
        println(io, spitree)
    end

    # Bootstrap (takes ~5hrs with 100 bootstraps on 5 processes in the background of my MacBook Pro)
    if args.nboot > 0
        @info "Starting Bootstrap with $(args.nboot) across $(nworkers()) processes"
        @time begin
        boottrees = @distributed (vcat) for i in 1:args.nboot
            nfeats = size(M,2)
            cols = sample(1:nfeats, nfeats, replace=true)
            distributed_calc_spi_tree(M[:, cols], adata.obs_names.vals)
        end
        end # time
        @info "Writing out Bootstrap trees"
        ## write out SPI boot trees
        open(joinpath(args.outputdir, name * "-boottrees.nw"), "w") do io
            for btree in boottrees
                println(io, btree)
            end
        end

        @info "using Booster to compute support values"
        ## calculate support
        run(pipeline(`$(gotree()) compute support tbe --silent \
            -i $(joinpath(args.outputdir, name * "-tree.nw")) \
            -b $(joinpath(args.outputdir, name * "-boottrees.nw")) \
            -o $(joinpath(args.outputdir, name * "-supporttree.nw"))`,
            stderr=joinpath(args.outputdir, "booster.log")))
    end 
    @info "Finishing run"
    return 0
end

@everywhere function distributed_calc_spi_tree(M, ids)
    nfeats = size(M,2)
    vals, vecs = eigen(Matrix(M*M'))
    dij = calc_spi_mtx(vecs, sqrt.(max.(vals, zero(eltype(vals)))))
    dij = dij ./ nfeats
    hc = hclust(dij, linkage=:average, branchorder=:optimal)
    nwstr(hc, ids; labelinternalnodes=false)
end

julia_main()