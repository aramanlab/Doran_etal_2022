using DrWatson
@quickactivate "Doran_etal_2022"
using ArgMacros
using TimerOutputs
using SPI
using StatsBase: sample
using Gotree_jll
using Logging
include(joinpath(srcdir(), "parsephylip.jl"))

@structarguments false Args begin
    @argumentrequired String inputfile "-i" "--inputfile"
    @argtest inputfile isfile "Couldn't find the input file."
    @argtest inputfile (f)->(split(f, ".")[end] ∈ ["txt", "phy", "phylip"]) "extension not recognized; should be .phy .phylip or .txt"
    @argumentrequired String outputdir "-o" "--outputdir"
    @argumentdefault Int 100 nboot "-b" "--nboot"
    @argtest nboot n->0≤n "nboot must be positive"
    @argumentdefault Int 0 loglevel "-l" "--loglevel"
end

function julia_main()::Cint
    # parse arguments
    args = Args()
    logger = ConsoleLogger(stdout, LogLevel(args.loglevel))
    global_logger(logger)
    time = TimerOutput()
    @timeit time "total" begin

    @info "Starting SPI inference"
    @info "Setting up workspace"
    # setup output dir
    mkpath(args.outputdir) != nothing || 
        ErrorException("Could not create outputdir: $(args.outputdir)")
    
    @info "Running SPI" 
    @timeit time "running SPI" begin
        phydf = readphylip(args.inputfile)
        M = onehotencode(phydf.seqs)
        spitree = SPI.calc_spi_tree(Matrix(M), phydf.ids; labelinternalnodes=false)
    end #time
    
    @info "Writing out SPI Tree"
    open(joinpath(args.outputdir, "SPI-tree.nw"), "w") do io
        write(io, spitree * "\n")
    end

    # Bootstrap
    if args.nboot > 0
        @info "Starting Bootstrap with $(args.nboot)"
        @timeit time "running bootstrap SPI" begin
            boottrees = Vector()
            chardf = _stringcolumntochardf(phydf.seqs)
            Threads.@threads for i in 1:args.nboot
                nchars = size(chardf, 2)
                colsmps = sample(1:nchars, nchars, replace=true)
                tmpM = DataFrame(Dict(string(n)=>chardf[:,c] for (n,c) in enumerate(colsmps)))
                tmpM = onehotencode(tmpM)
                push!(boottrees, SPI.calc_spi_tree(Matrix(tmpM), phydf.ids; labelinternalnodes=false))
            end
        end # time
        @info "Writing out Bootstrap trees"
        ## write out SPI boot trees
        open(joinpath(args.outputdir, "SPI-boottrees.nw"), "w") do io
            for btree in boottrees
                write(io, btree * "\n")
            end
        end

        @info "using Booster to compute support values"
        ## calculate support
        run(pipeline(`$(gotree()) compute support tbe --silent \
            -i $(joinpath(args.outputdir, "SPI-tree.nw")) \
            -b $(joinpath(args.outputdir, "SPI-boottrees.nw")) \
            -o $(joinpath(args.outputdir, "SPI-supporttree.nw"))`,
            stderr=joinpath(args.outputdir, "booster.log")))
    end 
    end # function timeit
    @info "Finishing run"
    @info "\ntiming" show(time) println("")
    return 0
end

julia_main()