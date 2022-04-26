using DrWatson
@quickactivate "Doran_etal_2022"
using ArgMacros
using SPI
using StatsBase: sample
using Gotree_jll
using Logging
include(joinpath(srcdir(), "parsephylip.jl"))

@structarguments false Args begin
    @argumentrequired String inputfile "-i" "--inputfile"
    @argtest inputfile isfile "Couldn't find the input file."
    @argtest inputfile (f)->split(f, ".")[end] ∈ ["txt", "phy", "phylip"] \
        "extension not recognized; should be .phy .phylip or .txt"
    @argumentrequired String outputdir "-o" "--outputdir"
    @argumentdefault Int 100 nboot "-b" "--nboot"
    @argtest nboot n->0≤n "nboot must be positive"
    @argumentdefault Int 0 loglevel "-l" "--loglevel"
end

function julia_main()::Cint
    # parse arguments
    args = Args()
    Logging.LogLevel(args.loglevel)
    @info "Starting SPI inference"
    @info "Setting up workspace"
    # setup output dir
    mkpath(args.outputdir) || 
        ErrorException("Could not create outputdir: $(args.outputdir)")
    
    @info "Running SPI" 
    @time begin
        phydf = readphylip(args.inputfile) ||
            ErrorException("""
            Could not parse inputfile: $(args.inputfile)
            Please check formatting should be in relaxed phylip format
            """)
        M = onehotencode(phydf.seqs)
        spitree = calc_spi_tree(M, phydf.ids)
    end #time
    
    @info "Writing out SPI Tree"
    open(joinpath(args.outputdir, "SPI-tree.nw"), "w") do io
        write(io, spitree * "\n")
    end

    # Bootstrap
    if args.nboot > 0
        @info "Starting Bootstrap with $(args.nboots)"
        @time begin
            boottrees = Vector()
            chardf = _stringcolumntochardf(phydf.seqs)
            Threads.@threads for i in 1:args.nboot
                nchars = size(chardf, 2)
                tmpM = chardf[:,sample(1:nchars, nchars, replace=true)]
                push!(boottrees, calc_spi_tree(tmpM, phydf.ids))
            end
        end # time
        @info "Writing out Bootstrap trees"
        ## write out SPI boot trees
        open(joinpath(args.outputdir, "SPI-boottrees.nw"), "w") do io
            for btree in boottrees
                write(io, btree * "\n")
            end
        end

        ## calculate support
        run(`$(gotree()) compute support booster \
            -i $(joinpath(args.outputdir, "SPI-tree.nw")) \
            -b $(joinpath(args.outputdir, "SPI-boottrees.nw")) \
            -o $(joinpath(args.outputdir, "SPI-supporttree.nw"))`)
    end 
    return 0
end

@time julia_main()