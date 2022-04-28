using DrWatson
@quickactivate "Doran_etal_2022"
using FastME_jll, Gotree_jll
using ArgMacros, TimerOutputs
using Logging

@structarguments false Args begin
    @argumentrequired String inputfile "-i" "--inputfile"
    @argtest inputfile isfile "Couldn't find the input file."
    @argtest inputfile (f)->(split(f, ".")[end] ∈ ["txt", "phy", "phylip"]) "input extension should be .phy .phylip or .txt"
    @argumentrequired String outputdir "-o" "--outputdir"
    @argumentrequired Symbol model "-m" "--model"
    @argtest model (m)->(m ∈ [:JC69, :WAG]) "model must be JC69 or WAG"
    @argumentdefault Int 0 loglevel "-l" "--loglevel"
    @argumentdefault Int 100 nboot "-b" "--nboot"
    @argtest nboot n->0≤n "nboot must be positive"
end

function julia_main()::Cint
    # parse arguments
    args = Args()
    logger = ConsoleLogger(stdout, LogLevel(args.loglevel))
    global_logger(logger)
    time = TimerOutput()
    @timeit time "total" begin
    
    mkpath(args.outputdir) != nothing || 
        ErrorException("Could not create outputdir: $(args.outputdir)")
    name = first(split(basename(args.inputfile), "."))
    modelparam = if args.model == :JC69
        "-dJC69"
    elseif args.model == :WAG
        "-pWAG"
    end
    
    @info "Starting FastME on $name"
    # protein WAG, general JTT
    @timeit time "FastME" begin
        run(pipeline(`$(fastmeMP()) \
            -i $(args.inputfile) \
            $(modelparam) \
            -b $(args.nboot) \
            --method=BIONJ \
            --spr \
            --gamma=1.0 \
            --output_tree=$(joinpath(args.outputdir, name * "-tree.nw")) \
            --output_matrix=$(joinpath(args.outputdir, name * "-distancematrix.txt")) \
            --output_info=$(joinpath(args.outputdir, name * "-outputinfo.txt")) \
            --output_boot=$(joinpath(args.outputdir, name * "-boottrees.nw")) \
            --verbose=3`,
        stdout=joinpath(args.outputdir, "fastme.stdout")))
    end

    @info "using Booster to compute support values"
    @timeit time "booster" begin 
        ## calculate support
        run(pipeline(`$(gotree()) compute support tbe --silent \
            -i $(joinpath(args.outputdir, name * "-tree.nw")) \
            -b $(joinpath(args.outputdir, name * "-boottrees.nw")) \
            -o $(joinpath(args.outputdir, name * "-supporttree.nw"))`,
            stderr=joinpath(args.outputdir, "booster.log")))
    end

    end # timeit
    @info "stopping run"
    @info "timing" show(time) println("")
    return 0
end

julia_main()
