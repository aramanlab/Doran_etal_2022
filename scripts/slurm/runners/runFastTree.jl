using DrWatson
@quickactivate "Doran_etal_2022"
using FastTree_jll, Gotree_jll, Goalign_jll
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
        "-nt" # by default JC69+CAT
    elseif args.model == :WAG
        "-wag"
    end
    
    @info "Converting to fasta format"
    @timeit time "fasta conversion" begin
        run(pipeline(`$(goalign()) reformat fasta --auto-detect -i $(args.inputfile)`,
            stdout=joinpath(args.outputdir, name * ".fasta")))
    end

    @info "Starting FastTree on $name"
    # protein WAG, general JTT
    @timeit time "fasttree" begin
        run(pipeline(`$(fasttreeMP()) \
        $(modelparam) \
        -gamma \
        -boot $(args.nboot)
        -log $(joinpath(args.outputdir, name * ".log"))
        `,
        stdin=joinpath(args.outputdir, name * ".fasta"),
        stderr=joinpath(args.outputdir, name * "_fasttree.out"),
        stdout=joinpath(args.outputdir, name * "-supporttree.nw")))
    end

    end # timeit
    @info "stopping run"
    @info "timing" show(time) println("")
    return 0
end

julia_main()
