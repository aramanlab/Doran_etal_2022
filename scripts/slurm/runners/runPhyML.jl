using DrWatson
@quickactivate "Doran_etal_2022"
using PhyML_jll #, Gotree_jll
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
        ["-dnt", "-mJC69", "-fe"]
    elseif args.model == :WAG
        [ "-daa", "-mWAG", "-fe"]
    end

    cp(args.inputfile, joinpath(args.outputdir, basename(args.inputfile)), force=true)
    
    @info "Starting PhyML on $name"
    # protein WAG, general JTT
    @timeit time "PhyML" begin
        run(pipeline(`$(phyml()) \
            $modelparam \
            -i $(joinpath(args.outputdir, basename(args.inputfile))) \
            -o tlr \
            --search SPR \
            --r_seed 123456 \
            --rand_start \
            --n_rand_starts 5 \
            --bootstrap -1`,
        stdout=joinpath(args.outputdir, name * "_phyml.out")))
    end
    mv(joinpath(args.outputdir, basename(args.inputfile) * "_phyml_tree.txt"),
       joinpath(args.outputdir, basename(args.inputfile) * "-supporttree.txt")
    )
    end # timeit
    @info "stopping run"
    @info "timing" show(time) println("")
    return 0
end

julia_main()
