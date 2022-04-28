using DrWatson
@quickactivate "Doran_etal_2022"
using RAxML_jll, Gotree_jll
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
    @argumentdefault Int 5 nthreads "-t" "--nthreads"
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
        "-mGTRGAMMAX"
    elseif args.model == :WAG
        "-mPROTGAMMAWAG"
    end

    cp(args.inputfile, joinpath(args.outputdir, basename(args.inputfile)), force=true)

    @info "changing dir to: $(args.outputdir)"
    rundir = pwd()
    cd(args.outputdir)
    
    @timeit time "raxml" begin
        run(pipeline(`$(raxmlHPC_PTHREADS_AVX()) \
            -f a \
            -T $(args.nthreads) \
            -s $(basename(args.inputfile)) \
            $modelparam \
            -N $(args.nboot) \
            -p 123456 \
            -x 123456 \
            -n T01`,
            stdout="raxml.out"
        ))
    end
    cp("RAxML_bipartitions.T01", "RAxML-supporttree.nw")
    
    # move back to run dir
    cd(rundir)
    @info "changing dir to: $(rundir)"

    end # timeit
    @info "stopping run"
    @info "timing" show(time) println("")
    return 0
end

julia_main()
