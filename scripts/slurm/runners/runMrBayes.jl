using DrWatson
@quickactivate "Doran_etal_2022"
using MrBayes_jll, Gotree_jll, Goalign_jll
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
    @argumentdefault Int 5 nproc "-p" "--nproc" 
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
        "Revmatpr=dirichlet(1,1,1,1,1,1) statefreqpr=fixed(empirical)"
    elseif args.model == :WAG
        "Aamodelpr=fixed(wag) statefreqpr=fixed(empirical)"
    end
    
    mrbayesinputfile = joinpath(args.outputdir, name * ".nex")
    @info "Converting input to nexus format"
    @timeit time "convert input to nexus" begin
        run(pipeline(`$(goalign()) reformat nexus --auto-detect -i $(args.inputfile)`,
        stdout=mrbayesinputfile))
    end
    
    @info "changing dir to: $(args.outputdir)"
    rundir = pwd()
    cd(args.outputdir)

    # format MrBayes commands
    mrbayescommands = """
    #NEXUS
    begin mrbayes;
        set autoclose=yes nowarnings=yes;
        set Seed=123456;
        set Swapseed=5;
        execute $(basename(mrbayesinputfile));
        prset $(modelparam);
        lset nst=6;
        lset rates=equal;
        mcmcp ngen=1000000
            nrun=1
            nchain=$(args.nproc)
            checkfreq=1000
            samplefreq=1000
            printfreq=1000
            stoprule=no
            burninfrac=0.25;
        mcmcp append=no;
        mcmc;
        sump;
        sumt;
        plot;
    end;
    """
    open(joinpath(pwd(), "commands.nex"), "w") do io
        println(io, mrbayescommands)
    end

    @timeit time "mrbayes" begin
        mb_MPI() do exe
            run(pipeline(`mpirun -np $(args.nproc) $(exe) commands.nex`, stdout="mb.out"))
        end
    end # timeit

    @timeit time "convert tree file" begin
        run(pipeline(`$(gotree()) reformat newick -f nexus -i $(name).nex.con.tre`,
            stdout=name * "-supporttree.nw"
        ))
    end

    # move back to run dir
    cd(rundir)
    @info "changing dir to: $(rundir)"

    end # timeit
    @info "stopping run"
    @info "timing" show(time) println("")
    return 0
end

julia_main()
