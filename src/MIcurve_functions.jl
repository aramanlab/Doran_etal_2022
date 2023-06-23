using NewickTree
using NewickTreeTools

function clusters_per_cutlevel(distfun::Function, tree::Node, ncuts::Number)
    minmax = extrema(mapinternalnodes(distfun, tree, tree))
    cuts = range(0, minmax[2], length=ncuts)
    clusts = [cuttree(distfun, tree, cut) for cut in cuts]
    clustmappings = map(c->getleafnames.(c), clusts)
    clustersmps = [vcat(clustmapping...) for clustmapping in clustmappings]
    clusterids = [Int.(vcat([zeros(length(c)) .+ j for (j, c) in enumerate(clustmapping)]...)) for clustmapping in clustmappings];
    return clusterids, clustersmps
end

function pairedMIagainstmetacolumn(metacolumns, IDS, clusterids, clustersmps; bootstrap=false, mask=nothing)
    tstat_MI = zeros(length(clusterids), size(metacolumns, 2))
    for (i, mcol) in enumerate(eachcol(metacolumns))
        # cat = levelorder(categorical(mcol))
        # pcat = cat .== cat'
        pcat = mcol .== permutedims(replace(mcol, ""=>"missing"))
        tstat_MI[:, i] .= collectMI_across_treedepth(clusterids, clustersmps, IDS, pcat; bootstrap, mask)
    end
    DataFrame(tstat_MI, names(metacolumns)) |> stack |> df->rename!(df,["taxaID","MI"]);
end

function pairedMIagainstnumericcolumn(metacolumns, IDS, clusterids, clustersmps; bootstrap=false, mask=nothing)
    tstat_MI = zeros(length(clusterids), size(metacolumns, 2))
    for (i, mcol) in enumerate(eachcol(metacolumns))
        # cat = levelorder(categorical(mcol))
        # pcat = cat .== cat'
        pcat = abs.(mcol .- permutedims(mcol))
        tstat_MI[:, i] .= collectMI_across_treedepth(clusterids, clustersmps, IDS, pcat; bootstrap, mask)
    end
    DataFrame(tstat_MI, names(metacolumns)) |> stack |> df->rename!(df,["taxaID","MI"]);
end

function collectMI_across_treedepth(clusterids, clustersmps, IDS, ptax; bootstrap=false, mask=nothing)
    uppertriangle = triu(trues(size(ptax)), 1);
    uppertriangle = isnothing(mask) ? uppertriangle : uppertriangle[mask, mask]
    # ptax = if isnothing(mask) ptax[uppertriangle] else ptax[mask, mask][uppertriangle] end
    map(clusterids, clustersmps) do cids, smps; ptax, IDS, mask
        clustorder = indexin(IDS, smps)
        pcids = cids[clustorder] .== cids[clustorder]'
        pcids = if isnothing(mask) pcids else pcids[mask, mask] end
        wptax = if isnothing(mask) ptax else ptax[mask, mask] end
        pcids = pcids[uppertriangle]
        wptax = wptax[uppertriangle]
        if bootstrap
            vals_idx = sample(axes(pcids, 1), length(pcids), replace=true)
            pcids = pcids[vals_idx]
            wptax = wptax[vals_idx]
        end
        empiricalMI(wptax, pcids)
    end
end