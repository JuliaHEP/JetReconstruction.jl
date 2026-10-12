"""
Reusable compact state for one hadron-collider N2Plain reconstruction.

The vectors have exactly the active event length. Their backing storage may be
retained by Julia when they are resized for a smaller event.
"""
mutable struct PPPlainScratch
    kt2::Vector{Float64}
    phi::Vector{Float64}
    rapidity::Vector{Float64}
    nn::Vector{Int}
    nndist::Vector{Float64}
    nndij::Vector{Float64}
    clusterseq_index::Vector{Int}
end

function PPPlainScratch(N::Integer = 0)
    N >= 0 || throw(ArgumentError("N must be non-negative, got $N"))

    return PPPlainScratch(Vector{Float64}(undef, N),
                          Vector{Float64}(undef, N),
                          Vector{Float64}(undef, N),
                          Vector{Int}(undef, N),
                          Vector{Float64}(undef, N),
                          Vector{Float64}(undef, N),
                          Vector{Int}(undef, N))
end

"""
Resize every compact pp reconstruction vector to the active event size `N`.
"""
function ensure_length!(scratch::PPPlainScratch, N::Int)
    N >= 0 || throw(ArgumentError("N must be non-negative, got $N"))

    resize!(scratch.kt2, N)
    resize!(scratch.phi, N)
    resize!(scratch.rapidity, N)
    resize!(scratch.nn, N)
    resize!(scratch.nndist, N)
    resize!(scratch.nndij, N)
    resize!(scratch.clusterseq_index, N)

    return nothing
end

"""Resize electron-positron compact reconstruction storage to length `N`."""
function ensure_length!(scratch::StructArray{EERecoJet}, N::Int)
    N >= 0 || throw(ArgumentError("N must be non-negative, got $N"))

    # Resizing the StructArray also gives every component vector the exact
    # active event length.
    resize!(scratch, N)
    return nothing
end

"""
    PlainWorkspace(algorithm::JetAlgorithm.Algorithm)

Reusable full-semantics storage for one N2Plain or electron-positron
reconstruction worker.

`algorithm` selects the internal jet and scratch storage. Input particles may
be any type supported by the corresponding reconstruction algorithm.

The `ClusterSequence` produced inside [`with_n2plain_reconstruction`](@ref)
or [`with_ee_reconstruction`](@ref) borrows the workspace's jets and history
and is overwritten by its next reconstruction. Copy values that must outlive
that operation. A workspace must not be used concurrently or reentrantly.
"""
mutable struct PlainWorkspace{A, J, S}
    scratch::S
    jets::Vector{J}
    history::Vector{HistoryElement}
end

function PlainWorkspace(algorithm::JetAlgorithm.Algorithm)
    if is_pp(algorithm)
        return PlainWorkspace{algorithm, PseudoJet, PPPlainScratch}(PPPlainScratch(),
                                                                    PseudoJet[],
                                                                    HistoryElement[])
    elseif is_ee(algorithm)
        scratch = StructArray{EERecoJet}(undef, 0)
        return PlainWorkspace{algorithm, EEJet, typeof(scratch)}(scratch,
                                                                 EEJet[],
                                                                 HistoryElement[])
    end
    throw(ArgumentError("Unsupported jet algorithm: $algorithm"))
end

function _prepare_plain_recombination_jets!(jets::Vector{J},
                                            particles::AbstractVector{T};
                                            preprocess = preprocess_escheme) where {
                                                                                    J <:
                                                                                    Union{PseudoJet,
                                                                                          EEJet},
                                                                                    T}
    Base.mightalias(jets, particles) &&
        throw(ArgumentError("reusable jet storage must not alias the input particle vector"))

    N = length(particles)
    empty!(jets)
    _sizehint_for_reuse!(jets, 2 * N)

    if isnothing(preprocess)
        if T == J
            append!(jets, particles)
        else
            for (i, particle) in enumerate(particles)
                push!(jets, J(particle; cluster_hist_index = i))
            end
        end
    else
        for (i, particle) in enumerate(particles)
            push!(jets, preprocess(particle, J; cluster_hist_index = i))
        end
    end

    return jets
end

function prepare_plain_recombination_jets!(workspace::PlainWorkspace,
                                           particles::AbstractVector;
                                           preprocess = preprocess_escheme)
    return _prepare_plain_recombination_jets!(workspace.jets,
                                              particles;
                                              preprocess = preprocess)
end

function _release_capacity!(scratch::PPPlainScratch)
    scratch.kt2 = Float64[]
    scratch.phi = Float64[]
    scratch.rapidity = Float64[]
    scratch.nn = Int[]
    scratch.nndist = Float64[]
    scratch.nndij = Float64[]
    scratch.clusterseq_index = Int[]
    return scratch
end

function _release_capacity!(::StructArray{EERecoJet})
    return StructArray{EERecoJet}(undef, 0)
end

"""
    release_plain_workspace_capacity!(workspace)

Release all event-size-dependent capacity retained by `workspace`.

The caller must ensure that `workspace` is not in use by another task.
"""
function release_plain_workspace_capacity!(workspace::PlainWorkspace{A, J}) where {A,
                                                                                   J}
    workspace.scratch = _release_capacity!(workspace.scratch)
    workspace.jets = J[]
    workspace.history = HistoryElement[]

    return nothing
end

"""
Reusable state for one full-semantics N2Tiled reconstruction worker.

The workspace owns reusable:

- tiled scratch;
- reconstructed jets;
- complete clustering history.

The ClusterSequence produced from this workspace is borrowed. It is overwritten
by the next reconstruction using the same workspace. Post-processing functions
such as `inclusive_jets` return independently owned results.

A workspace must never be used concurrently or reentrantly.
"""
mutable struct TiledWorkspace
    scratch::TiledScratch
    jets::Vector{PseudoJet}
    history::Vector{HistoryElement}
end

function TiledWorkspace()
    return TiledWorkspace(TiledScratch(),
                          PseudoJet[],
                          HistoryElement[])
end

"""
    release_tiled_workspace_capacity!(workspace)

Release storage retained by `workspace`.

The small fixed-size `tile_union` buffer is retained. All event-size-dependent
vectors and cached tiling matrices become eligible for garbage collection.

The caller must ensure that `workspace` is not in use by another task.
"""
function release_tiled_workspace_capacity!(workspace::TiledWorkspace)
    workspace.jets = PseudoJet[]
    workspace.history = HistoryElement[]

    scratch = workspace.scratch

    scratch.eta = Float64[]
    scratch.tiledjets = TiledJet[]
    scratch.NNs = TiledJet[]
    scratch.dij = Float64[]

    scratch.tiling_cache = Dict{Tuple{Int, Int}, TilingArrays}()

    return nothing
end
