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
    N2PlainWorkspace(algorithm::JetAlgorithm.Algorithm)

Reusable full-semantics storage for one N2Plain or electron-positron
reconstruction worker.

`algorithm` selects the internal jet and scratch storage. Input particles may
be any type supported by the corresponding reconstruction algorithm.

The `ClusterSequence` produced inside [`with_n2plain_reconstruction`](@ref)
or [`with_ee_reconstruction`](@ref) borrows the workspace's jets and history
and is overwritten by its next reconstruction. Copy values that must outlive
that operation. A workspace must not be used concurrently or reentrantly.
"""
mutable struct N2PlainWorkspace{A, J, S}
    scratch::S
    jets::Vector{J}
    history::Vector{HistoryElement}
end

function N2PlainWorkspace(algorithm::JetAlgorithm.Algorithm)
    if is_pp(algorithm)
        return N2PlainWorkspace{algorithm, PseudoJet, PPPlainScratch}(PPPlainScratch(),
                                                                      PseudoJet[],
                                                                      HistoryElement[])
    elseif is_ee(algorithm)
        scratch = StructArray{EERecoJet}(undef, 0)
        return N2PlainWorkspace{algorithm, EEJet, typeof(scratch)}(scratch,
                                                                   EEJet[],
                                                                   HistoryElement[])
    end
    throw(ArgumentError("Unsupported jet algorithm: $algorithm"))
end

function _prepare_n2plain_recombination_jets!(jets::Vector{J},
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

function prepare_n2plain_recombination_jets!(workspace::N2PlainWorkspace,
                                             particles::AbstractVector;
                                             preprocess = preprocess_escheme)
    return _prepare_n2plain_recombination_jets!(workspace.jets,
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

"""
    with_n2plain_reconstruction(
        f,
        workspace,
        particles;
        p = nothing,
        R = nothing,
        recombine = addjets_escheme,
        preprocess = preprocess_escheme,
    )

Run pp N2Plain reconstruction with workspace-owned storage and call `f` with
the borrowed [`ClusterSequence`](@ref). Construct the workspace with the pp
algorithm to use, for example `N2PlainWorkspace(JetAlgorithm.AntiKt)`.

When `R` is omitted, reconstruction uses `R = 1.0`, matching
[`plain_jet_reconstruct`](@ref).

The sequence, its jets, and its history are valid only until the workspace's
next reconstruction. Copy values that must escape the callback. Each
concurrently executing worker must own a distinct workspace, and the caller
must prevent concurrent or reentrant use of the same workspace.
"""
function with_n2plain_reconstruction(f::F,
                                     workspace::N2PlainWorkspace{A, PseudoJet,
                                                                 PPPlainScratch},
                                     particles::AbstractVector;
                                     p::Union{Real, Nothing} = nothing,
                                     R = nothing,
                                     recombine = addjets_escheme,
                                     preprocess = preprocess_escheme) where {F, A}
    clusterseq = _n2plain_reconstruct_with_workspace!(workspace,
                                                      particles;
                                                      p = p,
                                                      R = R,
                                                      recombine = recombine,
                                                      preprocess = preprocess)

    return f(clusterseq)
end

"""
    with_ee_reconstruction(
        f,
        workspace,
        particles;
        p = nothing,
        R = nothing,
        recombine = addjets_escheme,
        preprocess = preprocess_escheme,
        γ = nothing,
        β = nothing,
    )

Run electron-positron reconstruction with workspace-owned storage and call
`f` with the borrowed [`ClusterSequence`](@ref). Construct the workspace
with the e⁺e⁻ algorithm to use, for example
`N2PlainWorkspace(JetAlgorithm.Durham)`. Input particles may have any
supported four-momentum type.

When `R` is omitted, reconstruction uses `R = 4.0`, matching
[`ee_genkt_algorithm`](@ref). Durham always uses its conventional nominal
value `R = 4.0`. The sequence, its jets, and its history are overwritten by
the next reconstruction using this workspace. Copy values that must outlive
that operation, and do not use one workspace concurrently or reentrantly.
"""
function with_ee_reconstruction(f::F,
                                workspace::N2PlainWorkspace{A, EEJet},
                                particles::AbstractVector;
                                p::Union{Real, Nothing} = nothing,
                                R = nothing,
                                recombine = addjets_escheme,
                                preprocess = preprocess_escheme,
                                γ::Union{Real, Nothing} = nothing,
                                β::Union{Real, Nothing} = nothing) where {F, A}
    clusterseq = _n2plain_reconstruct_with_workspace!(workspace,
                                                      particles;
                                                      p = p,
                                                      R = R,
                                                      recombine = recombine,
                                                      preprocess = preprocess,
                                                      γ = γ,
                                                      β = β)

    return f(clusterseq)
end

function _release_capacity!(::StructArray{EERecoJet})
    return StructArray{EERecoJet}(undef, 0)
end

"""
    release_n2plain_workspace_capacity!(workspace)

Release all event-size-dependent capacity retained by `workspace`.

The caller must ensure that `workspace` is not in use by another task.
"""
function release_n2plain_workspace_capacity!(workspace::N2PlainWorkspace{A, J}) where {A,
                                                                                       J}
    workspace.scratch = _release_capacity!(workspace.scratch)
    workspace.jets = J[]
    workspace.history = HistoryElement[]

    return nothing
end
