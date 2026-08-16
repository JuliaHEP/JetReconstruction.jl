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
    N2PlainWorkspace([::Type{J} = PseudoJet])

Reusable full-semantics storage for one N2Plain reconstruction worker.

`J` selects the reconstruction family: use `PseudoJet` for pp algorithms and
`EEJet` for electron-positron algorithms.

The `ClusterSequence` produced inside [`with_n2plain_reconstruction`](@ref)
borrows the workspace's jets and history and is overwritten by its next
reconstruction. Callers may provide their own reusable output vectors to
post-processing functions such as [`inclusive_jets!`](@ref). Copy values that
must outlive those operations. A workspace must not be used concurrently or
reentrantly.
"""
mutable struct N2PlainWorkspace{J, S}
    scratch::S
    jets::Vector{J}
    history::Vector{HistoryElement}
end

N2PlainWorkspace() = N2PlainWorkspace(PseudoJet)

function N2PlainWorkspace(::Type{PseudoJet})
    return N2PlainWorkspace{PseudoJet, PPPlainScratch}(PPPlainScratch(),
                                                       PseudoJet[],
                                                       HistoryElement[])
end

function N2PlainWorkspace(::Type{EEJet})
    scratch = StructArray{EERecoJet}(undef, 0)
    return N2PlainWorkspace{EEJet, typeof(scratch)}(scratch,
                                                    EEJet[],
                                                    HistoryElement[])
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
        algorithm,
        p = nothing,
        R = nothing,
        recombine = addjets_escheme,
        preprocess = preprocess_escheme,
        γ = nothing,
        β = nothing,
    )

Run N2Plain reconstruction with workspace-owned storage and call `f` with the
borrowed [`ClusterSequence`](@ref).

When `R` is omitted, pp reconstruction uses `R = 1.0` and electron-positron
reconstruction uses `R = 4.0`, matching [`plain_jet_reconstruct`](@ref) and
[`ee_genkt_algorithm`](@ref), respectively. Durham always uses its conventional
nominal value `R = 4.0`.

The sequence, its jets, and its history are valid only until the workspace's
next reconstruction. Copy values that must escape the callback. Each
concurrently executing worker must own a distinct workspace, and the caller
must prevent concurrent or reentrant use of the same workspace.
"""
function with_n2plain_reconstruction(f::F,
                                     workspace::N2PlainWorkspace,
                                     particles::AbstractVector;
                                     algorithm::JetAlgorithm.Algorithm,
                                     p::Union{Real, Nothing} = nothing,
                                     R = nothing,
                                     recombine = addjets_escheme,
                                     preprocess = preprocess_escheme,
                                     γ::Union{Real, Nothing} = nothing,
                                     β::Union{Real, Nothing} = nothing) where {F}
    clusterseq = _n2plain_reconstruct_with_workspace!(workspace,
                                                      particles;
                                                      algorithm = algorithm,
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
function release_n2plain_workspace_capacity!(workspace::N2PlainWorkspace{J}) where {J}
    workspace.scratch = _release_capacity!(workspace.scratch)
    workspace.jets = J[]
    workspace.history = HistoryElement[]

    return nothing
end
