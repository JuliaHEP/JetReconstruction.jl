"""
Minimal C bindings for JetReconstruction.jl
"""
module C_JetReconstruction

using ..JetReconstruction: JetAlgorithm, RecoStrategy, PseudoJet, ClusterSequence,
                           HistoryElement, jet_reconstruct,
                           exclusive_jets, inclusive_jets

"""
    unsafe_wrap_c_array(ptr::Ptr{T}, array_length::Csize_t) where {T}

Wraps a C array into a Julia `Vector` for both bits and non-bits types.

# Arguments
- `ptr::Ptr{T}`: A pointer to the C array.
- `array_length::Csize_t`: The length of the C array.

# Returns
- A Julia `Vector{T}` containing the elements of the C array.

# Safety
This function use 'unsafe' methods and has undefined behaviour
if pointer isn't valid or length isn't correct.
"""
function unsafe_wrap_c_array(ptr::Ptr{T}, array_length::Csize_t) where {T}
    if isbitstype(T)
        return unsafe_wrap(Vector{T}, ptr, array_length)
    end

    vec = Vector{T}(undef, array_length)
    for i in eachindex(vec)
        @inbounds vec[i] = unsafe_load(ptr, i)
    end
    return vec
end

"""
    make_c_array(v::Vector{T}) where {T}

Helper function for converting a Julia vector to a C-style array.
A C-style array is dynamically allocated and contents of input vector copied to it.

# Arguments
- `v::Vector{T}`: A Julia vector of type `T`.

# Returns
- `ptr::Ptr{T}`: A pointer to the allocated C-style array.
- `len::Int`: The length of the vector.

# Throws
- `OutOfMemoryError`: If memory allocation fails.

# Notes
- The caller is responsible for freeing the allocated memory using `Libc.free(ptr)`.
"""
function make_c_array(v::Vector{T}) where {T}
    len = length(v)
    ptr = Ptr{T}(Libc.malloc(len * sizeof(T)))
    if ptr == C_NULL
        throw(OutOfMemoryError("Libc.malloc failed to allocate memory"))
    end
    for i in 1:len
        @inbounds unsafe_store!(ptr, v[i], i)
    end
    return ptr, len
end

Base.@ccallable function jetreconstruction_PseudoJet_init(ptr::Ptr{PseudoJet}, px::Cdouble,
                                                          py::Cdouble, pz::Cdouble,
                                                          E::Cdouble)::Cint
    pseudojet = PseudoJet(px, py, pz, E)
    unsafe_store!(ptr, pseudojet)
    return 0
end

struct C_ClusterSequence{T}
    algorithm::JetAlgorithm.Algorithm
    power::Cdouble
    R::Cdouble
    strategy::RecoStrategy.Strategy
    jets::Ptr{T}
    jets_length::Csize_t
    n_initial_jets::Clong
    history::Ptr{HistoryElement}
    history_length::Csize_t
    Qtot::Cdouble
end

function free_members(ptr::Ptr{C_ClusterSequence{T}}) where {T}
    if ptr != C_NULL
        clusterseq = unsafe_load(ptr)
        Libc.free(clusterseq.jets)
        Libc.free(clusterseq.history)
        # Struct is immutable so pointers can't be assigned NULL and lengths updated to zero (without making a copy)
    end
end
Base.@ccallable function jetreconstruction_ClusterSequence_free_members_(ptr::Ptr{C_ClusterSequence{PseudoJet}})::Cvoid
    free_members(ptr)
    return nothing
end

function ClusterSequence{T}(c::C_ClusterSequence{T}) where {T}
    jets_v = unsafe_wrap_c_array(c.jets, c.jets_length)
    history_v = unsafe_wrap_c_array(c.history, c.history_length)
    return ClusterSequence{T}(c.algorithm, c.power, c.R, c.strategy, jets_v,
                              c.n_initial_jets,
                              history_v, c.Qtot)
end

function C_ClusterSequence{T}(clustersequence::ClusterSequence{T}) where {T}
    jets_ptr, jets_length = make_c_array(clustersequence.jets)
    history_ptr, history_length = make_c_array(clustersequence.history)
    return C_ClusterSequence{T}(clustersequence.algorithm, clustersequence.power,
                                clustersequence.R, clustersequence.strategy, jets_ptr,
                                jets_length, clustersequence.n_initial_jets, history_ptr,
                                history_length, clustersequence.Qtot)
end

function c_jet_reconstruct(particles::Ptr{T},
                           particles_length::Csize_t,
                           algorithm::JetAlgorithm.Algorithm,
                           R::Cdouble,
                           strategy::RecoStrategy.Strategy,
                           result::Ptr{C_ClusterSequence{U}}) where {T, U}
    particles_v = unsafe_wrap_c_array(particles, particles_length)
    clusterseq = jet_reconstruct(particles_v; p = nothing, algorithm = algorithm, R = R,
                                 strategy = strategy)
    c_clusterseq = C_ClusterSequence{U}(clusterseq)
    unsafe_store!(result, c_clusterseq)
    return 0
end

Base.@ccallable function jetreconstruction_jet_reconstruct(particles::Ptr{PseudoJet},
                                                           particles_length::Csize_t,
                                                           algorithm::JetAlgorithm.Algorithm,
                                                           R::Cdouble,
                                                           strategy::RecoStrategy.Strategy,
                                                           result::Ptr{C_ClusterSequence{PseudoJet}})::Cint
    c_jet_reconstruct(particles, particles_length, algorithm, R, strategy, result)
    return 0
end

struct C_JetsResult{T}
    data::Ptr{T}
    length::Csize_t
end

function free_members(ptr::Ptr{C_JetsResult{T}}) where {T}
    if ptr != C_NULL
        result = unsafe_load(ptr)
        Libc.free(result.data)
        # Struct is immutable so pointer can't be assigned NULL and length updated to zero (without making a copy)
    end
end

Base.@ccallable function jetreconstruction_JetsResult_free_members_(ptr::Ptr{C_JetsResult{PseudoJet}})::Cvoid
    free_members(ptr)
    return nothing
end

function jets_selection(selector, clustersequence::Ptr{C_ClusterSequence{T}},
                        result::Ptr{C_JetsResult{U}}; kwargs...)::Cint where {T, U}
    c_clusterseq = unsafe_load(clustersequence)
    clusterseq = ClusterSequence{T}(c_clusterseq)
    jets_result = selector(clusterseq; kwargs...)
    println(jets_result)
    c_results = C_JetsResult{U}(C_NULL, 0) # TODO convert and write to result
    unsafe_store!(result, c_results)
    return 0
end

Base.@ccallable function jetreconstruction_exclusive_jets_dcut(clustersequence::Ptr{C_ClusterSequence{PseudoJet}},
                                                               dcut::Cdouble,
                                                               result::Ptr{C_JetsResult{PseudoJet}})::Cint
    return jets_selection(exclusive_jets, clustersequence, result; dcut = dcut)
end

Base.@ccallable function jetreconstruction_exclusive_jets_njets(clustersequence::Ptr{C_ClusterSequence{PseudoJet}},
                                                                njets::Csize_t,
                                                                result::Ptr{C_JetsResult{PseudoJet}})::Cint
    return jets_selection(exclusive_jets, clustersequence, result; njets = njets)
end

Base.@ccallable function jetreconstruction_inclusive_jets(clustersequence::Ptr{C_ClusterSequence{PseudoJet}},
                                                          ptmin::Cdouble,
                                                          result::Ptr{C_JetsResult{PseudoJet}})::Cint
    return jets_selection(inclusive_jets, clustersequence, result; ptmin = ptmin)
end

end # module C_JetReconstruction
