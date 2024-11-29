"""
Minimal C bindings for JetReconstruction.jl
"""
module C_JetReconstruction

using ..JetReconstruction: JetAlgorithm, RecoStrategy, PseudoJet, ClusterSequence,
                           HistoryElement, jet_reconstruct,
                           exclusive_jets, inclusive_jets
using EnumX

"""
    @enumx StatusCode

An enumeration of common status codes used in the JetReconstruction C-bindings.

## Fields
- `OK` - The operation succeeded.
- `GenericException`- An unspecified error, not covered by other status codes occurred.
- The rest of the status codes have corresponding standard Julia exception.
"""
@enumx StatusCode OK=0 GenericException=1 ArgumentError=2 BoundsError=3 CompositeException=4 DimensionMismatch=5 DivideError=6 DomainError=7 EOFError=8 ErrorException=9 InexactError=10 InitError=11 InterruptException=12 InvalidStateException=13 KeyError=14 LoadError=15 OutOfMemoryError=16 ReadOnlyMemoryError=17 RemoteException=18 MethodError=19 OverflowError=20 ParseError=21 SystemError=22 TypeError=23 UndefRefError=24 UndefVarError=25 StringIndexError=26

function handle_exception(exception)::Cint
    @error exception
    return Cint(exception_to_enum(exception))
end

"""
    exception_to_enum(::Any)::Cint

Helper function matching Julia exception with C-style status code.
# Returns
- A numerical representation of the corresponding status code from the `StatusCode.T`.
"""
exception_to_enum(::Any) = Cint(StatusCode.GenericException)
exception_to_enum(::ArgumentError) = Cint(StatusCode.ArgumentError)
exception_to_enum(::BoundsError) = Cint(StatusCode.BoundsError)
exception_to_enum(::CompositeException) = Cint(StatusCode.CompositeException)
exception_to_enum(::DimensionMismatch) = Cint(StatusCode.DimensionMismatch)
exception_to_enum(::DivideError) = Cint(StatusCode.DivideError)
exception_to_enum(::DomainError) = Cint(StatusCode.DomainError)
exception_to_enum(::EOFError) = Cint(StatusCode.EOFError)
exception_to_enum(::ErrorException) = Cint(StatusCode.ErrorException)
exception_to_enum(::InexactError) = Cint(StatusCode.InexactError)
exception_to_enum(::InitError) = Cint(StatusCode.InitError)
exception_to_enum(::InterruptException) = Cint(StatusCode.InterruptException)
exception_to_enum(::InvalidStateException) = Cint(StatusCode.InvalidStateException)
exception_to_enum(::KeyError) = Cint(StatusCode.KeyError)
exception_to_enum(::LoadError) = Cint(StatusCode.LoadError)
exception_to_enum(::OutOfMemoryError) = Cint(StatusCode.OutOfMemoryError)
exception_to_enum(::ReadOnlyMemoryError) = Cint(StatusCode.ReadOnlyMemoryError)
# exception_to_enum(::Distributed.RemoteException) = Cint(StatusCode.RemoteException)
exception_to_enum(::MethodError) = Cint(StatusCode.MethodError)
exception_to_enum(::OverflowError) = Cint(StatusCode.OverflowError)
# exception_to_enum(::Meta.ParseError) = Cint(StatusCode.ParseError)
exception_to_enum(::SystemError) = Cint(StatusCode.SystemError)
exception_to_enum(::TypeError) = Cint(StatusCode.TypeError)
exception_to_enum(::UndefRefError) = Cint(StatusCode.UndefRefError)
exception_to_enum(::UndefVarError) = Cint(StatusCode.UndefVarError)
exception_to_enum(::StringIndexError) = Cint(StatusCode.StringIndexError)

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

"""
    jetreconstruction_PseudoJet_init(ptr::Ptr{PseudoJet}, px::Cdouble, py::Cdouble, pz::Cdouble, E::Cdouble) -> Cint

C-binding for `PseudoJet` initialization.

# Arguments
- `ptr::Ptr{PseudoJet}`: A pointer to the memory location where the `PseudoJet` object will be stored.
- `px::Cdouble`: The x-component of the momentum.
- `py::Cdouble`: The y-component of the momentum.
- `pz::Cdouble`: The z-component of the momentum.
- `E::Cdouble`: The energy of the jet.

# Returns
- `Cint`: An integer status code indicating the success or failure.

"""
Base.@ccallable function jetreconstruction_PseudoJet_init(ptr::Ptr{PseudoJet}, px::Cdouble,
                                                          py::Cdouble, pz::Cdouble,
                                                          E::Cdouble)::Cint
    try
        pseudojet = PseudoJet(px, py, pz, E)
        unsafe_store!(ptr, pseudojet)
    catch e
        return handle_exception(e)
    end
    return Cint(StatusCode.OK)
end

"""
    struct C_ClusterSequence{T}

A C-compatible struct corresponding to `ClusterSequence`
"""
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

"""
    free_members(ptr::Ptr{C_ClusterSequence{T}}) where {T}

Internal function for freeing dynamically allocated `C_ClusterSequence` members memory.

# Arguments
- `ptr::Ptr{C_ClusterSequence{T}}`: A pointer to a `C_ClusterSequence` structure.
"""
function free_members(ptr::Ptr{C_ClusterSequence{T}}) where {T}
    if ptr != C_NULL
        clusterseq = unsafe_load(ptr)
        Libc.free(clusterseq.jets)
        Libc.free(clusterseq.history)
        # Struct is immutable so pointers can't be assigned NULL and lengths updated to zero (without making a copy)
    end
end

"""
    jetreconstruction_ClusterSequence_free_members_(ptr::Ptr{C_ClusterSequence{PseudoJet}}) -> Cvoid

C-binding for freeing the members of a `C_ClusterSequence` object pointed to by `ptr`.

# Arguments
- `ptr::Ptr{C_ClusterSequence{PseudoJet}}`: A pointer to the `C_ClusterSequence` object whose members are to be freed.

# Returns
- `Cvoid`: This function does not return a value.
"""
Base.@ccallable function jetreconstruction_ClusterSequence_free_members_(ptr::Ptr{C_ClusterSequence{PseudoJet}})::Cvoid
    free_members(ptr)
    return nothing
end

"""
    ClusterSequence(c::C_ClusterSequence{T}) where {T}

Convert a `C_ClusterSequence` object to a `ClusterSequence` object.


# Arguments
- `c::C_ClusterSequence{T}`: The `C_ClusterSequence` object to be converted.

# Returns
- `ClusterSequence{T}`: The converted `ClusterSequence` object.

# Notes
- The array members of output are wrapping the array members of input.
- The input object must remain valid while the output object is used.
"""
function ClusterSequence{T}(c::C_ClusterSequence{T}) where {T}
    jets_v = unsafe_wrap_c_array(c.jets, c.jets_length)
    history_v = unsafe_wrap_c_array(c.history, c.history_length)
    return ClusterSequence{T}(c.algorithm, c.power, c.R, c.strategy, jets_v,
                              c.n_initial_jets,
                              history_v, c.Qtot)
end

"""
    C_ClusterSequence(clustersequence::ClusterSequence{T}) where {T}

Convert a `ClusterSequence` object to a `C_ClusterSequence` object.

# Arguments
- `clustersequence::ClusterSequence{T}`: The `ClusterSequence` object to be converted.

# Returns
- `C_ClusterSequence{T}`: The converted `C_ClusterSequence` object.

# Notes
- The array members of input are deep-copied to output object, it's safe to use output even if input doesn't exist anymore.
- The memory allocated for output data-object members should be freed with  the `free_members` function or equivalent.
"""
function C_ClusterSequence{T}(clustersequence::ClusterSequence{T}) where {T}
    jets_ptr, jets_length = make_c_array(clustersequence.jets)
    history_ptr, history_length = make_c_array(clustersequence.history)
    return C_ClusterSequence{T}(clustersequence.algorithm, clustersequence.power,
                                clustersequence.R, clustersequence.strategy, jets_ptr,
                                jets_length, clustersequence.n_initial_jets, history_ptr,
                                history_length, clustersequence.Qtot)
end

"""
    c_jet_reconstruct(particles::Ptr{T},
                      particles_length::Csize_t,
                      algorithm::JetAlgorithm.Algorithm,
                      R::Cdouble,
                      strategy::RecoStrategy.Strategy,
                      result::Ptr{C_ClusterSequence{U}})::Cint where {T, U}

Internal helper functions for calling `jet_reconstruct` with C-compatible data-structers.

# Arguments
- `particles::Ptr{T}`: Pointer to an array of pseudojet objects used for jet reconstruction.
- `particles_length::Csize_t`: The length of `particles`` array.
- `algorithm::JetAlgorithm.Algorithm`: The algorithm to use for jet reconstruction.
- `R::Cdouble`: The jet radius parameter..
- `strategy::RecoStrategy.Strategy`: The jet reconstruction strategy to use.
- `result::Ptr{C_ClusterSequence{U}}`: A pointer to which a cluster sequence containing the reconstructed jets and the merging history will be stored.

# Returns
- `Cint`: An integer status code indicating the success or failure.

# Notes
- To avoid memory leaks the memory allocated for members of `result` should be freed with `free_members` function or equivalent.
"""
function c_jet_reconstruct(particles::Ptr{T},
                           particles_length::Csize_t,
                           algorithm::JetAlgorithm.Algorithm,
                           R::Cdouble,
                           strategy::RecoStrategy.Strategy,
                           result::Ptr{C_ClusterSequence{U}}) where {T, U}
    try
        particles_v = unsafe_wrap_c_array(particles, particles_length)
        clusterseq = jet_reconstruct(particles_v; p = nothing, algorithm = algorithm, R = R,
                                     strategy = strategy)
        c_clusterseq = C_ClusterSequence{U}(clusterseq)
        unsafe_store!(result, c_clusterseq)
    catch e
        return handle_exception(e)
    end
    return Cint(StatusCode.OK)
end

"""
    jetreconstruction_jet_reconstruct(particles::Ptr{PseudoJet},
                                      particles_length::Csize_t,
                                      algorithm::JetAlgorithm.Algorithm,
                                      R::Cdouble,
                                      strategy::RecoStrategy.Strategy,
                                      result::Ptr{C_ClusterSequence{PseudoJet}})::Cint

C-binding for `jet_reconstruct`.

# Arguments
- `particles::Ptr{PseudoJet}`: Pointer to an array of pseudojet objects used for jet reconstruction.
- `particles_length::Csize_t`: The length of `particles`` array.
- `algorithm::JetAlgorithm.Algorithm`: The algorithm to use for jet reconstruction.
- `R::Cdouble`: The jet radius parameter..
- `strategy::RecoStrategy.Strategy`: The jet reconstruction strategy to use.
- `result::Ptr{C_ClusterSequence{PseudoJet}}`: A pointer to which a cluster sequence containing the reconstructed jets and the merging history will be stored.

# Returns
- `Cint`: An integer status code indicating the success or failure.

# Notes
- To avoid memory leaks the memory allocated for members of `result` should be freed with `jetreconstruction_ClusterSequence_free_members_`.
"""
Base.@ccallable function jetreconstruction_jet_reconstruct(particles::Ptr{PseudoJet},
                                                           particles_length::Csize_t,
                                                           algorithm::JetAlgorithm.Algorithm,
                                                           R::Cdouble,
                                                           strategy::RecoStrategy.Strategy,
                                                           result::Ptr{C_ClusterSequence{PseudoJet}})::Cint
    return c_jet_reconstruct(particles, particles_length, algorithm, R, strategy, result)
end

"""
    struct C_JetsResult{T}

A C-compatible wrapper for array-like data

# Fields
- `data::Ptr{T}`: A pointer to the data of type `T`.
- `length::Csize_t`: The length of the data.
"""
struct C_JetsResult{T}
    data::Ptr{T}
    length::Csize_t
end

"""
    free_members(ptr::Ptr{C_JetsResult{T}}) where {T}

Internal function for freeing dynamically allocated `C_JetsResult` members memory.

# Arguments
- `ptr::Ptr{C_JetsResult{T}}`: A pointer to a `C_JetsResult` structure.
"""
function free_members(ptr::Ptr{C_JetsResult{T}}) where {T}
    if ptr != C_NULL
        result = unsafe_load(ptr)
        Libc.free(result.data)
        # Struct is immutable so pointer can't be assigned NULL and length updated to zero (without making a copy)
    end
end

"""
    jetreconstruction_JetsResult_free_members_(ptr::Ptr{C_JetsResult{PseudoJet}})::Cvoid

C-binding for freeing the members of a `C_JetsResult{PseudoJet}` object pointed to by `ptr`.

# Arguments
- `ptr::Ptr{C_JetsResult{PseudoJet}}`: A pointer to the `C_JetsResult{PseudoJet}` structure whose members are to be freed.

# Returns
- `Cvoid`: This function does not return any value.
"""
Base.@ccallable function jetreconstruction_JetsResult_free_members_(ptr::Ptr{C_JetsResult{PseudoJet}})::Cvoid
    free_members(ptr)
    return nothing
end

"""
    jets_selection(selector, clustersequence::Ptr{C_ClusterSequence{T}},
                   result::Ptr{C_JetsResult{U}}; kwargs...)::Cint where {T, U}

An internal helper function for calling calling functions selecting jets from a given cluster sequence.

# Arguments
- `selector`: A function that takes a `ClusterSequence{T}` and returns a list of jets.
- `clustersequence::Ptr{C_ClusterSequence{T}}`: A pointer to a C-compatible cluster sequence.
- `result::Ptr{C_JetsResult{U}}`: A pointer to a C-compatible structure where the result will be stored.
- `kwargs...`: Additional keyword arguments to be passed to the selector function.

# Returns
- `Cint`: An integer status code indicating the success or failure.
"""
function jets_selection(selector, clustersequence::Ptr{C_ClusterSequence{T}},
                        result::Ptr{C_JetsResult{U}}; kwargs...)::Cint where {T, U}
    try
        c_clusterseq = unsafe_load(clustersequence)
        clusterseq = ClusterSequence{T}(c_clusterseq)
        jets_result = selector(clusterseq; T = U, kwargs...)
        c_results = C_JetsResult{U}(make_c_array(jets_result)...)
        unsafe_store!(result, c_results)
    catch e
        return handle_exception(e)
    end
    return Cint(StatusCode.OK)
end

"""
    jetreconstruction_exclusive_jets_dcut(clustersequence::Ptr{C_ClusterSequence{PseudoJet}},
                                          dcut::Cdouble,
                                          result::Ptr{C_JetsResult{PseudoJet}}) -> Cint

C-binding for `exclusive_jets` with a cut on the maximum distance parameter.

# Arguments
- `clustersequence::Ptr{C_ClusterSequence{PseudoJet}}`: A pointer to the cluster sequence object containing the clustering history and jets.
- `dcut::Cdouble`: The distance parameter used to define the exclusive jets.
- `result::Ptr{C_JetsResult{PseudoJet}}`: A pointer to the results.

# Returns
- `Cint`: An integer status code indicating the success or failure.

# Notes
- To avoid memory leaks the memory allocated for members of `result` should be freed with `jetreconstruction_JetsResult_free_members_`.
"""
Base.@ccallable function jetreconstruction_exclusive_jets_dcut(clustersequence::Ptr{C_ClusterSequence{PseudoJet}},
                                                               dcut::Cdouble,
                                                               result::Ptr{C_JetsResult{PseudoJet}})::Cint
    return jets_selection(exclusive_jets, clustersequence, result; dcut = dcut)
end

"""
    jetreconstruction_exclusive_jets_njets(clustersequence::Ptr{C_ClusterSequence{PseudoJet}},
                                           njets::Csize_t,
                                           result::Ptr{C_JetsResult{PseudoJet}}) -> Cint

C-binding for `exclusive_jets` with a specific number of jets.

# Arguments
- `clustersequence::Ptr{C_ClusterSequence{PseudoJet}}`: A pointer to the cluster sequence object containing the clustering history and jets.
- `njets::Csize_t`: The number of exclusive jets to be calculated.
- `result::Ptr{C_JetsResult{PseudoJet}}`: A pointer to the results.

# Returns
- `Cint`: An integer status code indicating the success or failure.

# Notes
- To avoid memory leaks the memory allocated for members of `result` should be freed with `jetreconstruction_JetsResult_free_members_`.
"""
Base.@ccallable function jetreconstruction_exclusive_jets_njets(clustersequence::Ptr{C_ClusterSequence{PseudoJet}},
                                                                njets::Csize_t,
                                                                result::Ptr{C_JetsResult{PseudoJet}})::Cint
    return jets_selection(exclusive_jets, clustersequence, result; njets = njets)
end

"""
    jetreconstruction_inclusive_jets(clustersequence::Ptr{C_ClusterSequence{PseudoJet}},
                                     ptmin::Cdouble,
                                     result::Ptr{C_JetsResult{PseudoJet}}) -> Cint

C-binding for `inclusive_jets`.

# Arguments
- `clustersequence::Ptr{C_ClusterSequence{PseudoJet}}`: A pointer to the cluster sequence object containing the clustering history and jets.
- `ptmin::Cdouble`: The minimum transverse momentum (pt) threshold for the inclusive jets.
- `result::Ptr{C_JetsResult{PseudoJet}}`: A pointer to the results.

# Returns
- `Cint`: An integer status code indicating the success or failure.

# Notes
- To avoid memory leaks the memory allocated for members of `result` should be freed with `jetreconstruction_JetsResult_free_members_`.
"""
Base.@ccallable function jetreconstruction_inclusive_jets(clustersequence::Ptr{C_ClusterSequence{PseudoJet}},
                                                          ptmin::Cdouble,
                                                          result::Ptr{C_JetsResult{PseudoJet}})::Cint
    return jets_selection(inclusive_jets, clustersequence, result; ptmin = ptmin)
end

end # module C_JetReconstruction
