/*
 * SPDX-FileCopyrightText: 2022-2024 JetReconstruction.jl authors, CERN
 * SPDX-License-Identifier: MIT
 */

/**
 * @file JetReconstruction.h
 * @brief Header file for the C-bindings of JetReconstruction.jl.
 *
 * This file contains the definitions of data structures and functions used in
 * JetReconstruction.jl. The library provides functionality for efficient jet
 * reconstruction.
 */

#include <stddef.h>

/*Special states in a history.*/

/**Special history state: Invalid child for this jet, meaning it did not
 * recombine further */
#define JETRECONSTRUCTION_INVALID -3
/**Special history state: Original cluster so it has no parent */
#define JETRECONSTRUCTION_NONEXISTENTPARENT -2
/**Special history state: Cluster recombined with beam" */
#define JETRECONSTRUCTION_BEAMJET -1

/**
 * @enum jetreconstruction_JetAlgorithm
 * @brief Enumeration representing different jet algorithms used in
 * the JetReconstruction.
 */
typedef enum {
  JETRECONSTRUCTION_JETALGORITHM_ANTIKT = 0, /**< The Anti-Kt algorithm. */
  JETRECONSTRUCTION_JETALGORITHM_CA = 1, /**< The Cambridge/Aachen algorithm. */
  JETRECONSTRUCTION_JETALGORITHM_KT = 2, /**< The Inclusive-Kt algorithm. */
  JETRECONSTRUCTION_JETALGORITHM_GENKT =
      3, /**< The Generalised Kt algorithm (with arbitrary power). */
  JETRECONSTRUCTION_JETALGORITHM_EEKT =
      4, /**< The Generalised e+e- kt algorithm. */
  JETRECONSTRUCTION_JETALGORITHM_DURHAM =
      5 /**< The e+e- kt algorithm, aka Durham. */
} jetreconstruction_JetAlgorithm;

/**
 * @enum jetreconstruction_RecoStrategy
 * @brief Enumeration representing the different strategies for jet
 * reconstruction.
 */
typedef enum {
  JETRECONSTRUCTION_RECOSTRATEGY_BEST = 0,    /**< The best strategy. */
  JETRECONSTRUCTION_RECOSTRATEGY_N2PLAIN = 1, /**< The plain N2 strategy. */
  JETRECONSTRUCTION_RECOSTRATEGY_N2TILTED = 2 /**< The tiled N2 strategy. */
} jetreconstruction_RecoStrategy;

/**
 * @struct jetreconstruction_PseudoJet
 * @brief A struct representing a pseudojet, a four-momentum object used in
 * jet reconstruction algorithms. Additional information for the link back into
 * the history of the clustering is stored in the _cluster_hist_index field.
 * There is caching of the more expensive calculations for rapidity and
 * azimuthal angle.
 */
typedef struct {
  double px;                /**< The x-component of the momentum. */
  double py;                /**< The y-component of the momentum. */
  double pz;                /**< The z-component of the momentum. */
  double E;                 /**< The energy component of the momentum. */
  long _cluster_hist_index; /**< The index of the cluster history. */
  double _pt2;              /**< The squared transverse momentum. */
  double _inv_pt2;          /**< The inverse squared transverse momentum. */
  double _rap;              /**< The rapidity. */
  double _phi;              /**< The azimuthal angle. */
} jetreconstruction_PseudoJet;

/**
 * @brief Initializes a jetreconstruction_PseudoJet object with given momentum.
 *
 * @param[in] ptr Pointer to the jetreconstruction_PseudoJet object to be
 * initialized.
 * @param[in] px The x-component of the momentum.
 * @param[in] py The y-component of the momentum.
 * @param[in] pz The z-component of the momentum.
 * @param[in] E The energy component of the momentum.
 * @return An integer status code indicating the success or failure.
 */
int jetreconstruction_PseudoJet_init(jetreconstruction_PseudoJet *ptr,
                                     double px, double py, double pz, double E);

/**
 * @struct jetreconstruction_HistoryElement
 * @brief A struct holding a record of jet mergers and finalisations.
 */
typedef struct {
  long parent1; /**< Index in history where first parent of this jet was
                  created (@ref JETRECONSTRUCTION_NONEXISTENTPARENT if this jet is
                  an original    particle) */
  long parent2; /**< Index in history where second parent of this jet was
                  created (@ref JETRECONSTRUCTION_NONEXISTENTPARENT if this jet is
                  an original    particle); @ref JETRECONSTRUCTION_BEAMJET if this
                  history entry just labels the    fact that the jet has
                  recombined    with the beam */
  long child;   /**< Index in history where the current jet is recombined with
                                  another jet to form its child. It is Invalid
                  if this jet does   not further recombine. */
  long jetp_index; /**< Index in the jets vector where we will find the
                     PseudoJet object corresponding to this jet (i.e. the
                     jet created at this entry of the history). NB: if this
                     element of the history corresponds to a beam
                     recombination, then @p jetp_index = @ref JETRECONSTRUCTION_INVALID.
                   */
  double dij;      /**< The distance corresponding to the recombination at this
                      stage of the clustering. */
  double max_dij_so_far; /**< The largest recombination distance seen so far
                            in the clustering history */
} jetreconstruction_HistoryElement;

/**
 * @struct jetreconstruction_ClusterSequence
 * @brief A struct holding the full history of a jet clustering sequence,
 * including the final jets.
 */
typedef struct {
  jetreconstruction_JetAlgorithm
      algorithm; /**< The algorithm used for clustering. */
  double power;  /**< The power value used for the clustering algorithm */
  double R;      /**< The R parameter used for the clustering algorithm. */
  jetreconstruction_RecoStrategy
      strategy; /**< The strategy used for clustering. */
  jetreconstruction_PseudoJet
      *jets;           /**< A pointer to an array of jetreconstruction_PseudoJet
                          containing the actual jets in the cluster sequence. */
  size_t jets_length;  /**< The length of @ref jets array. */
  long n_initial_jets; /**< The initial number of particles used for exclusive
                         jets. */
  jetreconstruction_HistoryElement
      *history; /**< A pointer to an array of jetreconstruction_HistoryElement
                   containing the branching history of the cluster sequence.
                   Each stage in the history indicates where to look in the jets
                   vector to get the physical PseudoJet. */
  size_t history_length; /**< The length of @ref history array. */
  double Qtot;           /**< The total energy of the event. */
} jetreconstruction_ClusterSequence;

/** @private */
void jetreconstruction_ClusterSequence_free_members_(
    jetreconstruction_ClusterSequence *ptr);

/**
 * @brief Frees the members of a jetreconstruction_ClusterSequence structure.
 *
 * This function releases any resources held by the members of the given
 * jetreconstruction_ClusterSequence structure.
 *
 * @param[in,out] ptr Pointer to the jetreconstruction_ClusterSequence
 * structure.
 */
static inline void jetreconstruction_ClusterSequence_free_members(
    jetreconstruction_ClusterSequence *ptr) {
  jetreconstruction_ClusterSequence_free_members_(ptr);
  ptr->jets = NULL;
  ptr->jets_length = 0;
  ptr->history = NULL;
  ptr->history_length = 0;
}

/**
 * @brief Reconstructs jets from a collection of particles using a specified
 * algorithm and strategy.
 *
 * This function reconstructs jets from a collection of particles using a
 * specified algorithm and strategy.
 *
 * @note The memory allocated for members of @p result must be freed using
 * the @ref jetreconstruction_ClusterSequence_free_members function after the
 * @p result is no longer needed to prevent memory leaks.
 *
 * @param[in] particles Pointer to an array of pseudojet objects used for jet
 * reconstruction.
 * @param[in] particles_length The length of @p particles array.
 * @param[in] algorithm The algorithm to use for jet reconstruction.
 * @param[in] R The jet radius parameter.
 * @param[in] strategy The jet reconstruction strategy to use.
 * @param[out] result A pointer to which a cluster sequence containing the
 * reconstructed jets and the merging history will be stored.
 * @return An integer status code indicating the success or failure.
 */
int jetreconstruction_jet_reconstruct(
    const jetreconstruction_PseudoJet *particles, size_t particles_length,
    jetreconstruction_JetAlgorithm algorithm, double R,
    jetreconstruction_RecoStrategy strategy,
    jetreconstruction_ClusterSequence *result);

/**
 * @struct jetreconstruction_JetsResult
 * @brief A structure to hold the inclusive or exclusive jets.
 */
typedef struct {
  jetreconstruction_PseudoJet *data; /**< Pointer to an array of jetreconstruction_PseudoJet. */
  size_t length;                     /**< The length of the @ref data array. */
} jetreconstruction_JetsResult;

/** @private */
void jetreconstruction_JetsResult_free_members_(
    jetreconstruction_JetsResult *ptr);

/**
 * @brief Frees the members of a jetreconstruction_JetsResult structure.
 *
 * This function releases any resources held by the members of the given
 * jetreconstruction_JetsResult structure.
 *
 * @param[in,out] ptr A pointer to the jetreconstruction_JetsResult structure.
 */
static inline void
jetreconstruction_JetsResult_free_members(jetreconstruction_JetsResult *ptr) {
  jetreconstruction_JetsResult_free_members_(ptr);
  ptr->data = NULL;
  ptr->length = 0;
}

/**
 * @brief Return all exclusive jets of a jetreconstruction_ClusterSequence with
 * a cut on the maximum distance parameter.
 *
 * This function computes the exclusive jets from a given
 * jetreconstruction_ClusterSequence with a cut on the maximum distance
 * parameter.
 *
 * @note The memory allocated for members of @p result must be freed using
 * the @ref jetreconstruction_JetsResult_free_members function after the
 * @p result is no longer needed to prevent memory leaks.
 *
 * @param[in] clustersequence A pointer to the jetreconstruction_ClusterSequence
 * object containing the clustering history and jets.
 * @param[in] dcut The distance parameter used to define the exclusive jets.
 * @param[out] result A pointer to the jetreconstruction_JetsResult object where
 * the resulting jets will be stored.
 * @return An integer status code indicating the success or failure.
 */
int jetreconstruction_exclusive_jets_dcut(
    const jetreconstruction_ClusterSequence *clustersequence, double dcut,
    jetreconstruction_JetsResult *result);

/**
 * @brief Return all exclusive jets of a jetreconstruction_ClusterSequence with
 * a specific number of jets.
 *
 * This function computes a specific number of exclusive jets from a given
 * jetreconstruction_ClusterSequence.
 *
 * @note The memory allocated for members of @p result must be freed using
 * the @ref jetreconstruction_JetsResult_free_members function after the
 * @p result is no longer needed to prevent memory leaks.
 *
 * @param[in] clustersequence A pointer to the jetreconstruction_ClusterSequence
 * object containing the clustering history and jets.
 * @param[in] njets The number of exclusive jets to be calculated.
 * @param[out] result A pointer to the jetreconstruction_JetsResult object where
 * the resulting jets will be stored.
 * @return An integer status code indicating the success or failure.
 */
int jetreconstruction_exclusive_jets_njets(
    const jetreconstruction_ClusterSequence *clustersequence, size_t njets,
    jetreconstruction_JetsResult *result);

/**
 * @brief Return all inclusive jets of a jetreconstruction_ClusterSequence with
 * pt > @p ptmin.
 *
 * This function computes the inclusive jets from a given
 * jetreconstruction_ClusterSequence. It iterates over the clustering history
 * and checks the transverse momentum of each parent jet. If the transverse
 * momentum is greater than or equal to @p ptmin, the jet is added to results
 * structure.
 *
 * @note The memory allocated for members of @p result must be freed using
 * the @ref jetreconstruction_JetsResult_free_members function after the
 * @p result is no longer needed to prevent memory leaks.
 *
 * @param[in] clustersequence A pointer to the jetreconstruction_ClusterSequence
 * object containing the clustering history and jets.
 * @param[in] ptmin The minimum transverse momentum (pt) threshold for the
 *              inclusive jets.
 * @param[out] result A pointer to the jetreconstruction_JetsResult object where
 * the resulting jets will be stored.
 * @return An integer status code indicating the success or failure.
 */
int jetreconstruction_inclusive_jets(
    const jetreconstruction_ClusterSequence *clustersequence, double ptmin,
    jetreconstruction_JetsResult *result);
