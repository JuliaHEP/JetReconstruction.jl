/*
 * SPDX-FileCopyrightText: 2022-2024 JetReconstruction.jl authors, CERN
 * SPDX-License-Identifier: MIT
 */

#include <stddef.h>

/*
Enumeration representing different jet algorithms used in
the JetReconstruction module.
*/
typedef enum {
  JETRECONSTRUCTION_JETALGORITHM_ANTIKT = 0, /* The Anti-Kt algorithm. */
  JETRECONSTRUCTION_JETALGORITHM_CA = 1, /* The Cambridge/Aachen algorithm. */
  JETRECONSTRUCTION_JETALGORITHM_KT = 2, /* The Inclusive-Kt algorithm. */
  JETRECONSTRUCTION_JETALGORITHM_GENKT =
      3, /* The Generalised Kt algorithm (with arbitrary power). */
  JETRECONSTRUCTION_JETALGORITHM_EEKT =
      4, /* The Generalised e+e- kt algorithm. */
  JETRECONSTRUCTION_JETALGORITHM_DURHAM =
      5 /* The e+e- kt algorithm, aka Durham. */
} jetreconstruction_JetAlgorithm;

/*
Scoped enumeration (using EnumX) representing the different strategies for jet
reconstruction.
*/
typedef enum {
  JETRECONSTRUCTION_RECOSTRATEGY_BEST = 0,    /* The best strategy. */
  JETRECONSTRUCTION_RECOSTRATEGY_N2PLAIN = 1, /* The plain N2 strategy. */
  JETRECONSTRUCTION_RECOSTRATEGY_N2TILTED = 2 /* The tiled N2 strategy. */
} jetreconstruction_RecoStrategy;

/* The `PseudoJet` struct represents a pseudojet, a four-momentum object used in
jet reconstruction algorithms. Additional information for the link back into the
history of the clustering is stored in the `_cluster_hist_index` field. There is
caching of the more expensive calculations for rapidity and azimuthal angle.*/
typedef struct {
  double px;                /* The x-component of the momentum. */
  double py;                /* The y-component of the momentum. */
  double pz;                /* The z-component of the momentum. */
  double E;                 /* The energy component of the momentum. */
  long _cluster_hist_index; /* The index of the cluster history. */
  double _pt2;              /* The squared transverse momentum. */
  double _inv_pt2;          /* The inverse squared transverse momentum. */
  double _rap;              /* The rapidity. */
  double _phi;              /* The azimuthal angle. */
} jetreconstruction_PseudoJet;

int jetreconstruction_PseudoJet_init(jetreconstruction_PseudoJet *ptr,
                                     double px, double py, double pz, double E);
/*
A struct holding a record of jet mergers and finalisations
*/
typedef struct {
  long parent1;    /* Index in history where first parent of this jet was
                     created (NonexistentParent if this jet is an original
                     particle) */
  long parent2;    /* Index in history where second parent of this jet was
                     created (NonexistentParent if this jet is an original
                     particle); BeamJet if this history entry just labels the
                     fact that the jet has recombined with the beam */
  long child;      /* Index in history where the current jet is recombined with
                                     another jet to form its child. It is Invalid
                     if this jet does   not further recombine. */
  long jetp_index; /* Index in the jets vector where we will find the
                     PseudoJet object corresponding to this jet (i.e. the
                     jet created at this entry of the history). NB: if this
                     element of the history corresponds to a beam
                     recombination, then `jetp_index=Invalid`. */
  double dij;      /* The distance corresponding to the recombination at this
                      stage     of the clustering. */
  double max_dij_so_far; /* The largest recombination distance seen so far
                            in the clustering history */
} jetreconstruction_HistoryElement;

/*
A struct holding the full history of a jet clustering sequence, including the
final jets.
*/
typedef struct {
  jetreconstruction_JetAlgorithm
      algorithm; /* The algorithm used for clustering. */
  double power;  /* The power value used for the clustering algorithm (note that
                               this value is always stored as a Float64 to be
                    type stable) */
  double R;      /* The R parameter used for the clustering algorithm. */
  jetreconstruction_RecoStrategy
      strategy; /* The strategy used for clustering. */
  jetreconstruction_PseudoJet
      *jets; /* The actual jets in the cluster sequence, which are of type `T
  <: FourMomentum`. */
  size_t jets_length;  /* Length of jets. */
  long n_initial_jets; /* The initial number of particles used for exclusive
                         jets. */
  jetreconstruction_HistoryElement
      *history;          /* The branching history of the cluster sequence. Each
                               stage in the history indicates where to look in the
                               jets vector to get the physical PseudoJet. */
  size_t history_length; /* Length of history. */
  double Qtot;           /* The total energy of the event. */
} jetreconstruction_ClusterSequence;

void jetreconstruction_ClusterSequence_free_members_(
    jetreconstruction_ClusterSequence *ptr);
static inline void jetreconstruction_ClusterSequence_free_members(
    jetreconstruction_ClusterSequence *ptr) {
  jetreconstruction_ClusterSequence_free_members_(ptr);
  ptr->jets = NULL;
  ptr->jets_length = 0;
  ptr->history = NULL;
  ptr->history_length = 0;
}

int jetreconstruction_jet_reconstruct(
    const jetreconstruction_PseudoJet *particles, size_t particles_length,
    jetreconstruction_JetAlgorithm algorithm, double R,
    jetreconstruction_RecoStrategy strategy,
    jetreconstruction_ClusterSequence *result);

typedef struct {
  jetreconstruction_PseudoJet *data;
  size_t length;
} jetreconstruction_JetsResult;

void jetreconstruction_JetsResult_free_members_(
    jetreconstruction_JetsResult *ptr);
static inline void
jetreconstruction_JetsResult_free_members(jetreconstruction_JetsResult *ptr) {
  jetreconstruction_JetsResult_free_members_(ptr);
  ptr->data = NULL;
  ptr->length = 0;
}

int jetreconstruction_exclusive_jets_dcut(
    const jetreconstruction_ClusterSequence *clustersequence, double dcut,
    jetreconstruction_JetsResult *result);

int jetreconstruction_exclusive_jets_njets(
    const jetreconstruction_ClusterSequence *clustersequence, size_t njets,
    jetreconstruction_JetsResult *result);

int jetreconstruction_inclusive_jets(
    const jetreconstruction_ClusterSequence *clustersequence, double ptmin,
    jetreconstruction_JetsResult *result);
