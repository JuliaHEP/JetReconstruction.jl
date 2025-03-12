#include "JetReconstruction.h"
#ifdef JETRECONSTRUCTION_COMPILER_PACKAGECOMPILER
#include "julia_init.h"
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printPseudoJet(const jetreconstruction_PseudoJet *jet) {
  assert(jet != NULL);
  printf("PseudoJet(%f %f %f %f %ld %f %f %f %f)\n", jet->px, jet->py, jet->pz,
         jet->E, jet->_cluster_hist_index, jet->_pt2, jet->_inv_pt2, jet->_rap,
         jet->_phi);
}

void printHistoryElement(const jetreconstruction_HistoryElement *history) {
  assert(history != NULL);
  printf("HistoryElement(%ld %ld %ld %ld %lf %lf)\n", history->parent1,
         history->parent2, history->child, history->jetp_index, history->dij,
         history->max_dij_so_far);
}

void printClusterSequence(const jetreconstruction_ClusterSequence *sequence) {
  printf("Cluster Sequence Information:\n"
         "Algorithm: %d\n"
         "Power: %f\n"
         "R parameter: %f\n"
         "Strategy: %d\n"
         "Initial number of jets: %ld\n"
         "Total event energy (Qtot): %f\n",
         sequence->algorithm, sequence->power, sequence->R, sequence->strategy,
         sequence->n_initial_jets, sequence->Qtot);
  printf("Number of jets: %zu\n", sequence->jets_length);
  if (sequence->jets != NULL) {
    for (size_t i = 0; i < sequence->jets_length; i++) {
      printPseudoJet(sequence->jets + i);
    }
  }
  printf("History length: %zu\n", sequence->history_length);
  if (sequence->history != NULL) {
    for (size_t i = 0; i < sequence->history_length; i++) {
      printHistoryElement(sequence->history + i);
    }
  }
}

void printJetsResult(const jetreconstruction_JetsResult *results) {
  assert(results != NULL);
  for (size_t i = 0; i < results->length; ++i) {
    printPseudoJet(results->data + i);
  }
}

int main(int argc, char *argv[]) {
  clock_t start_time = clock();
  int sc = 0;
#ifdef JETRECONSTRUCTION_COMPILER_PACKAGECOMPILER
  init_julia(0, NULL);
#endif
  size_t len = 2;
  jetreconstruction_PseudoJet particles[2];
  sc = jetreconstruction_PseudoJet_init(&particles[0], 0.0, 1.0, 2.0, 3.0);
  assert(sc == JETRECONSTRUCTION_STATUSCODE_OK);
  sc = jetreconstruction_PseudoJet_init(&particles[1], 1.0, 2.0, 3.0, 4.0);
  assert(sc == JETRECONSTRUCTION_STATUSCODE_OK);

  jetreconstruction_JetAlgorithm algorithm = JETRECONSTRUCTION_JETALGORITHM_CA;
  double R = 3.0;
  jetreconstruction_RecoStrategy strategy = JETRECONSTRUCTION_RECOSTRATEGY_BEST;

  jetreconstruction_ClusterSequence cluster_seq;
  sc = jetreconstruction_jet_reconstruct(particles, len, algorithm, R, strategy,
                                         &cluster_seq);
  assert(sc == JETRECONSTRUCTION_STATUSCODE_OK);

  printClusterSequence(&cluster_seq);
  jetreconstruction_JetsResult result;
  sc = jetreconstruction_exclusive_jets_njets(&cluster_seq, 2, &result);
  assert(sc == JETRECONSTRUCTION_STATUSCODE_OK);
  printJetsResult(&result);

  jetreconstruction_JetsResult_free_members(&result);
  jetreconstruction_ClusterSequence_free_members(&cluster_seq);
#ifdef JETRECONSTRUCTION_COMPILER_PACKAGECOMPILER
  shutdown_julia(0);
#endif

  clock_t end_time = clock();
  double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  printf("Execution time: %f seconds\n", time_spent);
  return 0;
}
