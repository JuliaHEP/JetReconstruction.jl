#include <fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <iomanip>
#include <iostream>
#include <vector>

void write_momentum(const fastjet::PseudoJet &jet) {
    std::cout << '[' << jet.px() << ',' << jet.py() << ',' << jet.pz() << ',' << jet.E() << ']';
}

void write_particles(const std::vector<fastjet::PseudoJet> &particles) {
    std::cout << '[';
    bool first = true;
    for (const auto &particle : particles) {
        if (!first) std::cout << ',';
        first = false;
        write_momentum(particle);
    }
    std::cout << ']';
}

int main() {
    using namespace fastjet;
    ClusterSequence::set_fastjet_banner_stream(&std::cerr);
    std::cout << std::setprecision(17);
    std::vector<PseudoJet> particles;
    for (const auto &coordinates : std::vector<std::vector<double>>{
             {100.0, 0.2, 0.03}, {15.0, 0.25, 6.2},
             {80.0, -0.3, 3.0}, {8.0, 1.3, 1.5}}) {
        PseudoJet particle;
        particle.reset_PtYPhiM(coordinates[0], coordinates[1], coordinates[2]);
        particles.push_back(particle);
    }
    GhostedAreaSpec spec(1.5, 1, 0.08, 1.0, 0.1, 1e-100);
    spec.set_random_status({12345, 67890});
    std::vector<PseudoJet> ghosts;
    spec.add_ghosts(ghosts);
    std::cout << "{\"fastjet_version\":\"" << FASTJET_VERSION << "\",\"R\":0.4,\"ghost_area\":"
              << spec.actual_ghost_area() << ",\"particles\":";
    write_particles(particles);
    std::cout << ",\"ghosts\":";
    write_particles(ghosts);
    std::cout << ",\"results\":[";
    bool first_result = true;
    for (const auto algorithm : {antikt_algorithm, cambridge_algorithm, kt_algorithm, genkt_algorithm}) {
        for (const auto scheme : {E_scheme, pt_scheme, pt2_scheme}) {
            JetDefinition definition = algorithm == genkt_algorithm
                ? JetDefinition(algorithm, 0.4, 1.5, scheme, N2Plain)
                : JetDefinition(algorithm, 0.4, scheme, N2Plain);
            ClusterSequenceActiveAreaExplicitGhosts cs(particles, definition, ghosts, spec.actual_ghost_area());
            if (!first_result) std::cout << ',';
            first_result = false;
            std::cout << "{\"algorithm\":\""
                      << (algorithm == antikt_algorithm ? "AntiKt" : algorithm == cambridge_algorithm ? "CA"
                          : algorithm == kt_algorithm ? "Kt" : "GenKt")
                      << "\",\"scheme\":" << static_cast<int>(scheme) << ",\"jets\":[";
            bool first_jet = true;
            for (const auto &jet : sorted_by_pt(cs.inclusive_jets(1.0))) {
                if (!first_jet) std::cout << ',';
                first_jet = false;
                std::cout << "{\"momentum\":";
                write_momentum(jet);
                std::cout << ",\"area\":" << cs.area(jet) << ",\"area_4vector\":";
                write_momentum(cs.area_4vector(jet));
                std::cout << '}';
            }
            std::cout << "]}";
        }
    }
    std::cout << "]}\n";
}
