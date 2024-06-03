
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <kmc.H>

#include <thread>
#include <random>

using namespace amrex;

void SerialReactionLoop(Vector<Vector<Vector<Real>>>& ResultTensor, Vector<long> ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, const Real Volume, Real save_point, const Real runtime, const Real save_step, unsigned int init_seed, long i_max, size_t start, size_t end) {

    for (size_t i = start; i < end; i++ ) {
        // if (i % modulator == 0) {Print() << "Iteration " << i << std::endl; };
        ReactionLoop (ResultTensor[i], ReactantQuantity, amu, a0, Reactions, ReactionRates, Volume, save_point, runtime, save_step, init_seed + i, i_max);
    }
}


void ParallelReactionLoop(Vector<Vector<Vector<Real>>>& ResultTensor, Vector<long> ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, const Real Volume, Real save_point, const Real runtime, const Real save_step, unsigned int init_seed, long i_max, size_t num_threads) {
    std::vector<std::thread> threads;
    size_t chunk_size = ResultTensor.size() / num_threads;
    size_t start = 0;
    size_t end =0;

    for (size_t i = 0; i < num_threads - 1; ++i) {
        end = start + chunk_size;
        threads.emplace_back(SerialReactionLoop, std::ref(ResultTensor), ReactantQuantity, amu, a0, std::cref(Reactions), std::cref(ReactionRates), Volume, save_point, runtime, save_step, init_seed, i_max, start, end);
        start = end;
    }

    // Last thread handles any remaining entries
    threads.emplace_back(SerialReactionLoop, std::ref(ResultTensor), ReactantQuantity, amu, a0, std::cref(Reactions), std::cref(ReactionRates), Volume, save_point, runtime, save_step, init_seed, i_max, start, size_t(ResultTensor.size()));

    // Join all threads
    for (auto& thread : threads) {
        thread.join();
    }
}



void ReactionLoop (Vector<Vector<Real>>& ResultMatrix, Vector<long> ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, const Real Volume, Real save_point, const Real runtime, const Real save_step, unsigned int init_seed, long i_max) {

    // Vector<std::string> ReactantNames;
    std::seed_seq seed{init_seed, init_seed, init_seed, init_seed};
    // std::size_t seed_2 = size_t(init_seed) << 32 + size_t(init_seed);
    std::mt19937_64 generator (seed);
    std::uniform_real_distribution<Real> distribution (0,1);

    // System contains N reacting spicies and M Reaction Paths
    int N = ReactantQuantity.size();
    // int M = Reactions.size();

    Vector<Real> ResultVector(N+2, 0);

    int ii = 0;
    int k = 0;
    Real t = 0;
    
    Real tau, r1, r2;
    int mu;

    // Main Reaction Loop
    while (t<runtime) {

        if (a0*runtime*100 < 1.0) {
            break;
        }

        ResultVector[0] = a0;
        ResultVector[1] = std::accumulate(ReactantQuantity.begin(), ReactantQuantity.end(), 0.0);
        std::copy(ReactantQuantity.begin(), ReactantQuantity.end(), ResultVector.begin()+2);

        // Generate random variables
        if (ii >= i_max) {Print() << "ERROR: Reached end of simulation loop with insufficient time!!" << std::endl; break;}
        r1 = distribution(generator);
        r2 = distribution(generator);

        // generating tau and mu
        tau = 1/a0 * std::log(1/r1);
        auto mu0 = std::upper_bound(amu.begin(), amu.end(), r2);
        mu = std::distance(amu.begin(), mu0);

        // Performing reaction
        FacilitateReaction(ReactantQuantity, Reactions[mu]);
        Compute_Reaction_Schema(amu, a0, Volume, ReactantQuantity, Reactions, ReactionRates);

        // Advancing time
        ii++;
        t += tau;

        while (t >= save_point && save_point <= runtime) {
            k++;
            Real lerpt = (save_point - t + tau) / tau; // notice that t > savepoint, hence the tau term in the enumerator.
            AMREX_ASSERT_WITH_MESSAGE(lerpt >= 0.0 && lerpt <= 1.0, "Interpolation factor is out of bounds!");

            ResultMatrix[k][0] = ResultVector[0] + lerpt * (a0 - ResultVector[0]);
            ResultMatrix[k][1] = ResultVector[1] + lerpt * (std::accumulate(ReactantQuantity.begin(), ReactantQuantity.end(), 0.0) - ResultVector[1]);
            for (int j = 0; j < N; j++) {
                ResultMatrix[k][j+2] = ResultVector[j+2] + lerpt * (ReactantQuantity[j] - ResultVector[j+2]);
            }

            save_point *= save_step; // account for large periods of inactivity
            
            
        }
    }
    // Post-simlation save
    // Print() << "Simulation completed successfully in " << ii << " steps!" << std::endl;
    while (runtime >= save_point) {
            k++;
            Real lerpt = (save_point - t + tau) / tau;
            AMREX_ASSERT_WITH_MESSAGE(lerpt >= 0.0 && lerpt <= 1.0, "Interpolation factor is out of bounds!");

            ResultMatrix[k][0] = a0;
            ResultMatrix[k][1] = std::accumulate(ReactantQuantity.begin(), ReactantQuantity.end(), 0.0);
            for (int j = 0; j < N; j++) {
                ResultMatrix[k][j+2] = ReactantQuantity[j];
            }

            save_point *= save_step; // account for large periods of inactivity
            
        }

    // Print() << "Simulation completed successfully in " << i << " steps!" << std::endl;
}
