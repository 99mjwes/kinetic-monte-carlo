#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <kmc.H>

#include <thread>
#include <random>

using namespace amrex;

void SerialReactionLoop(Vector<Vector<Vector<Real>>>& ResultTensor, Vector<unsigned long> ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, const Real Volume, Real save_point, const Real runtime, const Real save_step, unsigned int init_seed, unsigned long i_max, unsigned long start, unsigned long end) {

    for (unsigned long i = start; i < end; i++ ) {
        ReactionLoop (ResultTensor[i], ReactantQuantity, amu, a0, Reactions, ReactionRates, Volume, save_point, runtime, save_step, init_seed + i, i_max);
    }
}


void ParallelReactionLoop(Vector<Vector<Vector<Real>>>& ResultTensor, Vector<unsigned long> ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, const Real Volume, Real save_point, const Real runtime, const Real save_step, unsigned int init_seed, unsigned long i_max, unsigned long num_workers) {
    std::vector<std::thread> workers;

    unsigned long mod = ResultTensor.size() / num_workers;
    unsigned long rem = ResultTensor.size() % num_workers;

    unsigned long start = 0;
    unsigned long end = 0;

    for (unsigned long i = 0; i < num_workers; ++i) {
        end = start + mod + (i < rem);
        workers.emplace_back(SerialReactionLoop, std::ref(ResultTensor), ReactantQuantity, amu, a0, std::cref(Reactions), std::cref(ReactionRates), Volume, save_point, runtime, save_step, init_seed, i_max, start, end);
        Print() << "Worker " << i << " will run sim " << start + 1 << " to " << end << std::endl;
        start = end;
    }

    // Join all workers
    for (auto& thread : workers) {
        thread.join();
    }
}



void ReactionLoop (Vector<Vector<Real>>& ResultMatrix, Vector<unsigned long> ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, const Real Volume, Real save_point, const Real runtime, const Real save_step, unsigned int init_seed, unsigned long i_max) {
    
    // Random number generator
    std::seed_seq seed{init_seed, init_seed, init_seed, init_seed};
    std::mt19937_64 generator (seed);
    std::uniform_real_distribution<Real> distribution (0,1);

    // System contains N reacting spicies and M Reaction Paths
    int N = ReactantQuantity.size();
    Vector<Real> ResultVector(N+2, 0);

    // Initialize the Reaction
    Real ii = 0;
    int k = 0;
    Real t = 0;
    
    Real tau, r1, r2;
    int mu;

    // Main Reaction Loop
    while (t<runtime) {

        if ((a0 * runtime * 100.0) < 1.0) {
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
        tau =  std::log(1/r1) / a0;
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
    if (ii >= i_max) { // if the simulation reached the maximum number of iterations
        Print() << "    Runtime:        " << runtime << std::endl;
        Print() << "    Save Point:     " << save_point << std::endl;
        Print() << "    Final time:     " << t << std::endl;
        Print() << "    Final a0:       " << a0 << std::endl;
        Print() << "    Last reaction:  " << ii << std::endl;
        Print() << "    Last Time Step: " << tau << std::endl;
        return;
    }
    while (runtime >= save_point) {
            k++;
            ResultMatrix[k][0] = a0;
            ResultMatrix[k][1] = std::accumulate(ReactantQuantity.begin(), ReactantQuantity.end(), 0.0);
            for (int j = 0; j < N; j++) {
                ResultMatrix[k][j+2] = ReactantQuantity[j];
            }

            save_point *= save_step; // account for large periods of inactivity
            
        }
        // Print() << " Completed simulation with " << ii << " reactions" << std::endl;
}
