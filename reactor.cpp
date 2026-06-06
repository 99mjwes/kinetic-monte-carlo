#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <kmc.H>

#include <thread>
#include <random>

using namespace amrex;

void SerialReactionLoop(Vector<Vector<ULong>>& ReactantQuantity, const Vector<Vector<ActiveSpecies>>& SparseReactions, const Vector<Vector<ActiveSpecies>>& SparseStateChange, const Vector<Real>& ReactionRates, Vector<Vector<ULong>>& ReactionTracker, Vector<SimulationData>& simdata, ULong start, ULong end) {

    for (size_t i = start; i < end; i++ ) {
        ReactionLoop (ReactantQuantity[i], SparseReactions, SparseStateChange, ReactionRates, ReactionTracker[i], simdata[i]);
    }
}


void ParallelReactionLoop(Vector<Vector<ULong>>& ReactantQuantity, const Vector<Vector<ActiveSpecies>>& SparseReactions, const Vector<Vector<ActiveSpecies>>& SparseStateChange, const Vector<Real>& ReactionRates, Vector<Vector<ULong>>& ReactionTracker, Vector<SimulationData>& simdata, ULong num_workers) {
    std::vector<std::thread> workers;

    size_t mod = ReactantQuantity.size() / num_workers;
    size_t rem = ReactantQuantity.size() % num_workers;

    size_t start = 0;
    size_t end = 0;

    for (size_t i = 0; i < num_workers; ++i) {
        end = start + mod + (i < rem);
        workers.emplace_back(SerialReactionLoop, std::ref(ReactantQuantity), std::cref(SparseReactions), std::cref(SparseStateChange), std::cref(ReactionRates), std::ref(ReactionTracker), std::ref(simdata), start, end);
        // Print() << "Worker " << i << " will run sim " << start + 1 << " to " << end << std::endl;
        start = end;
    }

    // Join all workers
    for (auto& thread : workers) {
        thread.join();
    }
}



void ReactionLoop (Vector<ULong>& ReactantQuantity, const Vector<Vector<ActiveSpecies>>& SparseReactions, const Vector<Vector<ActiveSpecies>>& SparseStateChange, const Vector<Real>& ReactionRates, Vector<ULong>& ReactionTracker, SimulationData& simdata) {
    
    // Read simulation parameters from simdata
    std::mt19937_64 generator = simdata.generator;
    std::uniform_real_distribution<Real> distribution = simdata.distribution;
    Real t = simdata.save_point;
    Real Volume = simdata.Volume;
    Real runtime = simdata.runtime;
    ULong i_max = simdata.i_max;
    ULong iteration = simdata.iteration;

    // Precompute reaction schema for the initial state
    size_t M = SparseReactions.size();
    Vector<Real> amu(M);
    Vector<Real> BaseSchema(M);
    Vector<Real> CummulativeSchema(M);
    Real a0;
    Compute_Reaction_Schema(amu, a0, Volume, ReactantQuantity, SparseReactions, ReactionRates, BaseSchema, CummulativeSchema); // Reusing amu and cummulative schema vectors to save memory allocations

    
    Real tau, r1, r2;
    size_t mu;

    // Main Reaction Loop
    while (t < runtime) {

        if ((a0 * (runtime - t) * 100.0) < 1.0) {
            break;
        }

        // Generate random variables
        if (iteration >= i_max) {Print() << "ERROR: Reached end of simulation loop with insufficient time!!" << std::endl; break;}
        r1 = distribution(generator);
        r2 = distribution(generator);

        // generating tau and mu
        tau =  std::log(1/r1) / a0;
        if (t + tau > runtime) { break; } // Check if the next reaction would exceed the runtime, if so, exit the loop
        auto mu0 = std::upper_bound(amu.begin(), amu.end(), r2);
        mu = std::distance(amu.begin(), mu0);
        if (mu >= M) { mu = M - 1; }; // Handle edge case where r2 is very close to 1, which can cause mu to be out of bounds

        // Performing reaction
        FacilitateReaction(ReactantQuantity, SparseStateChange[mu]);
        ReactionTracker[mu]++; // Track the reaction count for this reaction path
        Compute_Reaction_Schema(amu, a0, Volume, ReactantQuantity, SparseReactions, ReactionRates, BaseSchema, CummulativeSchema); // Reusing amu and cummulative schema vectors to save memory allocations

        // Advancing time
        iteration++;
        t += tau;

        // while (t >= save_point && save_point <= runtime) {
        //     k++;
        //     Real lerpt = (save_point - t + tau) / tau; // notice that t > savepoint, hence the tau term in the enumerator.
        //     AMREX_ASSERT_WITH_MESSAGE(lerpt >= 0.0 && lerpt <= 1.0, "Interpolation factor is out of bounds!");

        //     ResultMatrix[k][0] = ResultVector[0] + lerpt * (a0 - ResultVector[0]);
        //     ResultMatrix[k][1] = ResultVector[1] + lerpt * (std::accumulate(ReactantQuantity.begin(), ReactantQuantity.end(), 0.0) - ResultVector[1]);
        //     for (int j = 0; j < N; j++) {
        //         ResultMatrix[k][j+2] = ResultVector[j+2] + lerpt * (ReactantQuantity[j] - ResultVector[j+2]);
        //     }

        //     save_point *= save_step; // account for large periods of inactivity
            
            
        // }
    }

    // Update simdata
    simdata.save_point = runtime;   // Set save point to runtime to indicate completion
    simdata.last_point = t;         // Record the actual final time reached in the simulation
    simdata.a0 = a0;                // Record the final a0 value
    simdata.iteration = iteration;  // Record the final iteration count

    // Post-simlation save
    // if (ii >= i_max) { // if the simulation reached the maximum number of iterations
    //     Print() << "    Runtime:        " << runtime << std::endl;
    //     Print() << "    Save Point:     " << save_point << std::endl;
    //     Print() << "    Final time:     " << t << std::endl;
    //     Print() << "    Final a0:       " << a0 << std::endl;
    //     Print() << "    Last reaction:  " << ii << std::endl;
    //     Print() << "    Last Time Step: " << tau << std::endl;
    //     return;
    // }
    // while (runtime >= save_point) {
    //         k++;
    //         ResultMatrix[k][0] = a0;
    //         ResultMatrix[k][1] = std::accumulate(ReactantQuantity.begin(), ReactantQuantity.end(), 0.0);
    //         for (int j = 0; j < N; j++) {
    //             ResultMatrix[k][j+2] = ReactantQuantity[j];
    //         }

    //         save_point *= save_step; // account for large periods of inactivity
            
    //     }
        // Print() << " Completed simulation with " << ii << " reactions" << std::endl;
}
