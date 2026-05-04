#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <kmc.H>

#include <chrono>
#include <execution>
#include <thread>
#include <random>

using namespace amrex;

void ReactionPrint(const Vector<int> &reaction, const Vector<std::string> &ReactantNames, Real ReactionRate)
{

    int M = reaction.size();
    if (2 * ReactantNames.size() != M)
    {
        Print() << "Invalid reaction!!" << std::endl;
        return;
    }

    std::stringstream input;
    std::stringstream output;

    bool in1 = false, out1 = false;
    int div = M / 2;

    for (int i = 0; i < div; i++)
    {
        if (reaction[i])
        {
            if (in1)
            {
                input << " +";
            }
            if (reaction[i] > 1)
            {
                input << " " << reaction[i];
            }
            input << " " << ReactantNames[i];
            in1 = true;
        }

        if (reaction[i + div])
        {
            if (out1)
            {
                output << " +";
            }
            if (reaction[i + div] > 1)
            {
                output << " " << reaction[i + div];
            }
            output << " " << ReactantNames[i];
            out1 = true;
        }
    }

    Print() << input.str() << " -->" << output.str() << "  (" << ReactionRate << ")" << std::endl;
}

int main(int argc, char *argv[])
{
    amrex::Initialize(argc, argv);
    // {

    // Constants
    bool Electrons_Are_Ideal_Gas = false;

    Real avogadros_number = 6.02214076e23;
    Real R = 8.314;

    Vector<std::string> ReactantNames;
    Vector<ULong> ReactantQuantity;
    Vector<Vector<int>> Reactions;
    Vector<Real> ReactionRates;
    Vector<ParserExecutor<4>> ReactionRateExecutors;
    Vector<Parser> ReactionRateParsers;

    // System contains N reacting spicies and M Reaction Paths
    // Each reaction path is a vector of size 2N
    // Each reaction path has a rate constant
    // Each reaction path has a reaction rate

    Real T = 300; // Kelvin, assume T = 300
    Real P = 1e5; // Pascal
    Real E = 0.0; // Townsend, assume E = 0
    Real runtime = 1e-2;
    Real first_save = 1e-8;
    Real ne = 0;
    Real V = 0;
    int n_saves = 50;
    int niter = 1;
    int modulator = 1;
    int num_workers = 1;
    int iseed = 0;
    Real i_maxx = 0;
    bool export_coeffs = true;
    Real n_amb = 0;

    std::string filename = "reaction.tsv";
    std::string savename = "results.csv";

    // Parameters
    ParmParse pp;
    pp.query("T", T);
    pp.query("P", P);
    pp.query("E", E);
    pp.query("Electrons_Are_Ideal_Gas", Electrons_Are_Ideal_Gas);
    pp.query("volume", V);
    pp.query("filename", filename);
    pp.query("savename", savename);
    pp.query("runtime", runtime);
    pp.query("first_save", first_save);
    pp.query("n_saves", n_saves);
    pp.query("max_steps", i_maxx);
    pp.query("n_iter", niter);
    pp.query("ne", ne);
    pp.query("modulator", modulator);
    pp.query("num_workers", num_workers);
    pp.query("seed", iseed);
    pp.query("export_coeffs", export_coeffs);
    pp.query("n_amb", n_amb);

    ULong i_max = static_cast<ULong>(i_maxx);
    ULong seed = static_cast<ULong>(iseed);

    Real log_save_step = std::exp((std::log(runtime) - std::log(first_save)) / (n_saves));
    AMREX_ASSERT_WITH_MESSAGE(log_save_step > 1, "Save step improperly sized!"); // Assert that test_save_step is larger than 1
    Print() << "Logarithmic save step: " << log_save_step << std::endl;
    Print() << "First save: " << first_save << std::endl;
    Print() << "Runtime: " << runtime << std::endl;
    Print() << "Number of saves: " << n_saves << std::endl;

    int errorlevel = LoadData(filename, ReactantNames, ReactantQuantity, Reactions, ReactionRateExecutors, ReactionRateParsers);
    Real Ee = E < 65 ? 2.31 * E : 3 * std::pow(E / 65, 2.6) / (1 + std::pow(E / 65, 2));
    Real TeEN = 2 * Ee / 3;
    Real Te = T + TeEN;
    Print() << "Electron Temperature: " << Te << std::endl;

    // Overwrite the reaction rates with the parsed values
    for (int i = 0; i < ReactionRateExecutors.size(); i++)
    {
        ReactionRates.push_back((ReactionRateExecutors[i](T, E, Te, ne * 1e-6)));
    }

    if (errorlevel != 0)
    {
        Print() << "Exit return code: " << errorlevel << std::endl;
        return errorlevel;
    }

    // Simulation definitions
    ULong n = VectorSum(ReactantQuantity);
    if (n == 0)
    {
        Print() << "Error: Total Reactant Quantity is zero!" << std::endl;
        return 6;
    }

    int N = ReactantNames.size();
    int M = Reactions.size();
    int e_pos = -1;

    for (int i = 0; i < N; i++)
    { // Check if initial electron quantity is defined
        if (ReactantNames[i] == "e" || ReactantNames[i] == "*e")
        {
            if (!Electrons_Are_Ideal_Gas)
            {
                n = n - ReactantQuantity[i];
            }
            e_pos = i;
            break;
        }
    }

    Print() << "Electron Position: " << e_pos << std::endl;

    Print() << "Reactants: ";
    VectorPrint(ReactantNames);
    Print() << "Quantity: ";
    Vector<ULong>(N, 0);
    VectorPrint(ReactantQuantity);

    for (int i = 0; i < Reactions.size(); i++)
    {
        Print() << "Reaction R" << i + 1 << ": ";
        ReactionPrint(Reactions[i], ReactantNames, ReactionRates[i]);
    }

    ne *= 1e-6; // Convert to cm^-3
    Real Volume;
    if (V == 0.0)
    {
        Volume = 1e6 * n * R * T / (P * avogadros_number); // assume ideal gas in cm^3
    }
    else
    {
        Volume = V;
    }
    Real n0 = 1e6 * n / Volume;
    Real ef = n0 * 1e-21 * 1e-5 * E; // Electric field in kV/cm
    ULong n_electron = static_cast<ULong>(ne * Volume);

    if ((e_pos > -1) && (ne != 0.0) && (ReactantQuantity[e_pos] == 0))
    {
        ReactantQuantity[e_pos] = n_electron;
    }
    else
    {
        n_electron = ReactantQuantity[e_pos];
    }

    Print() << "Loaded following Parameters:" << std::endl;
    Print() << "Pressure is " << P << "Pa" << std::endl;
    Print() << "Temperature is " << T << "K" << std::endl;
    Print() << "Simulation Volume: " << Volume << " cm^3" << std::endl;
    Print() << "Particle Density: " << n0 << " m^-3" << std::endl;
    Print() << "Electric Field Strength: " << E << " Td  (" << ef << " kV/cm)" << std::endl;
    Print() << "Inital Electron Count: " << n_electron << std::endl;
    Print() << "Quantity of reactants: " << (n + n_electron) << " (" << (n + n_electron) / avogadros_number << " moles)" << std::endl;
    if (i_max == 0)
    {
        i_max = static_cast<ULong>(Reactions.size()) * n;
        Print() << "Maximum iterations set to: " << i_max << std::endl;
    }

    // Export the coefficients to a file
    if (export_coeffs)
    {
        std::ofstream coeffs;
        coeffs.open("coeffs.csv");
        coeffs << "T,E,Te,ne,V,P,";
        for (int i = 0; i < ReactionRates.size(); i++)
        {
            coeffs << "k" << i << ",";
        }
        coeffs << "\n";
        coeffs << T << "," << E << "," << Te << "," << ne << "," << Volume << "," << P << ",";
        for (int i = 0; i < ReactionRates.size(); i++)
        {
            coeffs << ReactionRates[i] << ",";
        }
        coeffs << "\n";
        coeffs.close();
    }

    // Set thread count
    int maxnum_workers = std::thread::hardware_concurrency();
    if (num_workers < 0)
    {
        num_workers = maxnum_workers;
    }
    if (num_workers > maxnum_workers)
    {
        Print() << "Warning: Requested worker count exceeds hardware capabilities!" << std::endl;
        num_workers = maxnum_workers;
    }
    if (num_workers > niter)
    {
        num_workers = niter;
    }
    if (num_workers == 0)
    {
        num_workers = 1;
    }
    Print() << "Number of workers: " << num_workers << "/" << maxnum_workers << std::endl;

    // Allocate memory for the results
    Vector<Real> ResultVector(N + 2, 0.0); // Vector of size N + 2 to store results of each save point
    Vector<Real> OutputVector(N + 2, 0.0); // Vector of size n_saves + 2 to store averaged results of each save point
    Vector<Real> ErrorVector(N + 2, 0.0);  // Vector of size n_saves + 2 to store standard deviation of results of each save point

    // Vector<Vector<Real> > ResultMatrix(niter, ResultVector);     // Matrix of size n_saves + 2 x N + 2 to store results of each iteration
    Vector<Real> ReactionCountVector(M, 0.0);                                // Vector of size M to store reaction counts of each save point
    Vector<Real> ReactionErrorVector(M, 0.0);                                // Vector of size M to store reaction count standard deviations of each save point
    Vector<Vector<ULong>> ReactionCountMatrix(niter, Vector<ULong>(M, 0)); // Matrix of size n_saves + 2 x M to store averaged reaction counts of each save
    Vector<Vector<ULong>> ReactantMatrix(niter, ReactantQuantity);         // Tensor of size niter x N to track the quantity of each reactant at each save point in each iteration

    Vector<Real> amu(M);

    // Compute the reaction initial conditions
    Real a0 = 1.0;
    Compute_Reaction_Schema(amu, a0, Volume, ReactantQuantity, Reactions, ReactionRates);
    Print() << "a0 initial: " << a0 << std::endl;

    if (seed == 0)
    {
        std::random_device device;
        seed = device();
        Print() << "Seed: " << seed << std::endl;
    }
    Vector<SimulationData> simdata(niter);
    for (int i = 0; i < niter; i++)
    {
        simdata[i].Volume = Volume;
        simdata[i].runtime = first_save;
        simdata[i].save_point = 0.0;
        simdata[i].last_point = 0.0;
        simdata[i].a0 = a0;
        simdata[i].generator = std::mt19937_64(seed + i);
        simdata[i].distribution = std::uniform_real_distribution<Real>(0.0, 1.0);
        simdata[i].iteration = 0;
        simdata[i].i_max = i_max;
    }

    std::ofstream results;
    std::ofstream reaction_counts;
    results.open(savename);
    results << "i,t,a0,,N,";
    for (int k = 0; k < N; k++)
    {
        if (ReactantNames[k][0] != '*')
        {
            results << "," << ReactantNames[k] << ",";
        }
    }
    results << "\n";

    reaction_counts.open("reaction_counts.csv");
    for (int i = 0; i < M; i++)
    {
        reaction_counts << i << ",";
    }
    reaction_counts << "\n";


    Real inverse_iter_count = 1.0 / niter; // Scale factor for averaging across iterations
    Real ReactionMass = VectorSum(ReactantQuantity);
    if (n_amb > 0)
    {
        ReactionMass = n_amb * Volume * 1e-6; // Account for ambient species, converting from cm^-3 to m^-3
        Print() << "Accounting for ambient species, total particle count is " << ReactionMass << std::endl;
    }



    auto start_time = std::chrono::high_resolution_clock::now();
    Print() << "Starting KMC model..." << std::endl;
    // Run the simulation
    for (int i = 0; i < n_saves + 1; i++)
    {
        Print() << "Running iteration " << i + 1 << " from t = " << simdata[0].save_point << " to " << simdata[0].runtime << std::endl;
        ParallelReactionLoop(ReactantMatrix, amu, a0, Reactions, ReactionRates, ReactionCountMatrix, simdata, num_workers);

        // Scale results by the number of iterations.
        for (int j = 0; j < niter; j++)
        {
            ReactionMass = n_amb > 0 ? n_amb * Volume * 1e-6 : std::accumulate(ReactantMatrix[j].begin(), ReactantMatrix[j].end(), 0.0); // Recalculate total particle count for current iteration, accounting for ambient species if specified
            OutputVector[0] += simdata[j].a0 * inverse_iter_count;
            OutputVector[1] += ReactionMass * inverse_iter_count;
            for (int k = 2; k < N + 2; k++)
            {
                OutputVector[k] += ReactantMatrix[j][k] * inverse_iter_count;
            }
            for (int k = 0; k < M; k++)
            {
                ReactionCountVector[k] += ReactionCountMatrix[j][k] * inverse_iter_count;
            }
        }

        // Compute Standard Deviations
        for (int j = 0; j < niter; j++)
        {
            ReactionMass = n_amb > 0 ? n_amb * Volume * 1e-6 : std::accumulate(ReactantMatrix[j].begin(), ReactantMatrix[j].end(), 0.0);
            ErrorVector[0] += std::pow(simdata[j].a0 - OutputVector[0], 2) * inverse_iter_count;
            ErrorVector[1] += std::pow(ReactionMass - OutputVector[1], 2) * inverse_iter_count;
            for (int k = 2; k < N + 2; k++)
            {
                ErrorVector[k] += std::pow(ReactantMatrix[j][k] - OutputVector[k], 2) * inverse_iter_count;
            }
            for (int k = 0; k < M; k++)
            {
                ReactionErrorVector[k] += std::pow(ReactionCountMatrix[j][k] - ReactionCountVector[k], 2) * inverse_iter_count;
                ReactionCountMatrix[j][k] = 0; // reset for next iteration
            }
        }


        // Write results to file
        results << i << "," << simdata[0].runtime << "," << OutputVector[0] << "," << std::sqrt(ErrorVector[0]) << "," << OutputVector[1] << "," << std::sqrt(ErrorVector[1]);
        for (int j = 0; j < N; j++)
        {
            if (ReactantNames[j][0] != '*')
            {
                results << "," << OutputVector[j + 2] << "," << std::sqrt(ErrorVector[j + 2]);
            }
        }
        results << "\n";
        for (int j = 0; j < M; j++)
        {
            reaction_counts << ReactionCountVector[j] << "," << std::sqrt(ReactionErrorVector[j]) << ",";
        }
        reaction_counts << "\n";


        // Update simdata for next save point
        for (int j = 0; j < niter; j++)
        {
            simdata[j].runtime *= log_save_step;
        }

        // Reset output and error vectors for next iteration
        std::fill(OutputVector.begin(), OutputVector.end(), 0.0);
        std::fill(ErrorVector.begin(), ErrorVector.end(), 0.0);
        std::fill(ReactionCountVector.begin(), ReactionCountVector.end(), 0.0);
        std::fill(ReactionErrorVector.begin(), ReactionErrorVector.end(), 0.0);
    }

    // if (num_workers > niter || num_workers == 1) {
    //     SerialReactionLoop(ReactantMatrix, amu, a0, Reactions, ReactionRates, ReactionTracker, simdata, 0, niter);
    // } else {
    //     ParallelReactionLoop(ReactantMatrix, amu, a0, Reactions, ReactionRates, ReactionTracker, simdata, num_workers);
    // }

    // Print() << "Scaling results by " << invn << std::endl;
    // for (int i = 0; i < niter; i++) {
    //     for (int j = 0; j < n_saves + 2; j++) {

    //         for (int k = 0; k < N+2; k++) {
    //             OutputMatrix[j][k] += ResultTensor[i][j][k] * invn;
    //         }
    //     }
    // }

    // if (niter > 1) {
    //     Print() << "Computing standard deviations.. " << std::endl;
    //     for (int i = 0; i < niter; i++) {
    //         for (int j = 0; j < n_saves + 2; j++) {
    //             for (int k = 0; k < N+2; k++) {
    //                 ErrorMatrix[j][k] += std::pow(ResultTensor[i][j][k] - OutputMatrix[j][k], 2) /(niter * std::pow(OutputMatrix[j][k], 2));
    //             }
    //         }
    //     }
    // }

    // Print() << "Writing results to file" << std::endl;
    // Real stave;
    // for (int i = 0; i < n_saves + 2; i++) {
    //     stave = (i) ? first_save : 0.0;
    //     results << i << "," << stave << "," << OutputMatrix[i][0] << "," << std::sqrt(ErrorMatrix[i][0]) << "," << OutputMatrix[i][1] << "," << std::sqrt(ErrorMatrix[i][1]);
    //     // errors << std::sqrt(ErrorMatrix[i][0]) << "," << std::sqrt(ErrorMatrix[i][1]);

    //     for (int j = 0; j < N; j++) {
    //         if (ReactantNames[j][0] != '*') {
    //             results << "," << OutputMatrix[i][j+2] << "," << std::sqrt(ErrorMatrix[i][j+2]);
    //         }
    //     }
    //     if (i) {first_save *= log_save_step;} // account for large periods of inactivity
    //     results << "\n";
    //     // errors << "\n";
    // }

    results.close();
    reaction_counts.close();
    // errors.close();
    auto stop_time = std::chrono::high_resolution_clock::now();
    Print() << "Simulation runtime was " << std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count() << "ms" << std::endl;
    amrex::Finalize();

    // }

    return 0;
}
