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


void ReactionPrint(const Vector<int>& reaction, const Vector<std::string>& ReactantNames, Real ReactionRate) {

    int M = reaction.size();
    if (2 * ReactantNames.size() != M) {Print() << "Invalid reaction!!" << std::endl; return;}
    
    std::stringstream input;
    std::stringstream output;

    bool in1 = false, out1  = false;
    int div = M/2;

    for (int i = 0; i < div; i++) {
        if(reaction[i]) {
            if (in1) {input << " +"; }
            if (reaction[i] > 1) {input << " " << reaction[i];}
            input << " " << ReactantNames[i];
            in1 = true;
        }

        if(reaction[i + div]) {
            if (out1) {output << " +"; }
            if (reaction[i + div] > 1) {output << " " << reaction[i + div];}
            output << " " << ReactantNames[i];
            out1 = true;
        }
    }

    Print() << input.str() << " -->" << output.str() << "  (" << ReactionRate << ")" << std::endl;

}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    // {

        // Constants
        bool Electrons_Are_Ideal_Gas = false;

        Real avogadros_number = 6.02214076e23;
        Real R = 8.314;

        Vector<std::string> ReactantNames;
        Vector<unsigned long> ReactantQuantity;
        Vector<Vector<int>> Reactions;
        Vector<Real> ReactionRates;
        Vector<ParserExecutor <4>> ReactionRateExecutors;
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
        int i_maxx = 0;
        bool export_coeffs = true;

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

        unsigned long i_max = static_cast<unsigned long>(i_maxx);
        unsigned long seed = static_cast<unsigned long>(iseed);



        Real log_save_step =  std::exp( (std::log(runtime) - std::log(first_save)) / (n_saves) );
        AMREX_ASSERT_WITH_MESSAGE(log_save_step > 1, "Save step improperly sized!"); // Assert that test_save_step is larger than 1
        Print() << "Logarithmic save step: " << log_save_step << std::endl;
        Print() << "First save: " << first_save << std::endl;
        Print() << "Runtime: " << runtime << std::endl;
        Print() << "Number of saves: " << n_saves << std::endl;

        int errorlevel = LoadData(filename,  ReactantNames,  ReactantQuantity,  Reactions, ReactionRateExecutors, ReactionRateParsers);
        Real Ee = E < 65 ? 2.31*E : 3* std::pow(E/65, 2.6) /(1 + std::pow(E/65,2));
        Real TeEN = 2 * Ee / 3;
        Real Te = T + TeEN;
        Print() << "Electron Temperature: " << Te << std::endl;

        // Overwrite the reaction rates with the parsed values
        for (int i = 0; i < ReactionRateExecutors.size(); i++) {
            ReactionRates.push_back((ReactionRateExecutors[i](T, E, Te, ne*1e-6)));
        }

        if (errorlevel != 0) {
            Print() << "Exit return code: " << errorlevel << std::endl;
            return errorlevel;
        }

        // Simulation definitions
        unsigned long n = VectorSum(ReactantQuantity);
        if (n == 0) {            
            Print() << "Error: Total Reactant Quantity is zero!" << std::endl;
            return 6;
        }

        int N = ReactantNames.size();
        int M = Reactions.size();
        int e_pos = -1;

        for(int i = 0; i < N; i++) { // Check if initial electron quantity is defined
            if (ReactantNames[i] == "e" || ReactantNames[i] == "*e") {
                if (!Electrons_Are_Ideal_Gas) {n = n - ReactantQuantity[i];} 
                e_pos = i;
                break;
            }
        }

        Print() << "Electron Position: " << e_pos << std::endl;

        Print() << "Reactants: ";
        VectorPrint(ReactantNames);
        Print() << "Quantity: ";
        VectorPrint(ReactantQuantity);

        for (int i = 0; i < Reactions.size(); i++) {
            Print() << "Reaction R" << i + 1 << ": ";
            ReactionPrint(Reactions[i], ReactantNames, ReactionRates[i]);
            }

        ne *= 1e-6; // Convert to cm^-3
        Real Volume;
        if (V == 0.0) {
            Volume = 1e6 * n * R * T / (P * avogadros_number); // assume ideal gas in cm^3
        } else {
            Volume = V;
        }
        Real n0 = 1e6*n/Volume;
        Real ef = n0 * 1e-21 * 1e-5 * E; // Electric field in kV/cm
        unsigned long n_electron = static_cast<unsigned long>(ne * Volume);

        if ((e_pos > -1) && (ne != 0.0) && (ReactantQuantity[e_pos] == 0)) {
            ReactantQuantity[e_pos] = n_electron;
        } else {
            n_electron = ReactantQuantity[e_pos];
        }

        Print() << "Loaded following Parameters:" << std::endl;
        Print() << "Pressure is " << P << "Pa" << std::endl;
        Print() << "Temperature is " << T << "K" << std::endl; 
        Print() << "Simulation Volume: " << Volume << " cm^3" << std::endl;
        Print() << "Particle Density: " << n0 << " m^-3" << std::endl;
        Print() << "Electric Field Strength: " << E << " Td  (" << ef << " kV/cm)" << std::endl;
        Print() << "Inital Electron Count: " << n_electron << std::endl;
        Print() << "Quantity of reactants: " << (n + n_electron) << " (" << (n + n_electron)/avogadros_number << " moles)" << std::endl;
        if (i_max == 0) {
            i_max = static_cast<unsigned long>(Reactions.size()) * n;
            Print() << "Maximum iterations set to: " << i_max << std::endl;
        } 
        

        // Export the coefficients to a file
        if (export_coeffs) {
            std::ofstream coeffs;
            coeffs.open ("coeffs.csv");
            coeffs << "T,E,Te,ne,V,P,";
            for (int i = 0; i < ReactionRates.size(); i++) {
                coeffs << "k" << i << ",";
            }
            coeffs << "\n";
            coeffs << T << "," << E << "," << Te << "," << ne << "," << Volume << "," << P << ",";
            for (int i = 0; i < ReactionRates.size(); i++) {
                coeffs << ReactionRates[i] << ",";
            }
            coeffs << "\n";
            coeffs.close();
        }

    
        // Set thread count
        int maxnum_workers = std::thread::hardware_concurrency();
        if (num_workers < 0) {num_workers = maxnum_workers;}
        if (num_workers > maxnum_workers) {
            Print() << "Warning: Requested worker count exceeds hardware capabilities!" << std::endl;
            num_workers = maxnum_workers;
            }
        if (num_workers > niter) {num_workers = niter;}
        if (num_workers == 0) {num_workers = 1;}
        Print() << "Number of workers: " << num_workers << "/" << maxnum_workers << std::endl;
        
        
        // Allocate memory for the results
        Vector<Real> ResultVector(N + 2, 0.0);                              // Vector of size N + 2 to store results of each save point
        Vector<Vector<Real>> OutputMatrix(n_saves + 2, ResultVector);       // Matrix of size n_saves + 2 x N + 2 to store results of each iteration
        Vector<Vector<Real>> ErrorMatrix(n_saves + 2, ResultVector);        // Matrix of size n_saves + 2 x N + 2 to store relative errors of each iteration
        Vector<Vector<Vector<Real>>> ResultTensor(niter, OutputMatrix);     // Tensor of size niter x n_saves + 2 x N + 2 to store and combine all results

        Vector<Real> amu(M);

        // Compute the reaction initial conditions
        Real a0 = 1.0;
        Compute_Reaction_Schema( amu, a0, Volume, ReactantQuantity, Reactions, ReactionRates);
        Print() << "a0 initial: " << a0 << std::endl;

        std::ofstream results;
        results.open (savename);
        results << "i,t,a0,,N,";
        for (int k = 0; k < N; k++) {
            if (ReactantNames[k][0] != '*') {
            results << "," << ReactantNames[k] << "," ;
            } 
        }
        results << "\n";
        
        auto start_time = std::chrono::high_resolution_clock::now();
        Print() << "Starting KMC model..." << std::endl;

        Real invn = 1.0 / niter; // Scale factor, the inverse of the number of iterations for average computation 
        Real ReactionMass = VectorSum(ReactantQuantity);

        if (seed == 0){        
            std::random_device device;
            seed = device();
            Print() << "Seed: " << seed << std::endl;}

        for (int i = 0; i < niter; i++) {

            // Preload the tensor with the initial conditions
            ResultTensor[i][0][0] = a0;
            ResultTensor[i][0][1] = ReactionMass;
            for (int j = 0; j < N; j++) {
                ResultTensor[i][0][j+2] = static_cast<Real>(ReactantQuantity[j]);
            }  
        };
        
        // Run the simulation
        if (num_workers > niter || num_workers == 1) {
            SerialReactionLoop(ResultTensor, ReactantQuantity, amu, a0, Reactions, ReactionRates, Volume, first_save, runtime, log_save_step, seed, i_max, 0, niter);
        } else {
            ParallelReactionLoop(ResultTensor, ReactantQuantity, amu, a0, Reactions, ReactionRates, Volume, first_save, runtime, log_save_step, seed, i_max, num_workers);
        }



        Print() << "Scaling results by " << invn << std::endl;
        for (int i = 0; i < niter; i++) {
            for (int j = 0; j < n_saves + 2; j++) {
                
                for (int k = 0; k < N+2; k++) {
                    OutputMatrix[j][k] += ResultTensor[i][j][k] * invn;
                }
            }
        }

        if (niter > 1) {
            Print() << "Computing standard deviations.. " << std::endl;
            for (int i = 0; i < niter; i++) {
                for (int j = 0; j < n_saves + 2; j++) {
                    for (int k = 0; k < N+2; k++) {
                        ErrorMatrix[j][k] += std::pow(ResultTensor[i][j][k] - OutputMatrix[j][k], 2) /(niter * std::pow(OutputMatrix[j][k], 2));
                    }
                }
            } 
        }

        Print() << "Writing results to file" << std::endl;
        Real stave; 
        for (int i = 0; i < n_saves + 2; i++) {
            stave = (i) ? first_save : 0.0;
            results << i << "," << stave << "," << OutputMatrix[i][0] << "," << std::sqrt(ErrorMatrix[i][0]) << "," << OutputMatrix[i][1] << "," << std::sqrt(ErrorMatrix[i][1]);
            // errors << std::sqrt(ErrorMatrix[i][0]) << "," << std::sqrt(ErrorMatrix[i][1]);
    
            for (int j = 0; j < N; j++) {
                if (ReactantNames[j][0] != '*') {
                    results << "," << OutputMatrix[i][j+2] << "," << std::sqrt(ErrorMatrix[i][j+2]);
                }
            }
            if (i) {first_save *= log_save_step;} // account for large periods of inactivity
            results << "\n";
            // errors << "\n";
        }

        results.close();
        // errors.close();
        auto stop_time = std::chrono::high_resolution_clock::now();
        Print() << "Simulation runtime was " << std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count() << "ms" << std::endl;
        amrex::Finalize();

    // }
    
    return 0;
}

