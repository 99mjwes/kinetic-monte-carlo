// #include <AMReX.H>
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <kmc.H>

#include <random>


using namespace amrex;

Real power(long base, int exponent) {
    if (exponent == 0) {return 1.0;}      // x^0 = 1, special case 0^0 = 1!!
    if (base < exponent) {return 0.0;}    // x*(x-1)*(x-2)*...*0*... = 0
    Real p = 1.0;
    for (int i = 0; i < exponent; i++) {
        p *= Real(base - i);
    }
    return p;
}

// Real randnum ()
// {
//         static std::random_device device;
//         static std::seed_seq seed{device(), device(), device(), device()};
//         // static std::default_random_engine generator(seed);
//         static  std::mt19937_64 generator (seed);
//         std::uniform_real_distribution<Real> distribution (0,1);
//         return distribution(generator);

// }


void Compute_Reaction_Schema( Vector<Real>& ReactionSchema, Real& SchemaSum, const Real Volume, const Vector<long>& ReactantQuantity, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates)
{
    int N = ReactantQuantity.size();
    int M = Reactions.size();

    AMREX_ASSERT_WITH_MESSAGE(Reactions[0].size() == 2*N, "ERROR: Vector mismatch in Reactions");

    Vector<int> reaction(N);
    Vector<Real> BaseSchema(M);
    Vector<Real> CummulativeSchema(M);

    long Reactant_sum;
    Real tracker;
    Real ssinv;

   
    for (int i = 0; i < M; i++) {
        std::copy(Reactions[i].begin(), Reactions[i].end() - N, reaction.begin());
        tracker = 1;

        Reactant_sum = std::accumulate(reaction.begin(), reaction.end() , long(0));
        AMREX_ASSERT_WITH_MESSAGE(Reactant_sum == Reactant_sum, "ERROR: NaN Reactant sum"); 

        for (int j = 0; j < N; j++) {
            tracker *= power(ReactantQuantity[j], reaction[j]);
            AMREX_ASSERT_WITH_MESSAGE(tracker == tracker, "ERROR: NaN Tracker");
        }

        BaseSchema[i] = tracker * ReactionRates[i] * std::pow(Volume, 1 - Reactant_sum);
        AMREX_ASSERT_WITH_MESSAGE(BaseSchema[i] == BaseSchema[i], "ERROR: NaN Base Schema");
    }

    std::partial_sum(BaseSchema.begin(), BaseSchema.end(), CummulativeSchema.begin());
    SchemaSum = std::accumulate(BaseSchema.begin(), BaseSchema.end(), Real(0.0));

    ssinv = 1/SchemaSum; // Inverse Schema Sum

    for (int i = 0; i < M; i++) {
        ReactionSchema[i] = CummulativeSchema[i]*ssinv;
    }

    AMREX_ASSERT_WITH_MESSAGE(std::abs(ReactionSchema[M - 1] - Real(1.0)) < 1e-15, "Cummulative Density Function is not proper!!");
    ReactionSchema[M - 1] = 1.0; // Ensuring that the last element is exactly 1

   
}


void FacilitateReaction(Vector<long>& ReactantQuantity, const Vector<int>& reaction){
    int N = ReactantQuantity.size();

    AMREX_ASSERT_WITH_MESSAGE(reaction.size() == 2*N, "ERROR: Size mismatch in reaction");

    for (int i = 0; i < N; i++) {
        ReactantQuantity[i] += reaction[N+i] - reaction[i];
    }
}

void interpret_line(std::string Header, Vector<std::string>& content, char delimiter)
{
    long oldpos = 0;
    // char delimiter = ';';
    long newpos = Header.find(delimiter);

    content.clear( );
    
    while(newpos != long(-1))
    {
        std::string substring = Header.substr(oldpos,newpos - oldpos);
        oldpos = newpos + 1;
        content.push_back(substring);
        newpos = Header.find(delimiter, oldpos);
    }

    std::string substring = Header.substr(oldpos,newpos - oldpos);
    content.push_back(substring);

}

int LoadData(std::string filename, Vector<std::string>& ReactantNames, Vector<long>& ReactantQuantity, Vector<Vector<int>>& Reactions, Vector<ParserExecutor <4>>& ReactionRateExecutors, Vector<Parser>& ReactionRateParsers) {

    std::ifstream Data;
    char delimiter = '\t';
    Data.open(filename);
    if (!Data.is_open()) {Print() << "Error: Could not open file " << filename << std::endl; return 1;}
    std::string Header;
    std::getline(Data, Header);

    // Skipping initial lines - may be desciption!
    std::size_t first_B = 0;
    std::string substring = "";

    std::string subcheck = "A";
    subcheck += delimiter;

    while(substring != subcheck) {
        if (!std::getline(Data, Header)) { Print() << "Error: Could not read data in file " << filename << std::endl << "Exected: " << subcheck << std::endl; return 2;}
        first_B = Header.find("B");
        substring = Header.substr(0,first_B);
    }


    // Reactant Name phase
    int n_params = 4;
    Vector<std::string> Entries;
    interpret_line(Header, Entries, delimiter);
    if (Entries.size() <= n_params) { Print() << "Error: Simulation in " << filename << " contains no reactants!" << std::endl; return 3;}

    for (int i = n_params; i < Entries.size(); i++) {
        substring = Entries[i];
        if (substring.empty()) {break;}
        ReactantNames.push_back( substring );
    }

    int N = ReactantNames.size();

    // Reactant Quantity phase 
    std::getline(Data, Header);
    interpret_line(Header, Entries, delimiter);
    if (Entries.size() < N + n_params) {
        if (Entries.size() > n_params) {
        Print() << "Error: Simulation in " << filename << " lacks initial concentration of reactant " << ReactantNames[Entries.size() - n_params] << std::endl; return 4;
        } else {
            Print() << "Error: Simulation in " << filename << " lacks initial concentration data!" << std::endl; return 5;
        }
    }

    for (int i = n_params; i < n_params + N; i++) {
        ReactantQuantity.push_back( long( std::atof(Entries[i].c_str())) );
    }
   
    Vector<int> reaction (2*N);
    Real A, B, C;
    int R = 0;
    std::string rr_expr;
    std::string rr_base_expr = "A * (T/300)^B *exp(C/T)";

    while(std::getline(Data, Header))
    {
        interpret_line(Header, Entries, delimiter);
        if (Entries[0].empty()) {continue;}
        if (Entries.size() != n_params + 2*N) {Print() << "Error: Simulation in " << filename << " has invalid reaction path in R" << R + 1 << std::endl << "Expected: " << n_params + 2*N << " entries --- got: " << Entries.size() << std::endl; return 7;}

        // Interpretation of reaction rate expression
        A = std::atof(Entries[0].c_str());
        B = std::atof(Entries[1].c_str());
        C = std::atof(Entries[2].c_str());
        rr_expr = Entries[3];
        if (rr_expr.empty() || rr_expr == "default") {rr_expr = rr_base_expr;} else {Print() << "Notice: Custom reaction rate expression in R" << R + 1 << " - " << rr_expr << std::endl;}
        
        ReactionRateParsers.push_back(Parser());
        ReactionRateParsers[R].define(rr_expr);
        ReactionRateParsers[R].setConstant("A", A);
        ReactionRateParsers[R].setConstant("B", B);
        ReactionRateParsers[R].setConstant("C", C);
        ReactionRateParsers[R].registerVariables({"T","E","Te","ne"});

        ReactionRateExecutors.push_back(ReactionRateParsers[R].compile<4>());


        for (int i = 0; i < 2*N; i++) {
            reaction[i] = std::atoi(Entries[i+n_params].c_str());
        }

        Reactions.push_back(reaction);
        R++;
    }

    check_unique_reations(Reactions);

    Data.close();

    return 0;
}

void check_unique_reations(const Vector<Vector<int>>& Reactions) {
    int M = Reactions.size();
    for (int i = 0; i < M; i++) {
        for (int j = i + 1; j < M; j++) {
            if (Reactions[i] == Reactions[j]) {
                Print() << "Warning: Reactions " << i + 1 << " and " << j + 1 << " are identical!" << std::endl;
            }
        }
    }
}
