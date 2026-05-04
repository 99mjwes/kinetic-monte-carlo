#ifndef KMC_H_
#define KMC_H_

#include <AMReX.H>
#include <random>

using namespace amrex;

template <class T>
void VectorPrint(const Vector<T>& thevector){
    int I = thevector.size();
    if (I == 0) {
        return;
        }
    Print() << thevector[0];
    for (int i = 1; i < I; i++) {
        Print() << ", " << thevector[i];
    }
    Print() << std::endl;
}

// template <class T>
// T VectorSum(const Vector<T>& thevector){
//     T sum = 0.0;
//     T c = 0.0; // Compensation
//     for (auto& val : thevector) { //  Kahan summation
//         T y = val - c;
//         T t = sum + y;
//         c = (t - sum) - y;
//         sum = t;
//     }
//     return sum;
// }
    
template <typename T>
T VectorSum(const std::vector<T>& vec, uint skip = 0) {
    T sum = 0;
    size_t i = 0;
    T y = 0;
    T c = 0;
    T t = 0;
    int n = vec.size() - static_cast<int>(skip);
    if (n <= 0) return sum;
    for (; i + 3 < n; i += 4) {//  Kahan summation with unrolled loop
        y = vec[i] + vec[i + 1] + vec[i + 2] + vec[i + 3] - c;
        t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    for (; i < n; i++) {
        y = vec[i];
        t = sum + y - c;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}


struct SimulationData {
    Real Volume;
    Real runtime;
    Real save_point;
    Real last_point;
    Real a0;
    std::uniform_real_distribution<Real> distribution;
    std::mt19937_64 generator;
    ULong iteration;
    ULong i_max;
};


Real power(ULong base, int exponent);

// Real randnum();

void Compute_Reaction_Schema( Vector<Real>& ReactionSchema, Real& SchemaSum, const Real Volume, const Vector<ULong>& ReactantQuantity, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates);

void FacilitateReaction(Vector<ULong>& ReactantQuantity, const Vector<int>& reaction);

void interpret_line(std::string Header, Vector<std::string>& content, char delimiter);

int LoadData(std::string filename, Vector<std::string>& ReactantNames,  Vector<ULong>& ReactantQuantity, Vector<Vector<int>>& Reactions, Vector<ParserExecutor <4>>& ReactionRateExecutors, Vector<Parser>& ReactionRateParsers);

void ReactionPrint(const Vector<int>& reaction, const Vector<std::string>& ReactantNames, Real ReactionRate);

void check_unique_reations(const Vector<Vector<int>>& Reactions);

void ReactionLoop (Vector<ULong>& ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, Vector<ULong>& ReactionTracker, SimulationData& simdata);

void SerialReactionLoop(Vector<Vector<ULong>>& ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, Vector<Vector<ULong>>& ReactionTracker, Vector<SimulationData>& simdata, ULong start, ULong end);

void ParallelReactionLoop(Vector<Vector<ULong>>& ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, Vector<Vector<ULong>>& ReactionTracker, Vector<SimulationData>& simdata, ULong num_workers);

#endif
