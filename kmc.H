#ifndef KMC_H_
#define KMC_H_

#include <AMReX.H>

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

Real power(long base, int exponent);

// Real randnum();

void Compute_Reaction_Schema( Vector<Real>& ReactionSchema, Real& SchemaSum, const Real Volume, const Vector<long>& ReactantQuantity, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates);

void FacilitateReaction(Vector<long>& ReactantQuantity, const Vector<int>& reaction);

void interpret_line(std::string Header, Vector<std::string>& content, char delimiter);

int LoadData(std::string filename, Vector<std::string>& ReactantNames,  Vector<long>& ReactantQuantity, Vector<Vector<int>>& Reactions, Vector<ParserExecutor <4>>& ReactionRateExecutors, Vector<Parser>& ReactionRateParsers);

void ReactionPrint(const Vector<int>& reaction, const Vector<std::string>& ReactantNames, Real ReactionRate);

void check_unique_reations(const Vector<Vector<int>>& Reactions);

void ReactionLoop (Vector<Vector<Real>>& ResultMatrix, Vector<long> ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, const Real Volume, Real save_point, const Real runtime, const Real save_step, unsigned int init_seed, long i_max);

void SerialReactionLoop(Vector<Vector<Vector<Real>>>& ResultTensor, Vector<long> ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, const Real Volume, Real save_point, const Real runtime, const Real save_step, unsigned int init_seed, long i_max, size_t start, size_t end);

void ParallelReactionLoop(Vector<Vector<Vector<Real>>>& ResultTensor, Vector<long> ReactantQuantity, Vector<Real> amu, Real a0, const Vector<Vector<int>>& Reactions, const Vector<Real>& ReactionRates, const Real Volume, Real save_point, const Real runtime, const Real save_step, unsigned int init_seed, long i_max, size_t num_workers);

#endif
