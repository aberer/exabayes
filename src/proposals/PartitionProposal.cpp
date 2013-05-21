#include "PartitionProposal.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"

template<> double PartitionProposal<MultiplierProposal,RateHetParameter>::relativeWeight = 0.5;
template<> double PartitionProposal<SlidingProposal,RateHetParameter>::relativeWeight = 0 ;

template<> double PartitionProposal<DirichletProposal,FrequencyParameter>::relativeWeight = 0.5;
template<> double PartitionProposal<SlidingProposal,FrequencyParameter>::relativeWeight = 0.5 ;

template<> double PartitionProposal<DirichletProposal,RevMatParameter>::relativeWeight = 0.5;
template<> double PartitionProposal<SlidingProposal,RevMatParameter>::relativeWeight = 0.5 ;




