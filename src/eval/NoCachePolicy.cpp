#include "NoCachePolicy.hpp"

#include <cassert>

NoCachePolicy::NoCachePolicy(
    const TreeAln&traln)
{}


auto                    NoCachePolicy::clone() const
    ->ArrayPolicy::UPtr
{
    return std::make_unique<NoCachePolicy>(*this);
}


void                    NoCachePolicy::
    accountForRejectionPolicy(
    TreeAln&                traln,
    const std::vector<bool>&partitions,
    const std::vector<nat>& invalidNodes,
    ArrayOrientation&       arrayOrient,
    ArrayReservoir&         res)
{
    for (nat i = 0; i < partitions.size(); ++i)
    {
        if (not partitions[i])
            continue;

        for (auto&elem : invalidNodes)
        {
            nat id = elem - traln.getNumberOfTaxa()  - 1;
            arrayOrient.setInvalid(i, id);
        }
    }
}
