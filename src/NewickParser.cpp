#include "NewickParser.hpp"
#include "Topology.hpp"

template<> class NewickParser<Topology>;
template<> class NewickParser<BareTopology>; 
