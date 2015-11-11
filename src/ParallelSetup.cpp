#include "ParallelSetup.hpp"

// just delegates 

#if HAVE_PLL == 0
#include "ParallelSetupExamlCode.hpp"
#else 

static void emptyFunToAvoidWarning()
{
}


#endif 
