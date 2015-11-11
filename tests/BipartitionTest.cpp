#include "../src/Bipartition.hpp"

TEST(BipartitionTest, general)
{
  Bipartition bp;
  bp.reserve(23); 

  bp.set(3); 
  bp.set(18); 
  bp.set(5); 
  bp.set(1); 
  bp.set(13); 
  bp.set(22); 

  Bipartition bp2 ; 
  bp2.reserve(23); 

  for(nat i = 0; i < 23; ++i)
    {
      if(bp.isSet(i))
	bp2.set(i); 
    }

  ASSERT_TRUE(bp == bp2); 
}
