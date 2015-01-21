
#ifndef BITVECTOR_TEST_H
#define BITVECTOR_TEST_H


#include "bitvector.hpp"


TEST(bitvectortest, standard)
{
  auto b = bitvector{};
  b.set(3);

  ASSERT_TRUE(not b.get(0)); 
  ASSERT_TRUE(not b.get(1)); 
  ASSERT_TRUE(not b.get(2)); 
  ASSERT_TRUE(b.get(3));
  ASSERT_EQ(b.count(), 1);
  b.unset(3);
  ASSERT_EQ(b.count(), 0);
  b.set(31);
  ASSERT_EQ(b.count(), 1);
  b.set(100);
  ASSERT_EQ(b.count(), 2);
}


TEST(bitvectortest, friendOperators)
{
  auto b = bitvector{};
  b.set(1);
  b.set(60);
  ASSERT_EQ(b.count(), 2); 

  auto a = bitvector();
  a.set(0);
  a.set(60);
  a.set(100);
  ASSERT_EQ(a.count(), 3); 

  auto c = a | b ;
  ASSERT_EQ(c.count() , 4);

  // std::cout << "====" << std::endl;
  // std::cout << a << std::endl;
  // std::cout << b << std::endl;
  
  auto d = a & b;
  ASSERT_EQ(d.count(), 1); 
}



TEST(bitvectortest, smallerComp )
{
  auto a = bitvector{ 0, 1, 2 ,8};
  auto b = bitvector{ 31};
  auto c = (a & (~ b) );
  ASSERT_EQ(c.count() ,  4);
}



TEST(bitvectortest, counting)
{
  auto bv = bitvector();
  
  for(int i = 1; i < 60; ++i)
    {
      bv.set(i);
      ASSERT_EQ(bv.count(), i ); 
    }
}


TEST(bitvectortest, complement)
{
  auto a = bitvector{  0, 10, 63 , 12};
  auto b = bitvector{ 0, 1 ,10, 11, 50, 51, 63};
  
  ASSERT_EQ(a.count(), 4);
  ASSERT_EQ(b.count(), 7);

  ASSERT_EQ( (~b).count(), 57  );

  a.resize(b.size()); 
  
  auto c = (a & (~b));
  ASSERT_EQ(c.count(),1);
  
  c = ( ~ a & b ); 

  ASSERT_EQ( c.count(), 4  );

  c = (a | ~ b );
  auto d = ~(~ a &  b) ; 

  ASSERT_TRUE( c == d ); 
}


TEST(bitvectortest, plusMinus)
{
  auto a = bitvector{ 0,1,2};
  auto b = bitvector{ 2,3,4};

  ASSERT_EQ(  bitvector(a - b).count() , 2   ); 
}


TEST(bitvectortest, symDiff)
{
  auto a = bitvector{ 0,1,2};
  auto b = bitvector{ 2,3,4};
  auto c = a.symmetricDifference(b);
  ASSERT_EQ( c.count(), 4); 
}



#endif /* BITVECTOR_TEST_H */



