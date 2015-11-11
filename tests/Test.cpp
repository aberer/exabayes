
TEST(general, bla)
{
  auto bla = std::vector<uint8_t>{1,2,3,4,5,6}; 
  
  auto cpy = std::vector<uint32_t>(6,0); 
  std::copy(begin(bla), end(bla), begin(cpy)); 

  for(int i = 0; i < 6; ++i)
    ASSERT_EQ(cpy[i], i+1); 
}
