#include <limits>
#include <algorithm>
#include <thread>

#include "comm/threads/MessageQueue.hpp"


typedef uint8_t byte; 

#define ITER 10000


void doSomeWork(MessageQueue &q, int myId)
{
  for(int i = 0; i < ITER ; ++i)
    {
      auto gotMessage = false; 
      while(not gotMessage)
	{
	  auto msg = std::vector<byte>{}; 
	  std::tie(gotMessage, msg) = q.consume<byte>(myId); 
	  if(gotMessage)
	    {
	      auto res = std::vector<int>(reinterpret_cast<int*>(msg.data()), reinterpret_cast<int*>(msg.data() + msg.size() / sizeof(int)));

	      auto total = std::accumulate(begin(res), end(res), 0); 
	      ASSERT_EQ(10 * (i+1) , total); 
	    }
	}
    }
}


#define NUM_THREADS 4


TEST(MessageQueueTest, test)
{
  auto&& queue = MessageQueue(NUM_THREADS );   

  auto threads = std::vector<std::thread>{}; 
  for(int i = 0; i < NUM_THREADS ; ++i)
    {
      auto fun = std::bind(doSomeWork, std::ref(queue), i);
      threads.emplace_back(fun );
    }

  auto data = std::vector<int>{1 ,2 ,3 ,4 };

  for(int i = 0; i < ITER ; ++i)
    {
      auto dataNow = data; 
      for(auto &elem : dataNow)
	elem *= (i+1); 
      
      auto msg = std::vector<byte>(reinterpret_cast<byte*>(dataNow.data()), reinterpret_cast<byte*>(dataNow.data() + sizeof(int) * dataNow.size()) );

      auto readers = std::vector<int> (NUM_THREADS,1);
      queue.produce(msg, readers ); 
    }

  for(auto &t : threads)
    t.join();
}
