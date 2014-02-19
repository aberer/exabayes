
template<typename T>
void MessageQueue::post(const std::vector<T> &message, int numRead, int tag) 
{
  auto &&lock = std::lock_guard<std::mutex>{_mtx} ;
  assert(_asyncMessages.find(tag) == _asyncMessages.end()); 

  auto&& elem = ConsumableMessage(message, numRead); 
  _asyncMessages.insert(std::make_pair(tag, elem)); 
}

template<typename T>
std::tuple<bool,std::vector<T>> MessageQueue::consumeOne(int tag)
{
  _mtx.lock();
  auto found = _asyncMessages.find(tag); 
  if(found ==  end(_asyncMessages))
    {
      _mtx.unlock();
      return make_tuple(false, std::vector<T>{});
    }
  else
    _mtx.unlock();

  auto result = std::vector<T>{}; 
  bool readsLeft = false ;  
  std::tie(result, readsLeft) = found->second.consume<T>();
  
  if(readsLeft == 0 )
    {
      auto &&lock = std::lock_guard<std::mutex>(_mtx); 
      _asyncMessages.erase(found); 
    }

  return make_tuple(true, result);
}




