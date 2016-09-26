#ifndef _THREAD_DEFS_HPP
#define _THREAD_DEFS_HPP

#include <thread>

#define MY_TID std::this_thread::get_id()

using tid_t =  std::thread::id;

#endif
