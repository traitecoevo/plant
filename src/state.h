// -*-c++-*-
#ifndef TREE_STATE_H_
#define TREE_STATE_H_

namespace util {

// Little utility for loading state in Plants, etc (used in
// interface.cpp).  This might more naturally belong in plant.cpp if
// it's only used there (but perhaps specialised just on Plant*), or
// in util.h if it's generally useful.

template <typename T>
void r_set_state(T* obj, typename T::state x) {
  util::check_length(x.size(), obj->state_size());
  obj->set_state(x.begin());
}

template <typename T>
typename T::state r_get_state(T* obj) {
  typename T::state ret(obj->state_size());
  obj->get_state(ret.begin());
  return ret;
}

}

#endif
