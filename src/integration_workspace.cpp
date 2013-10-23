#include "integration_workspace.h"

namespace integration {
namespace internal {
void workspace::clear() {
  return points.clear();
}

workspace::point workspace::worst_point() const {
  return points.front();
}

void workspace::update(point el1, point el2) {
  drop_worst();
  if (el1.error > el2.error) {
    insert_forward(el1);
    insert_backward(el2);
  } else {
    insert_forward(el2);
    insert_backward(el1);
  }
}

double workspace::total_area() const {
  double tot = 0.0;
  for (points_const_iter p = points.begin(); p != points.end(); ++p)
    tot += p->area;
  return tot;
}

double workspace::total_error() const {
  double tot = 0.0;
  for (points_const_iter p = points.begin(); p != points.end(); ++p)
    tot += p->error;
  return tot;
}

void workspace::insert_forward(point el) {
  points.insert(search_forward(el.error), el);
}

void workspace::insert_backward(point el) {
  points.insert(search_backward(el.error), el);
}

void workspace::drop_worst() {
  points.pop_front();
}

// Search forward, from the beginning of the list, until we find an
// error estimate *smaller* than error.  Return the iterator pointing
// at that position, so that std::insert will insert a new value
// *before* this point.
workspace::points_iter workspace::search_forward(double error) {
  points_iter p = points.begin();
  while (p != points.end() && p->error > error) {
    ++p;
  }
  return p;
}

// Search backwards
workspace::points_iter workspace::search_backward(double error) {
  points_iter p = points.end();
  do {
    --p;
  } while (p != points.begin() && p->error < error);
  return ++p;
}
}
}
