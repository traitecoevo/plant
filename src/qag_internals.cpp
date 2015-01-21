#include <tree2/qag_internals.h>
#include <tree2/util.h>

#include <algorithm> // min_element, max_element

namespace tree2 {
namespace quadrature {
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
  for (points_const_iter p = points.begin(); p != points.end(); ++p) {
    tot += p->area;
  }
  return tot;
}

double workspace::total_error() const {
  double tot = 0.0;
  for (points_const_iter p = points.begin(); p != points.end(); ++p) {
    tot += p->error;
  }
  return tot;
}

void workspace::push_back(point el) {
  points.push_back(el);
}

void workspace::insert_forward(point el) {
  points.insert(search_forward(el.error), el);
}

void workspace::insert_backward(point el) {
  points.insert(search_backward(el.error), el);
}

intervals_type workspace::get_intervals() const {
  std::vector<double> a, b;
  for (points_const_iter p = points.begin(); p != points.end(); ++p) {
    a.push_back(p->a);
    b.push_back(p->b);
  }
  intervals_type ret;
  ret.push_back(a);
  ret.push_back(b);
  return ret;
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

intervals_type rescale_intervals(intervals_type x,
                                 double min, double max) {
  std::vector<double> a = x[0], b = x[1];
  const double
    min_old = *std::min_element(a.begin(), a.end()),
    max_old = *std::max_element(b.begin(), b.end());

  util::rescale(a.begin(), a.end(), min_old, max_old, min, max);
  util::rescale(b.begin(), b.end(), min_old, max_old, min, max);

  intervals_type ret;
  ret.push_back(a);
  ret.push_back(b);
  return ret;
}

}
}
}
