/*
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
 * https://github.com/ttk592/spline/
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */


#ifndef _tk_spline_h
#define _tk_spline_h

/* Defining NDEBUG before cassert is included disables assert functions in 
 * tk_spline.cpp
 * (see http://en.cppreference.com/w/cpp/error/assert for details)
 * Assert functions are not allowed in R extensions because 
 * 
 */ 

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

namespace tk {

// band matrix solver
class band_matrix {
private:
   std::vector< std::vector<double> > m_upper;  // upper band
   std::vector< std::vector<double> > m_lower;  // lower band
public:
   band_matrix() {};                             // constructor
   band_matrix(int dim, int n_u, int n_l);       // constructor
   ~band_matrix() {};                            // destructor
   void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
   int dim() const;                             // matrix dimension
   int num_upper() const {
      return m_upper.size()-1;
   }
   int num_lower() const {
      return m_lower.size()-1;
   }
   // access operator
   double & operator () (int i, int j);            // write
   double   operator () (int i, int j) const;      // read
   // we can store an additional diogonal (in m_lower)
   double& saved_diag(int i);
   double  saved_diag(int i) const;
   void lu_decompose();
   std::vector<double> r_solve(const std::vector<double>& b) const;
   std::vector<double> l_solve(const std::vector<double>& b) const;
   std::vector<double> lu_solve(const std::vector<double>& b,
                                bool is_lu_decomposed=false);

};


// spline interpolation
class spline {
private:
   std::vector<double> m_x,m_y;           // x,y coordinates of points
   // interpolation parameters
   // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
   std::vector<double> m_a,m_b,m_c,m_d;
   // Fast O(1) index lookup for (near-)equidistant x-grids. When the knots are
   // evenly spaced the std::lower_bound() binary search in operator() can be
   // replaced by direct arithmetic, which profiling showed to be a large share
   // of runtime in the TF24 hydraulics. See traitecoevo/plant#435.
   bool   m_uniform = false;  // set in set_points() if the x-grid is equidistant
   double m_x0 = 0.0;         // first knot (m_x[0])
   double m_inv_dx = 0.0;     // 1 / mean knot spacing
public:
   void set_points(const std::vector<double>& x,
                   const std::vector<double>& y, bool cubic_spline=true);
   // Defined inline in the header (not tk_spline.cpp) so it can inline into the
   // hot assimilation quadrature loop. Without LTO an out-of-line definition
   // forces an un-inlinable cross-TU call at every quadrature point; this is the
   // same lever used for the FF16 competition helpers. (traitecoevo/plant#435)
   double operator() (double x) const;
};

inline double spline::operator() (double x) const {
   size_t n=m_x.size();
   // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
   int idx;
   if(m_uniform) {
      // O(1) index from the uniform spacing, then nudge by at most a step or
      // two so the result is bit-identical to std::lower_bound (covers knot
      // rounding from grid construction and the exact-knot edge case).
      // (traitecoevo/plant#435)
      idx=static_cast<int>((x-m_x0)*m_inv_dx);
      // clamp to [0, n-1]: idx == n-1 is the right-extrapolation case, where
      // lower_bound() also returns n-1 (so h = x - m_x[n-1]).
      const int last=static_cast<int>(n)-1;
      if(idx<0) idx=0;
      else if(idx>last) idx=last;
      while(idx>0 && m_x[idx]>=x) --idx;
      while(idx<last && m_x[idx+1]<x) ++idx;
   } else {
      // Non-uniform grid: seed idx with an O(1) proportional guess based on the
      // average spacing, then nudge to the exact segment. The nudge loops
      // converge to the same idx as std::lower_bound (m_x[idx] < x <= m_x[idx+1],
      // clamped to [0,last]), so this is bit-identical; on a smoothly graded
      // grid the guess lands within a step or two, replacing the O(log n) binary
      // search with O(1)-amortised stepping. (traitecoevo/plant#435)
      const int last=static_cast<int>(n)-1;
      idx=static_cast<int>((x-m_x[0])*(last/(m_x[last]-m_x[0])));
      if(idx<0) idx=0; else if(idx>last) idx=last;
      while(idx>0 && m_x[idx]>=x) --idx;
      while(idx<last && m_x[idx+1]<x) ++idx;
   }

   double h=x-m_x[idx];
   double interpol;
   if(x<m_x[0]) {
      // extrapolation to the left
      interpol=((m_b[0])*h + m_c[0])*h + m_y[0];
   } else if(x>m_x[n-1]) {
      // extrapolation to the right
      interpol=((m_b[n-1])*h + m_c[n-1])*h + m_y[n-1];
   } else {
      // interpolation
      interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
   }
   return interpol;
}

} // namespace tk

#endif /* _tk_spline_h */
