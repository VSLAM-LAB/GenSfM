/*
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014, 2016, 2021 Tino Kluge (ttk448 at gmail.com)
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


#ifndef TK_SPLINE_H
#define TK_SPLINE_H

#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ceres/jet.h>
#ifdef HAVE_SSTREAM
#include <sstream>
#include <string>
#endif // HAVE_SSTREAM

 // not ideal but disable unused-function warnings
 // (we get them because we have implementations in the header file,
 // and this is because we want to be able to quickly separate them
 // into a cpp file if necessary)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"

// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files
namespace
{

    namespace tk
    {



        // spline interpolation
        template <typename T>
        class spline
        {
        public:
            // spline types
            enum spline_type {
                linear = 10,            // linear interpolation
                cspline = 30,           // cubic splines (classical C^2)
                cspline_hermite = 31    // cubic hermite splines (local, only C^1)
            };

            // boundary condition type for the spline end-points
            enum bd_type {
                first_deriv = 1,
                second_deriv = 2
            };

        protected:
            std::vector<T> m_x, m_y;            // x,y coordinates of points
            // interpolation parameters
            // f(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
            // where a_i = y_i, or else it won't go through grid points
            std::vector<T> m_b, m_c, m_d;        // spline coefficients
            T m_c0;                            // for left extrapolation
            spline_type m_type;
            bd_type m_left, m_right;
            T  m_left_value, m_right_value;
            bool m_made_monotonic;
            void set_coeffs_from_b();               // calculate c_i, d_i from b_i
            size_t find_closest(T x) const;    // closest idx so that m_x[idx]<=x

        public:
            // default constructor: set boundary condition to be zero curvature
            // at both ends, i.e. natural splines
            spline() : m_type(cspline),
                m_left(second_deriv), m_right(second_deriv),
                m_left_value(T(0.0)), m_right_value(T(0.0)), m_made_monotonic(false)
            {
                ;
            }
            spline(const std::vector<T>& X, const std::vector<T>& Y,
                spline_type type = cspline,
                bool make_monotonic = false,
                bd_type left = second_deriv, T left_value = T(0),
                bd_type right = second_deriv, T right_value = T(0)
            ) :
                m_type(type),
                m_left(left), m_right(right),
                m_left_value(left_value), m_right_value(right_value),
                m_made_monotonic(false) // false correct here: make_monotonic() sets it
            {
                this->set_points(X, Y, m_type);
                if (make_monotonic) {
                    this->make_monotonic();
                }
            }


            // modify boundary conditions: if called it must be before set_points()
            void set_boundary(bd_type left, T left_value,
                bd_type right, T right_value);

            // set all data points (cubic_spline=false means linear interpolation)
            void set_points(const std::vector<T>& x,
                const std::vector<T>& y,
                spline_type type = cspline);

            // adjust coefficients so that the spline becomes piecewise monotonic
            // where possible
            //   this is done by adjusting slopes at grid points by a non-negative
            //   factor and this will break C^2
            //   this can also break boundary conditions if adjustments need to
            //   be made at the boundary points
            // returns false if no adjustments have been made, true otherwise
            bool make_monotonic();

            // evaluates the spline at point x
            T operator() (T x) const;
            T deriv(int order, T x) const;

            // returns the input data points
            std::vector<T> get_x() const { return m_x; }
            std::vector<T> get_y() const { return m_y; }
            T get_x_min() const { assert(!m_x.empty()); return m_x.front(); }
            T get_x_max() const { assert(!m_x.empty()); return m_x.back(); }
            T get_m_left_value() const { return m_left_value; }
            T get_m_right_value() const { return m_right_value; }

#ifdef HAVE_SSTREAM
            // spline info string, i.e. spline type, boundary conditions etc.
            std::string info() const;
#endif // HAVE_SSTREAM

        };



        namespace internal
        {

            // band matrix solver
            template <typename T>
            class band_matrix
            {
            private:
                std::vector< std::vector<T> > m_upper;  // upper band
                std::vector< std::vector<T> > m_lower;  // lower band
            public:
                band_matrix() {};                             // constructor
                band_matrix(int dim, int n_u, int n_l);       // constructor
                ~band_matrix() {};                            // destructor
                void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
                int dim() const;                             // matrix dimension
                int num_upper() const
                {
                    return (int)m_upper.size() - 1;
                }
                int num_lower() const
                {
                    return (int)m_lower.size() - 1;
                }
                // access operator
                T& operator () (int i, int j);            // write
                T   operator () (int i, int j) const;      // read
                // we can store an additional diagonal (in m_lower)
                T& saved_diag(int i);
                T  saved_diag(int i) const;
                void lu_decompose();
                std::vector<T> r_solve(const std::vector<T>& b) const;
                std::vector<T> l_solve(const std::vector<T>& b) const;
                std::vector<T> lu_solve(const std::vector<T>& b,
                    bool is_lu_decomposed = false);

            };

        } // namespace internal




        // ---------------------------------------------------------------------
        // implementation part, which could be separated into a cpp file
        // ---------------------------------------------------------------------

        // spline implementation
        // -----------------------
        template <typename T>
        void spline<T>::set_boundary(spline<T>::bd_type left, T left_value,
            spline<T>::bd_type right, T right_value)
        {
            assert(m_x.size() == 0);          // set_points() must not have happened yet
            m_left = left;
            m_right = right;
            m_left_value = left_value;
            m_right_value = right_value;
        }

        template <typename T>
        void spline<T>::set_coeffs_from_b()
        {
            assert(m_x.size() == m_y.size());
            assert(m_x.size() == m_b.size());
            assert(m_x.size() > 2);
            size_t n = m_b.size();
            if (m_c.size() != n)
                m_c.resize(n);
            if (m_d.size() != n)
                m_d.resize(n);

            for (size_t i = 0; i < n - 1; i++) {
                const T h = m_x[i + 1] - m_x[i];
                // from continuity and differentiability condition
                m_c[i] = (T(3.0) * (m_y[i + 1] - m_y[i]) / h - (T(2.0) * m_b[i] + m_b[i + 1])) / h;
                // from differentiability condition
                m_d[i] = ((m_b[i + 1] - m_b[i]) / (T(3.0) * h) - T(2.0) / T(3.0) * m_c[i]) / h;
            }

            // for left extrapolation coefficients
            m_c0 = (m_left == first_deriv) ? T(0) : m_c[0];
        }
        template <typename T>
        void spline<T>::set_points(const std::vector<T>& x,
            const std::vector<T>& y,
            spline_type type)
        {
            assert(x.size() == y.size());
            assert(x.size() > 2);
            m_type = type;
            m_made_monotonic = false;
            m_x = x;
            m_y = y;
            int n = (int)x.size();
            // check strict monotonicity of input vector x
            for (int i = 0; i < n - 1; i++) {
                assert(m_x[i] < m_x[i + 1]);
            }


            if (type == spline<T>::linear) {
                // linear interpolation
                m_d.resize(n);
                m_c.resize(n);
                m_b.resize(n);
                for (int i = 0; i < n - 1; i++) {
                    m_d[i] = T(0);
                    m_c[i] = T(0);
                    m_b[i] = (m_y[i + 1] - m_y[i]) / (m_x[i + 1] - m_x[i]);
                }
                // ignore boundary conditions, set slope equal to the last segment
                m_b[n - 1] = m_b[n - 2];
                m_c[n - 1] = T(0);
                m_d[n - 1] = T(0);
            }
            else if (type == spline<T>::cspline) {
                // classical cubic splines which are C^2 (twice cont differentiable)
                // this requires solving an equation system

                // setting up the matrix and right hand side of the equation system
                // for the parameters b[]
                internal::band_matrix<T> A(n, 1, 1);
                std::vector<T>  rhs(n);
                for (int i = 1; i < n - 1; i++) {
                    A(i, i - 1) = T(1.0) / T(3.0) * (x[i] - x[i - 1]);
                    A(i, i) = T(2.0) / T(3.0) * (x[i + 1] - x[i - 1]);
                    A(i, i + 1) = T(1.0) / T(3.0) * (x[i + 1] - x[i]);
                    rhs[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
                }
                // boundary conditions
                if (m_left == spline<T>::second_deriv) {
                    // 2*c[0] = f''
                    A(0, 0) = T(2.0);
                    A(0, 1) = T(0.0);
                    rhs[0] = m_left_value;
                }
                else if (m_left == spline<T>::first_deriv) {
                    // b[0] = f', needs to be re-expressed in terms of c:
                    // (2c[0]+c[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
                    A(0, 0) = T(2.0) * (x[1] - x[0]);
                    A(0, 1) = T(1.0) * (x[1] - x[0]);
                    rhs[0] = T(3.0) * ((y[1] - y[0]) / (x[1] - x[0]) - m_left_value);
                }
                else {
                    assert(false);
                }
                if (m_right == spline<T>::second_deriv) {
                    // 2*c[n-1] = f''
                    A(n - 1, n - 1) = T(2.0);
                    A(n - 1, n - 2) = T(0.0);
                    rhs[n - 1] = m_right_value;
                }
                else if (m_right == spline<T>::first_deriv) {
                    // b[n-1] = f', needs to be re-expressed in terms of c:
                    // (c[n-2]+2c[n-1])(x[n-1]-x[n-2])
                    // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
                    A(n - 1, n - 1) = T(2.0) * (x[n - 1] - x[n - 2]);
                    A(n - 1, n - 2) = T(1.0) * (x[n - 1] - x[n - 2]);
                    rhs[n - 1] = T(3.0) * (m_right_value - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
                }
                else {
                    assert(false);
                }

                // solve the equation system to obtain the parameters c[]
                m_c = A.lu_solve(rhs);

                // calculate parameters b[] and d[] based on c[]
                m_d.resize(n);
                m_b.resize(n);
                for (int i = 0; i < n - 1; i++) {
                    m_d[i] = T(1.0) / T(3.0) * (m_c[i + 1] - m_c[i]) / (x[i + 1] - x[i]);
                    m_b[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
                        - T(1.0) / T(3.0) * (T(2.0) * m_c[i] + m_c[i + 1]) * (x[i + 1] - x[i]);
                }
                // for the right extrapolation coefficients (zero cubic term)
                // f_{n-1}(x) = y_{n-1} + b*(x-x_{n-1}) + c*(x-x_{n-1})^2
                T h = x[n - 1] - x[n - 2];
                // m_c[n-1] is determined by the boundary condition
                m_d[n - 1] = T(0.0);
                m_b[n - 1] = T(3.0) * m_d[n - 2] * h * h + T(2.0) * m_c[n - 2] * h + m_b[n - 2];   // = f'_{n-2}(x_{n-1})
                if (m_right == first_deriv)
                    m_c[n - 1] = T(0.0);   // force linear extrapolation

            }
            else if (type == spline<T>::cspline_hermite) {
                // hermite cubic splines which are C^1 (cont. differentiable)
                // and derivatives are specified on each grid point
                // (here we use 3-point finite differences)
                m_b.resize(n);
                m_c.resize(n);
                m_d.resize(n);
                // set b to match 1st order derivative finite difference
                for (int i = 1; i < n - 1; i++) {
                    const T h = m_x[i + 1] - m_x[i];
                    const T hl = m_x[i] - m_x[i - 1];
                    m_b[i] = -h / (hl * (hl + h)) * m_y[i - 1] + (h - hl) / (hl * h) * m_y[i]
                        + hl / (h * (hl + h)) * m_y[i + 1];
                }
                // boundary conditions determine b[0] and b[n-1]
                if (m_left == first_deriv) {
                    m_b[0] = m_left_value;
                }
                else if (m_left == second_deriv) {
                    const T h = m_x[1] - m_x[0];
                    m_b[0] = T(0.5) * (-m_b[1] - T(0.5) * m_left_value * h + T(3.0) * (m_y[1] - m_y[0]) / h);
                }
                else {
                    assert(false);
                }
                if (m_right == first_deriv) {
                    m_b[n - 1] = m_right_value;
                    m_c[n - 1] = T(0.0);
                }
                else if (m_right == second_deriv) {
                    const T h = m_x[n - 1] - m_x[n - 2];
                    m_b[n - 1] = T(0.5) * (-m_b[n - 2] + T(0.5) * m_right_value * h + T(3.0) * (m_y[n - 1] - m_y[n - 2]) / h);
                    m_c[n - 1] = T(0.5) * m_right_value;
                }
                else {
                    assert(false);
                }
                m_d[n - 1] = T(0.0);

                // parameters c and d are determined by continuity and differentiability
                set_coeffs_from_b();

            }
            else {
                assert(false);
            }

            // for left extrapolation coefficients
            m_c0 = (m_left == spline<T>::first_deriv) ? T(0.0) : m_c[0];
        }

        template <typename T>
        bool spline<T>::make_monotonic()
        {
            assert(m_x.size() == m_y.size());
            assert(m_x.size() == m_b.size());
            assert(m_x.size() > 2);
            bool modified = false;
            const int n = (int)m_x.size();
            // make sure: input data monotonic increasing --> b_i>=0
            //            input data monotonic decreasing --> b_i<=0
            for (int i = 0; i < n; i++) {
                int im1 = std::max(i - 1, 0);
                int ip1 = std::min(i + 1, n - 1);
                if (((m_y[im1] <= m_y[i]) && (m_y[i] <= m_y[ip1]) && m_b[i] < T(0.0)) ||
                    ((m_y[im1] >= m_y[i]) && (m_y[i] >= m_y[ip1]) && m_b[i] > T(0.0))) {
                    modified = true;
                    m_b[i] = T(0.0);
                }
            }
            // if input data is monotonic (b[i], b[i+1], avg have all the same sign)
            // ensure a sufficient criteria for monotonicity is satisfied:
            //     sqrt(b[i]^2+b[i+1]^2) <= 3 |avg|, with avg=(y[i+1]-y[i])/h,
            for (int i = 0; i < n - 1; i++) {
                T h = m_x[i + 1] - m_x[i];
                T avg = (m_y[i + 1] - m_y[i]) / h;
                if (avg == T(0.0) && (m_b[i] != T(0.0) || m_b[i + 1] != T(0.0))) {
                    modified = true;
                    m_b[i] = T(0.0);
                    m_b[i + 1] = T(0.0);
                }
                else if ((m_b[i] >= T(0.0) && m_b[i + 1] >= T(0.0) && avg > T(0.0)) ||
                    (m_b[i] <= T(0.0) && m_b[i + 1] <= T(0.0) && avg < T(0.0))) {
                    // input data is monotonic
                    T r = ceres::sqrt(m_b[i] * m_b[i] + m_b[i + 1] * m_b[i + 1]) / T(ceres::abs(avg));
                    if (r > T(3.0)) {
                        // sufficient criteria for monotonicity: r<=3
                        // adjust b[i] and b[i+1]
                        modified = true;
                        m_b[i] *= (T(3.0) / r);
                        m_b[i + 1] *= (T(3.0) / r);
                    }
                }
            }

            if (modified == true) {
                set_coeffs_from_b();
                m_made_monotonic = true;
            }

            return modified;
        }

        template <typename T>
        // return the closest idx so that m_x[idx] <= x (return 0 if x<m_x[0])
        size_t spline<T>::find_closest(T x) const
        {
            typename std::vector<T>::const_iterator it;
            it = std::upper_bound(m_x.begin(), m_x.end(), x);       // *it > x
            size_t idx = std::max<int>((it - m_x.begin()) - 1, 0);    // m_x[idx] <= x
            return idx;
        }

        template <typename T>

        T spline<T>::operator() (T x) const
        {
            // polynomial evaluation using Horner's scheme
            // TODO: consider more numerically accurate algorithms, e.g.:
            //   - Clenshaw
            //   - Even-Odd method by A.C.R. Newbery
            //   - Compensated Horner Scheme
            size_t n = m_x.size();
            size_t idx = find_closest(x);

            T h = x - m_x[idx];
            T interpol;
            if (x < m_x[0]) {
                // extrapolation to the left
                interpol = (m_c0 * h + m_b[0]) * h + m_y[0];
            }
            else if (x > m_x[n - 1]) {
                // extrapolation to the right
                interpol = (m_c[n - 1] * h + m_b[n - 1]) * h + m_y[n - 1];
            }
            else {
                // interpolation
                interpol = ((m_d[idx] * h + m_c[idx]) * h + m_b[idx]) * h + m_y[idx];
            }
            return interpol;
        }

        template <typename T>
        T spline<T>::deriv(int order, T x) const
        {
            assert(order > 0);
            size_t n = m_x.size();
            size_t idx = find_closest(x);

            T h = x - m_x[idx];
            T interpol;
            if (x < m_x[0]) {
                // extrapolation to the left
                switch (order) {
                case 1:
                    interpol = T(2.0) * m_c0 * h + m_b[0];
                    break;
                case 2:
                    interpol = T(2.0) * m_c0;
                    break;
                default:
                    interpol = T(0.0);
                    break;
                }
            }
            else if (x > m_x[n - 1]) {
                // extrapolation to the right
                switch (order) {
                case 1:
                    interpol = T(2.0) * m_c[n - 1] * h + m_b[n - 1];
                    break;
                case 2:
                    interpol = T(2.0) * m_c[n - 1];
                    break;
                default:
                    interpol = T(0.0);
                    break;
                }
            }
            else {
                // interpolation
                switch (order) {
                case 1:
                    interpol = (T(3.0) * m_d[idx] * h + T(2.0) * m_c[idx]) * h + m_b[idx];
                    break;
                case 2:
                    interpol = T(6.0) * m_d[idx] * h + T(2.0) * m_c[idx];
                    break;
                case 3:
                    interpol = T(6.0) * m_d[idx];
                    break;
                default:
                    interpol = T(0.0);
                    break;
                }
            }
            return interpol;
        }

#ifdef HAVE_SSTREAM
        std::string spline::info() const
        {
            std::stringstream ss;
            ss << "type " << m_type << ", left boundary deriv " << m_left << " = ";
            ss << m_left_value << ", right boundary deriv " << m_right << " = ";
            ss << m_right_value << std::endl;
            if (m_made_monotonic) {
                ss << "(spline has been adjusted for piece-wise monotonicity)";
            }
            return ss.str();
        }
#endif // HAVE_SSTREAM


        namespace internal
        {

            // band_matrix implementation
            // -------------------------
            template <typename T>
            band_matrix<T>::band_matrix(int dim, int n_u, int n_l)
            {
                resize(dim, n_u, n_l);
            }
            template <typename T>
            void band_matrix<T>::resize(int dim, int n_u, int n_l)
            {
                assert(dim > 0);
                assert(n_u >= 0);
                assert(n_l >= 0);
                m_upper.resize(n_u + 1);
                m_lower.resize(n_l + 1);
                for (size_t i = 0; i < m_upper.size(); i++) {
                    m_upper[i].resize(dim);
                }
                for (size_t i = 0; i < m_lower.size(); i++) {
                    m_lower[i].resize(dim);
                }
            }
            template <typename T>
            int band_matrix<T>::dim() const
            {
                if (m_upper.size() > 0) {
                    return m_upper[0].size();
                }
                else {
                    return 0;
                }
            }


            // defines the new operator (), so that we can access the elements
            // by A(i,j), index going from i=0,...,dim()-1
            template <typename T>
            T& band_matrix<T>::operator () (int i, int j)
            {
                int k = j - i;       // what band is the entry
                assert((i >= 0) && (i < dim()) && (j >= 0) && (j < dim()));
                assert((-num_lower() <= k) && (k <= num_upper()));
                // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
                if (k >= 0)    return m_upper[k][i];
                else        return m_lower[-k][i];
            }

            template <typename T>
            T band_matrix<T>::operator () (int i, int j) const
            {
                int k = j - i;       // what band is the entry
                assert((i >= 0) && (i < dim()) && (j >= 0) && (j < dim()));
                assert((-num_lower() <= k) && (k <= num_upper()));
                // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
                if (k >= 0)    return m_upper[k][i];
                else        return m_lower[-k][i];
            }
            // second diag (used in LU decomposition), saved in m_lower

            template <typename T>
            T band_matrix<T>::saved_diag(int i) const
            {
                assert((i >= 0) && (i < dim()));
                return m_lower[0][i];
            }
            template <typename T>
            T& band_matrix<T>::saved_diag(int i)
            {
                assert((i >= 0) && (i < dim()));
                return m_lower[0][i];
            }

            // LR-Decomposition of a band matrix
            template <typename T>
            void band_matrix<T>::lu_decompose()
            {
                int  i_max, j_max;
                int  j_min;
                T x;

                // preconditioning
                // normalize column i so that a_ii=1
                for (int i = 0; i < this->dim(); i++) {
                    assert(this->operator()(i, i) != T(0.0));
                    this->saved_diag(i) = T(1.0) / this->operator()(i, i);
                    j_min = std::max(0, i - this->num_lower());
                    j_max = std::min(this->dim() - 1, i + this->num_upper());
                    for (int j = j_min; j <= j_max; j++) {
                        this->operator()(i, j) *= this->saved_diag(i);
                    }
                    this->operator()(i, i) = T(1.0);          // prevents rounding errors
                }

                // Gauss LR-Decomposition
                for (int k = 0; k < this->dim(); k++) {
                    i_max = std::min(this->dim() - 1, k + this->num_lower());  // num_lower not a mistake!
                    for (int i = k + 1; i <= i_max; i++) {
                        assert(this->operator()(k, k) != T(0.0));
                        x = -this->operator()(i, k) / this->operator()(k, k);
                        this->operator()(i, k) = -x;                         // assembly part of L
                        j_max = std::min(this->dim() - 1, k + this->num_upper());
                        for (int j = k + 1; j <= j_max; j++) {
                            // assembly part of R
                            this->operator()(i, j) = this->operator()(i, j) + x * this->operator()(k, j);
                        }
                    }
                }
            }
            // solves Ly=b
            template <typename T>
            std::vector<T> band_matrix<T>::l_solve(const std::vector<T>& b) const
            {
                assert(this->dim() == (int)b.size());
                std::vector<T> x(this->dim());
                int j_start;
                T sum;
                for (int i = 0; i < this->dim(); i++) {
                    sum = T(0);
                    j_start = std::max(0, i - this->num_lower());
                    for (int j = j_start; j < i; j++) sum += this->operator()(i, j) * x[j];
                    x[i] = (b[i] * this->saved_diag(i)) - sum;
                }
                return x;
            }
            // solves Rx=y
            template <typename T>
            std::vector<T> band_matrix<T>::r_solve(const std::vector<T>& b) const
            {
                assert(this->dim() == (int)b.size());
                std::vector<T> x(this->dim());
                int j_stop;
                T sum;
                for (int i = this->dim() - 1; i >= 0; i--) {
                    sum = T(0);
                    j_stop = std::min(this->dim() - 1, i + this->num_upper());
                    for (int j = i + 1; j <= j_stop; j++) sum += this->operator()(i, j) * x[j];
                    x[i] = (b[i] - sum) / this->operator()(i, i);
                }
                return x;
            }

            template <typename T>
            std::vector<T> band_matrix<T>::lu_solve(const std::vector<T>& b,
                bool is_lu_decomposed)
            {
                assert(this->dim() == (int)b.size());
                std::vector<T>  x, y;
                if (is_lu_decomposed == false) {
                    this->lu_decompose();
                }
                y = this->l_solve(b);
                x = this->r_solve(y);
                return x;
            }

        } // namespace internal


    } // namespace tk


} // namespace

#pragma GCC diagnostic pop

#endif /* TK_SPLINE_H */