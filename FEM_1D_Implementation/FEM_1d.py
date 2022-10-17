#*****************************************************************************80
#
## barebones() defines a "bare bones" finite element solution scheme.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    03 March 2021
#
#  Author:
#
#    John Burkardt
#
# This is a modified version by Andres Rubiano

import numpy as np
import scipy.linalg as la
import time
import matplotlib.pyplot as plt

def barebones ( n, a, c, f, x ):

    quad_num = 2
    abscissa = np.array ( [ -0.577350269189625764509148780502, \
                            +0.577350269189625764509148780502 ] )
    weight = np.array ( [ 1.0, 1.0 ] )
    
    A = np.zeros ( [ n, n ] )
    b = np.zeros ( n )
    
    e_num = n - 1
    
    for e in range ( 0, e_num ):
    
        l = e
        xl = x[l]
        r = e + 1
        xr = x[r]
      
        for q in range ( 0, quad_num ):
      
            xq = ( ( 1.0 - abscissa[q] ) * xl   \
                 + ( 1.0 + abscissa[q] ) * xr ) \
                 / 2.0
            wq = weight[q] * ( xr - xl ) / 2.0
        
            vl = ( xr - xq ) / ( xr - xl )
            vlp =     - 1.0  / ( xr - xl )
        
            vr = ( xq - xl ) / ( xr - xl )
            vrp = +1.0       / ( xr - xl )
        
            axq = a ( xq )
            cxq = c ( xq )
            fxq = f ( xq )
        
            A[l,l] = A[l,l] + wq * ( vlp * axq * vlp + vl * cxq * vl )
            A[l,r] = A[l,r] + wq * ( vlp * axq * vrp + vl * cxq * vr )
            b[l]   = b[l]   + wq * ( vl * fxq )
        
            A[r,l] = A[r,l] + wq * ( vrp * axq * vlp + vr * cxq * vl )
            A[r,r] = A[r,r] + wq * ( vrp * axq * vrp + vr * cxq * vr )
            b[r]   = b[r]   + wq * ( vr * fxq )
    
    for j in range ( 0, n ):
        A[0,j] = 0.0
    A[0,0] = 1.0
    b[0] = 0.0
    
    for j in range ( 0, n ):
        A[n-1,j] = 0.0
    A[n-1,n-1] = 1.0
    b[n-1] = 0.0
    
    u = la.solve ( A, b )
    
    return u

def fem1d_bvp_linear ( n, a, c, f, x ):

#*****************************************************************************80
#
## fem1d_bvp_linear() solves a two point boundary value problem.
#
#  Discussion:
#
#    The program uses the finite element method, with piecewise linear basis
#    functions to solve a boundary value problem in one dimension.
#
#    The problem is defined on the region 0 <= x <= 1.
#
#    The following differential equation is imposed between 0 and 1:
#
#      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
#
#    where a(x), c(x), and f(x) are given functions.
#
#    At the boundaries, the following conditions are applied:
#
#      u(0.0) = 0.0
#      u(1.0) = 0.0
#
#    A set of N equally spaced nodes is defined on this
#    interval, with 0 = X(1) < X(2) < ... < X(N) = 1.0.
#
#    At each node I, we associate a piecewise linear basis function V(I,X),
#    which is 0 at all nodes except node I.  This implies that V(I,X) is
#    everywhere 0 except that
#
#    for X(I-1) <= X <= X(I):
#
#      V(I,X) = ( X - X(I-1) ) / ( X(I) - X(I-1) ) 
#
#    for X(I) <= X <= X(I+1):
#
#      V(I,X) = ( X(I+1) - X ) / ( X(I+1) - X(I) )
#
#    We now assume that the solution U(X) can be written as a linear
#    sum of these basis functions:
#
#      U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X)
#
#    where U(X) on the left is the function of X, but on the right,
#    is meant to indicate the coefficients of the basis functions.
#
#    To determine the coefficient U(J), we multiply the original
#    differential equation by the basis function V(J,X), and use
#    integration by parts, to arrive at the I-th finite element equation:
#
#        Integral A(X) * U'(X) * V'(I,X) + C(X) * U(X) * V(I,X) dx 
#      = Integral F(X) * V(I,X) dx
#
#    We note that the functions U(X) and U'(X) can be replaced by
#    the finite element form involving the linear sum of basis functions,
#    but we also note that the resulting integrand will only be nonzero
#    for terms where J = I - 1, I, or I + 1.
#
#    By writing this equation for basis functions I = 2 through N - 1,
#    and using the boundary conditions, we have N linear equations
#    for the N unknown coefficients U(1) through U(N), which can
#    be easily solved.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 January 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of nodes.
#
#    function A ( X ), evaluates a(x);
#
#    function C ( X ), evaluates c(x);
#
#    function F ( X ), evaluates f(x);
#
#    real X(N), the mesh points.
#
#  Output:
#
#    real U(N), the finite element coefficients, which are also
#    the value of the computed solution at the mesh points.
#
#
#  Define a 2 point Gauss-Legendre quadrature rule on [-1,+1].
#
    quad_num = 2
    abscissa = np.array ( [ -0.577350269189625764509148780502, \
                            +0.577350269189625764509148780502 ] )
    weight = np.array ( [ 1.0, 1.0 ] )
#
#  Make room for the matrix A and right hand side b.
#
    A = np.zeros ( [ n, n ] )
    b = np.zeros ( n )
#
#  We assemble the finite element matrix by looking at each element,
#  that is, each interval [ X(L), X(R) ].
#
#  There are only two nonzero piecewise linear basis functions in this
#  element, and we can call them VL and VR.  So the only integrals we
#  need to compute involve products of:
#
#    VL * VL   VR * VL    F * VL
#    VR * VL   VR * VR    F * VR
#
#  and
#
#    VL' * VL'   VR' * VL'
#    VR' * VL'   VR' * VR'
#
    e_num = n - 1
    
    for e in range ( 0, e_num ):

        l = e
        xl = x[l]
    
        r = e + 1
        xr = x[r]

        for q in range ( 0, quad_num ):
    
            xq = ( ( 1.0 - abscissa[q] ) * xl   \
                 + ( 1.0 + abscissa[q] ) * xr ) \
                 / 2.0
        
            wq = weight[q] * ( xr - xl ) / 2.0
        
            vl = ( xr - xq ) / ( xr - xl )
            vlp =     - 1.0  / ( xr - xl )
        
            vr = ( xq - xl ) / ( xr - xl )
            vrp = +1.0       / ( xr - xl )
        
            axq = a ( xq )
            cxq = c ( xq )
            fxq = f ( xq )
        
            A[l,l] = A[l,l] + wq * ( vlp * axq * vlp + vl * cxq * vl )
            A[l,r] = A[l,r] + wq * ( vlp * axq * vrp + vl * cxq * vr )
            b[l]   = b[l]   + wq * ( vl * fxq )
        
            A[r,l] = A[r,l] + wq * ( vrp * axq * vlp + vr * cxq * vl )
            A[r,r] = A[r,r] + wq * ( vrp * axq * vrp + vr * cxq * vr )
            b[r]   = b[r]   + wq * ( vr * fxq )
#
#  Equation 0 is the left boundary condition, U(0) = 0.0;
#
    for j in range ( 0, n ):
        A[0,j] = 0.0
    A[0,0] = 1.0
    b[0] = 1.0
#
#  We can keep the matrix symmetric by using the boundary condition
#  to eliminate U(0) from all equations but #0.
#
    for i in range ( 1, n ):
        b[i] = b[i] - A[i,0] * b[0]
        A[i,0] = 0.0
#
#  Equation N-1 is the right boundary condition, U(N-1) = 0.0;
#
    for j in range ( 0, n ):
        A[n-1,j] = 0.0
    A[n-1,n-1] = 1.0
    b[n-1] = 1.0
#
#  We can keep the matrix symmetric by using the boundary condition
#  to eliminate U(N-1) from all equations but #N-1.
#
    for i in range ( 0, n - 1 ):
        b[i] = b[i] - A[i,n-1] * b[n-1]
        A[i,n-1] = 0.0
#
#  Solve the linear system for the finite element coefficients U.
#
    u = la.solve ( A, b )

    return u

def fem1d_bvp_linear_test00 (n):

#*****************************************************************************80
#
## fem1d_bvp_linear_test00() tests fem1d_bvp_linear().
#
#  Discussion:
#
#    - uxx + u = x  for 0 < x < 1
#    u(0) = u(1) = 0
#
#    exact  = x - sinh(x) / sinh(1)
#    exact' = 1 - cosh(x) / sinh(1)
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 January 2015
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Dianne O'Leary,
#    Scientific Computing with Case Studies,
#    SIAM, 2008,
#    ISBN13: 978-0-898716-66-5,
#    LC: QA401.O44.
#
#
#  Geometry definitions.
#
    x_lo = 0.0
    x_hi = 1.0
    x = np.linspace ( x_lo, x_hi, n )
    
    u = fem1d_bvp_linear ( n, a00, c00, f00, x )
    
    g = np.zeros ( n )
    for i in range ( 0, n ):
        g[i] = exact00 ( x[i] )

#  Plot the computed solution.
#
    fig = plt.figure ( )
    plt.plot ( x, u, 'bo-' )
    plt.xlabel ( '<---x--->' )
    plt.ylabel ( '<---u(x)--->' )
    plt.title ( 'fem1d_bvp_linear_test00 Solution n='+str(n) )
    plt.savefig ( 'fem1d_bvp_linear_test00.png' )
    plt.show ( block = False )
    plt.close ( )
    
    return u

def a00 ( x ):

#*****************************************************************************80
#
## a00() evaluates the A coefficient.
#
    value = 1.0
    
    return value

def c00 ( x ):

#*****************************************************************************80
#
## c00() evaluates the C coefficient.
#
    value = 0
    
    return value

def exact00 ( x ):

#*****************************************************************************80
#
## exact00() evaluates the exact solution.
#

    value = ((x**2)/2) - ((1/2) *x) + 1
    
    return value

def exactp00 ( x ):

#*****************************************************************************80
#
## exactp00() evaluates the exact derivative of the solution.
#

    value = (x) - (1/2)
    
    return value

def f00 ( x ):

#*****************************************************************************80
#
## f00() evaluates the right hand side function.
#
    value = -1
    
    return value


def fem1d_bvp_linear_test ( ):

    print ( '' )
    print ( '  Solve -u\'\'(x) = f(x)' )
    print ( '  for 0 < x < 1, with u(0) = u(1) = 1.0' )
    print ( '  f(x)  = -1' )
    print ( '' )
    print ( ' The analytic solution is: ' )
    print ( '  u(x) = ((x**2)/2) - ((1/2) *x) + 1' )
#
#  Library functions.
#
    listOfTests = [4,7,13,25,49,97,193,385]
  
    print ( '' )
    print ( '  h1s_error_linear computes the H1 seminorm approximation error' )
    print ( '  between the exact derivative of a function and the derivative' )
    print ( '  of a piecewise linear approximation to the function,' )
    print ( '  associated with n mesh points x().'  )
    print ( '' )
    print ( '   N    H1S_Error' )
    print ( '' )

    for t in listOfTests:
        u = fem1d_bvp_linear_test00 (t)
        h1s_error_linear_test (u,t)


def h1s_error_linear ( n, x, u, exact_ux ):

#*****************************************************************************80
#
## h1s_error_linear(): seminorm error of a finite element solution.
#
#  Discussion:
#
#    We assume the finite element method has been used, over an interval [A,B]
#    involving N nodes, with piecewise linear elements used for the basis.
#    The finite element solution U(x) has been computed, and a formula for the
#    exact derivative V'(x) is known.
#
#    This function estimates the H1 seminorm of the error:
#
#      H1S = sqrt ( integral ( A <= x <= B ) ( U'(x) - V'(x) )^2 dx
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    12 July 2015
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    integer N, the number of nodes.
#
#    real X(N), the mesh points.
#
#    real U(N), the finite element coefficients.
#
#    function EQ = EXACT_UX ( X ), returns the value of the exact
#    derivative at the point X.
#
#  Output:
#
#    real H1S, the estimated seminorm of the error.
#

    h1s = 0.0
#
#  Define a 2 point Gauss-Legendre quadrature rule on [-1,+1].
#
    quad_num = 2
    abscissa = np.array ( [ -0.577350269189625764509148780502, \
                            +0.577350269189625764509148780502 ] )
    weight = np.array ( [ 1.0, 1.0 ] )
#
#  Integrate over each interval.
#
    e_num = n - 1

    for e in range ( 0, e_num ):

        l = e
        xl = x[l]
        ul = u[l]
        
        r = e + 1
        xr = x[r]
        ur = u[r]

        for q in range ( 0, quad_num ):
        
            xq = ( ( 1.0 - abscissa[q] ) * xl   \
                 + ( 1.0 + abscissa[q] ) * xr ) \
                 /   2.0
          
            wq = weight[q] * ( xr - xl ) / 2.0
#
#  The piecewise linear derivative is a constant in the interval.
#
        uxq = ( ur - ul ) / ( xr - xl )
    
        exq = exact_ux ( xq )
    
        h1s = h1s + wq * ( uxq - exq ) ** 2

    h1s = np.sqrt ( h1s )

    return h1s

def h1s_error_linear_test (u,t):

    x_n = t
    
    
    x_lo = 0.0
    x_hi = 1.0
    x = np.linspace ( x_lo, x_hi, x_n )
   
    #
    #  Compare the derivative of the piecewise interpolant of U
    #  to the actual derivative
    #
    e1 = h1s_error_linear ( x_n, x, u, exactp00 )
   
    print ( '  %2d  %g' % ( x_n, e1 ) )


    return


def timestamp ( ):

  t = time.time ( )
  print ( time.ctime ( t ) )

  return None

if ( __name__ == '__main__' ):
  timestamp ( )
  fem1d_bvp_linear_test ( )
  timestamp ( )
