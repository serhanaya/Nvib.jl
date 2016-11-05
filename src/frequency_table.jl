## Finding Roots.
#
#
#

""" find_all_roots(f::Function; step_size=0.5, start=0, finish=10e10, tol=1e-6, max_iter=10000, root_num_max=5)

Find all roots of a given function in a given range using Bisection method.

Parameters
---------
f             : Continuous function in a given interval
step_size     : Step size for root finding
start         : Start value (Beginning of the interval)
finish        : End value (End of the interval)
tol           : Tolerance (stop criteria_1 for this function)
max_iter      : Maximum allowed iteration (stop criteria_2 for this function)
root_num_max  : How many roots to find in a given interval.

Returns
-------
roots         : Vector containing the roots of the continuous function [in given range].

Example
-------
using PyPlot

f(x) = sin(x)
x = linspace(-2pi, 2pi, 100)
r_arr = find_all_roots(f, step_size=0.1, start=-2pi, finish=2pi, tol=1e-10, max_iter=10000, root_num_max=2)

plot(x, f(x))
hlines(y=0, xmin=-2pi, xmax=2pi)
scatter(r_arr, f(r_arr))
"""
function find_all_roots(f::Function; step_size=0.5, start=0, finish=10e10, tol=1e-6, max_iter=10000, root_num_max=5)

  roots = zeros(root_num_max)

  a = start
  for no_of_root in 1:root_num_max
    b = a + step_size
    c = (a + b) / 2

    if b > finish
      print("Interval exceeded")
      return roots
    end

    if f(a) == 0
      zero = a
    elseif f(b) == 0
      zero = b
    end

    while f(a) * f(b) > 0
      a = b
      b = a + step_size
      c = (a + b) * 0.5
    end

    niter = 0
    interval = (b - a) * 0.5

    while (interval >= tol && niter < max_iter)
      niter += 1
      if f(a) * f(c) < 0
        b = c
        c = a + (b - a) * 0.5
        interval = (b - a) * 0.5
      elseif f(c) * f(b) < 0
        a = c
        c = a + (b - a) * 0.5
        interval = (b - a) * 0.5
      end
    end

    if (niter == max_iter && interval > tol)
      println("""Bisection stopped without converging
                 to the desired tolerance because the
                 maximum number of iterations was reached.""")
    end
    zero = c
    roots[no_of_root] = zero
    a = b
  end

  return roots
end  # end of function find_all_roots.

""" tab(atype::Beam; min_lambda=10, max_lambda=150, no_lambda=10, min_rad0=0.2,
    max_rad0=2, no_rad0=10, min_phit=pi/18, max_phit=5pi/6, no_phit=10, min_step=0,
    max_step=10e10, size_step=10, tolerance=1e-6, iter_max=10000, no_of_roots=5)

A parametric approach for root finding. Use find_all_roots method and return
a table containing mode frequencies for a given analysis type.

Parameters
----------
min_lambda          : Minimum value for slenderness ratio.
max_lambda          : Maximum value for slenderness ratio.
no_lambda           : This determines the step size (s_size=(max_l-min_l)/no_lambda).
min_rad0            : Minimum value for initial radius R0.
max_rad0            : Maximum value for initial radius R0.
no_rad0             : This determines the step size (s_size=(max_r0-min_r0)/no_rad0).
min_phit            : Minimum value for opening angle.
max_phit            : Maximum value for angle.
no_phit             : This determines the step size (s_size=(max_p-min_p)/no_phit).
min_step            : Start point for searching roots (see find_all_roots - start).
max_step            : End point for searching roots (see find_all_roots - finish).
size_step           : Step size for searching roots (see find_all_roots - step).
tolerance           : (see find_all_roots - tol)
iter_max            : (see find_all_roots - max_iter)
no_of_roots         : How many roots are required? (see find_all_roots - root_num_max)

Returns
-------
result              :

"""
function freqTab(atype::Beam; min_lambda=10, max_lambda=150, no_lambda=10, min_rad0=0.2,
    max_rad0=2, no_rad0=10, min_phit=pi/18, max_phit=5pi/6, no_phit=10, min_step=0,
    max_step=10e10, size_step=10, tolerance=1e-6, iter_max=10000, no_of_roots=5)

    start_time = time()

    l = linspace(min_lambda, max_lambda, no_lambda)
    r = linspace(min_rad0, max_rad0, no_rad0)
    p = linspace(min_phit, max_phit, no_phit)

    detfunc(om) = detrmt(atype, om)
    rootfunc(f::Function) = find_all_roots(f, step_size=size_step, start=min_step, finish=max_step,
                                           tol=tolerance, max_iter=iter_max, root_num_max = no_of_roots)
    result = zeros(no_lambda, no_rad0, no_phit, no_of_roots)
    n = length(l) * length(r) * length(p) * length(rootfunc(detfunc))
    pr = Progress(n, 1)

    for (i, lam) in enumerate(l)
        atype.lam = lam
        for (j, rad0) in enumerate(r)
            atype.rad0 = rad0
            for (k, phit) in enumerate(p)
                atype.phit = phit
                for (m, root) in enumerate(rootfunc(detfunc))
                    h = sqrt(12) * rad0 * phit/lam  # h : height = variable
                    b = h * 2/3  # b : width = constant
                    c = root * (rad0 * phit)^2 * sqrt(atype.ro * b * h / (atype.el * b * h^3 / 12))
                    result[i, j, k, m] = c  # sonuc: boyutsuz frekans elde edildi.
                    next!(pr)
                end
            end
        end
    end
    
    return result

end  # end of function tab.
