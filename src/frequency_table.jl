""" find_all_roots(f::Function; step_size=0.5, start=0, finish=10e10, tol=1e-6, max_iter=10000, root_num_max=5, wrn=true)

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
wrn           : Warns when the maximum iteration is reached.

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
function find_all_roots(f::Function; step_size=0.5, start=0, finish=10e10, tol=1e-6, max_iter=10000, root_num_max=5,
                        wrn=true)

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

    if (wrn == true && (niter == max_iter && interval > tol))
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

""" freqTab(atype::Beam; min_lambda=10, max_lambda=150, no_lambda=10, min_rad0=1.56,
    max_rad0=15.6, no_rad0=10, min_thetat=pi/18, max_thetat=5pi/6, no_thetat=10, min_step=0,
    max_step=10e10, size_step=10, tolerance=1e-6, iter_max=10000, no_of_roots=5, wrn=false)

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
min_thetat          : Minimum value for opening angle.
max_thetat          : Maximum value for angle.
no_thetat           : This determines the step size (s_size=(max_p-min_p)/no_thetat).
min_step            : Start point for searching roots (see find_all_roots - start).
max_step            : End point for searching roots (see find_all_roots - finish).
size_step           : Step size for searching roots (see find_all_roots - step).
tolerance           : (see find_all_roots - tol)
iter_max            : (see find_all_roots - max_iter)
no_of_roots         : How many roots are required? (see find_all_roots - root_num_max)
wrn                 : Warning for maximum iteration while root searching.

Returns
-------
result              : result[i,j,k,m] is the nondimensional frequency. Here, i stands for the
                      slenderness ratio, j is the radius which determines the ratio R/γ (the
                      small scale parameter), k is the value of the opening angle of the beam,
                      and lastly, m is the mode number.

Example
-------
This function is created for a parametric analysis in mind, however we can calculate nondimensional
frequencies for specific conditions:

If we want to determine the nondimensional frequency for the first and second mode of a nonlocal
beam (in-plane vibration analysis) having slenderness ratio of Λ=150, opening angle θT=120 and small
scale parameter R0/γ=1 where γ=1.56 nm we will follow the given steps:

julia> Pkg.clone("git@github.com:serhanaya/Nvib.jl")
julia> analysis1 = IPBeam()
julia> nondimfrq = freqTab(analysis1, min_lambda=150, max_lambda=150, no_lambda=1, min_rad0=1.56,
    max_rad0=1.56, no_rad0=1, min_thetat=2pi/3, max_thetat=2pi/3, no_thetat=1, min_step=0,
    max_step=10e10, size_step=2, tolerance=1e-6, iter_max=10000, no_of_roots=2)
julia> nondimfrq[1, 1, 1, :]  # gives an array consist of the nondimensional frequencies for the
first and second mode.
"""
function freqTab(atype::Beam; min_lambda=10, max_lambda=150, no_lambda=10, min_rad0=1.56,
    max_rad0=15.6, no_rad0=10, min_thetat=pi/18, max_thetat=5pi/6, no_thetat=10, min_step=0,
    max_step=10e10, size_step=10, tolerance=1e-6, iter_max=10000, no_of_roots=5, wrn=false)

    l = linspace(min_lambda, max_lambda, no_lambda)
    r = linspace(min_rad0, max_rad0, no_rad0)
    p = linspace(min_thetat, max_thetat, no_thetat)

    detfunc(om) = detrmt(atype, om)
    rootfunc(f::Function) = find_all_roots(f, step_size=size_step, start=min_step, finish=max_step,
                                           tol=tolerance, max_iter=iter_max, root_num_max = no_of_roots,
                                           wrn=wrn)
    result = zeros(no_lambda, no_rad0, no_thetat, no_of_roots)
    n = length(l) * length(r) * length(p) * length(rootfunc(detfunc))
    pr = Progress(n, 1)

    for (i, lam) in enumerate(l)
        atype.lam = lam
        for (j, rad0) in enumerate(r)
            atype.rad0 = rad0
            for (k, thetat) in enumerate(p)
                atype.thetat = thetat
                for (m, root) in enumerate(rootfunc(detfunc))
                    h = sqrt(12) * rad0 * thetat/lam  # h : height = variable
                    b = h * 2/3  # b : width = constant
                    area_moment(analysis_type::IPBeam) = b * h^3 / 12
                    area_moment(analysis_type::OPBeam) = h * b^3 / 12
                    c = root * (rad0 * thetat)^2 * sqrt(atype.ro * b * h / (atype.el * area_moment(atype)))
                    result[i, j, k, m] = c  # result[i,j,k,m] is the nondimensional frequency
                    next!(pr)
                end
            end
        end
    end

    return result

end  # end of function freqTab.
