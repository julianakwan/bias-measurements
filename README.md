Scripts to measure cosmological observables from N-body simulations

fit_bias.py fits for linear bias, b1, and its next order perturbation,
b2, from halo and dark matter density fields measured from N-body
simulations. I apply a simple non-linear least squares fit (scipy
curve_fit) to determine b1, b2 from delta_h = b1 * delta_m + b2 *
delta_m**2 + eps.

scatter_plot.py makes nice contour plots of the dark matter and halo
density fields for papers, talks etc.

xi.py measures the 2-point correlation function using Corrfunc, a
python package by Manodeep Singha
(http://corrfunc.readthedocs.io/en/master/index.html)

convert_xi.py takes the 2-point functions measured in xi.py and fits a
smooths them using a Savitzky-Golay filter. 
