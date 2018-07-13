Scripts to measure b1, b2 from simulations

These fit for linear bias, b1, and its next order perturbation, b2, from
halo and dark matter density fields measured from N-body simulations. I apply
a simple non-linear least squares fit (scipy curve_fit) to determine b1,b2
from delta_h = b1*delta_m + b2*delta_m**2 + eps. 