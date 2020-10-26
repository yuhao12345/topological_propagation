# Tutorial 2.3.1. Matrix structure of on-site and hopping elements
# ================================================================
#
# Physics background
# ------------------
#  Gaps in quantum wires with spin-orbit coupling and Zeeman splititng,
#  as theoretically predicted in
#   http://prl.aps.org/abstract/PRL/v90/i25/e256601
#  and (supposedly) experimentally oberved in
#   http://www.nature.com/nphys/journal/v6/n5/abs/nphys1626.html
#
# Kwant features highlighted
# --------------------------
#  - Numpy matrices as values in Builder

import kwant

# For plotting
from matplotlib import pyplot

# For matrix support
import tinyarray

# define Pauli-matrices for convenience
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])


def make_system(a=1, t=1.0, alpha=0.5, e_z=0.08, W=10, L=30):
    # Start with an empty tight-binding system and a single square lattice.
    # `a` is the lattice constant (by default set to 1 for simplicity).
    lat = kwant.lattice.square(a)

    sys = kwant.Builder()

    #### Define the scattering region. ####
    sys[(lat(x, y) for x in range(L) for y in range(W))] = \
        4 * t * sigma_0 + e_z * sigma_z
    # hoppings in x-direction
    sys[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_y
    # hoppings in y-directions
    sys[kwant.builder.HoppingKind((0, 1), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_x

    #### Define the left lead. ####
    lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))

    lead[(lat(0, j) for j in xrange(W))] = 4 * t * sigma_0 + e_z * sigma_z
    # hoppings in x-direction
    lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_y
    # hoppings in y-directions
    lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_x

    #### Attach the leads and return the finalized system. ####
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys


def plot_conductance(sys, energies):
    # Compute conductance
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(sys, energy)
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()


def main():
    sys = make_system()

    # Check that the system looks as intended.
    kwant.plot(sys)

    # Finalize the system.
    sys = sys.finalized()

    # We should see non-monotonic conductance steps.
    plot_conductance(sys, energies=[0.01 * i - 0.3 for i in xrange(100)])


# Call the main function if the script gets executed (as opposed to imported).
# See <http://docs.python.org/library/__main__.html>.
if __name__ == '__main__':
    main()
