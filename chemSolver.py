import numpy as np; from python.my_header import *
from numpy import power as npow; from numpy import loadtxt
from numpy import exp as nexp; from numpy import sqrt as nsqrt
from scipy.integrate import ode
#
# Cython
#
import pyximport; pyximport.install(setup_args={
    "include_dirs":np.get_include()},
    reload_support = True)
from Solvers import chemSolve1, getelectrons, getH2
#
# Read in the species
#
def readSpecs(infile='species.dat'):
    """
    Read the species.dat
    return it
    """
    dum = loadtxt(infile, dtype={'names': ['id', 'species'],
                            'formats': ['i', 'S5']}, usecols=(0,1,))
    specs = np.array([a for a in dum])
    return specs
""""""
#
# Indexer
#
def getIndx(a, specs):
    if a == 'PHOTON':
        return 'PH'
    elif a == 'CRP':
        return 'CRP'
    elif a == 'CRPHOTON':
        return 'CRPH'
    elif a == 'ACC': # accrete
        return 'ACC'
    else:
        try:
            return (a== specs['species']).nonzero()[0][0]
        except IndexError:
            return -99
    """"""
""""""
#
# Read in the rates: 2 body reaction only
#
def readRates(specs, infile='rates.dat'):
    """
    Read the rates.dat file and create
    """
    fin     = open('rates.dat','r')
    rates   = []
    #
    # This has to be 7 species
    #
    for line in fin:
        columns = line.strip().split()
        ncol    = len(columns)
        rate    = {}
        if ncol == 11:
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = getIndx(columns[3], specs)
            rate['P1']      = getIndx(columns[4], specs)
            rate['P2']      = getIndx(columns[5], specs)
            rate['P3']      = getIndx(columns[6], specs)
            rate['P4']      = getIndx(columns[7], specs)
            rate['alp']     = float(columns[8])
            rate['beta']    = float(columns[9])
            rate['gamma']   = float(columns[10])
        elif ncol == 10:
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = getIndx(columns[3], specs)
            rate['P1']      = getIndx(columns[4], specs)
            rate['P2']      = getIndx(columns[5], specs)
            rate['P3']      = getIndx(columns[6], specs)
            rate['P4']      = -99
            rate['alp']     = float(columns[7])
            rate['beta']    = float(columns[8])
            rate['gamma']   = float(columns[9])
        elif ncol == 9:
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = -99
            rate['P1']      = getIndx(columns[3], specs)
            rate['P2']      = getIndx(columns[4], specs)
            rate['P3']      = getIndx(columns[5], specs)
            rate['P4']      = -99
            rate['alp']     = float(columns[6])
            rate['beta']    = float(columns[7])
            rate['gamma']   = float(columns[8])
        elif ncol == 8:
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = -99
            rate['P1']      = getIndx(columns[3], specs)
            rate['P2']      = getIndx(columns[4], specs)
            rate['P3']      = -99
            rate['P4']      = -99
            rate['alp']     = float(columns[5])
            rate['beta']    = float(columns[6])
            rate['gamma']   = float(columns[7])
        elif ncol == 7: # ICE
            rate['R1']      = getIndx(columns[1], specs)
            rate['R2']      = getIndx(columns[2], specs)
            rate['R3']      = -99
            rate['P1']      = getIndx(columns[3], specs)
            rate['P2']      = -99
            rate['P3']      = -99
            rate['P4']      = -99
            rate['alp']     = float(columns[4])
            rate['beta']    = float(columns[5])
            rate['gamma']   = float(columns[6])
        else:
            if columns[0] == '9999':
                break
            else:
                print 'UNKNOWN text format! %d \n'%(ncol)
                print line, columns[0]
                print 'Check input file\n'
                raise SystemExit
        """"""
        rates.append(rate)
    """"""
    fin.close()
    del rate
    #
    # Get the indices and rates directly
    #
    rates1 = np.zeros((len(rates),11))
    for irate in range(len(rates)):
        """
        Options for rate reactions
        """
        rates1[irate,0]         = 0.0
        if rates[irate]['R2'] == 'CRP':
            rates1[irate, 0]    = 1.
            rates1[irate, 2]    = -99
        elif rates[irate]['R2'] == 'CRPH':
            rates1[irate, 0]    = 2.
            rates1[irate, 2]    = -99
        elif rates[irate]['R2'] == 'PH':
            rates1[irate, 0]    = 3.
            rates1[irate, 2]    = -99
        elif rates[irate]['R2'] == 'ACC':
            rates1[irate, 0]    = 4.
            rates1[irate, 2]    = -99
        else:
            rates1[irate, 2]    = rates[irate]['R2']
        """"""
        rates1[irate, 1]    = rates[irate]['R1']
        rates1[irate, 3]    = rates[irate]['R3']
        rates1[irate, 4]    = rates[irate]['P1']
        try:
            rates1[irate, 5]    = rates[irate]['P2']
        except ValueError:
            rates1[irate, 5]    = -99
        rates1[irate, 6]    = rates[irate]['P3']
        rates1[irate, 7]    = rates[irate]['P4']
        rates1[irate, 8]    = rates[irate]['alp']
        rates1[irate, 9]    = rates[irate]['beta']
        rates1[irate, 10]   = rates[irate]['gamma']
    """"""
    rates = []
    del rates
    return np.array(rates1)
""""""
#
# Init gas species
#
def initGas(specs, H2=1e8):
    """
    Initialize the gas species
    """
    ngas = np.zeros(specs.shape[0])
    raise SystemExit
""""""
#
# Plot Chemistry
#
def PlotChem(wsol, specs, arrs):
    """
    Plot it
    """
    for ispec in xrange(len(specs)):
        print '%s \t %2.3e  -- %2.3e'%(specs['species'][ispec], wsol[0,ispec+1],
            wsol[-1, ispec+1])
    """"""
    inorm   = ('H' == specs['species']).nonzero()[0][0]
    norm    = wsol[:,inorm+1]
    #
    # Find where these species are
    #
    arrs = [(a == specs['species']).nonzero()[0][0] for a in arrs]
    #
    #
    #
    fig, ax = subplots(1,1,figsize=[aaonecol, aaonecol])
    subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95)
    #
    #
    #
    x   = wsol[:,0]
    for ar in arrs:
        ax.loglog(x, wsol[:,ar+1], '-', alpha=0.8,
            label=r'\boldmath${\rm %s}$'%(specs[ar]) )
    #
    #
    #
    leg = ax.legend(loc='upper right', fontsize=7)
    leg.draw_frame(False)
    #
    #
    #
    ax = fig_labs(ax, xlim=[1e1, 1e8], ylim=[1e-12, 3e-4],
        xlab=r'\boldmath${\rm Time \ [yrs]}$',
        ylab=r'\boldmath${N/N_{\rm H}}$',
        xform=1, yform=1)
    show()
    close()
""""""
#
# Grab for test
#
def grab(x, specs):
    arrs    = [0, 1, 2, 3, 4, 5, 6, 7, 8 , 9, 20, 22, 65, 29]
    for ar in arrs:
        print '%s %2.4e'%(specs['species'][ar], x[ar+1])
    print
""""""
#
# Main
#
def main():
    specs   = readSpecs()
    rates   = readRates(specs)
    #
    # initialize number densities
    #
    ngas    = np.zeros(len(specs), dtype=np.float) + 1e-25
    #
    # initialize fractions
    #
    density         = 2e5
    temperature     = 10.0
    rad             = 1e2
    avmag           = 10.0
    #
    #
    #
    eFrac           = 0.0 # 1e-7
    H2Frac          = 0.5
    HeFrac          = 0.14
    CFrac           = 7.3e-5
    NFrac           = 2.14e-5
    RCO             = 7.3/17.6 # 1.2 or 7.3/17.6
    OFrac           = CFrac/RCO
    MgFrac          = 3.0e-9
    #
    # Set the abundances
    #
    FH      = 1e-6
    TOTAL   = np.zeros(3)
    TOTAL[0]    = eFrac*density
    TOTAL[1]    = H2Frac*density
    TOTAL[2]    = HeFrac*density
    ngas[specs['species'] == 'H']   = FH * density
    ngas[specs['species'] == 'N']   = NFrac * density
    ngas[specs['species'] == 'C+']  = CFrac * density
    ngas[specs['species'] == 'O']   = OFrac * density
    ngas[specs['species'] == 'Mg']  = MgFrac * density
    #
    # Correct electrons, H2 and He
    #
    xe      = (0.0 + TOTAL[0] + getelectrons(ngas))
    xh      = (0.0 + TOTAL[1] + (-0.5* getH2(ngas)))
    xhe     = (0.0 + TOTAL[2] - ngas[4])
    ngas[69]    = xhe
    ngas[68]    = xh
    ngas[67]    = xe
    #
    # Solve this
    #
    tstart  = 1e1*np.pi*1e7
    tend    = 1e10*np.pi*1e7
    nt      = 1000
    print
    print 'START: ---'
    print 'He: %2.4e, N: %2.4e, Mg: %2.4e, C: %2.4e, O: %2.4e'%(
        ngas[69]/density,
        ngas[11]/density,
        ngas[35]/density,
        ngas[5]/density,
        ngas[22]/density)
    #
    # Vode init
    #
    wsol    = []
    w0      = ngas
    p       = [rates, temperature, rad, avmag, TOTAL, density]
    for it in xrange(nt-1):
        if it == 0:
            dt      = tstart
            tnext   = tstart
            vode = ode(chemSolve1).set_integrator('vode',
                atol=1e-30, rtol=1e-8, method='bdf', nsteps=1e4,
                first_step=1e-10, with_jacobian=True, order=5)
            vode.set_initial_value(w0, 0.0).set_f_params(p)
        else:
            tnext   = npow(10.0, np.log10(tstart)+0.05)
            dt      = tnext - tstart
            vode = ode(chemSolve1).set_integrator('vode',
                atol=1e-30, rtol=1e-8, method='bdf', nsteps=1e4,
                first_step=tstart*1e-10, with_jacobian=True, order=5)
            vode.set_initial_value(w0, tstart).set_f_params(p)
        """"""
        if (tstart > tend):
            break
        elif (tstart == tend):
            break
        elif (tnext > tend):
            dt  = tend - tstart
        w1  = [vode.t/(np.pi*1e7)] + [a/density for a in w0]
        #        print it, vode.t/(np.pi*1e7)
        #        grab(w1, specs)
        #        if it == 5: raise SystemExit
        #
        # Add to the solution
        #
        wsol.append(w1)
        #
        # Integrate
        #
        vode.integrate(vode.t+dt)
        #
        # Set new inputs for the next step
        #
        w0 = vode.y
        for iw in xrange(len(w0)):
            if w0[iw] < 1e-90: w0[iw] = 1e-90
        tstart = vode.t
        if vode.successful():
            pass
        else:
            print vode.t, it
            print 'SOLVER FAILS'
            raise SystemExit
    """"""
    #
    # Plot results
    #
    PlotChem(np.array(wsol), specs, ['CN', 'HCN', 'HN2+'])
""""""
#
# call function
#
if __name__ == "__main__":
    main()
""""""