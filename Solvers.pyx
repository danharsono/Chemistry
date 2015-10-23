"""
Cython version of the chemical solver
"""
from __future__ import division
from numpy cimport ndarray, dtype
import numpy as np; cimport numpy as np
cimport cython; from cpython cimport bool
#
# Types
#
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
#
# Constants here and externs from math
#
cdef extern from "math.h":
    double sqrt(double m)
    double pow(double m, double n)
    double exp(double a)
    double fabs(double a)
#
# Functions here
#
cpdef double getelectrons(np.ndarray[DTYPE_t, ndim=1] x):
    cdef double dum
    cdef int iiter
    cdef list bi    = [0, 2, 3, 4, 5, 8, 10, 12, 15, 16, 17, 18,
        20, 23, 24, 27, 28, 30, 31, 32, 34, 36, 38, 39, 45,
        46, 48, 49, 51, 52, 53, 55, 56, 58, 61, 62, 63, 64]
    dum = 0.0
    for iiter in xrange(len(bi)):
        dum += x[<unsigned int> bi[iiter]]
    return dum
""""""
cpdef double getH2(np.ndarray[DTYPE_t, ndim=1] x):
    cdef double dum
    cdef int iiter
    dum     = 0.0
    dum     += x[0]+x[1] + 2.*x[2] + 3.*x[3]+x[7]+x[8]+2.*x[9]
    dum     += 2.*x[10] +x[13] + 3.*x[14]+3.*x[15]+x[16]+4.*x[17]+2.*x[18]
    dum     += 2.*x[19]
    dum     += 4.*x[21]+x[23]+3.*x[24]+3.*x[25]+x[26]+5.*x[27]+4.*x[28]+2.*x[29]
    dum     += 2.*x[30]+3.*x[31]+x[36]+x[37]+2.*x[38]+2.*x[41]+x[42]+x[43]
    dum     += 3.*x[44]+3.*x[45]+x[46]+2.*x[51]+4.*x[52]+x[53]+x[54]+x[55]
    dum     += 2.*x[56]+2.*x[59]+x[60]+3.*x[61]+x[62]+2.*x[64]
    return dum
""""""
#
# Ices
#
@cython.boundscheck(False) # turns off bounds-checking for functions
cdef iceform(double T, double eff = 1.0):
    return <double> (4.55e-18*sqrt(T/28.)*eff*2e5)
""""""
@cython.boundscheck(False) # turns off bounds-checking for functions
cdef chemSolve(double t, np.ndarray[DTYPE_t, ndim=1] x, np.ndarray[DTYPE_t, ndim=2] rate, double Temp, double rad, double avmag, np.ndarray[DTYPE_t, ndim=1] tot, double ndens):
    """
    Solve chemistry
    """
    #
    # define variables
    #
    cdef int irate, ir, ip1, ip2, ip3, ip4, ip5, ip6, ip7
    cdef double k, d1, d2, d3, d4, d5, d6, d7
    cdef double xe, xh, xhe
    #
    # set up the change of number densities
    #
    cdef np.ndarray[DTYPE_t, ndim=1] dndt    = (np.zeros(x.shape[0],
        dtype=DTYPE))
    #
    # Correct electrons, H2 and He
    #
    xe      = (0.0 + tot[0] + getelectrons(x))
    xh      = (0.0 + tot[1] + (-0.5* getH2(x)))
    xhe     = (0.0 + tot[2] - x[4])
    x[69]   = xhe
    x[68]   = xh
    x[67]   = xe
    #
    # loop for every species and get the rates
    # alpha - 8, beta - 9, gamma - 10
    #
    cdef double albedo  = 0.5
    cdef double zeta    = 1.0
    for irate in range(rate.shape[0]):
        """
        Options for rate reactions
        """
        k = 0.0
        if <unsigned int> rate[irate,0] == 1:
            k   += zeta*rate[irate, 8] # alpha
        elif <unsigned int> rate[irate,0] == 2:
            k   += (zeta*rate[irate, 8]*pow( Temp / 300.0,
                rate[irate,9]) * rate[irate,10]/ (1.0-albedo))
        elif <unsigned int> rate[irate,0] == 3:
            if (-rate[irate,9]*avmag > -30.0):
                k   += (rad * rate[irate,8] * exp(-rate[irate, 10] *
                    avmag) )
            else: k += 0.0
        elif <unsigned int> rate[irate,0] == 4:
            k   += iceform(Temp)
        else:
            if (-rate[irate,10]/Temp > -30.0):
                k   += (rate[irate, 8]*pow(Temp / 300.0, rate[irate, 9])*
                    exp(-rate[irate, 10] / Temp) )
            else: k += 0.0
        """"""
        #
        # Check rates
        #
        if irate == 53: k=0.0
        if irate == 78: k=0.0
        if irate == 111: k=0.0
        if irate == 164: k=0.0
        if irate == 172: k=0.0
        if irate == 178: k=0.0
        if irate == 179: k=0.0
        if irate == 961: k=0.0
        if irate == 962: k=0.0
        if irate == 964: k=0.0
        if irate == 965: k=0.0
        if irate == 980: k=0.0
        if irate == 981: k=0.0
        if irate == 983: k=0.0
        if irate == 984: k=0.0
        if irate == 986: k=0.0
        if irate == 987: k=0.0
        if irate == 1006: k=0.0
        if irate == 1007: k=0.0
        if irate == 1008: k=0.0
        if irate == 1009: k=0.0
        if irate == 1010: k=0.0
        if irate == 1011: k=0.0
        if irate == 1040: k=0.0
        if irate == 1064: k=0.0
        if k > 1e10:
            print irate, rate[irate,0], rate[irate, 8:]
            print rate[irate,8] * pow(Temp/300.0, rate[irate,9])
            print exp(-rate[irate,10]/Temp)
            raise SystemExit
        #
        #
        #
        if (k > 1e-90) and (~np.isinf(k)):
            """
            #
            # Get the seven indices of products and reactants
            # include this rate
            #
            """
            d1  = 0.0
            d2  = 0.0
            d3  = d4 = d5 = d6 = d7 = 0.0
            ip1 =   <unsigned int> rate[irate, 1]
            ip2 =   <unsigned int> rate[irate, 2]
            ip3 =   <unsigned int> rate[irate, 3]
            ip4 =   <unsigned int> rate[irate, 4]
            ip5 =   <unsigned int> rate[irate, 5]
            ip6 =   <unsigned int> rate[irate, 6]
            ip7 =   <unsigned int> rate[irate, 7]
            #
            # destruction
            #
            if ip2 != -99:
                d1          += (-k*x[ip2] if ip3 == -99 else -k*x[ip2]*x[ip3])
                dndt[ip1]   += d1*x[ip1]
                d2          += (-k*x[ip1] if ip3 == -99 else -k*x[ip1]*x[ip3])
                dndt[ip2]   += d2*x[ip2]
                if ip3 > -99:
                    d3          += (-k*x[ip1]*x[ip2])
                    dndt[ip3]   += d3*x[ip3]
            else:
                d1          += -k
                dndt[ip1]   += d1*x[ip1]
            #
            # formation
            #
            if ip2 != -99:
                d4          += (k*x[ip1]*x[ip2] if ip3 == -99
                    else k*x[ip1]*x[ip2]*x[ip3])
                dndt[ip4]   += d4
            else:
                d4          += k*x[ip1]
                dndt[ip4]   += d4
            if ip5 > -99:
                if ip2 != -99:
                    d5          += (k*x[ip1]*x[ip2] if ip3 == -99
                        else k*x[ip1]*x[ip2]*x[ip3])
                    dndt[ip5]   += d5
                else:
                    d5          += k*x[ip1]
                    dndt[ip5]   += d5
            if ip6 > -99:
                if ip2 != -99:
                    d6          +=  (k*x[ip1]*x[ip2] if ip3 == -99
                        else k*x[ip1]*x[ip2]*x[ip3])
                    dndt[ip6]   += d6
                else:
                    d6          += k*x[ip1]
                    dndt[ip6]   += d6
            if ip7 > -99:
                d7          +=  (k*x[ip1]*x[ip2] if ip3 == -99
                    else k*x[ip1]*x[ip2]*x[ip3])
                dndt[ip7]   += d7
            """"""
    """"""
    #
    # These rates of changes are in term of seconds
    # cm^3 per seconds
    #
    #print '%2.5e %2.5e %2.5e %2.5e'%(getelectrons(x), tot[0], tot[1], tot[2])
    #print '%2.5e %2.5e %2.5e'%(xe, xh, xhe)
    #for ir in xrange(x.shape[0]):
    #if fabs(dndt[ir]) < 1e-90: dndt[ir] = 0.0
    #print '%d  %2.5e %2.5e'%(ir, x[ir], dndt[ir])
    #
    # Hydrogen is tricky
    #
    cdef double GR  = 5.2e-17*sqrt(Temp * 3.33e-3)
    dndt[1]     -=  (GR *ndens)*x[1]
    #
    # e, H2, He
    #
    dndt[67]    = 0.0
    dndt[68]    = 0.0
    dndt[69]    = 0.0
    #    raise SystemExit
    return dndt
""""""
#
# Function to call from outside
#
cpdef chemSolve1(double t, np.ndarray[DTYPE_t, ndim=1] x, list p):
    """
    Calls the more pure c code
    """
    return chemSolve(t, x, p[0], p[1], p[2], p[3], p[4], p[5])
""""""