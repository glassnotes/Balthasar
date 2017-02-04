#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# lbmle.py: A class implementing least-biased MLE tomography
# according to Rehacek et. al, PRA 92, 052303 (2015)
#                                                                                  
# © 2017 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Balthasar.                                      
# Licensed under BSD-3-Clause                                                      
# 


import numpy as np
from pynitefields import *

from balthasar.mubs import MUBs

from math import log

class LBMLE():
    """ Class to hold routines for least-biased MLE tomography. 

        Args:
            mubs (MUBs): The set of MUBs in which the system is measured.
                         
        Class MUBs contains the following public attributes.

        Attributes:
            mubs (MUBs): The MUBs. 
            p (int): A prime number, the dimension of a single particle
            n (int): The number of particles
            dim (int): The dimension of the system dim = :math:`p^n`
            projectors (list): The projectors corresponding to the MUB vectors.

    """

    def __init__(self, mubs):
        """ Initialize a least-bias MLE container.

            Arguments:
                mubs: The set of MUBs in which the system is measured.
        """
        self.mubs = mubs
        self.p = mubs.p
        self.n = mubs.n
        self.dim = mubs.dim
        
        # Construct the projector versions of the MUB vectors.
        self.projectors = generate_projectors(self.mubs) 
    

    def estimate(self, bases, freqs, **kwargs):
        """ Perform least-bias MLE tomography. 
            
            Arguments:
                bases (list): The list of indices (as powers of the primitive element, 
                              or -1 for infinite slope) of the bases in which measurements
                              were done. 
                freqs (list): Frequency counts of the measurements done w.r.t the 
                              specified bases. 

            Keyword Args:
                mu (double): A parameter to tune the estimation. Defaults to 1e-4.
                eps (double): A parameter to tune the estimation. Defaults to 1e-4.
                
            Returns:
                An estimate for a density matrix that best fits the provided data. 
        """

        # Go through the keyword arguments and set mu and eps if required.
        mu = 1e-4
        eps = 1e-4

        if "mu" in kwargs:
            mu = kwargs["mu"]
        if "eps" in kwargs:
            eps = kwargs["eps"]

        # Separate the bases out into measured and unmeasured
        meas_bs = [projectors[x] if x in bases for x in range(self.dim + 1)]
        unmeas_bs [projectors[x] if x not in bases for x in range(self.dim + 1)]

        # Begin with the initial state, the maximally mixed state
        rho_0 = (1.0 / self.dim) * np.eye(self.dim)
        rho_n = rho_0
     
        n = 1

        # Iterate
        while (n < 1000):
            ########################################################
            #                    Compute W(rho)
            # I might eventually put this in a separate method, but
            # for now I'm going to leave it here to avoid having to
            # repeatedly pass the same (large chunk of) information 
            # to some helper function.
            ########################################################
            term_1 = np.zeros((self.dim, self.dim))
            term_2 = np.zeros((self.dim, self.dim))

            # Compute the first sum, which contains the measurement 
            # frequencies and the measured bases
            for basis_idx in range(len(meas_bs)):
                for proj_idx in range(len(meas_bs[basis_idx])):
                    this_projector = meas_bs[basis_idx][proj_idx]

                    p_num = freqs[basis_idx][proj_idx] 
                    p_denom = np.trace(np.dot(rho_n, this_projector))
                    prefactor = p_num / p_denom

                    term_1 = term_1 + (prefactor * this_projector)
            
            # Compute the second sum, which is over all the unmeasured bases.
            for basis_idx in range(len(unmeas_bs)):
                for proj_idx in range(len(unmeas_bs[basis_idx])):
                    this_projector = unmeas_bs[basis_idx][proj_idx]

                    prefactor = log(np.trace(np.dot(rho_n, this_projector)))

                    term_2 = term_2 + (prefactor * this_projector)
            
    
            # Finally, compute W(rho)
            W_rho_n = term1 - mu * term2
            ########################################################


            # Check if we've got a good estimate. If the desired accuracy 
            # is satisfied by the most recent rho_n, then we're done. 
            # Return the estimator and the number of steps.
            # If not, increment n and keep going.
            if check_accuracy(W_rho_n, rho_n):
                return rho_n, n 
            else:
                n += 1

            # Compute the next term in the series. It's a big ugly expression,
            # so I've separated out a term 'clump', and also the num/denom
            clump = W_rho_n - np.trace(np.dot(W_rho_n, rho_n)) * np.eye(self.dim)
            
            numerator = np.dot(np.eye(self.dim) + eps * clump, \
                             np.dot(rho_n, np.eye(self.dim) + eps * clump))
            denominator = 1 + (eps ** 2) * np.trace(np.dot(np.dot(clump, clump), rho_n))
            rho_np1 = numerator / denominator

            rho_n = rho_np1


    def check_accuracy(self, W_rho, rho):
        """ Determines if a state estimator is *good enough*.

        """
        tol = 1e-5
        LHS = np.dot(W_rho, rho)
        RHS = np.trace(np.dot(W_rho, rho)) * rho
        
        return all(np.isclose(LHS, RHS, atol = tol)):
