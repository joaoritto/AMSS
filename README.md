# AMSS
Replication of AMSS (2002): Optimal Taxation without State-Contingent Debt

This code solves for a numeric example of Optimal Taxation without State-Contingent Debt. Production is given by y=zl. z takes two values z_l or z_h exogenously
(independent draw each period). Government spending also takes two value g_l and g_h exogenously. Labor choice is endogenous and government chooses the evolution of government debt
(non-state-contingent) to maximize welfare.

There are files for the solution using parallelization and without this (code is still pretty quick).
The files starting with IM solve for the incomplete markets (non-state-contingent debt) and CM for the case where debt is state contingent. The solution for the incomplete markets
case starts by solving the complete markets and uses this value function as initial guess.
