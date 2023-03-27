# ApocalypseFreeLowRankOptimization

This repository contains the Matlab code used to realize the numerical experiment presented in [OGA23, section 6.2.7], which compares PGD [JKMW22, Algorithm 3.1], P2GD [SU15, Algorithm 3], P2GDR [OGA23, Algorithms 5.3 & 6.3], and [LKB22, Algorithm 1] on the problem introduced in [LKB22, section 2.2].

(*) Files related to the numerical experiment. The script PGDvsP2GDvsP2GDR.m calls the functions PGDinfo.m and P2GDRinfo.m, and produces the file PGDvsP2GDvsP2GDR_LKB22.mat. The script LKB22Algo1_LKB22.m calls the function LKB22Algo1info.m, and produces the file LKB22Algo1_LKB22.mat. The script PGDvsP2GDvsP2GDRvsLKB22Algo1.m produces the figures.

(*) Files related to the check of the Hessian. The script CheckHessianLKB22.m calls the functions CheckHessian.m and CheckHessianLift.m.

[ApocalypseFreeLowRankOptimization.zip](https://github.com/golikier/ApocalypseFreeLowRankOptimization/files/11082550/ApocalypseFreeLowRankOptimization.zip)
