# ApocalypseFreeLowRankOptimization

This repository contains the Matlab code used to realize the numerical experiment presented in [OGA23, section 4.2.5], which compares PGD [JKMW22, Algorithm 3.1], P2GD [SU15, Algorithm 3], P2GDR [OGA23, Algorithms 3 & 6], and [LKB22, Algorithm 1] on the problem introduced in [LKB22, section 2.2].

(*) Files related to the numerical experiment.
The script PGDvsP2GDvsP2GDR.m calls the functions PGDinfo.m and P2GDRinfo.m, and produces the file PGDvsP2GDvsP2GDR_LKB22.mat.
The script LKB22Algo1.m calls the function LKB22Algo1info.m, and produces the file LKB22Algo1_LKB22.mat.

(*) Files related to the check of the Hessian.
The script CheckHessianLKB22.m calls the functions CheckHessian.m and CheckHessianLift.m.

[ApocalypseFreeLowRankOptimization.zip](https://github.com/golikier/ApocalypseFreeLowRankOptimization/files/10406833/ApocalypseFreeLowRankOptimization.zip)
