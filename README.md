# BouligandStationarityLowRankOptimization

This repository contains the Matlab code used to realize the numerical experiment presented in [OGA24, section 6.2], which compares PGD [JKMW23, Algorithm 3.1], P2GD [SU15, Algorithm 3], P2GDR [OGA24, Algorithm 3], and [LKB23, Algorithm 1] on the problem introduced in [LKB23, section 2.2].

(*) Files related to the numerical experiment. The script PGDvsP2GDvsP2GDR.m calls the functions PGDinfo.m and P2GDRinfo.m, and produces the file PGDvsP2GDvsP2GDR_LKB23.mat. The script LKB23Algo1_LKB23.m calls the function LKB23Algo1info.m, and produces the file LKB23Algo1_LKB23.mat. The script PGDvsP2GDvsP2GDRvsLKB23Algo1.m produces the figures.

(*) Files related to the check of the Hessian. The script CheckHessianLKB23.m calls the functions CheckHessian.m and CheckHessianLift.m.

[BouligandStationarityLowRankOptimization.zip](https://github.com/golikier/BouligandStationarityLowRankOptimization/files/15075646/BouligandStationarityLowRankOptimization.zip)
