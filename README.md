# BouligandStationarityLowRankOptimization

This repository contains the Matlab code used to realize the numerical experiment presented in [OGA24, section 8.2], which compares monotone PGD [JKMW23, Algorithm 3.1 with m = 0], P2GD [SU15, Algorithm 3], P2GDR [OGA24, Algorithm 6.2], and HRTR [LKB23, Algorithm 1] on the problem introduced in [LKB23, section 2.2].

(*) Files related to the numerical experiment. The script PGDvsP2GDvsP2GDR.m calls the functions PGDinfo.m and P2GDRinfo.m, and produces the file PGDvsP2GDvsP2GDR_LKB23.mat. The script HRTR_LKB23.m calls the function HRTRinfo.m, and produces the file HRTR_LKB23.mat. The script PGDvsP2GDvsP2GDRvsHRTR.m produces the figures.

(*) Files related to the check of the Hessian. The script CheckHessianLKB23.m calls the functions CheckHessian.m and CheckHessianLift.m.

[BouligandStationarityLowRankOptimization.zip](https://github.com/user-attachments/files/15937877/BouligandStationarityLowRankOptimization.zip)
