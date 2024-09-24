# BouligandStationarityLowRankOptimization

This repository contains two zip files related to the problem of finding B-stationary points in low-rank optimization.

The zip file [BouligandStationarityLowRankOptimization.zip](https://github.com/user-attachments/files/15937940/BouligandStationarityLowRankOptimization.zip) contains the Matlab code used to realize the numerical experiment presented in [OGA24, section 8.2], which compares monotone PGD [JKMW23, Algorithm 3.1 with m = 0], P2GD [SU15, Algorithm 3], P2GDR [OGA24, Algorithm 6.2], and HRTR [LKB23, Algorithm 1] on the problem introduced in [LKB23, section 2.2].

(*) Files related to the numerical experiment. The script PGDvsP2GDvsP2GDR.m calls the functions PGDinfo.m and P2GDRinfo.m, and produces the file PGDvsP2GDvsP2GDR_LKB23.mat. The script HRTR_LKB23.m calls the function HRTRinfo.m, and produces the file HRTR_LKB23.mat. The script PGDvsP2GDvsP2GDRvsHRTR.m produces the figures.

(*) Files related to the check of the Hessian. The script CheckHessianLKB23.m calls the functions CheckHessian.m and CheckHessianLift.m.

The zip file [ERFDR.zip](https://github.com/user-attachments/files/17119095/ERFDR.zip) contains Octave implementations of four subclasses of ERFDR [OA24, Algorithm 4.2]: RFDR and CRFDR with each of the three cones from [OA24, Table 6.1].

The script ERFDR.m calls the function ERFDRinfo.m on the problem introduced in [LKB23, section 2.2]. The script ERFDRplot.m produces the figures.

[OGA24] https://arxiv.org/abs/2201.03962v2

[OA24] https://arxiv.org/abs/2409.12298
