The files for this project are organized into three R scripts:

PlotSR.zip contains the following files:
PlotSR.R is the main script for plotting the stock and recruitment time series and stock-recruitment model fits
Read.ASAP3.rep.file.R is a generic script for reading ASAP report files:
RUN14.rep (2022 assessment)
2019_HKW_UNIT_MODEL_BASE.rep
2017_HKW_UNIT_MOD.rep
2015_HKW_UNIT_MOD.rep
dynBH.rep contains the results of the dynamic Beverton-Holt model
the image of this script is saved in WhiteHake.RData

CalculateReferencePoints.zip contains the following files:
CalculateReferencePoints.R is the main script. It compiles the necessary life-history information and saves is as WHT_HAKE_NE.RData
Read.ASAP3.dat.file.R is a generic script for reading ASAP dat files
2022_HKW_UNIT_FINAL_ASAP_REV.dat
WHT_HAKE_NE.RData contains the life-history parameters for calulating in the reference points

POPSIM.R runs stochastic projections.
It reads from WHT_HAKE.NE.RData and RUN14.rep