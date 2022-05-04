setup:
	julia scripts/00_downloaddata.jl
	julia scripts/01_preprocessdata.jl
	julia scripts/02_makesimulationdata.jl
	Rscript scripts/03_calculatehomoplasyonMSAs.R