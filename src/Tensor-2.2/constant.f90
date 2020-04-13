module constant
	implicit none
	real*8,parameter::const_eVJ=1.602176565*1d-19!1eV=1.602176565*1d-19 J
	real*8,parameter::const_muN_eV=3.1524512605*1d-8!muN=3.1524512605*1d-8 eV,nuclear magneton
	real*8,parameter::const_h_bar_eV=6.58211928*1d-16!eV s,Planck constant over 2 pi in eV s
	
	real*8,parameter::const_h_bar=1.054571726*1d-34!J s, Planck constant over 2 pi in J s
	real*8,parameter::const_mu0=12.566370614*1d-7!N A^-2,magnetic constant in SI units
	real*8,parameter::const_muB=927.400968*1d-26!J T^-1,Bohr magneton
	real*8,parameter::const_muN=5.05078353*1d-27!J T^-1,nuclear magneton
	real*8,parameter::const_muI=1.439*const_muN!J T^-1,nuclear magnetic moment
	real*8,parameter::const_a0=0.52917721092*1d-10 !m,Bohr radius
	real*8,parameter::const_ge=-2.00231930436153!the g factor of the electron
	real*8,parameter::const_gI=5.585694713!the g factor of the nuclear
	real*8,parameter::const_pi=3.14159265358979323846
	real*8,parameter::const_gamma_n=1.83247179*1d8 !s^-1 T^-1 neutron gyromagnetic ratio 
	real*8,parameter::const_gamma_p=2.675222005*1d8!s^-1 T^-1 proton gyromagnetic ratio 
	real*8,parameter::const_gamma_e=1.760859708*1d11 !s^-1 T^-1 electron gyromagnetic ratio 
end module
