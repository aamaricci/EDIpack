MODULE ED_INPUT_VARS
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str
  USE ED_VERSION
  implicit none

  !input variables
  !=========================================================
  integer              :: Nbath               !Nbath=# of bath sites (per orbital or not depending on bath_type)
  integer              :: Norb                !Norb =# of impurity orbitals
  integer              :: Nspin               !Nspin=# spin degeneracy (max 2)

  integer              :: nloop               !max dmft loop variables
  integer              :: Nph                 !max number of phonons allowed (cut off)
  real(8),dimension(5) :: Uloc                !local interactions
  real(8)              :: Ust                 !intra-orbitals interactions
  real(8)              :: Jh                  !J_Hund: Hunds' coupling constant 
  real(8)              :: Jx                  !J_X: coupling constant for the spin-eXchange interaction term
  real(8)              :: Jp                  !J_P: coupling constant for the Pair-hopping interaction term 
  real(8)              :: xmu                 !chemical potential
  real(8)              :: beta                !inverse temperature
  !
  integer              :: ph_type             !shape of the e part of the e-ph interaction: 1=orbital occupation, 2=orbital hybridization
  real(8),dimension(10):: g_ph                !g_ph: electron-phonon coupling constant
  real(8)              :: w0_ph               !w0_ph: phonon frequency (constant)
  !
  integer              :: Nsuccess            !# of repeated success to fall below convergence threshold  
  real(8)              :: dmft_error          !dmft convergence threshold
  real(8)              :: eps                 !broadening
  real(8)              :: wini,wfin           !frequency range
  real(8)              :: xmin,xmax           !x-range for the local lattice probability distribution function (phonons)
  logical              :: HFmode              !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)              :: cutoff              !cutoff for spectral summation
  real(8)              :: gs_threshold        !Energy threshold for ground state degeneracy loop up
  real(8)              :: sb_field            !symmetry breaking field
  !
  logical              :: chispin_flag        !evaluate spin susceptibility
  logical              :: chidens_flag        !evaluate dens susceptibility
  logical              :: chipair_flag        !evaluate pair susceptibility
  logical              :: chiexct_flag        !evaluate excitonic susceptibility
  !
  logical              :: ed_finite_temp      !flag to select finite temperature method. note that if T then lanc_nstates_total must be > 1 
  logical              :: ed_sparse_H         !flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--
  logical              :: ed_total_ud         !flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw
  logical              :: ed_solve_offdiag_gf !flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal 
  logical              :: ed_print_Sigma      !flag to print impurity Self-energies
  logical              :: ed_print_G          !flag to print impurity Green`s functions
  logical              :: ed_print_G0         !flag to print impurity non-interacting Green`s functions
  logical              :: ed_all_G            !flag to evaluate all the components of the impurity Green`s functions irrespective of the symmetries
  logical              :: ed_twin             !flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.
  logical              :: ed_sectors          !flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift
  integer              :: ed_sectors_shift    !shift to the ed_sectors scan
  integer              :: ed_verbose          !
  real(8)              :: ed_offset_bath      !half-bandwidth for the bath initialization: flat in -hwband:hwband
  real(8)              :: ed_hw_bath          !half-bandwidth for the bath initialization: flat in -hwband:hwband

  !
  character(len=12)    :: lanc_method         !select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only) 
  real(8)              :: lanc_tolerance      !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos. 
  integer              :: lanc_niter          !Max number of Lanczos iterations
  integer              :: lanc_ngfiter        !Max number of iteration in resolvant tri-diagonalization
  integer              :: lanc_ncv_factor     !Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_ncv_add        !Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_nstates_sector !Max number of required eigenvalues per sector
  integer              :: lanc_nstates_total  !Max number of states hold in the finite T calculation
  integer              :: lanc_nstates_step   !Number of states added at each step to determine the optimal spectrum size at finite T
  integer              :: lanc_dim_threshold  !Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.
  !
  character(len=5)     :: cg_Scheme           !fit scheme: delta (default), weiss for G0
  integer              :: cg_method           !fit routine type:0=CGnr (default), 1=minimize (old f77)
  integer              :: cg_grad             !gradient evaluation: 0=analytic, 1=numeric
  integer              :: cg_Niter            !Max number of iteration in the fit
  real(8)              :: cg_Ftol             !Tolerance in the cg fit
  integer              :: cg_stop             !fit stop condition:0-3, 0=C1.AND.C2, 1=C1, 2=C2 with C1=|F_n-1 -F_n|<tol*(1+F_n), C2=||x_n-1 -x_n||<tol*(1+||x_n||).
  integer              :: cg_Weight           !CGfit mode 0=1, 1=1/n , 2=1/w_n weight
  integer              :: cg_pow              !fit power for the calculation of the Chi distance function as |G0 - G0and|**cg_pow
  logical              :: cg_minimize_ver     !flag to pick old (Krauth) or new (Lichtenstein) version of the minimize CG routine
  real(8)              :: cg_minimize_hh      !unknown parameter used in the CG minimize procedure.  
  !
  logical              :: finiteT             !flag for finite temperature calculation
  character(len=7)     :: bath_type           !flag to set bath type: normal (1bath/imp), hybrid(1bath)
  !
  real(8)              :: nread               !fixed density. if 0.d0 fixed chemical potential calculation.
  real(8)              :: nerr                !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)              :: ndelta              !initial chemical potential step
  real(8)              :: ncoeff              !multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff)
  integer              :: niter               !


  !Some parameters for function dimension:
  !=========================================================
  integer              :: Lmats
  integer              :: Lreal
  integer              :: Lfit
  integer              :: Ltau
  integer              :: Lpos

  !LOG AND Hamiltonian UNITS
  !=========================================================
  character(len=100)   :: Hfile,HLOCfile,SectorFile
  integer,save         :: LOGfile

  !THIS IS JUST A RELOCATED GLOBAL VARIABLE
  character(len=200)                                 :: ed_input_file=""


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine ed_read_input(INPUTunit,comm)
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    integer,optional :: comm
    logical          :: master=.true.
    integer          :: i,rank=0
#ifdef _MPI
    if(present(comm))then
       master=get_Master_MPI(comm)
       rank  =get_Rank_MPI(comm)
    endif
#endif
    !
    !Store the name of the input file:
    ed_input_file=str(INPUTunit)
    !
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of impurity orbitals (max 5).")
    call parse_input_variable(Nbath,"NBATH",INPUTunit,default=6,comment="Number of bath sites:(normal=>Nbath per orb)(hybrid=>Nbath total)(replica=>Nbath=Nreplica)")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy (max 2)")
    call parse_input_variable(Nph,"NPH",INPUTunit,default=0,comment="Max number of phonons allowed (cut off)")
    call parse_input_variable(bath_type,"BATH_TYPE",INPUTunit,default='normal',comment="flag to set bath type: normal (1bath/imp), hybrid(1bath), replica(1replica/imp)")
    !
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[2d0,0d0,0d0,0d0,0d0],comment="Values of the local interaction per orbital (max 5)")
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="S-E coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="P-H coupling")
    !
    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0,comment="Inverse temperature, at T=0 is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")
    !
    call parse_input_variable(ph_type,"PH_TYPE",INPUTunit,default=1,comment="Shape e-ph interaction: 1=orbital density, 2=orbital hybridization")
    call parse_input_variable(g_ph,"G_PH",INPUTunit,default=[0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0],comment="Electron-phonon coupling constant")
    call parse_input_variable(w0_ph,"W0_PH",INPUTunit,default=0.d0,comment="Phonon frequency")
    call parse_input_variable(nloop,"NLOOP",INPUTunit,default=100,comment="Max number of DMFT iterations.")
    call parse_input_variable(dmft_error,"DMFT_ERROR",INPUTunit,default=0.00001d0,comment="Error threshold for DMFT convergence")
    call parse_input_variable(nsuccess,"NSUCCESS",INPUTunit,default=1,comment="Number of successive iterations below threshold for convergence")
    call parse_input_variable(sb_field,"SB_FIELD",INPUTunit,default=0.1d0,comment="Value of a symmetry breaking field for magnetic solutions.")
    !
    call parse_input_variable(ed_finite_temp,"ED_FINITE_TEMP",INPUTunit,default=.false.,comment="flag to select finite temperature method. note that if T then lanc_nstates_total must be > 1")
    call parse_input_variable(ed_twin,"ED_TWIN",INPUTunit,default=.false.,comment="flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.")
    call parse_input_variable(ed_sectors,"ED_SECTORS",INPUTunit,default=.false.,comment="flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift.")
    call parse_input_variable(ed_sectors_shift,"ED_SECTORS_SHIFT",INPUTunit,1,comment="shift to ed_sectors")
    call parse_input_variable(ed_sparse_H,"ED_SPARSE_H",INPUTunit,default=.true.,comment="flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE ")
    call parse_input_variable(ed_total_ud,"ED_TOTAL_UD",INPUTunit,default=.true.,comment="flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw")
    call parse_input_variable(ed_solve_offdiag_gf,"ED_SOLVE_OFFDIAG_GF",INPUTunit,default=.false.,comment="flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal") 
    call parse_input_variable(ed_print_Sigma,"ED_PRINT_SIGMA",INPUTunit,default=.true.,comment="flag to print impurity Self-energies")
    call parse_input_variable(ed_print_G,"ED_PRINT_G",INPUTunit,default=.true.,comment="flag to print impurity Greens function")
    call parse_input_variable(ed_print_G0,"ED_PRINT_G0",INPUTunit,default=.true.,comment="flag to print non-interacting impurity Greens function")
    call parse_input_variable(ed_verbose,"ED_VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 5:all. Really: all")
    call parse_input_variable(ed_hw_bath,"ed_hw_bath",INPUTunit,default=2d0,comment="half-bandwidth for the bath initialization: flat in -hwband:hwband")
    call parse_input_variable(ed_offset_bath,"ed_offset_bath",INPUTunit,default=1d-1,comment="offset for the initialization of diagonal terms in replica bath: -offset:offset")

    !
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=4096,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of real-axis frequencies.")
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1024,comment="Number of imaginary time points.")
    call parse_input_variable(Lfit,"LFIT",INPUTunit,default=1000,comment="Number of Matsubara frequencies used in the \Chi2 fit.")
    call parse_input_variable(Lpos,"LPOS",INPUTunit,default=100,comment="Number of points for the lattice PDF.")
    !
    call parse_input_variable(nread,"NREAD",INPUTunit,default=0.d0,comment="Objective density for fixed density calculations.")
    call parse_input_variable(nerr,"NERR",INPUTunit,default=1.d-4,comment="Error threshold for fixed density calculations.")
    call parse_input_variable(ndelta,"NDELTA",INPUTunit,default=0.1d0,comment="Initial step for fixed density calculations.")
    call parse_input_variable(ncoeff,"NCOEFF",INPUTunit,default=1d0,comment="multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff).")
    !
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    call parse_input_variable(xmin,"XMIN",INPUTunit,default=-3.d0,comment="Smallest position for the lattice PDF")
    call parse_input_variable(xmax,"XMAX",INPUTunit,default=3.d0,comment="Largest position for the lattice PDF")
    !
    call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(chidens_flag,"CHIDENS_FLAG",INPUTunit,default=.false.,comment="Flag to activate density susceptibility calculation.")
    call parse_input_variable(chipair_flag,"CHIPAIR_FLAG",INPUTunit,default=.false.,comment="Flag to activate pair susceptibility calculation.")
    call parse_input_variable(chiexct_flag,"CHIEXCT_FLAG",INPUTunit,default=.false.,comment="Flag to activate excitonis susceptibility calculation.")
    !
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    !    
    call parse_input_variable(lanc_method,"LANC_METHOD",INPUTunit,default="arpack",comment="select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only), DVDSON (no MPI)")
    call parse_input_variable(lanc_nstates_sector,"LANC_NSTATES_SECTOR",INPUTunit,default=2,comment="Initial number of states per sector to be determined.")
    call parse_input_variable(lanc_nstates_total,"LANC_NSTATES_TOTAL",INPUTunit,default=1,comment="Initial number of total states to be determined.")
    call parse_input_variable(lanc_nstates_step,"LANC_NSTATES_STEP",INPUTunit,default=2,comment="Number of states added to the spectrum at each step.")
    call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",INPUTunit,default=10,comment="Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",INPUTunit,default=0,comment="Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512,comment="Number of Lanczos iteration in spectrum determination.")
    call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200,comment="Number of Lanczos iteration in GF determination. Number of momenta.")
    call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",INPUTunit,default=1d-18,comment="Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.")
    call parse_input_variable(lanc_dim_threshold,"LANC_DIM_THRESHOLD",INPUTunit,default=1024,comment="Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.")
    !
    call parse_input_variable(cg_method,"CG_METHOD",INPUTunit,default=0,comment="Conjugate-Gradient method: 0=NR, 1=minimize.")
    call parse_input_variable(cg_grad,"CG_GRAD",INPUTunit,default=0,comment="Gradient evaluation method: 0=analytic (default), 1=numeric.")
    call parse_input_variable(cg_ftol,"CG_FTOL",INPUTunit,default=0.00001d0,comment="Conjugate-Gradient tolerance.")
    call parse_input_variable(cg_stop,"CG_STOP",INPUTunit,default=0,comment="Conjugate-Gradient stopping condition: 0-3, 0=C1.AND.C2, 1=C1, 2=C2 with C1=|F_n-1 -F_n|<tol*(1+F_n), C2=||x_n-1 -x_n||<tol*(1+||x_n||).")
    call parse_input_variable(cg_niter,"CG_NITER",INPUTunit,default=500,comment="Max. number of Conjugate-Gradient iterations.")
    call parse_input_variable(cg_weight,"CG_WEIGHT",INPUTunit,default=1,comment="Conjugate-Gradient weight form: 1=1.0, 2=1/n , 3=1/w_n.")
    call parse_input_variable(cg_scheme,"CG_SCHEME",INPUTunit,default='weiss',comment="Conjugate-Gradient fit scheme: delta or weiss.")
    call parse_input_variable(cg_pow,"CG_POW",INPUTunit,default=2,comment="Fit power for the calculation of the Chi distance function as 1/L*|G0 - G0and|**cg_pow ")
    call parse_input_variable(cg_minimize_ver,"CG_MINIMIZE_VER",INPUTunit,default=.false.,comment="Flag to pick old/.false. (Krauth) or new/.true. (Lichtenstein) version of the minimize CG routine")
    call parse_input_variable(cg_minimize_hh,"CG_MINIMIZE_HH",INPUTunit,default=1d-4,comment="Unknown parameter used in the CG minimize procedure.")
    !
    call parse_input_variable(SectorFile,"SectorFile",INPUTunit,default="sectors",comment="File where to retrieve/store the sectors contributing to the spectrum.")
    call parse_input_variable(Hfile,"Hfile",INPUTunit,default="hamiltonian",comment="File where to retrieve/store the bath parameters.")
    call parse_input_variable(HLOCfile,"HLOCfile",INPUTunit,default="inputHLOC.in",comment="File read the input local H.")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")


#ifdef _MPI
    if(present(comm))then
       if(.not.master)then
          LOGfile=1000-rank
          open(LOGfile,file="stdOUT.rank"//str(rank)//".ed")
          do i=1,get_Size_MPI(comm)
             if(i==rank)write(*,"(A,I0,A,I0)")"Rank ",rank," writing to unit: ",LOGfile
          enddo
       endif
    endif
#endif
    !
    !
    Ltau=max(int(beta),Ltau)
    if(master)then
       call print_input()
       call save_input(INPUTunit)
       call scifor_version()
       call code_version(version)
    endif
    !Act on the input variable only after printing.
    !In the new parser variables are hard-linked into the list:
    !any change to the variable is immediately copied into the list... (if you delete .ed it won't be printed out)
    call substring_delete(Hfile,".restart")
    call substring_delete(Hfile,".ed")
  end subroutine ed_read_input




  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop


END MODULE ED_INPUT_VARS
