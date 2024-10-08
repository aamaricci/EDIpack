MODULE EDIPACK
  USE ED_INPUT_VARS  , only: &
       ed_read_input , &
       Norb          , &
       Nspin         , &
       Nloop         , &
       Nph           , &
       Uloc          , &
       Ust           , &
       Jh            , &
       Jx            , &
       Jp            , &
       xmu           , &
       beta          , &
       g_ph          , &
       w0_ph         , &
       eps           , &
       wini          , &
       wfin          , &
       xmin          , &
       xmax          , &
       Nsuccess      , &
       dmft_error    , &
       sb_field      , &
       cg_Scheme     , &
       nread         , &
       Lmats         , &
       Lreal         , &
       Lpos          , &
       Hfile         , &
       LOGfile       , &
       bath_type 


  USE ED_BATH, only:                                                    &
       ed_get_bath_dimension           => get_bath_dimension           ,&
       ed_get_bath_component_dimension => get_bath_component_dimension ,&
       ed_get_bath_component           => get_bath_component           ,&
       ed_set_bath_component           => set_bath_component           ,&
       ed_copy_bath_component          => copy_bath_component          ,&
       ed_spin_symmetrize_bath         => spin_symmetrize_bath         ,&
       ed_orb_symmetrize_bath          => orb_symmetrize_bath          ,&
       ed_orb_equality_bath            => orb_equality_bath            ,&
       ed_ph_symmetrize_bath           => ph_symmetrize_bath           ,&
       ed_ph_trans_bath                => ph_trans_bath                ,&
       ed_break_symmetry_bath          => break_symmetry_bath          ,&
       ed_set_Hreplica                 => set_Hreplica

  

  USE ED_AUX_FUNX, only:                        &
       ed_set_suffix                                                     , &
       ed_reset_suffix                                                   , &
       ed_search_variable

  USE ED_IO, only:                              &
       ed_get_sigma_matsubara                 , &
       ed_get_sigma_realaxis                  , &
       ed_get_gimp_matsubara                  , &
       ed_get_gimp_realaxis                   , &
       ed_get_dens                            , &
       ed_get_mag                             , &
       ed_get_docc                            , &
       ed_get_eimp                            , &
       ed_get_doubles

  USE ED_BATH_FIT,  only: ed_chi2_fitgf

  USE ED_MAIN, only:                            &
       ed_init_solver                         , &
       ed_solve                               , &
       ed_finalize_solver


END MODULE EDIPACK

