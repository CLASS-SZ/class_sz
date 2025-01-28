/** @file szpowerspectrum.h Documented includes for sz module. */
#ifndef __SZ__
#define __SZ__

#include "common.h"
#include "lensing.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_sf_expint.h"
#include "gsl/gsl_sf_lambert.h"
# include <fftw3.h>
// # include <gsl_integration.h>
// #include "fft.h"


#define _pk_at_z_1h_ ((pclass_sz->has_pk_at_z_1h == _TRUE_) && (index_md == pclass_sz->index_md_pk_at_z_1h))
#define _pk_at_z_2h_ ((pclass_sz->has_pk_at_z_2h == _TRUE_) && (index_md == pclass_sz->index_md_pk_at_z_2h))
#define _pk_gg_at_z_1h_ ((pclass_sz->has_pk_gg_at_z_1h == _TRUE_) && (index_md == pclass_sz->index_md_pk_gg_at_z_1h))
#define _pk_gg_at_z_2h_ ((pclass_sz->has_pk_gg_at_z_2h == _TRUE_) && (index_md == pclass_sz->index_md_pk_gg_at_z_2h))
#define _pk_bb_at_z_1h_ ((pclass_sz->has_pk_bb_at_z_1h == _TRUE_) && (index_md == pclass_sz->index_md_pk_bb_at_z_1h))
#define _pk_bb_at_z_2h_ ((pclass_sz->has_pk_bb_at_z_2h == _TRUE_) && (index_md == pclass_sz->index_md_pk_bb_at_z_2h))
#define _pk_b_at_z_2h_ ((pclass_sz->has_pk_b_at_z_2h == _TRUE_) && (index_md == pclass_sz->index_md_pk_b_at_z_2h))
#define _pk_em_at_z_1h_ ((pclass_sz->has_pk_em_at_z_1h == _TRUE_) && (index_md == pclass_sz->index_md_pk_em_at_z_1h))
#define _pk_em_at_z_2h_ ((pclass_sz->has_pk_em_at_z_2h == _TRUE_) && (index_md == pclass_sz->index_md_pk_em_at_z_2h))
#define _pk_HI_at_z_1h_ ((pclass_sz->has_pk_HI_at_z_1h == _TRUE_) && (index_md == pclass_sz->index_md_pk_HI_at_z_1h))
#define _pk_HI_at_z_2h_ ((pclass_sz->has_pk_HI_at_z_2h == _TRUE_) && (index_md == pclass_sz->index_md_pk_HI_at_z_2h))
#define _bk_at_z_1h_ ((pclass_sz->has_bk_at_z_1h == _TRUE_) && (index_md == pclass_sz->index_md_bk_at_z_1h))
#define _bk_at_z_2h_ ((pclass_sz->has_bk_at_z_2h == _TRUE_) && (index_md == pclass_sz->index_md_bk_at_z_2h))
#define _bk_at_z_3h_ ((pclass_sz->has_bk_at_z_3h == _TRUE_) && (index_md == pclass_sz->index_md_bk_at_z_3h))
#define _bk_ttg_at_z_1h_ ((pclass_sz->has_bk_ttg_at_z_1h == _TRUE_) && (index_md == pclass_sz->index_md_bk_ttg_at_z_1h))
#define _bk_ttg_at_z_2h_ ((pclass_sz->has_bk_ttg_at_z_2h == _TRUE_) && (index_md == pclass_sz->index_md_bk_ttg_at_z_2h))
#define _bk_ttg_at_z_3h_ ((pclass_sz->has_bk_ttg_at_z_3h == _TRUE_) && (index_md == pclass_sz->index_md_bk_ttg_at_z_3h))

//#define _bk_at_z_hf_ ((pclass_sz->has_bk_at_z_hf == _TRUE_) && (index_md == pclass_sz->index_md_bk_at_z_hf))
#define _mean_y_ ((pclass_sz->has_mean_y == _TRUE_) && (index_md == pclass_sz->index_md_mean_y))
#define _cib_monopole_ ((pclass_sz->has_cib_monopole == _TRUE_) && (index_md == pclass_sz->index_md_cib_monopole))
#define _cib_shotnoise_ ((pclass_sz->has_cib_shotnoise == _TRUE_) && (index_md == pclass_sz->index_md_cib_shotnoise))
#define _dcib0dz_ ((pclass_sz->has_dcib0dz == _TRUE_) && (index_md == pclass_sz->index_md_dcib0dz))
#define _dydz_ ((pclass_sz->has_dydz == _TRUE_) && (index_md == pclass_sz->index_md_dydz))
#define _hmf_ ((pclass_sz->has_hmf == _TRUE_) && (index_md == pclass_sz->index_md_hmf))
#define _tSZ_power_spectrum_ ((pclass_sz->has_sz_ps == _TRUE_) && (index_md == pclass_sz->index_md_sz_ps))
#define _trispectrum_ ((pclass_sz->has_sz_trispec == _TRUE_) && (index_md == pclass_sz->index_md_trispectrum))
#define _2halo_ ((pclass_sz->has_sz_2halo == _TRUE_) && (index_md == pclass_sz->index_md_2halo))
#define _te_y_y_ ((pclass_sz->has_sz_te_y_y == _TRUE_) && (index_md == pclass_sz->index_md_te_y_y))
#define _m_y_y_1h_ ((pclass_sz->has_sz_m_y_y_1h == _TRUE_) && (index_md == pclass_sz->index_md_m_y_y_1h))
#define _m_y_y_2h_ ((pclass_sz->has_sz_m_y_y_2h == _TRUE_) && (index_md == pclass_sz->index_md_m_y_y_2h))
#define _cov_Y_N_ ((pclass_sz->has_sz_cov_Y_N == _TRUE_) && (index_md == pclass_sz->index_md_cov_Y_N))
#define _cov_N_N_ ((pclass_sz->has_sz_cov_N_N == _TRUE_) && (index_md == pclass_sz->index_md_cov_N_N))
#define _cov_N_N_hsv_ ((pclass_sz->has_sz_cov_N_N_hsv == _TRUE_) && (index_md == pclass_sz->index_md_cov_N_N_hsv))
#define _cov_Y_Y_ssc_ ((pclass_sz->has_sz_cov_Y_Y_ssc == _TRUE_) && (index_md == pclass_sz->index_md_cov_Y_Y_ssc))
#define _cov_Y_N_next_order_ ((pclass_sz->has_sz_cov_Y_N_next_order == _TRUE_) && (index_md == pclass_sz->index_md_cov_Y_N_next_order))
#define _kSZ_kSZ_gal_1h_ ((pclass_sz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gal_1h))
#define _kSZ_kSZ_gal_1h_fft_ ((pclass_sz->has_kSZ_kSZ_gal_1h_fft == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gal_1h_fft))
#define _kSZ_kSZ_gal_2h_fft_ ((pclass_sz->has_kSZ_kSZ_gal_2h_fft == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gal_2h_fft))
#define _kSZ_kSZ_gal_3h_fft_ ((pclass_sz->has_kSZ_kSZ_gal_3h_fft == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gal_3h_fft))
#define _kSZ_kSZ_gallens_1h_fft_ ((pclass_sz->has_kSZ_kSZ_gallens_1h_fft == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gallens_1h_fft))
#define _kSZ_kSZ_gallens_2h_fft_ ((pclass_sz->has_kSZ_kSZ_gallens_2h_fft == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gallens_2h_fft))
#define _kSZ_kSZ_gallens_3h_fft_ ((pclass_sz->has_kSZ_kSZ_gallens_3h_fft == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gallens_3h_fft))
#define _kSZ_kSZ_lens_1h_fft_ ((pclass_sz->has_kSZ_kSZ_lens_1h_fft == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_lens_1h_fft))
#define _kSZ_kSZ_lens_2h_fft_ ((pclass_sz->has_kSZ_kSZ_lens_2h_fft == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_lens_2h_fft))
#define _kSZ_kSZ_lens_3h_fft_ ((pclass_sz->has_kSZ_kSZ_lens_3h_fft == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_lens_3h_fft))
#define _gal_gal_lens_1h_fft_ ((pclass_sz->has_gal_gal_lens_1h_fft == _TRUE_) && (index_md == pclass_sz->index_md_gal_gal_lens_1h_fft))
#define _gal_gal_lens_2h_fft_ ((pclass_sz->has_gal_gal_lens_2h_fft == _TRUE_) && (index_md == pclass_sz->index_md_gal_gal_lens_2h_fft))
#define _gal_gal_lens_3h_fft_ ((pclass_sz->has_gal_gal_lens_3h_fft == _TRUE_) && (index_md == pclass_sz->index_md_gal_gal_lens_3h_fft))
#define _kSZ_kSZ_gal_2h_ ((pclass_sz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gal_2h))
#define _kSZ_kSZ_gal_3h_ ((pclass_sz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gal_3h))
#define _kSZ_kSZ_gal_hf_ ((pclass_sz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gal_hf))
#define _kSZ_kSZ_gallens_hf_ ((pclass_sz->has_kSZ_kSZ_gallens_hf == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_gallens_hf))
#define _kSZ_kSZ_lens_hf_ ((pclass_sz->has_kSZ_kSZ_lens_hf == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_lens_hf))
#define _kSZ_kSZ_lensmag_1halo_ ((pclass_sz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_lensmag_1halo))
#define _gal_gal_1h_ ((pclass_sz->has_gal_gal_1h == _TRUE_) && (index_md == pclass_sz->index_md_gal_gal_1h))
#define _gal_gal_2h_ ((pclass_sz->has_gal_gal_2h == _TRUE_) && (index_md == pclass_sz->index_md_gal_gal_2h))
#define _n5k_ ((pclass_sz->has_n5k == _TRUE_) && (index_md == pclass_sz->index_md_n5k))
#define _gal_gal_hf_ ((pclass_sz->has_gal_gal_hf == _TRUE_) && (index_md == pclass_sz->index_md_gal_gal_hf))
#define _tau_gal_1h_ ((pclass_sz->has_tau_gal_1h == _TRUE_) && (index_md == pclass_sz->index_md_tau_gal_1h))
#define _tau_gal_2h_ ((pclass_sz->has_tau_gal_2h == _TRUE_) && (index_md == pclass_sz->index_md_tau_gal_2h))
#define _tau_tau_1h_ ((pclass_sz->has_tau_tau_1h == _TRUE_) && (index_md == pclass_sz->index_md_tau_tau_1h))
#define _tau_tau_2h_ ((pclass_sz->has_tau_tau_2h == _TRUE_) && (index_md == pclass_sz->index_md_tau_tau_2h))
#define _gal_lens_2h_ ((pclass_sz->has_gal_lens_2h == _TRUE_) && (index_md == pclass_sz->index_md_gal_lens_2h))
#define _gal_lens_hf_ ((pclass_sz->has_gal_lens_hf == _TRUE_) && (index_md == pclass_sz->index_md_gal_lens_hf))
#define _gal_lens_1h_ ((pclass_sz->has_gal_lens_1h == _TRUE_) && (index_md == pclass_sz->index_md_gal_lens_1h))
#define _gal_lensmag_2h_ ((pclass_sz->has_gal_lensmag_2h == _TRUE_) && (index_md == pclass_sz->index_md_gal_lensmag_2h))
#define _gal_lensmag_1h_ ((pclass_sz->has_gal_lensmag_1h == _TRUE_) && (index_md == pclass_sz->index_md_gal_lensmag_1h))
#define _gal_gallens_2h_ ((pclass_sz->has_gal_gallens_2h == _TRUE_) && (index_md == pclass_sz->index_md_gal_gallens_2h))
#define _gal_gallens_1h_ ((pclass_sz->has_gal_gallens_1h == _TRUE_) && (index_md == pclass_sz->index_md_gal_gallens_1h))
#define _gallens_gallens_2h_ ((pclass_sz->has_gallens_gallens_2h == _TRUE_) && (index_md == pclass_sz->index_md_gallens_gallens_2h))
#define _gallens_gallens_1h_ ((pclass_sz->has_gallens_gallens_1h == _TRUE_) && (index_md == pclass_sz->index_md_gallens_gallens_1h))
#define _gallens_lens_2h_ ((pclass_sz->has_gallens_lens_2h == _TRUE_) && (index_md == pclass_sz->index_md_gallens_lens_2h))
#define _gallens_lens_1h_ ((pclass_sz->has_gallens_lens_1h == _TRUE_) && (index_md == pclass_sz->index_md_gallens_lens_1h))
#define _gal_lensmag_hf_ ((pclass_sz->has_gal_lensmag_hf == _TRUE_) && (index_md == pclass_sz->index_md_gal_lensmag_hf))
#define _tSZ_lensmag_2h_ ((pclass_sz->has_tSZ_lensmag_2h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_lensmag_2h))
#define _tSZ_lensmag_1h_ ((pclass_sz->has_tSZ_lensmag_1h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_lensmag_1h))
#define _gallens_lensmag_2h_ ((pclass_sz->has_gallens_lensmag_2h == _TRUE_) && (index_md == pclass_sz->index_md_gallens_lensmag_2h))
#define _gallens_lensmag_1h_ ((pclass_sz->has_gallens_lensmag_1h == _TRUE_) && (index_md == pclass_sz->index_md_gallens_lensmag_1h))
#define _lensmag_lensmag_hf_ ((pclass_sz->has_lensmag_lensmag_hf == _TRUE_) && (index_md == pclass_sz->index_md_lensmag_lensmag_hf))
#define _lensmag_lensmag_2h_ ((pclass_sz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == pclass_sz->index_md_lensmag_lensmag_2h))
#define _lensmag_lensmag_1h_ ((pclass_sz->has_lensmag_lensmag_1h == _TRUE_) && (index_md == pclass_sz->index_md_lensmag_lensmag_1h))
#define _lens_lensmag_2h_ ((pclass_sz->has_lens_lensmag_2h == _TRUE_) && (index_md == pclass_sz->index_md_lens_lensmag_2h))
#define _lens_lensmag_1h_ ((pclass_sz->has_lens_lensmag_1h == _TRUE_) && (index_md == pclass_sz->index_md_lens_lensmag_1h))
#define _lens_lensmag_hf_ ((pclass_sz->has_lens_lensmag_hf == _TRUE_) && (index_md == pclass_sz->index_md_lens_lensmag_hf))
#define _lens_lens_1h_ ((pclass_sz->has_lens_lens_1h == _TRUE_) && (index_md == pclass_sz->index_md_lens_lens_1h))
#define _custom1_custom1_1h_ ((pclass_sz->has_custom1_custom1_1h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_custom1_1h))
#define _custom1_custom1_2h_ ((pclass_sz->has_custom1_custom1_2h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_custom1_2h))
#define _custom1_lens_1h_ ((pclass_sz->has_custom1_lens_1h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_lens_1h))
#define _custom1_lens_2h_ ((pclass_sz->has_custom1_lens_2h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_lens_2h))
#define _custom1_tSZ_1h_ ((pclass_sz->has_custom1_tSZ_1h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_tSZ_1h))
#define _custom1_tSZ_2h_ ((pclass_sz->has_custom1_tSZ_2h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_tSZ_2h))
#define _custom1_gal_1h_ ((pclass_sz->has_custom1_gal_1h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_gal_1h))
#define _custom1_gal_2h_ ((pclass_sz->has_custom1_gal_2h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_gal_2h))
#define _custom1_gallens_1h_ ((pclass_sz->has_custom1_gallens_1h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_gallens_1h))
#define _custom1_gallens_2h_ ((pclass_sz->has_custom1_gallens_2h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_gallens_2h))
#define _custom1_cib_1h_ ((pclass_sz->has_custom1_cib_1h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_cib_1h))
#define _custom1_cib_2h_ ((pclass_sz->has_custom1_cib_2h == _TRUE_) && (index_md == pclass_sz->index_md_custom1_cib_2h))
#define _lens_lens_2h_ ((pclass_sz->has_lens_lens_2h == _TRUE_) && (index_md == pclass_sz->index_md_lens_lens_2h))
#define _lens_lens_hf_ ((pclass_sz->has_lens_lens_hf == _TRUE_) && (index_md == pclass_sz->index_md_lens_lens_hf))
#define _tSZ_gal_1h_ ((pclass_sz->has_tSZ_gal_1h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_gal_1h))
#define _tSZ_gal_2h_ ((pclass_sz->has_tSZ_gal_2h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_gal_2h))
#define _IA_gal_2h_ ((pclass_sz->has_IA_gal_2h == _TRUE_) && (index_md == pclass_sz->index_md_IA_gal_2h))
#define _tSZ_gallens_1h_ ((pclass_sz->has_tSZ_gallens_1h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_gallens_1h))
#define _tSZ_gallens_2h_ ((pclass_sz->has_tSZ_gallens_2h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_gallens_2h))
#define _tSZ_cib_1h_ ((pclass_sz->has_tSZ_cib_1h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_cib_1h))
#define _tSZ_cib_2h_ ((pclass_sz->has_tSZ_cib_2h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_cib_2h))
#define _gal_cib_1h_ ((pclass_sz->has_gal_cib_1h == _TRUE_) && (index_md == pclass_sz->index_md_gal_cib_1h))
#define _gal_cib_2h_ ((pclass_sz->has_gal_cib_2h == _TRUE_) && (index_md == pclass_sz->index_md_gal_cib_2h))
#define _gallens_cib_1h_ ((pclass_sz->has_gallens_cib_1h == _TRUE_) && (index_md == pclass_sz->index_md_gallens_cib_1h))
#define _gallens_cib_2h_ ((pclass_sz->has_gallens_cib_2h == _TRUE_) && (index_md == pclass_sz->index_md_gallens_cib_2h))
#define _ngal_ngal_1h_ ((pclass_sz->has_ngal_ngal_1h == _TRUE_) && (index_md == pclass_sz->index_md_ngal_ngal_1h))
#define _ngal_ngal_2h_ ((pclass_sz->has_ngal_ngal_2h == _TRUE_) && (index_md == pclass_sz->index_md_ngal_ngal_2h))
#define _ngal_ngal_hf_ ((pclass_sz->has_ngal_ngal_hf == _TRUE_) && (index_md == pclass_sz->index_md_ngal_ngal_hf))
#define _ngal_lens_1h_ ((pclass_sz->has_ngal_lens_1h == _TRUE_) && (index_md == pclass_sz->index_md_ngal_lens_1h))
#define _ngal_lens_2h_ ((pclass_sz->has_ngal_lens_2h == _TRUE_) && (index_md == pclass_sz->index_md_ngal_lens_2h))
#define _ngal_lens_hf_ ((pclass_sz->has_ngal_lens_hf == _TRUE_) && (index_md == pclass_sz->index_md_ngal_lens_hf))
#define _ngal_tsz_1h_ ((pclass_sz->has_ngal_tsz_1h == _TRUE_) && (index_md == pclass_sz->index_md_ngal_tsz_1h))
#define _ngal_tsz_2h_ ((pclass_sz->has_ngal_tsz_2h == _TRUE_) && (index_md == pclass_sz->index_md_ngal_tsz_2h))
#define _ngal_gallens_1h_ ((pclass_sz->has_ngal_gallens_1h == _TRUE_) && (index_md == pclass_sz->index_md_ngal_gallens_1h))
#define _ngal_gallens_2h_ ((pclass_sz->has_ngal_gallens_2h == _TRUE_) && (index_md == pclass_sz->index_md_ngal_gallens_2h))
#define _ngal_IA_2h_ ((pclass_sz->has_ngal_IA_2h == _TRUE_) && (index_md == pclass_sz->index_md_ngal_IA_2h))
#define _nlensmag_gallens_1h_ ((pclass_sz->has_nlensmag_gallens_1h == _TRUE_) && (index_md == pclass_sz->index_md_nlensmag_gallens_1h))
#define _nlensmag_gallens_2h_ ((pclass_sz->has_nlensmag_gallens_2h == _TRUE_) && (index_md == pclass_sz->index_md_nlensmag_gallens_2h))
#define _nlensmag_tsz_1h_ ((pclass_sz->has_nlensmag_tsz_1h == _TRUE_) && (index_md == pclass_sz->index_md_nlensmag_tsz_1h))
#define _nlensmag_tsz_2h_ ((pclass_sz->has_nlensmag_tsz_2h == _TRUE_) && (index_md == pclass_sz->index_md_nlensmag_tsz_2h))
#define _cib_cib_1h_ ((pclass_sz->has_cib_cib_1h == _TRUE_) && (index_md == pclass_sz->index_md_cib_cib_1h))
#define _cib_cib_2h_ ((pclass_sz->has_cib_cib_2h == _TRUE_) && (index_md == pclass_sz->index_md_cib_cib_2h))
#define _lens_cib_1h_ ((pclass_sz->has_lens_cib_1h == _TRUE_) && (index_md == pclass_sz->index_md_lens_cib_1h))
#define _lens_cib_2h_ ((pclass_sz->has_lens_cib_2h == _TRUE_) && (index_md == pclass_sz->index_md_lens_cib_2h))
#define _tSZ_lens_1h_ ((pclass_sz->has_tSZ_lens_1h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_lens_1h))
#define _tSZ_lens_2h_ ((pclass_sz->has_tSZ_lens_2h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_lens_2h))
#define _isw_lens_ ((pclass_sz->has_isw_lens == _TRUE_) && (index_md == pclass_sz->index_md_isw_lens))
#define _isw_tsz_ ((pclass_sz->has_isw_tsz == _TRUE_) && (index_md == pclass_sz->index_md_isw_tsz))
#define _isw_auto_ ((pclass_sz->has_isw_auto == _TRUE_) && (index_md == pclass_sz->index_md_isw_auto))
#define _dndlnM_ ((pclass_sz->has_dndlnM == _TRUE_) && (index_md == pclass_sz->index_md_dndlnM))
#define _szrates_ ((pclass_sz->has_sz_rates == _TRUE_) && (index_md == pclass_sz->index_md_szrates))
#define _kSZ_kSZ_1h_ ((pclass_sz->has_kSZ_kSZ_1h == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_1h))
#define _kSZ_kSZ_2h_ ((pclass_sz->has_kSZ_kSZ_2h == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_2h))
#define _kSZ_kSZ_tSZ_1h_ ((pclass_sz->has_kSZ_kSZ_tSZ_1h == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_tSZ_1h))
#define _kSZ_kSZ_tSZ_2h_ ((pclass_sz->has_kSZ_kSZ_tSZ_2h == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_tSZ_2h))
#define _kSZ_kSZ_tSZ_3h_ ((pclass_sz->has_kSZ_kSZ_tSZ_3h == _TRUE_) && (index_md == pclass_sz->index_md_kSZ_kSZ_tSZ_3h))
#define _tSZ_tSZ_tSZ_1halo_ ((pclass_sz->has_tSZ_tSZ_tSZ_1halo == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_tSZ_tSZ_1halo))
#define _tSZ_tSZ_tSZ_2h_ ((pclass_sz->has_tSZ_tSZ_tSZ_2h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_tSZ_tSZ_2h))
#define _tSZ_tSZ_tSZ_3h_ ((pclass_sz->has_tSZ_tSZ_tSZ_3h == _TRUE_) && (index_md == pclass_sz->index_md_tSZ_tSZ_tSZ_3h))
//#define _tSZ_trispectrum_ ((pclass_sz->has_sz_trispec == _TRUE_))
//#define _tSZ_2halo_ ((pclass_sz->has_sz_2halo == _TRUE_))
//#define _tSZ_te_y_y_ ((pclass_sz->has_sz_te_y_y == _TRUE_))
// #define _cov_N_Cl_ ((pclass_sz->has_sz_cov_N_Cl == _TRUE_))





struct class_sz_structure {



  #define __ALLOCATE_TSZ_PARAMETER__
  #include "class_sz_precisions.h"
  #undef __ALLOCATE_TSZ_PARAMETER__


  fftw_plan forward_plan, reverse_plan;
  fftw_plan forward_plan_counts_fft, reverse_plan_counts_fft;
  // int N_samp_fftw;


  int use_analytical_truncated_nfw;
  int use_hod; // Eq. 15 or 16 of KA20
  int unwise_galaxy_sample_id;
  int galaxy_sample;


  int no_b2;
  //double unwise_m_min_cut;

  int cib_nu0_norm;
  double sn_cutoff;

  double f_sky;
  double fsky_from_skyfracs;
  double szunbinned_loglike;
  double Omega_survey;

  double chi_star; //comoving distance to the surface of last scattering [Mpc/h]

  double hmf_int;
  double y_monopole;
  double * cib_monopole;
  double * cib_shotnoise;
  double * pk_at_z_1h;
  double * pk_at_z_2h;
  double * pk_gg_at_z_1h;
  double * pk_gg_at_z_2h;
  double * pk_bb_at_z_1h;
  double * pk_bb_at_z_2h;
  double * pk_b_at_z_2h;
  double * pk_em_at_z_1h;
  double * pk_em_at_z_2h;
  double * pk_HI_at_z_1h;
  double * pk_HI_at_z_2h;
  double * bk_at_z_1h;
  double * bk_at_z_2h;
  double * bk_at_z_3h;
  double * bk_ttg_at_z_1h;
  double * bk_ttg_at_z_2h;
  double * bk_ttg_at_z_3h;
  double * cl_sz_1h;
  double * cl_gal_gal_1h;
  double * cl_gal_gal_2h;
  double * cl_gal_gal_hf;
  double * cl_gal_lens_hf;
  double * cl_lens_lens_hf;
  double * cl_tau_gal_2h;
  double * cl_tau_gal_1h;
  double * cl_tau_tau_2h;
  double * cl_tau_tau_1h;
  double * cl_gal_lens_2h;
  double * cl_gal_lens_1h;
  double * cl_gal_lensmag_hf;
  double * cl_gal_lensmag_2h;
  double * cl_gal_gallens_1h;
  double * cl_gal_gallens_2h;
  double * cl_gallens_gallens_1h;
  double * cl_gallens_gallens_2h;
  double * cl_gallens_lens_1h;
  double * cl_gallens_lens_2h;
  double * thetas_arcmin;
  double * gamma_gal_gallens_1h;
  double * gamma_gal_gallens_2h;
  double * cl_gal_lensmag_1h;
  double * cl_tSZ_lensmag_2h;
  double * cl_tSZ_lensmag_1h;
  double * cl_gallens_lensmag_2h;
  double * cl_gallens_lensmag_1h;
  double * cl_lensmag_lensmag_hf;
  double * cl_lensmag_lensmag_2h;
  double * cl_lensmag_lensmag_1h;
  double * cl_lens_lensmag_hf;
  double * cl_lens_lensmag_2h;
  double * cl_lens_lensmag_1h;
  double * cl_custom1_custom1_1h;
  double * cl_custom1_custom1_2h;
  double * cl_custom1_lens_1h;
  double * cl_custom1_lens_2h;
  double * cl_custom1_tSZ_1h;
  double * cl_custom1_tSZ_2h;
  double ** cl_custom1_cib_1h;
  double ** cl_custom1_cib_2h;
  double * cl_custom1_gal_1h;
  double * cl_custom1_gal_2h;
  double * cl_custom1_gallens_1h;
  double * cl_custom1_gallens_2h;
  double * cl_lens_lens_1h;
  double * cl_lens_lens_2h;
  double * cl_tSZ_gal_1h;
  double * cl_tSZ_gal_2h;
  double * cl_IA_gal_2h;
  double * cl_tSZ_gallens_1h;
  double * cl_tSZ_gallens_2h;
  double *** cl_ngal_ngal_1h;
  double *** cl_ngal_ngal_2h;
  double *** cl_ngal_ngal_hf;
  double ** cl_ngal_lens_1h;
  double ** cl_ngal_lens_2h;
  double ** cl_ngal_lens_hf;
  double ** cl_ngal_tsz_1h;
  double ** cl_ngal_tsz_2h;
  double ** cl_ngal_gallens_1h;
  double ** cl_ngal_gallens_2h;
  double ** cl_ngal_IA_2h;
  double ** cl_nlensmag_gallens_1h;
  double ** cl_nlensmag_gallens_2h;
  double ** cl_nlensmag_tsz_1h;
  double ** cl_nlensmag_tsz_2h;
  double *** cl_cib_cib_1h;
  double *** cl_cib_cib_2h;
  double ** cl_tSZ_cib_1h;
  double ** cl_tSZ_cib_2h;
  double ** cl_gal_cib_1h;
  double ** cl_gal_cib_2h;
  double ** cl_gallens_cib_1h;
  double ** cl_gallens_cib_2h;
  double ** cl_lens_cib_1h;
  double ** cl_lens_cib_2h;
  double * cl_tSZ_lens_1h;
  double * cl_tSZ_lens_2h;
  double * szrate;
  double * cl_isw_lens;
  double * cl_isw_tsz;
  double * cl_isw_auto;
  double * cov_ll_kSZ_kSZ_gal;
  double * cl_t2t2f;
  double * cl_kSZ_kSZ_gal_lensing_term;
  double * cl_kSZ_kSZ_gal_1h;
  double * cl_kSZ_kSZ_gal_1h_fft;
  double * cl_kSZ_kSZ_gal_2h_fft;
  double * cl_kSZ_kSZ_gal_3h_fft;
  double * cl_kSZ_kSZ_gallens_1h_fft;
  double * cl_kSZ_kSZ_gallens_2h_fft;
  double * cl_kSZ_kSZ_gallens_3h_fft;
  double * cl_kSZ_kSZ_gallens_hf;
  double * cov_ll_kSZ_kSZ_gallens;
  double * cl_kSZ_kSZ_gallens_lensing_term;
  double * cl_kSZ_kSZ_lens_1h_fft;
  double * cl_kSZ_kSZ_lens_2h_fft;
  double * cl_kSZ_kSZ_lens_3h_fft;
  double * cl_gal_gal_lens_1h_fft;
  double * cl_gal_gal_lens_2h_fft;
  double * cl_gal_gal_lens_3h_fft;
  double * cl_kSZ_kSZ_lens_hf;
  double * cov_ll_kSZ_kSZ_lens;
  double * cl_kSZ_kSZ_lens_lensing_term;
  double * cl_kSZ_kSZ_gal_2h;
  double * cl_kSZ_kSZ_gal_3h;
  double * cl_kSZ_kSZ_gal_hf;
  double * cl_kSZ_kSZ_lensmag_1h;
  double * b_tSZ_tSZ_tSZ_1halo;
  double * b_tSZ_tSZ_tSZ_2h;
  double * b_tSZ_tSZ_tSZ_3h;
  double * cl_kSZ_kSZ_1h;
  double * cl_kSZ_kSZ_2h;
  double * b_kSZ_kSZ_tSZ_1h;
  double * b_kSZ_kSZ_tSZ_2h;
  double * b_kSZ_kSZ_tSZ_3h;
  double * cl_te_y_y;
  double * m_y_y_1h;
  double * m_y_y_2h;
  double ** tllprime_sz;
  double ** dndlnM_at_z_and_M;
  double * cl_sz_2h;
  double ** cov_N_cl;
  double ** r_N_cl;
  double ** cov_Y_N;
  double ** cov_Y_N_next_order;
  double ** cov_Y_Y_ssc;
  double * cov_N_N;
  double ** cov_N_N_hsv;
  double ** r_Y_N;
  double ** r_cl_clp;
  double ** trispectrum_ref;
  double * cov_cl_cl;
  double * sig_cl_squared_binned;

  int profile_matter_density;
  double matter_nfw_power_law_index;


  int delta_def_galaxies;
  int delta_def_cib;
  int delta_def_matter_density;
  int delta_def_electron_pressure;
  int delta_def_electron_density;
  int delta_def_custom1;


  int delta_def_HI_pressure;
  int delta_def_HI_density;

  int bispec_conf_id;

  double M_min_ng_bar;
  double M_max_ng_bar;


  int index_d_tot;
  int index_phi;
  int index_psi;
  int number_of_titles;

  int need_m200c_to_mvir;
  int need_m200m_to_mvir;
  int need_m200m_to_m200c;
  int need_m200c_to_m200m;
  int need_m200m_to_m500c;
  int need_hmf;
  int need_sigma;
  int need_m200c_to_m500c;
  int need_m500c_to_m200c;

  int need_ng_bias;
  int nz_ng_bias;
  int nk_ng_bias;
  double * array_ln_1pz_ng_bias;
  double * array_ln_k_ng_bias;
  double * array_ln_ng_bias_at_z_and_k;

  double * array_ln_density_norm_at_z_and_m;

  double * array_ln_matter_density_norm_at_z_and_m;

  int need_ksz_template;
  int need_tt_noise;
  int need_lensing_noise;


  int integrate_wrt_mvir;
  int integrate_wrt_m500c;
  int integrate_wrt_m200m;
  int integrate_wrt_m200c;

  int has_electron_pressure;
  int has_electron_density;
  int has_HI_density;
  int has_galaxy;
  int has_matter_density;
  int has_lensing;
  int has_cib;
  int has_dcib0dz;
  int has_dydz;
  int has_isw;

  int has_vir;
  int has_500c;
  int has_200m;
  int has_200c;


  int index_integrate_wrt_mvir;
  int index_integrate_wrt_m500c;
  int index_integrate_wrt_m200m;

  int index_has_electron_pressure;
  int index_has_electron_density;
  int index_has_HI_density;
  int index_has_galaxy;
  int index_has_matter_density;
  int index_has_lensing;
  int index_has_cib;
  int index_has_isw;
  int index_has_custom1;

  int index_has_vir;
  int index_has_500c;
  int index_has_200m;
  int index_has_200c;

  int index_md;

  int has_sz_counts;
  int has_sz_counts_fft;

  int create_ref_trispectrum_for_cobaya;


  double alpha_break_pressure;
  double M_break_pressure;
  int use_broken_pressure;


  int use_m500c_in_ym_relation;
  int use_m200c_in_ym_relation;
  //int has_sz_te_y_y;
  int has_sz_cov_N_Cl;

  int has_sz_cov_Y_N;
  int index_md_cov_Y_N;
  int index_integrand_id_cov_Y_N_first;
  int index_integrand_id_cov_Y_N_last;

  int has_sz_cov_Y_N_next_order;
  int index_md_cov_Y_N_next_order;
  int index_integrand_id_cov_Y_N_next_order_first;
  int index_integrand_id_cov_Y_N_next_order_last;

  int has_sz_cov_Y_Y_ssc;
  int index_md_cov_Y_Y_ssc;
  int index_integrand_id_cov_Y_Y_ssc_first;
  int index_integrand_id_cov_Y_Y_ssc_last;

  int has_sz_cov_N_N;
  int index_md_cov_N_N;
  int index_integrand_id_cov_N_N_first;
  int index_integrand_id_cov_N_N_last;

  int has_sz_cov_N_N_hsv;
  int index_md_cov_N_N_hsv;
  int index_integrand_id_cov_N_N_hsv_first;
  int index_integrand_id_cov_N_N_hsv_last;


  int has_sz_rates;
  int index_md_szrates;
  int index_integrand_id_szrates_first;
  int index_integrand_id_szrates_last;

  int has_hmf;
  int index_md_hmf;
  int index_integrand_id_hmf;

  int has_pk_bb_at_z_1h;
  int index_md_pk_bb_at_z_1h;
  int index_integrand_id_pk_bb_at_z_1h_first;
  int index_integrand_id_pk_bb_at_z_1h_last;

  int has_pk_bb_at_z_2h;
  int index_md_pk_bb_at_z_2h;
  int index_integrand_id_pk_bb_at_z_2h_first;
  int index_integrand_id_pk_bb_at_z_2h_last;


  int has_gas_pressure_profile_2h;
  int has_gas_density_profile_2h;

  int has_pk_b_at_z_2h;
  int index_md_pk_b_at_z_2h;
  int index_integrand_id_pk_b_at_z_2h_first;
  int index_integrand_id_pk_b_at_z_2h_last;

  int has_pk_em_at_z_1h;
  int index_md_pk_em_at_z_1h;
  int index_integrand_id_pk_em_at_z_1h_first;
  int index_integrand_id_pk_em_at_z_1h_last;

  int has_pk_em_at_z_2h;
  int index_md_pk_em_at_z_2h;
  int index_integrand_id_pk_em_at_z_2h_first;
  int index_integrand_id_pk_em_at_z_2h_last;


  int has_pk_HI_at_z_1h;
  int index_md_pk_HI_at_z_1h;
  int index_integrand_id_pk_HI_at_z_1h_first;
  int index_integrand_id_pk_HI_at_z_1h_last;

  int has_pk_HI_at_z_2h;
  int index_md_pk_HI_at_z_2h;
  int index_integrand_id_pk_HI_at_z_2h_first;
  int index_integrand_id_pk_HI_at_z_2h_last;



  int has_pk_gg_at_z_1h;
  int index_md_pk_gg_at_z_1h;
  int index_integrand_id_pk_gg_at_z_1h_first;
  int index_integrand_id_pk_gg_at_z_1h_last;

  int has_pk_gg_at_z_2h;
  int index_md_pk_gg_at_z_2h;
  int index_integrand_id_pk_gg_at_z_2h_first;
  int index_integrand_id_pk_gg_at_z_2h_last;

  int has_pk_at_z_1h;
  int index_md_pk_at_z_1h;
  int index_integrand_id_pk_at_z_1h_first;
  int index_integrand_id_pk_at_z_1h_last;

  int has_pk_at_z_2h;
  int index_md_pk_at_z_2h;
  int index_integrand_id_pk_at_z_2h_first;
  int index_integrand_id_pk_at_z_2h_last;

  int has_bk_at_z_1h;
  int index_md_bk_at_z_1h;
  int index_integrand_id_bk_at_z_1h_first;
  int index_integrand_id_bk_at_z_1h_last;

  int has_bk_at_z_2h;
  int index_md_bk_at_z_2h;
  int index_integrand_id_bk_at_z_2h_first;
  int index_integrand_id_bk_at_z_2h_last;

  int has_bk_at_z_3h;
  int index_md_bk_at_z_3h;
  int index_integrand_id_bk_at_z_3h_first;
  int index_integrand_id_bk_at_z_3h_last;

  int has_bk_ttg_at_z_1h;
  int index_md_bk_ttg_at_z_1h;
  int index_integrand_id_bk_ttg_at_z_1h_first;
  int index_integrand_id_bk_ttg_at_z_1h_last;

  int has_bk_ttg_at_z_2h;
  int index_md_bk_ttg_at_z_2h;
  int index_integrand_id_bk_ttg_at_z_2h_first;
  int index_integrand_id_bk_ttg_at_z_2h_last;

  int has_bk_ttg_at_z_3h;
  int index_md_bk_ttg_at_z_3h;
  int index_integrand_id_bk_ttg_at_z_3h_first;
  int index_integrand_id_bk_ttg_at_z_3h_last;



  int has_bk_at_z_hf;
  int has_bk_ttg_at_z_hf;

  int has_mean_galaxy_bias;
  int has_ng_in_bh;
  double fNL;

  int has_mean_y;
  int index_md_mean_y;
  int index_integrand_id_mean_y;

  int has_sz_ps;
  int index_md_sz_ps;
  int index_integrand_id_sz_ps_first;
  int index_integrand_id_sz_ps_last;

  int has_sz_2halo;
  int index_md_2halo;
  int index_integrand_id_sz_ps_2halo_first;
  int index_integrand_id_sz_ps_2halo_last;

  int has_sz_te_y_y;
  int index_md_te_y_y;
  int index_integrand_id_sz_ps_te_y_y_first;
  int index_integrand_id_sz_ps_te_y_y_last;

  int has_sz_m_y_y_1h;
  int index_md_m_y_y_1h;
  int index_integrand_id_sz_ps_m_y_y_1h_first;
  int index_integrand_id_sz_ps_m_y_y_1h_last;

  int has_sz_m_y_y_2h;
  int index_md_m_y_y_2h;
  int index_integrand_id_sz_ps_m_y_y_2h_first;
  int index_integrand_id_sz_ps_m_y_y_2h_last;

  int has_kSZ_kSZ_gal_1h;
  int index_md_kSZ_kSZ_gal_1h;
  int index_integrand_id_kSZ_kSZ_gal_1h_first;
  int index_integrand_id_kSZ_kSZ_gal_1h_last;

  int has_kSZ_kSZ_gal_1h_fft;
  int index_md_kSZ_kSZ_gal_1h_fft;
  int index_integrand_id_kSZ_kSZ_gal_1h_fft_first;
  int index_integrand_id_kSZ_kSZ_gal_1h_fft_last;

  int has_kSZ_kSZ_gal_2h_fft;
  int index_md_kSZ_kSZ_gal_2h_fft;
  int index_integrand_id_kSZ_kSZ_gal_2h_fft_first;
  int index_integrand_id_kSZ_kSZ_gal_2h_fft_last;

  int has_kSZ_kSZ_gal_3h_fft;
  int index_md_kSZ_kSZ_gal_3h_fft;
  int index_integrand_id_kSZ_kSZ_gal_3h_fft_first;
  int index_integrand_id_kSZ_kSZ_gal_3h_fft_last;


  int has_kSZ_kSZ_gal_2h;
  int index_md_kSZ_kSZ_gal_2h;
  int index_integrand_id_kSZ_kSZ_gal_2h_first;
  int index_integrand_id_kSZ_kSZ_gal_2h_last;

  int has_kSZ_kSZ_gal_3h;
  int index_md_kSZ_kSZ_gal_3h;
  int index_integrand_id_kSZ_kSZ_gal_3h_first;
  int index_integrand_id_kSZ_kSZ_gal_3h_last;

  int has_kSZ_kSZ_gal_covmat;
  int has_kSZ_kSZ_gallens_covmat;
  int has_kSZ_kSZ_lens_covmat;
  int has_kSZ_kSZ_gal_lensing_term;
  int has_kSZ_kSZ_gallens_lensing_term;
  int has_kSZ_kSZ_lens_lensing_term;

  int has_kSZ_kSZ_gal_hf;
  int index_md_kSZ_kSZ_gal_hf;
  int index_integrand_id_kSZ_kSZ_gal_hf_first;
  int index_integrand_id_kSZ_kSZ_gal_hf_last;

  int has_kSZ_kSZ_gallens_1h_fft;
  int index_md_kSZ_kSZ_gallens_1h_fft;
  int index_integrand_id_kSZ_kSZ_gallens_1h_fft_first;
  int index_integrand_id_kSZ_kSZ_gallens_1h_fft_last;

  int has_kSZ_kSZ_gallens_2h_fft;
  int index_md_kSZ_kSZ_gallens_2h_fft;
  int index_integrand_id_kSZ_kSZ_gallens_2h_fft_first;
  int index_integrand_id_kSZ_kSZ_gallens_2h_fft_last;

  int has_kSZ_kSZ_gallens_3h_fft;
  int index_md_kSZ_kSZ_gallens_3h_fft;
  int index_integrand_id_kSZ_kSZ_gallens_3h_fft_first;
  int index_integrand_id_kSZ_kSZ_gallens_3h_fft_last;

  int has_kSZ_kSZ_gallens_hf;
  int index_md_kSZ_kSZ_gallens_hf;
  int index_integrand_id_kSZ_kSZ_gallens_hf_first;
  int index_integrand_id_kSZ_kSZ_gallens_hf_last;

  int has_kSZ_kSZ_lens_1h_fft;
  int index_md_kSZ_kSZ_lens_1h_fft;
  int index_integrand_id_kSZ_kSZ_lens_1h_fft_first;
  int index_integrand_id_kSZ_kSZ_lens_1h_fft_last;

  int has_kSZ_kSZ_lens_2h_fft;
  int index_md_kSZ_kSZ_lens_2h_fft;
  int index_integrand_id_kSZ_kSZ_lens_2h_fft_first;
  int index_integrand_id_kSZ_kSZ_lens_2h_fft_last;

  int has_kSZ_kSZ_lens_3h_fft;
  int index_md_kSZ_kSZ_lens_3h_fft;
  int index_integrand_id_kSZ_kSZ_lens_3h_fft_first;
  int index_integrand_id_kSZ_kSZ_lens_3h_fft_last;

  int has_gal_gal_lens_1h_fft;
  int index_md_gal_gal_lens_1h_fft;
  int index_integrand_id_gal_gal_lens_1h_fft_first;
  int index_integrand_id_gal_gal_lens_1h_fft_last;

  int has_gal_gal_lens_2h_fft;
  int index_md_gal_gal_lens_2h_fft;
  int index_integrand_id_gal_gal_lens_2h_fft_first;
  int index_integrand_id_gal_gal_lens_2h_fft_last;

  int has_gal_gal_lens_3h_fft;
  int index_md_gal_gal_lens_3h_fft;
  int index_integrand_id_gal_gal_lens_3h_fft_first;
  int index_integrand_id_gal_gal_lens_3h_fft_last;

  int has_kSZ_kSZ_lens_hf;
  int index_md_kSZ_kSZ_lens_hf;
  int index_integrand_id_kSZ_kSZ_lens_hf_first;
  int index_integrand_id_kSZ_kSZ_lens_hf_last;


  int has_kSZ_kSZ_lensmag_1halo;
  int index_md_kSZ_kSZ_lensmag_1halo;
  int index_integrand_id_kSZ_kSZ_lensmag_1halo_first;
  int index_integrand_id_kSZ_kSZ_lensmag_1halo_last;




  int has_kSZ_kSZ_1h;
  int index_md_kSZ_kSZ_1h;
  int index_integrand_id_kSZ_kSZ_1h_first;
  int index_integrand_id_kSZ_kSZ_1h_last;

  int has_kSZ_kSZ_2h;
  int index_md_kSZ_kSZ_2h;
  int index_integrand_id_kSZ_kSZ_2h_first;
  int index_integrand_id_kSZ_kSZ_2h_last;


  int has_kSZ_kSZ_tSZ_1h;
  int index_md_kSZ_kSZ_tSZ_1h;
  int index_integrand_id_kSZ_kSZ_tSZ_1h_first;
  int index_integrand_id_kSZ_kSZ_tSZ_1h_last;

  int has_kSZ_kSZ_tSZ_2h;
  int index_md_kSZ_kSZ_tSZ_2h;
  int index_integrand_id_kSZ_kSZ_tSZ_2h_first;
  int index_integrand_id_kSZ_kSZ_tSZ_2h_last;

  int has_kSZ_kSZ_tSZ_3h;
  int index_md_kSZ_kSZ_tSZ_3h;
  int index_integrand_id_kSZ_kSZ_tSZ_3h_first;
  int index_integrand_id_kSZ_kSZ_tSZ_3h_last;

  int has_tSZ_tSZ_tSZ_1halo;
  int index_md_tSZ_tSZ_tSZ_1halo;
  int index_integrand_id_tSZ_tSZ_tSZ_1halo_first;
  int index_integrand_id_tSZ_tSZ_tSZ_1halo_last;

  int has_tSZ_tSZ_tSZ_2h;
  int index_md_tSZ_tSZ_tSZ_2h;
  int index_integrand_id_tSZ_tSZ_tSZ_2h_first;
  int index_integrand_id_tSZ_tSZ_tSZ_2h_last;

  int has_tSZ_tSZ_tSZ_3h;
  int index_md_tSZ_tSZ_tSZ_3h;
  int index_integrand_id_tSZ_tSZ_tSZ_3h_first;
  int index_integrand_id_tSZ_tSZ_tSZ_3h_last;



  int has_tSZ_lens_1h;
  int index_md_tSZ_lens_1h;
  int index_integrand_id_tSZ_lens_1h_first;
  int index_integrand_id_tSZ_lens_1h_last;

  int has_gal_gal_1h;
  int index_md_gal_gal_1h;
  int index_integrand_id_gal_gal_1h_first;
  int index_integrand_id_gal_gal_1h_last;

  int has_gal_gal_2h;
  int index_md_gal_gal_2h;
  int index_integrand_id_gal_gal_2h_first;
  int index_integrand_id_gal_gal_2h_last;

  int has_gal_gal_hf;
  int index_md_gal_gal_hf;
  int index_integrand_id_gal_gal_hf_first;
  int index_integrand_id_gal_gal_hf_last;

  int has_n5k;
  int index_md_n5k;

  int has_gal_lens_hf;
  int index_md_gal_lens_hf;
  int index_integrand_id_gal_lens_hf_first;
  int index_integrand_id_gal_lens_hf_last;

  int has_lens_lens_hf;
  int index_md_lens_lens_hf;
  int index_integrand_id_lens_lens_hf_first;
  int index_integrand_id_lens_lens_hf_last;


  int has_tau_gal_2h;
  int index_md_tau_gal_2h;
  int index_integrand_id_tau_gal_2h_first;
  int index_integrand_id_tau_gal_2h_last;

  int has_tau_gal_1h;
  int index_md_tau_gal_1h;
  int index_integrand_id_tau_gal_1h_first;
  int index_integrand_id_tau_gal_1h_last;

  int has_tau_tau_2h;
  int index_md_tau_tau_2h;
  int index_integrand_id_tau_tau_2h_first;
  int index_integrand_id_tau_tau_2h_last;

  int has_tau_tau_1h;
  int index_md_tau_tau_1h;
  int index_integrand_id_tau_tau_1h_first;
  int index_integrand_id_tau_tau_1h_last;


  int has_gal_lens_2h;
  int index_md_gal_lens_2h;
  int index_integrand_id_gal_lens_2h_first;
  int index_integrand_id_gal_lens_2h_last;

  int has_gal_lens_1h;
  int index_md_gal_lens_1h;
  int index_integrand_id_gal_lens_1h_first;
  int index_integrand_id_gal_lens_1h_last;

  int has_gal_lensmag_hf;
  int index_md_gal_lensmag_hf;
  int index_integrand_id_gal_lensmag_hf_first;
  int index_integrand_id_gal_lensmag_hf_last;

  int has_gal_lensmag_2h;
  int index_md_gal_lensmag_2h;
  int index_integrand_id_gal_lensmag_2h_first;
  int index_integrand_id_gal_lensmag_2h_last;

  int has_gal_lensmag_1h;
  int index_md_gal_lensmag_1h;
  int index_integrand_id_gal_lensmag_1h_first;
  int index_integrand_id_gal_lensmag_1h_last;

  int convert_cls_to_gamma;
  int has_gal_gallens_2h;
  int index_md_gal_gallens_2h;
  int index_integrand_id_gal_gallens_2h_first;
  int index_integrand_id_gal_gallens_2h_last;

  int has_gal_gallens_1h;
  int index_md_gal_gallens_1h;
  int index_integrand_id_gal_gallens_1h_first;
  int index_integrand_id_gal_gallens_1h_last;


  int has_gallens_gallens_2h;
  int index_md_gallens_gallens_2h;
  int index_integrand_id_gallens_gallens_2h_first;
  int index_integrand_id_gallens_gallens_2h_last;

  int has_gallens_gallens_1h;
  int index_md_gallens_gallens_1h;
  int index_integrand_id_gallens_gallens_1h_first;
  int index_integrand_id_gallens_gallens_1h_last;


  int has_gallens_lens_2h;
  int index_md_gallens_lens_2h;
  int index_integrand_id_gallens_lens_2h_first;
  int index_integrand_id_gallens_lens_2h_last;

  int has_gallens_lens_1h;
  int index_md_gallens_lens_1h;
  int index_integrand_id_gallens_lens_1h_first;
  int index_integrand_id_gallens_lens_1h_last;

  int has_gallens_lensmag_2h;
  int index_md_gallens_lensmag_2h;
  int index_integrand_id_gallens_lensmag_2h_first;
  int index_integrand_id_gallens_lensmag_2h_last;

  int has_gallens_lensmag_1h;
  int index_md_gallens_lensmag_1h;
  int index_integrand_id_gallens_lensmag_1h_first;
  int index_integrand_id_gallens_lensmag_1h_last;

  int has_tSZ_lensmag_2h;
  int index_md_tSZ_lensmag_2h;
  int index_integrand_id_tSZ_lensmag_2h_first;
  int index_integrand_id_tSZ_lensmag_2h_last;

  int has_tSZ_lensmag_1h;
  int index_md_tSZ_lensmag_1h;
  int index_integrand_id_tSZ_lensmag_1h_first;
  int index_integrand_id_tSZ_lensmag_1h_last;

  int has_nlensmag_nlensmag_hf;
  int index_md_nlensmag_nlensmag_hf;
  int index_integrand_id_nlensmag_nlensmag_hf_first;
  int index_integrand_id_nlensmag_nlensmag_hf_last;


  int has_lensmag_lensmag_hf;
  int index_md_lensmag_lensmag_hf;
  int index_integrand_id_lensmag_lensmag_hf_first;
  int index_integrand_id_lensmag_lensmag_hf_last;

  int has_lensmag_lensmag_2h;
  int index_md_lensmag_lensmag_2h;
  int index_integrand_id_lensmag_lensmag_2h_first;
  int index_integrand_id_lensmag_lensmag_2h_last;

  int has_lensmag_lensmag_1h;
  int index_md_lensmag_lensmag_1h;
  int index_integrand_id_lensmag_lensmag_1h_first;
  int index_integrand_id_lensmag_lensmag_1h_last;


  int has_lens_lensmag_hf;
  int index_md_lens_lensmag_hf;
  int index_integrand_id_lens_lensmag_hf_first;
  int index_integrand_id_lens_lensmag_hf_last;


  int has_lens_nlensmag_hf;
  int index_md_lens_nlensmag_hf;
  int index_integrand_id_lens_nlensmag_hf_first;
  int index_integrand_id_lens_nlensmag_hf_last;


  int has_lens_lensmag_2h;
  int index_md_lens_lensmag_2h;
  int index_integrand_id_lens_lensmag_2h_first;
  int index_integrand_id_lens_lensmag_2h_last;

  int has_lens_lensmag_1h;
  int index_md_lens_lensmag_1h;
  int index_integrand_id_lens_lensmag_1h_first;
  int index_integrand_id_lens_lensmag_1h_last;

  int has_lens_lens_1h;
  int index_md_lens_lens_1h;
  int index_integrand_id_lens_lens_1h_first;
  int index_integrand_id_lens_lens_1h_last;

  int has_custom1;
  int has_b_custom1;

  int has_custom1_custom1_1h;
  int index_md_custom1_custom1_1h;
  int index_integrand_id_custom1_custom1_1h_first;
  int index_integrand_id_custom1_custom1_1h_last;

  int has_custom1_custom1_2h;
  int index_md_custom1_custom1_2h;
  int index_integrand_id_custom1_custom1_2h_first;
  int index_integrand_id_custom1_custom1_2h_last;

  int has_custom1_lens_1h;
  int index_md_custom1_lens_1h;
  int index_integrand_id_custom1_lens_1h_first;
  int index_integrand_id_custom1_lens_1h_last;

  int has_custom1_lens_2h;
  int index_md_custom1_lens_2h;
  int index_integrand_id_custom1_lens_2h_first;
  int index_integrand_id_custom1_lens_2h_last;

  int has_custom1_tSZ_1h;
  int index_md_custom1_tSZ_1h;
  int index_integrand_id_custom1_tSZ_1h_first;
  int index_integrand_id_custom1_tSZ_1h_last;

  int has_custom1_tSZ_2h;
  int index_md_custom1_tSZ_2h;
  int index_integrand_id_custom1_tSZ_2h_first;
  int index_integrand_id_custom1_tSZ_2h_last;

  int has_custom1_cib_1h;
  int index_md_custom1_cib_1h;
  int index_integrand_id_custom1_cib_1h_first;
  int index_integrand_id_custom1_cib_1h_last;

  int has_custom1_cib_2h;
  int index_md_custom1_cib_2h;
  int index_integrand_id_custom1_cib_2h_first;
  int index_integrand_id_custom1_cib_2h_last;

  int has_custom1_gal_1h;
  int index_md_custom1_gal_1h;
  int index_integrand_id_custom1_gal_1h_first;
  int index_integrand_id_custom1_gal_1h_last;

  int has_custom1_gal_2h;
  int index_md_custom1_gal_2h;
  int index_integrand_id_custom1_gal_2h_first;
  int index_integrand_id_custom1_gal_2h_last;

  int has_custom1_gallens_1h;
  int index_md_custom1_gallens_1h;
  int index_integrand_id_custom1_gallens_1h_first;
  int index_integrand_id_custom1_gallens_1h_last;

  int has_custom1_gallens_2h;
  int index_md_custom1_gallens_2h;
  int index_integrand_id_custom1_gallens_2h_first;
  int index_integrand_id_custom1_gallens_2h_last;


  int has_lens_lens_2h;
  int index_md_lens_lens_2h;
  int index_integrand_id_lens_lens_2h_first;
  int index_integrand_id_lens_lens_2h_last;

  int has_tSZ_gal_1h;
  int index_md_tSZ_gal_1h;
  int index_integrand_id_tSZ_gal_1h_first;
  int index_integrand_id_tSZ_gal_1h_last;

  int has_tSZ_gal_2h;
  int index_md_tSZ_gal_2h;
  int index_integrand_id_tSZ_gal_2h_first;
  int index_integrand_id_tSZ_gal_2h_last;

  int has_IA_gal_2h;
  int index_md_IA_gal_2h;
  int index_integrand_id_IA_gal_2h_first;
  int index_integrand_id_IA_gal_2h_last;

  int has_tSZ_gallens_1h;
  int index_md_tSZ_gallens_1h;
  int index_integrand_id_tSZ_gallens_1h_first;
  int index_integrand_id_tSZ_gallens_1h_last;

  int has_tSZ_gallens_2h;
  int index_md_tSZ_gallens_2h;
  int index_integrand_id_tSZ_gallens_2h_first;
  int index_integrand_id_tSZ_gallens_2h_last;

  int has_gallens_cib_1h;
  int index_md_gallens_cib_1h;
  int index_integrand_id_gallens_cib_1h_first;
  int index_integrand_id_gallens_cib_1h_last;

  int has_gallens_cib_2h;
  int index_md_gallens_cib_2h;
  int index_integrand_id_gallens_cib_2h_first;
  int index_integrand_id_gallens_cib_2h_last;

  int has_gal_cib_1h;
  int index_md_gal_cib_1h;
  int index_integrand_id_gal_cib_1h_first;
  int index_integrand_id_gal_cib_1h_last;

  int has_gal_cib_2h;
  int index_md_gal_cib_2h;
  int index_integrand_id_gal_cib_2h_first;
  int index_integrand_id_gal_cib_2h_last;

  int has_tSZ_cib_1h;
  int index_md_tSZ_cib_1h;
  int index_integrand_id_tSZ_cib_1h_first;
  int index_integrand_id_tSZ_cib_1h_last;

  int has_tSZ_cib_2h;
  int index_md_tSZ_cib_2h;
  int index_integrand_id_tSZ_cib_2h_first;
  int index_integrand_id_tSZ_cib_2h_last;

  int has_lens_cib_1h;
  int index_md_lens_cib_1h;
  int index_integrand_id_lens_cib_1h_first;
  int index_integrand_id_lens_cib_1h_last;

  int has_lens_cib_2h;
  int index_md_lens_cib_2h;
  int index_integrand_id_lens_cib_2h_first;
  int index_integrand_id_lens_cib_2h_last;

  int has_cib_monopole;
  int index_md_cib_monopole;
  int index_integrand_id_cib_monopole_first;
  int index_integrand_id_cib_monopole_last;

  int has_cib_shotnoise;
  int index_md_cib_shotnoise;
  int index_integrand_id_cib_shotnoise_first;
  int index_integrand_id_cib_shotnoise_last;

  int has_ngal_ngal_1h;
  int index_md_ngal_ngal_1h;
  int index_integrand_id_ngal_ngal_1h_first;
  int index_integrand_id_ngal_ngal_1h_last;

  int has_ngal_ngal_2h;
  int index_md_ngal_ngal_2h;
  int index_integrand_id_ngal_ngal_2h_first;
  int index_integrand_id_ngal_ngal_2h_last;

  int has_ngal_ngal_hf;
  int index_md_ngal_ngal_hf;
  int index_integrand_id_ngal_ngal_hf_first;
  int index_integrand_id_ngal_ngal_hf_last;

  int has_ngal_lens_1h;
  int index_md_ngal_lens_1h;
  int index_integrand_id_ngal_lens_1h_first;
  int index_integrand_id_ngal_lens_1h_last;

  int has_ngal_lens_2h;
  int index_md_ngal_lens_2h;
  int index_integrand_id_ngal_lens_2h_first;
  int index_integrand_id_ngal_lens_2h_last;

  int has_ngal_lens_hf;
  int index_md_ngal_lens_hf;
  int index_integrand_id_ngal_lens_hf_first;
  int index_integrand_id_ngal_lens_hf_last;


  int has_ngal_nlensmag_hf;
  int index_md_ngal_nlensmag_hf;
  int index_integrand_id_ngal_nlensmag_hf_first;
  int index_integrand_id_ngal_nlensmag_hf_last;

  int has_ngal_tsz_1h;
  int index_md_ngal_tsz_1h;
  int index_integrand_id_ngal_tsz_1h_first;
  int index_integrand_id_ngal_tsz_1h_last;

  int has_ngal_tsz_2h;
  int index_md_ngal_tsz_2h;
  int index_integrand_id_ngal_tsz_2h_first;
  int index_integrand_id_ngal_tsz_2h_last;

  int has_nlensmag_tsz_1h;
  int index_md_nlensmag_tsz_1h;
  int index_integrand_id_nlensmag_tsz_1h_first;
  int index_integrand_id_nlensmag_tsz_1h_last;

  int has_nlensmag_tsz_2h;
  int index_md_nlensmag_tsz_2h;
  int index_integrand_id_nlensmag_tsz_2h_first;
  int index_integrand_id_nlensmag_tsz_2h_last;

  int has_ngal_gallens_1h;
  int index_md_ngal_gallens_1h;
  int index_integrand_id_ngal_gallens_1h_first;
  int index_integrand_id_ngal_gallens_1h_last;

  int has_ngal_gallens_2h;
  int index_md_ngal_gallens_2h;
  int index_integrand_id_ngal_gallens_2h_first;
  int index_integrand_id_ngal_gallens_2h_last;

  int has_nlensmag_gallens_1h;
  int index_md_nlensmag_gallens_1h;
  int index_integrand_id_nlensmag_gallens_1h_first;
  int index_integrand_id_nlensmag_gallens_1h_last;

  int has_nlensmag_gallens_2h;
  int index_md_nlensmag_gallens_2h;
  int index_integrand_id_nlensmag_gallens_2h_first;
  int index_integrand_id_nlensmag_gallens_2h_last;

  int has_ngal_IA_2h;
  int index_md_ngal_IA_2h;
  int index_integrand_id_ngal_IA_2h_first;
  int index_integrand_id_ngal_IA_2h_last;

  int has_cib_cib_1h;
  int index_md_cib_cib_1h;
  int index_integrand_id_cib_cib_1h_first;
  int index_integrand_id_cib_cib_1h_last;

  int has_cib_cib_2h;
  int index_md_cib_cib_2h;
  int index_integrand_id_cib_cib_2h_first;
  int index_integrand_id_cib_cib_2h_last;

  int has_tSZ_lens_2h;
  int index_md_tSZ_lens_2h;
  int index_integrand_id_tSZ_lens_2h_first;
  int index_integrand_id_tSZ_lens_2h_last;

  int has_isw_lens;
  int index_md_isw_lens;
  int index_integrand_id_isw_lens_first;
  int index_integrand_id_isw_lens_last;

  int has_isw_tsz;
  int index_md_isw_tsz;
  int index_integrand_id_isw_tsz_first;
  int index_integrand_id_isw_tsz_last;

  int has_isw_auto;
  int index_md_isw_auto;
  int index_integrand_id_isw_auto_first;
  int index_integrand_id_isw_auto_last;


  int has_dndlnM;
  int index_md_dndlnM;
  int index_integrand_id_dndlnM_first;
  int index_integrand_id_dndlnM_last;

  int has_sz_trispec;
  //int index_md_sz_trispec;
  int index_integrand_id_trispectrum_first;
  int index_integrand_id_trispectrum_last;
  int index_md_trispectrum;


  int number_of_integrands;
  int index_integrand;
  int index_integrand_te_y_y;
  int index_integrand_2halo_term;

  int index_integrand_trispectrum_first; //for trispectrum
  int index_integrand_trispectrum_last;  //for trispectrum

  int index_integrand_cov_N_cl_first;
  int index_integrand_cov_N_cl_last;


  int index_integrand_N_for_cov_N_cl_first;
  int index_integrand_N_for_cov_N_cl_last;


  int index_integrand_id;

  int number_of_integrals_per_thread;

  int index_integrands_first;
  int index_integrands_last;

  int index_md_dcib0dz;
  int index_md_dydz;








  //double  pk;

  // FileName root; /**< root for all file names */
  // FileName path_to_class; /**< root for all file names */
  // FileName append_name_cobaya_ref;
  // FileName path_to_ref_trispectrum_for_cobaya;
  // FileName full_path_to_noise_curve_for_y_y;
  //FileName full_path_to_dndz_gal;

 /* vector of all SZ quantities function of redshift*/

  int  tsz_size;

  int  index_flag_cov_N_cl;
  int  index_Rho_crit;
  int  index_Delta_c;
  int  index_rVIR;
  int  index_cVIR;
  int  index_c200m;
  int  index_r200m;
  int  index_mVIR;
  int  index_m500;
  int  index_r500;
  int  index_l500;
  int  index_ls;
  int  index_rs;
  int  index_m200;
  int  index_m180m;
  int  index_m200m;
  int  index_m1600m;
  int  index_m500c;
  int  index_mass_for_hmf;
  int  index_mass_for_custom1;
  int  index_mass_for_galaxies;
  int  index_mass_for_cib;
  int  index_mass_for_matter_density;
  int  index_mass_for_electron_pressure;
  int  index_mass_for_electron_density;
  int  index_mass_for_HI_pressure;
  int  index_mass_for_HI_density;
  int  index_concentration_for_galaxies;
  int  index_concentration_for_custom1;
  int  index_concentration_for_cib;
  int  index_concentration_for_matter_density;
  int  index_concentration_for_electron_pressure;
  int  index_concentration_for_electron_density;
  int  index_concentration_for_HI_pressure;
  int  index_concentration_for_HI_density;
  int  index_radius_for_galaxies;
  int  index_radius_for_cib;
  int  index_radius_for_custom1;
  int  index_radius_for_matter_density;
  int  index_radius_for_electron_pressure;
  int  index_radius_for_electron_density;
  int  index_radius_for_HI_pressure;
  int  index_radius_for_HI_density;
  int  index_r500c;
  int  index_Rh;
  int  index_mf;
  int  index_dlognudlogRh;
  int  index_lognu;
  int  index_dlogSigma2dlogRh;
  int  index_dndlogRh;
  int  index_logSigma2;
  int  index_z;
  int  index_c200c;
  int  index_m200c;
  int  index_l200c;
  int  index_characteristic_multipole_for_nfw_profile;
  int  index_r200c;
  int  index_multipole;
  int  index_szrate;
  int  index_multipole_prime;
  int  index_mass_bin_1;
  int  index_mass_bin_2;
  int  index_multipole_1;
  int  index_multipole_2;
  int  index_multipole_3;
  int  index_redshift_for_dndlnM;
  int  index_mass_for_dndlnM;
  int  index_multipole_for_pressure_profile;
  int  index_pressure_profile;
  int  index_multipole_for_tau_profile;
  int  index_multipole_for_nfw_profile;
  int  index_tau_profile;
  int  index_lensing_profile;
  int  index_multipole_for_lensing_profile;
  int  index_completeness;
  int  index_te_of_m;
  int  index_volume;
  int  index_chi2; // in [Mpc/h]^2
  int  index_dgdz; // d(D/a)/dz = D(1-f)
  int  index_lensing_Sigma_crit;
  int  index_vrms2;
  int  index_pk_for_halo_bias;
  int  index_dlnMdeltadlnM;
  int  index_part_id_cov_hsv;

  int  index_mean_y;
  int  index_hmf;

  int index_sigma2_hsv;

  int  index_halo_bias;
  int  index_halo_bias_b2;
  int  index_k_value_for_halo_bias;

  int index_phi_galaxy_counts;
  int index_mean_galaxy_number_density;
  int index_mean_galaxy_bias;
  int index_c500c;
  int index_multipole_for_galaxy_profile;
  int index_multipole_for_truncated_nfw_profile;
  int index_galaxy_profile;

  int index_ngal_for_galaxy_profile;
  int index_ngal_prime_for_galaxy_profile;


  int index_multipole_for_cib_profile;
  int index_frequency_for_cib_profile;
  int index_frequency_prime_for_cib_profile;
  int index_cib_profile;

  int index_W_lensmag;

  int index_W_gallens_sources;

  int index_k_for_pk_hm;
  int index_density_profile;

  int index_multipole_for_pk;

  int index_ell_1;
  int index_ell_2;
  int index_ell_3;

  //////////////

  int index_integral;
  int index_integral_te_y_y;
  int index_integral_2halo_term;

  int index_integral_trispectrum_first;
  int index_integral_trispectrum_last;

  int index_integral_cov_N_cl_first;
  int index_integral_cov_N_cl_last;

  int index_integral_N_for_cov_N_cl_first;
  int index_integral_N_for_cov_N_cl_last;


  int  index_integrals_over_m_first;
  int  index_integrals_over_m_last;

  int  index_integrals_over_z_first;
  int  index_integrals_over_z_last;



  int  index_integral_over_m;
  int  index_integral_te_y_y_over_m;
  int  index_integral_2halo_term_over_m;
  int  index_integral_trispectrum_first_over_m;
  int  index_integral_trispectrum_last_over_m;
  int  index_integral_cov_N_cl_first_over_m;
  int  index_integral_cov_N_cl_last_over_m;
  int  index_integral_N_for_cov_N_cl_first_over_m;
  int  index_integral_N_for_cov_N_cl_last_over_m;




  //mass bins for covariance between cluster counts and power spectrum
  int nbins_M;
  double * M_bins;
  double dlogM;
  double * cov_Y_N_mass_bin_edges;

double * szcounts_fft_qobs;
double * szcounts_fft_z;
double * szcounts_fft_sigmayobs;
int ** szcounts_fft_index_zsig;
// double ** szcounts_fft_rates_at_z_sigy_qobs;
// double * szcounts_fft_nexpected_dndzdqgt;
double * szcounts_fft_dndzdq;
double * szcounts_fft_nexpected_qobs;
int ** szcounts_fft_index_zq;
int ** szcounts_fft_index_zq_final;
double ** szcounts_fft_qmconv_all_patches;

double szcounts_ntot;

  //HOD
  double M_min_HOD;
  double M_min_HOD_cib;
  double M0_HOD;
  double sigma_log10M_HOD;
  double alpha_s_HOD;
  double M1_prime_HOD;

  int centrals_only_HOD;
  int satellites_only_HOD;

  double * M_min_HOD_ngal;
  double * M0_HOD_ngal;
  double * sigma_log10M_HOD_ngal;
  double * alpha_s_HOD_ngal;
  double * M1_prime_HOD_ngal;

  int * centrals_only_ngal;
  int * satellites_only_ngal;

  double * dndz_shift_ngal;
  double * dndz_stretch_ngal;

  double rho_y_gal;

  int M0_Mmin_flag;

  int cosmo_model;

  int use_Amod;
  double Amod;

  int use_nl_bias;
  double bnl;

  double M_min_HOD_mass_factor_unwise;
  double x_out_truncated_nfw_profile;
  double x_out_truncated_nfw_profile_electrons;
  double x_out_truncated_density_profile;
  double M_min_HOD_satellite_mass_factor_unwise;
  double M1_prime_HOD_factor;
  double cvir_tau_profile_factor;


  double x_out_custom1;

  double effective_galaxy_bias;
  double * effective_galaxy_bias_ngal;
  double * effective_galaxy_bias_nl_ngal;

  int use_bg_eff_in_ksz2g_eff;

  int hm_consistency;
  int hm_consistency_ngbar;

  int use_class_sz_fast_mode;
  double * array_lnk;
  double * array_pknl_at_z_and_k;
  double * array_pkl_at_z_and_k;

  // int cszfast_pk_grid_nk;
  // int cszfast_pk_grid_nz;

  int check_consistency_conditions;

  // noise curve for cov(y,y)

  int include_noise_cov_y_y;

  //units for tSZ spectrum
  double exponent_unit;

  //completeness
  double theta_bin_min;
  double theta_bin_max;
  int nthetas;
  double * thetas;

  double *skyfracs;
  int nskyfracs;

  int Ny;
  int Nth;
  double * erfs_2d_to_1d_th_array;
  double * erfs_2d_to_1d_y_array;

  double ** ylims;
  double * sky_averaged_ylims;

  //SZ catalog
  double * szcat_z;
  double * szcat_snr;
  int  szcat_size;

  double shape_noise_siggamma2;
  double ns_gal_per_arcmin2;
  double cl_gal_gal_A_sn;

  int experiment;
  //SO completeness
  double * SO_thetas;
  double * SO_Qfit;
  int  SO_Q_size;

  double * SO_RMS;
  double * SO_skyfrac;
  int  SO_RMS_size;

  double csat_over_cdm;
  //INPUT PARAMETERS
  int nlSZ;
  int n_ell_independent_integrals;
  int n_frequencies_for_cib;

  double * l_unwise_filter;
  double * f_unwise_filter;
  int unwise_filter_size;

  double * M_min_of_z_z;
  double * M_min_of_z_M_min;
  int M_min_of_z_size;

  double * nl_lensing_noise;
  double * l_lensing_noise;
  int lensing_noise_size;

  double * l_ksz_template;
  double * cl_ksz_template;
  int ksz_template_size;

  int damping_1h_term;
  double kstar_damping_1h_term_Mpc; // inverse Mpc

  /*Redshift limits for the integration*/
  double z1SZ;
  double z2SZ;

  double z1SZ_dndlnM;
  double z2SZ_dndlnM;

  double M1SZ_dndlnM;
  double M2SZ_dndlnM;


  int include_y_counterterms_in_yk; // switch for counter term yk 2h calculation for y-part only.
  int include_g_counterterms_in_gk; // switch for counter term gk 2h calculation for g-part only.
  int include_k_counterterms_in_gk; // switch for counter term gk 2h calculation for k-part only.
  int include_gk_counterterms_in_gk; // switch for counter term gk 2h calculation for k-part only.

  double m_min_counter_terms;
  double m_max_counter_terms;


  double y_min;
  double y_max;


  double theta_ej_bcm;
  double delta_bcm;
  double gamma_bcm;
  double eta_star_bcm;
  double log10Mc_bcm;
  double mu_bcm;
  double nu_log10Mc_bcm;
  // double xxx_bcdm;

  // int n_y_y_to_m;
  // int n_z_y_to_m;
  // int n_mass_y_to_m;
  double * array_y_to_m_y;
  double * array_y_to_m_at_z_y;
  double * array_y_to_m_redshift;

  // double z1SZ_L_sat;
  // double z2SZ_L_sat;
  //
  // double M1SZ_L_sat;
  // double M2SZ_L_sat;
  //
  // double epsabs_L_sat;
  // double epsrel_L_sat;

  double mass_epsrel_cluster_counts;
  double mass_epsabs_cluster_counts;

  double redshift_epsrel_cluster_counts;
  double redshift_epsabs_cluster_counts;

  double dlnM_cluster_count_completeness_grid;
  double dz_cluster_count_completeness_grid_low_z;
  double dz_cluster_count_completeness_grid_mid_z;
  double dz_cluster_count_completeness_grid_high_z;

  double lnymin;
  double lnymax;
  double dlny;

  double cluster_count_completeness_grid_z_cutoff_low;
  double cluster_count_completeness_grid_z_cutoff_mid;

  int n_z_W_lensmag;
  int n_z_W_gallens_sources;

  /*Array size*/
  int ndim_redshifts;//number of z in the interpolation
  int ndim_redshifts_for_integral;//number of z in the integration



  int n_k;
  int n_z_dndlnM;
  int n_m_dndlnM;

  int compute_ksz2ksz2;

  // int n_z_L_sat;
  // int n_m_L_sat;
  // int n_nu_L_sat;

  int N_redshift_dndlnM;
  int N_mass_dndlnM;

  //mass limits: h^-1 Msun
  double M1SZ;
  double M2SZ;

  double delta_alpha;
  double alpha_p;

  double alpha_b;
  double Ap;
  int mass_dependent_bias;

  double A_IA;
  double eta_IA;
  double C1_IA;
  //Planck pressure profile
  double P0GNFW;
  double c500;
  double gammaGNFW;
  double alphaGNFW;
  double betaGNFW;

  double ln_x_size_for_pp;
  double * ln_x_for_pp;

  double x_size_for_pp;
  double * x_for_pp;

  int use_websky_m200m_to_m200c_conversion;
  //Battaglia pressure profile
  double alpha_B12;
  double gamma_B12;
  double P0_B12;
  double xc_B12;
  double beta_B12;

  double alpha_m_P0_B12;
  double alpha_m_xc_B12;
  double alpha_m_beta_B12;

  double alpha_z_P0_B12;
  double alpha_z_xc_B12;
  double alpha_z_beta_B12;


  // B.H.
  double mcut_B12;
  double c_B12;
  double c_B16;
  double alphap_m_P0_B12;
  double alphap_m_xc_B12;
  double alphap_m_beta_B12;

  double alpha_c_P0_B12;
  double alpha_c_xc_B12;
  double alpha_c_beta_B12;


    // B.H.
  double mcut;
  double alphap_m_rho0;
  double alphap_m_alpha;
  double alphap_m_beta;

  double alpha_c_rho0;
  double alpha_c_alpha;
  double alpha_c_beta;

  // Battaglia density profile:
  double A_rho0;
  double A_alpha;
  double A_beta;

  double alpha_m_rho0;
  double alpha_m_alpha;
  double alpha_m_beta;

  double alpha_z_rho0;
  double alpha_z_alpha;
  double alpha_z_beta;

  double gamma_B16;
  double xc_B16;


// JCH
//double precision :: f_free=0.85d0 !for kSZ calculations, fraction of free electrons w.r.t. total
//double precision :: mu_e=1.14d0 !mean molecular weight per electron, for primordial composition

  double f_free;
  double mu_e;
  double f_b_gas;

  /*Pressure profile is considered between x_in and x_out*/
  double x_inSZ;
  double x_outSZ;

  double HSEbias;

  /*For the computation of sigma2*/
  int  ndim_masses;
  double logR1SZ; // 0.0034Mpc/h, 1.8e4  solar mass
  double logR2SZ; // 54.9Mpc/h, 7.5e16 solar mass
  double delta_cSZ;



  /*Multplicity function Tinker 2010*/

  double alphaSZ;
  double beta0SZ;
  double gamma0SZ;

  double phi0SZ;
  double eta0SZ;
  int T10_alpha_fixed;


  /*Multplicity function Bocquet 2015*/

  double Ap0;
  double a0;
  double b0;
  double c0;

  int pk_nonlinear_for_vrms2;

  int MF;
  //1:Tinker 2010 (T10)
  //2:Bocquet 2015 (B15)
  //3:Jenkins 2001 (J01)
  //4:Tinker 2008 (T08)
  //5:Tinker 2008 interpolated @ M500 (T08@M500)
  int SHMF;

  //Precision Parameters For qromb_sz_integrand
  int K;
  double EPS;
  double JMAX;


  //Precision Parameters For qromb_sz_sigma
  int K_sigma;
  double EPS_sigma;
  double JMAX_sigma;

  ////////////////////////
  //integration method and parameters (mass)
  int integration_method_mass;

  double redshift_epsrel;
  double redshift_epsabs;

  double mass_epsrel;
  double mass_epsabs;

  double pressure_profile_epsabs;
  double pressure_profile_epsrel;
  double nu_y_dist_GHz;


  int * galaxy_samples_list;
  int galaxy_samples_list_num;
  int ngal_dim;

  int cib_frequency_list_num;
  int cib_dim;
  double * cib_frequency_list;
  double * cib_Snu_cutoff_list_in_mJy;

  int id_nu_cib_to_save;
  int id_nu_prime_cib_to_save;

  double ystar_ym;
  double alpha_ym;
  double sigmaM_ym;
  double beta_ym;
  double A_ym;
  double B_ym;
  double C_ym;
  double m_pivot_ym;

  double alpha_theta;
  int y_m_relation;
  double thetastar;

  int use_maniyar_cib_model;
  double maniyar_cib_tau;
  double maniyar_cib_zc;
  double maniyar_cib_etamax;
  double maniyar_cib_fsub;

  //BB: added for class_sz
  int ln_k_size_for_tSZ;
  double k_per_decade_for_tSZ;
  double k_min_for_pk_in_tSZ;
  double k_max_for_pk_in_tSZ;
  double * ln_k_for_tSZ;


  int ln_k_size_for_vrms2;
  double k_per_decade_for_vrms2;
  double k_min_for_pk_in_vrms2;
  double k_max_for_pk_in_vrms2;
  double * ln_k_for_vrms2;


int nsteps_m;
int nsteps_z;

double * steps_z;
double * steps_m;

  // Table 1  of MM20
  double alpha_cib; //redshift evolution of dust temperature
  double T0_cib; // dust temperature today
  double beta_cib; // emissivity index of sed
  double gamma_cib; // Power law index of SED at high frequency
  double delta_cib; // Redshift evolution of L − M normalisation
  double m_eff_cib; // Most efficient halo mass in Msun/h
  double L0_cib; // Normalisation of L − M relation
  double sigma2_LM_cib; // Size of of halo masses sourcing CIB emission
  int has_cib_flux_cut;
  double z_obs_cib;
  double z_plateau_cib;
  double M_min_subhalo_in_Msun;
  int use_nc_1_for_all_halos_cib_HOD;
  int use_redshift_dependent_M_min;

  double nfw_profile_epsabs;
  double nfw_profile_epsrel;

  int patterson_show_neval;

  int number_of_mass_bins; //for trapezoidal rule
  ////////////////////////

  ////////////////////////
  //integration method and parameters (pressure profile)
  int integration_method_pressure_profile;

  //Foreground parameters
  double A_cib, A_rs, A_ir, A_cn;

  double * k_for_pk_hm;
  double dlnk_for_pk_hm;
  double k_min_for_pk_hm;
  double k_max_for_pk_hm;
  int n_k_for_pk_hm;
  double z_for_pk_hm;

  //Cl spectrum
  double * ell;
  double * ell_plc;
  double * ell_plc_no_low_ell;
  double * ell_plc_low;
  double * ell_mock;
  double * ell_trispectrum;
  double * x_gauss;
  double * w_gauss;

  double * frequencies_for_cib;

  double * ell_kSZ2_gal_multipole_grid;
  int N_kSZ2_gal_multipole_grid;


  double * theta_kSZ2_gal_theta_grid;
  int N_kSZ2_gal_theta_grid;

  int photo_z_params;
  double dndz_shift_gal;
  double dndz_shift_source_gal;
  double dndz_stretch_gal;
  double dndz_stretch_source_gal;
  double shear_calibration_m;

  double dlogell;
  double dell;
  double ell_min_mock;
  double ell_max_mock;


  double freq_max;
  double freq_min;
  double dlogfreq;
  double dfreq;

  double Tcmb_gNU_at_150GHz;
  double Tcmb_gNU;

  double sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb;

  double Rho_crit_0;
  double D_0;
  double D_z1SZ;
  double Omega_m_0;
  double Omega_r_0;
  double Omega_ncdm_0;
  double Omega0_b;
  double Omega0_cdm;
  double bispectrum_lambda_k2;
  double bispectrum_lambda_k3;

  double Sigma8OmegaM_SZ;
  double sigma8_Pcb;

  short has_knl;
  short has_nl_index;
  short has_vrms2;
  short has_sigma2_hsv;

  short has_class_sz_structure;  //do we need tSZ spectrum? */
  short sz_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */
  short write_sz;  //do we need SZ quatitiies vs redshift? */

  int use_planck_binned_proba;
  double bin_z_min_cluster_counts;
  double bin_z_max_cluster_counts;
  double bin_dz_cluster_counts;
  int apply_relativistic_correction_to_y_m;
  double bin_dlog10_snr;
  double bin_dlog10_snr_last_bin;
  double log10_snr_min;
  double log10_snr_max;

  double x_out_truncated_nfw_profile_satellite_galaxies;
  double * x_out_truncated_nfw_profile_satellite_galaxies_ngal;

  double f_cen_HOD;
  double * f_cen_HOD_ngal;
  double Delta_z_lens;
  double Delta_z_source;

  double fEDE;
  double log10z_c;
  double thetai_scf;


  short has_completeness_for_ps_SZ;
  short has_completeness;
  short which_ps_sz;
  double H0_in_class_units;
  double sky_area_deg2;
  int ell_sz;
  // Figure 7 of KS02 -> KS02
  // Planck 2015 effective multipoles -> P15
  // SZFASTDKS -> DKS

  // halo occupation distribution
  int hod_model;


  int concentration_parameter;
  //Duffy 2008: D08
  //Seljak 2000: S00

  int pressure_profile;
  //Planck 2013 (P13)
  //Arnaud et al 2010 (A10)
  //Custom. GNFW

  int tau_profile;
  int tau_profile_mode;

  int HMF_prescription_NCDM;
  int effective_temperature;
  int temperature_mass_relation;
  int mean_y;

  double * PP_lnx;
  double * PP_lnI;
  double * PP_d2lnI;

  int PP_lnx_size;

  double * RNFW_lnx;
  double * RNFW_lnI;

  int RNFW_lnx_size;


  double * T10_ln1pz;
  double * T10_lnalpha;
  int T10_lnalpha_size;

  double * normalized_source_dndz_z;
  double * normalized_source_dndz_phig;

  int normalized_source_dndz_size;

  double * normalized_dndz_z;
  double * normalized_dndz_phig;

  double ** normalized_dndz_ngal_z;
  double ** normalized_dndz_ngal_phig;

  int normalized_dndz_size;
  int * normalized_dndz_ngal_size;

  double * normalized_fdndz_z;
  double * normalized_fdndz_phig;

  int normalized_fdndz_size;

  double * normalized_cosmos_dndz_z;
  double * normalized_cosmos_dndz_phig;

  int normalized_cosmos_dndz_size;

  double * unbinned_nl_yy_ell;
  double * unbinned_nl_yy_n_ell;
  int nl_yy_is_binned;
  int unbinned_nl_yy_size;

  double * unbinned_nl_tt_ell;
  double * unbinned_nl_tt_n_ell;
  int unbinned_nl_tt_size;

  int truncate_gas_pressure_wrt_rvir;

  int no_tt_noise_in_kSZ2X_cov;

  double * CM_redshift;
  double * CM_logM;

  int CM_redshift_size;
  int CM_logM_size;
  double * CM_logC;


  // double * array_profile_2h_ln_1pz;
  double * array_profile_ln_rho_2h_at_k_and_z;
  double * array_profile_rho_2h_at_r_and_z;


  double * array_m_m200m_to_m200c;
  double * array_ln_1pz_m200m_to_m200c;
  double * array_m200m_to_m200c_at_z_and_M;

  double * array_m_m200c_to_mvir;
  double * array_ln_1pz_m200c_to_mvir;
  double * array_m200c_to_mvir_at_z_and_M;

  double * array_m_m200m_to_mvir;
  double * array_ln_1pz_m200m_to_mvir;
  double * array_m200m_to_mvir_at_z_and_M;

  double * array_m_m200c_to_m200m;
  double * array_ln_1pz_m200c_to_m200m;
  double * array_m200c_to_m200m_at_z_and_M;


  double * array_m_m200m_to_m500c;
  double * array_ln_1pz_m200m_to_m500c;
  double * array_m200m_to_m500c_at_z_and_M;

  double * array_m_m200c_to_m500c;
  double * array_ln_1pz_m200c_to_m500c;
  double * array_m200c_to_m500c_at_z_and_M;


  double * array_m_m500c_to_m200c;
  double * array_ln_1pz_m500c_to_m200c;
  double * array_m500c_to_m200c_at_z_and_M;


  double * array_m_m500c_to_m200m;
  double * array_ln_1pz_m500c_to_m200m;
  double * array_m500c_to_m200m_at_z_and_M;

  double ** array_custom1_profile_u_at_lnk_lnm_ln1pz;
  double ** array_custom1_profile_u_at_lnx_lnm_ln1pz;
  double * array_custom1_profile_ln_k;
  double * array_custom1_profile_ln_x;
  double * array_custom1_profile_ln_m;
  double * array_custom1_profile_ln_1pz;



  // int array_b_custom1_n_z;
  double * array_b_custom1_ln1pz;
  double * array_b_custom1_bias;

  double ** array_pressure_profile_ln_p_at_lnk_lnm_z;
  double * array_pressure_profile_ln_p_at_lnk;
  double * array_pressure_profile_ln_k;
  double * array_pressure_profile_2h_ln_k;
  double * array_pressure_profile_ln_r;
  double * array_pressure_profile_ln_m;
  double * array_pressure_profile_ln_1pz;
  double * array_pressure_profile_ln_pressure_2h_at_k_and_z;
  double * array_pressure_profile_pressure_2h_at_r_and_z;

  double * array_pressure_profilel_ln_m;
  double * array_pressure_profilel_ln_1pz;
  double * array_pressure_profile_ln_l;
  double ** array_pressure_profile_ln_p_at_lnl_lnm_z;

  double ** array_profile_ln_rho_at_lnk_lnM_z;
  double * array_profile_ln_r;
  double * array_profile_ln_k;
  double * array_profile_ln_m;
  double * array_profile_ln_1pz;


  // int n_m_matter_density_profile;

  double * array_matter_density_profile_ln_k;
  double * array_matter_density_profile_ln_m;
  double * array_matter_density_profile_ln_1pz;
  double ** array_profile_ln_rho_matter_at_lnk;

  // int array_profile_ln_PgNFW_at_lnl_over_ls_size; defined in class_sz_precisions.h
  double * array_profile_ln_l_over_ls;
  double * array_profile_ln_PgNFW_at_lnl_over_ls;

  double * dndlnM_array_z;
  double * dndlnM_array_m;

  double * array_m_dndlnM;
  double * array_z_dndlnM;
  double * array_dndlnM_at_z_and_M;

  double * array_m_L_sat;
  double * array_z_L_sat;
  double * array_nu_L_sat;

  double ** array_L_sat_at_M_z_nu;
  double ** array_L_sat_at_z_and_M_at_nu;
  //double * array_L_sat_at_z_and_M_at_nu_prime;

  double ** array_z_W_nlensmag;
  double ** array_W_nlensmag;


  double * array_z_W_lensmag;
  double * array_W_lensmag;

  double * array_z_W_gallens_sources;
  double * array_W_gallens_sources;

  double * array_redshift;
  double * array_radius;
  // double * array_k;
  double * array_nl_index_at_z_and_k;
  double * array_nl_index_at_z_and_k_no_wiggles;
  double * array_sigma_at_z_and_R;
  double * array_dsigma2dR_at_z_and_R;

  double * array_knl_at_z;
  double * array_vrms2_at_z;
  double * array_sigma2_hsv_at_z;

  double * array_mean_galaxy_number_density;
  double ** array_mean_galaxy_number_density_ngal;


  double * array_custom1_redshift_kernel_W;
  double * array_custom1_redshift_kernel_ln1pz;

  double * array_mean_galaxy_bias;

  // int n_z_hmf_counter_terms;
  int hm_consistency_counter_terms_done;
  double * array_redshift_hmf_counter_terms;
  double * array_hmf_counter_terms_nmin;
  double * array_hmf_counter_terms_b1min;
  double * array_hmf_counter_terms_b2min;
  ErrorMsg error_message; /**< zone for writing error messages */


  double * array_n5k_F1_F;
  double * array_n5k_F1_k;
  int * array_n5k_F1_l;

  double * n5k_pk_z;
  double * n5k_pk_k;
  double * n5k_pk_pk;
  int n5k_pk_z_size;
  int n5k_pk_k_size;


  double * cib_Snu_z;
  double * cib_Snu_nu;
  double * cib_Snu_snu;
  int cib_Snu_z_size;
  int cib_Snu_nu_size;

  double * n5k_cl_K1_K1;
  double * n5k_cl_K1_chi;
  int n5k_cl_K1_size;

  double * n5k_z_of_chi_z;
  double * n5k_z_of_chi_chi;
  int n5k_z_of_chi_size;


  double * array_psi_b2t_redshift;
  double * array_psi_b2t_multipole;
  double * array_psi_b2t_psi;

  double * array_psi_b2g_redshift;
  double * array_psi_b2g_multipole;
  double * array_psi_b2g_psi;

  double * array_psi_b1g_redshift;
  double * array_psi_b1g_multipole;
  double * array_psi_b1g_psi;

  double * array_psi_b2kg_redshift;
  double * array_psi_b2kg_multipole;
  double * array_psi_b2kg_psi;

  double * array_psi_b1kg_redshift;
  double * array_psi_b1kg_multipole;
  double * array_psi_b1kg_psi;



  double * array_psi_b1t_redshift;
  double * array_psi_b1t_multipole;
  double * array_psi_b1t_psi;

  double * array_dcib0dz_nu;
  double * array_dcib0dz_redshift;
  double * array_dcib0dz_at_z_nu;


  double * array_m_to_xout_mass;
  double * array_m_to_xout_redshift;
  double * array_m_to_xout_at_z_m;


  double * array_dydz_redshift;
  double * array_dydz_at_z;


  double * array_psi_b1gt_redshift;
  double * array_psi_b1gt_multipole;
  double ** array_psi_b1gt_psi;

  double * array_psi_b1kgt_redshift;
  double * array_psi_b1kgt_multipole;
  double ** array_psi_b1kgt_psi;

  double * array_psi_b1kgg_redshift;
  double * array_psi_b1kgg_multipole;
  double ** array_psi_b1kgg_psi;
  // int n_z_psi_b1g;
  // int n_l_psi_b1g;


};

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


int class_sz_cosmo_init(struct background * pba,
                         struct thermo * pth,
                         struct perturbs * ppt,
                         struct nonlinear * pnl,
                         struct primordial * ppm,
                         struct spectra * psp,
                         struct lensing * ple,
                         struct class_sz_structure * pclass_sz,
                         struct precision * ppr);

int class_sz_tabulate_init(
                          struct background * pba,
                          struct thermo * pth,
                          struct perturbs * ppt,
                          struct nonlinear * pnl,
                          struct primordial * ppm,
                          struct spectra * psp,
                          struct lensing * ple,
                          struct class_sz_structure * pclass_sz,
                          struct precision * ppr
                        );

int class_sz_integrate_init(struct background * pba,
                         struct thermo * pth,
                         struct perturbs * ppt,
                         struct nonlinear * pnl,
                         struct primordial * ppm,
                         struct spectra * psp,
                         struct lensing * ple,
                         struct class_sz_structure * pclass_sz,
                         struct precision * ppr);


  int class_sz_free(struct class_sz_structure *pclass_sz);


  int compute_sz(struct background * pba,
                 struct nonlinear * pnl,
                 struct primordial * ppm,
                 struct perturbs * ppt,
                 struct class_sz_structure * pclass_sz,
                 double * pvecback,
                 double * Pvectsz);



  //This evaluates the integrand which will be integrated
  //over M and then over z
  double integrand_at_m_and_z(double logM ,
                              double * pvecback,
                              double * pvectsz,
                              struct background * pba,
                              struct primordial * ppm,
                              struct nonlinear * pnl,
                              struct perturbs * ppt,
                              struct class_sz_structure * pclass_sz);

  double delta_ell_lens_at_ell_and_z( double * pvecback,
                                  double * pvectsz,
                                  struct background * pba,
                                  struct primordial * ppm,
                                  struct nonlinear * pnl,
                                  struct class_sz_structure * pclass_sz);

  double delta_ell_isw_at_ell_and_z( double * pvecback,
                                  double * pvectsz,
                                  struct background * pba,
                                  struct primordial * ppm,
                                  struct nonlinear * pnl,
                                  struct class_sz_structure * pclass_sz);
  int do_mass_conversions(
                         double logM,
                         double z,
                         double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct class_sz_structure * pclass_sz);
  int evaluate_HMF_at_logM_and_z(
                   double logM ,
                   double z,
                   double * pvecback,
                   double * pvectsz,
                   struct background * pba,
                   struct nonlinear * pnl,
                   struct class_sz_structure * pclass_sz);

  int evaluate_completeness(double * pvecback,
                            double * pvectsz,
                            struct background * pba,
                            struct class_sz_structure * pclass_sz);

  double evaluate_pressure_profile(double kl,
                                double * pvecback,
                                double * pvectsz,
                                struct background * pba,
                                struct class_sz_structure * pclass_sz);

  // double evaluate_IA(double kl,
  //                    double * pvecback,
  //                    double * pvectsz,
  //                    struct background * pba,
  //                    struct class_sz_structure * pclass_sz);
double get_IA_of_z(double z,
                   struct background * pba,
                   struct class_sz_structure * pclass_sz);

double get_A_IA_of_z(double z,
                   struct background * pba,
                   struct class_sz_structure * pclass_sz);
  int evaluate_tau_profile(double k,
                           double * pvecback,
                           double * pvectsz,
                           struct background * pba,
                           struct class_sz_structure * pclass_sz);

  double evaluate_lensing_profile(double kl,
                                  double m_delta,
                                  double r_delta,
                                  double c_delta,
                                  double * pvecback,
                                  double * pvectsz,
                                  struct background * pba,
                                  struct class_sz_structure * pclass_sz);

double get_fstar_of_m_at_z(double m,
                           double z,
                           struct class_sz_structure * pclass_sz);

double get_tau_profile_at_z_m_l(double z,
                                double m,
                                double k,
                                struct class_sz_structure * pclass_sz,
                                struct background * pba);


double get_ksz_filter_at_l(double l,
                           struct class_sz_structure * pclass_sz);

double get_M_min_of_z(double l,
                      struct class_sz_structure * pclass_sz);

  int write_output_to_files_ell_indep_ints(struct nonlinear * pnl,
                                           struct background * pba,
                                           struct class_sz_structure * pclass_sz);

  int write_output_to_files_cl(struct nonlinear * pnl,
                               struct background * pba,
                               struct primordial * ppm,
                               struct class_sz_structure * pclass_sz);


  int show_preamble_messages(struct background * pba,
                             struct thermo * pth,
                             struct nonlinear * pnl,
                             struct primordial * ppm,
                             struct class_sz_structure * pclass_sz);

  double gnu_tsz_of_nu_in_ghz(double nu_in_ghz,double Tcmb);

  int show_results(struct background * pba,
                   struct nonlinear * pnl,
                   struct primordial * ppm,
                   struct class_sz_structure * pclass_sz);

  int select_multipole_array(struct class_sz_structure * pclass_sz);

  int evaluate_halo_bias(double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct primordial * ppm,
                         struct nonlinear * pnl,
                         struct perturbs * ppt,
                         struct class_sz_structure * pclass_sz);

  int evaluate_halo_bias_b2(double * pvecback,
                            double * pvectsz,
                            struct background * pba,
                            struct primordial * ppm,
                            struct nonlinear * pnl,
                            struct class_sz_structure * pclass_sz);

 int evaluate_effective_galaxy_bias_ngal(int index_g,
                                         double * pvecback,
                                         double * pvectsz,
                                         struct background * pba,
                                         struct primordial * ppm,
                                         struct nonlinear * pnl,
                                         struct class_sz_structure * pclass_sz);

 int evaluate_effective_galaxy_bias(double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct primordial * ppm,
                                    struct nonlinear * pnl,
                                    struct class_sz_structure * pclass_sz);
double get_pk_lin_at_k_and_z(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct class_sz_structure * pclass_sz);
double get_pk_lin_at_k_and_z_fast(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct class_sz_structure * pclass_sz);

double get2_pk_lin_at_k_and_z(//double * pvecback,//double * pvectsz,
  double * r,double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct class_sz_structure * pclass_sz);
double get_pk_nonlin_at_k_and_z(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct class_sz_structure * pclass_sz);

double get_pk_nonlin_at_k_and_z_fast(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct class_sz_structure * pclass_sz);

 int evaluate_pk_at_ell_plus_one_half_over_chi(double * pvecback,
                                              double * pvectsz,
                                              struct background * pba,
                                              struct primordial * ppm,
                                              struct nonlinear * pnl,
                                              struct class_sz_structure * pclass_sz);

 int evaluate_pk_at_ell_plus_one_half_over_chi_today(double * pvecback,
                                                      double * pvectsz,
                                                      struct background * pba,
                                                      struct primordial * ppm,
                                                      struct nonlinear * pnl,
                                                      struct class_sz_structure * pclass_sz);


double evaluate_pk_halofit_over_pk_linear_at_ell_plus_one_half_over_chi(double * pvecback,
                                                                     double * pvectsz,
                                                                     struct background * pba,
                                                                     struct primordial * ppm,
                                                                     struct nonlinear * pnl,
                                                                     struct class_sz_structure * pclass_sz);
int load_cl_ksz_template(struct class_sz_structure * pclass_sz);

int load_nl_lensing_noise(struct class_sz_structure * pclass_sz);


  int initialise_and_allocate_memory(struct class_sz_structure * pclass_sz);


  int evaluate_temperature_mass_relation(double * pvecback,
                                         double * pvectsz,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz);

double evaluate_dlnMdeltadlnM(double logM,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct nonlinear * pnl,
                             struct class_sz_structure * pclass_sz);

int evaluate_vrms2(double * pvecback,
                   double * pvectsz,
                   struct background * pba,
                   struct nonlinear * pnl,
                   struct class_sz_structure * pclass_sz);


int evaluate_sigma2_hsv(double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct nonlinear * pnl,
                         struct class_sz_structure * pclass_sz);

double integrand_mass(double xi, void *p);


int write_redshift_dependent_quantities(struct background * pba,
                                        struct class_sz_structure * pclass_sz);


// int evaluate_tau_profile(double * pvecback,
//                         double * pvectsz,
//                         struct background * pba,
//                         struct class_sz_structure * pclass_sz);

int tabulate_normalization_matter_density_profile(struct class_sz_structure * pclass_sz,struct background * pba);


int tabulate_normalization_gas_density_profile(struct class_sz_structure * pclass_sz,struct background * pba);

int tabulate_gas_pressure_profile_gNFW(struct background * pba,
                                       struct class_sz_structure * pclass_sz);

int tabulate_gas_pressure_profile_gNFW_fft(struct background * pba,
                                           struct class_sz_structure * pclass_sz);

int tabulate_gas_pressure_profile_B12(struct background * pba,
                                      struct class_sz_structure * pclass_sz);

int tabulate_gas_pressure_profile_B12_l(struct background * pba,
                                      struct class_sz_structure * pclass_sz);


int tabulate_gas_pressure_profile_B12_fft(struct background * pba,
                                          struct class_sz_structure * pclass_sz);

double evaluate_mean_galaxy_number_density_at_z(double z,
                                                struct class_sz_structure * pclass_sz);

double evaluate_mean_galaxy_number_density_at_z_ngal(double z,
                                                     int index_g,
                                                     struct class_sz_structure * pclass_sz);

double get_mean_galaxy_bias_at_z(double z,
                                 struct class_sz_structure * pclass_sz);



double get_dyldzdlnm_at_l_z_and_m(double l,
                                  double z,
                                  double m,
                                  struct background * pba,
                                  struct nonlinear * pnl,
                                  struct class_sz_structure * pclass_sz);

double get_dyl2dzdlnm_at_l_z_and_m(double l,
                                   double z,
                                   double m,
                                   struct background * pba,
                                   struct nonlinear * pnl,
                                   struct class_sz_structure * pclass_sz);

double get_dygldzdlnm_at_l_z_and_m(double l,
                                  double z,
                                  double m,
                                  struct background * pba,
                                  struct nonlinear * pnl,
                                  struct class_sz_structure * pclass_sz);


double HOD_mean_number_of_central_galaxies(double z,
                                           double M_halo,
                                           double M_min,
                                           double sigma_lnM,
                                           double f_cen,
                                           struct class_sz_structure * pclass_sz,
                                           struct background * pba);

double HOD_mean_number_of_satellite_galaxies(double z,
                                             double M_halo,
                                             double Nc_mean,
                                             double M_min,
                                             double alpha_s,
                                             double M1_prime,
                                             struct class_sz_structure * pclass_sz,
                                             struct background * pba);

double get_galaxy_profile_at_z_m_l_1h(double z,
                                      double m,
                                      double r_delta,
                                      double c_delta,
                                      double l,
                                      struct class_sz_structure * pclass_sz,
                                      struct background * pba);


int evaluate_galaxy_profile_1h(double kl,
                               double m_delta,
                               double r_delta,
                               double c_delta,
                               double * pvecback,
                               double * pvectsz,
                               struct background * pba,
                               struct class_sz_structure * pclass_sz);

int evaluate_galaxy_profile_ngal(double kl,
                               double m_delta,
                               double r_delta,
                               double c_delta,
                               double * pvecback,
                               double * pvectsz,
                               struct background * pba,
                               struct class_sz_structure * pclass_sz);


int evaluate_galaxy_profile_2h(double kl,
                               double m_delta,
                               double r_delta,
                               double c_delta,
                               double * pvecback,
                               double * pvectsz,
                               struct background * pba,
                               struct class_sz_structure * pclass_sz);

double get_truncated_nfw_profile_at_z_m_k_xout(//double * pvecback,
                                      double z,
                                      double m,
                                      double r_delta,
                                      double c_delta,
                                      double k,
                                      double xout,
                                      // double delta,
                                      struct background * pba,
                                      struct class_sz_structure * pclass_sz);

double get_nfw_with_power_law_profile_at_z_m_q(
                                               double z,
                                               double m,
                                               double q,
                                               struct class_sz_structure * pclass_sz);

double get_nfw_with_power_law_profile_at_x(double x,
                                           double n,
                                           double c);

double evaluate_truncated_nfw_profile(
                                   double z,
                                   double k,
                                   double r_delta,
                                   double c_delta,
                                   double xout);



int evaluate_c200m_D08(double * pvecback,
                        double * pvectsz,
                        struct background * pba,
                        struct class_sz_structure * pclass_sz);

double get_c200c_at_m_and_z_D08(double M,
                                double z);

double get_c200c_at_m_and_z(//double * pvecback,
                        double m,
                        double z,
                        struct background * pba,
                        struct class_sz_structure * pclass_sz);

double get_c500c_at_m_and_z(//double * pvecback,
                        double m,
                        double z,
                        struct background * pba,
                        struct class_sz_structure * pclass_sz);

double get_galaxy_number_counts(double z,
                                struct class_sz_structure * pclass_sz);


double get_f_of_sigma_at_m_and_z(double m,
                                 double z,
                                 struct background * pba,
                                 struct nonlinear * pnl,
                                 struct class_sz_structure * pclass_sz);


double get_source_galaxy_number_counts(double z,
                                struct class_sz_structure * pclass_sz);
double radial_kernel_W_galaxy_at_z( double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct class_sz_structure * pclass_sz);

double radial_kernel_W_galaxy_ngal_at_z(  int index_g,
                                          double * pvecback,
                                          double z,
                                          struct background * pba,
                                          struct class_sz_structure * pclass_sz);

double radial_kernel_W_galaxy_lensing_magnification_nlensmag_at_z( int index_g,
                                                          double * pvectsz, //double * pvecback,
                                                          double z,
                                                          // double z,
                                                          // double * pvectsz,
                                                          struct background * pba,
                                                          struct class_sz_structure * pclass_sz);
// used for the linear bias cases.
double radial_kernel_W_lensing_at_z(double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct primordial * ppm,
                                    struct nonlinear * pnl,
                                    struct class_sz_structure * pclass_sz);

double radial_kernel_W_lensing_magnification_at_z(double * pvecback,
                                                  double * pvectsz,
                                                  struct background * pba,
                                                  struct primordial * ppm,
                                                  struct nonlinear * pnl,
                                                  struct class_sz_structure * pclass_sz);

double radial_kernel_W_cmb_lensing_at_z(double z,
                                        double * pvectsz,
                                        struct background * pba,
                                        struct class_sz_structure * pclass_sz);



double radial_kernel_W_galaxy_lensing_at_z(double z,
                                           // double * pvectsz,
                                           // struct background * pba,
                                           struct class_sz_structure * pclass_sz);

double radial_kernel_W_galaxy_lensing_magnification_at_z(double z,
                                                         double * pvectsz,
                                                         struct background * pba,
                                                         struct class_sz_structure * pclass_sz);


double evaluate_galaxy_number_counts( double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct class_sz_structure * pclass_sz);

double evaluate_galaxy_number_counts_fdndz( double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct class_sz_structure * pclass_sz);
double evaluate_unwise_m_min_cut(double z,
                                 int sample_id,
                                 struct class_sz_structure * pclass_sz);


int evaluate_cib_profile(double m_delta,
                         double r_delta,
                         double c_delta,
                         double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct class_sz_structure * pclass_sz);

double Luminosity_of_central_galaxies(double z,
                                      double  M_halo,
                                      double nu,
                                      double * pvectsz,
                                      struct class_sz_structure * pclass_sz,
                                      struct background * pba);

double Luminosity_of_satellite_galaxies(double z,
                                        double  M_halo,
                                        double nu,
                                        struct class_sz_structure * pclass_sz,
                                        struct background * pba);

double maniyar_cib_Mdot(double M, double z, struct class_sz_structure * pclass_sz);
double evaluate_Sigma_cib(double M, struct class_sz_structure * pclass_sz);
double evaluate_phi_cib(double z, struct class_sz_structure * pclass_sz);
double evaluate_sed_cib(double z, double nu, struct class_sz_structure * pclass_sz);
double evaluate_dust_temperature(double z, struct class_sz_structure * pclass_sz);
double evaluate_galaxy_luminosity(double z, double M, double nu, struct class_sz_structure * pclass_sz);
double subhalo_hmf_dndlnMs(double M_host,double M_sub,struct class_sz_structure * pclass_sz);

double integrand_kSZ2_X_at_theta(double ell_prime, void *p);
double integrand_kSZ2_X(double theta, void *p);

double integrand_kSZ2_X_lensing_term_at_theta(double ell_prime, void *p);
double integrand_kSZ2_X_lensing_term(double theta, void *p);

double integrand_lensmag(double ln1pzs, void *p);
double integrand_nlensmag(double ln1pzs, void *p);

int evaluate_matter_density_profile(
                             double k,
                             double r_delta,
                             double c_delta,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct class_sz_structure * pclass_sz);
double get_matter_bispectrum_at_z_tree_level_PT(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct class_sz_structure * pclass_sz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double get_te_of_m500c_at_z_arnaud(double m, double z, struct background * pba,struct class_sz_structure * pclass_sz);
double get_te_of_m500c_at_z_lee(double m, double z, struct background * pba,struct class_sz_structure * pclass_sz);


int  evaluate_ttg_bispectrum_at_z_tree_level_PT(double * r,
                                                      double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct class_sz_structure * pclass_sz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);


int  evaluate_ttg_bispectrum_at_z_effective_approach(double * r,
                                                      double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct class_sz_structure * pclass_sz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double get_ttg_bispectrum_at_z_effective_approach(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct class_sz_structure * pclass_sz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);


double get_ttg_bispectrum_at_z_tree_level_PT(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct class_sz_structure * pclass_sz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double get_matter_bispectrum_at_z_effective_approach_smoothed(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct class_sz_structure * pclass_sz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double get_matter_bispectrum_at_z_effective_approach(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct class_sz_structure * pclass_sz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);
double get_matter_bispectrum_at_z_effective_approach_SC(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct class_sz_structure * pclass_sz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double bispectrum_f2_kernel(double k,
                            double k_prime,
                            double k_prime_prime);
double bispectrum_f2_kernel_eff_SC(double k1,
                            double k2,
                            double k3,
                            double n1,
                            double n2,
                            double sig8_at_z,
                            double knl);

double bispectrum_f2_kernel_eff_a_SC(double k1,double n1,double sig8_at_z,double knl);
double bispectrum_f2_kernel_eff_b_SC(double k1,double n1,double knl);
double bispectrum_f2_kernel_eff_c_SC(double k1,double n1,double knl);
double bispectrum_f2_kernel_eff(double k1,
                            double k2,
                            double k3,
                            double n1,
                            double n2,
                            double sig8_at_z,
                            double knl);

double bispectrum_f2_kernel_eff_a(double k1,double n1,double sig8_at_z,double knl);
double bispectrum_f2_kernel_eff_b(double k1,double n1,double knl);
double bispectrum_f2_kernel_eff_c(double k1,double n1,double knl);
double bispectrum_f2_kernel_eff_Q3(double n1);

double get_rho_crit_at_z(double z_asked,
                         struct background * pba,
                         struct class_sz_structure * pclass_sz);

double get_c200m_at_m_and_z(double M,
                            double z,
                            struct background * pba,
                            struct class_sz_structure * pclass_sz);

double get_c200m_at_m_and_z_D08(double M,
                                double z);

double get_c200m_at_m_and_z_B13(double M,
                                double z,
                                struct background * pba,
                                struct class_sz_structure * pclass_sz);

double get_c200c_at_m_and_z_B13(double M,
                                double z,
                                struct background * pba,
                                struct class_sz_structure * pclass_sz);
double get_gas_profile_at_x_M_z_nfw_200c(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz);

double get_rvir_of_m200c_at_z(
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz);

double get_mass_profile_at_x_M_z_nfw_200m(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz);

double get_gas_profile_at_x_M_z_nfw_200m(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz);

double get_gas_profile_at_x_M_z_bcm_200c(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz);


double get_gas_profile_at_x_M_z_b16_200c(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         double c,
                                         double A_rho,
                                         double A_alpha,
                                         double A_beta,
                                         double alpha_m_rho0,
                                         double alpha_m_alpha,
                                         double alpha_m_beta,
                                         double alpha_z_rho0,
                                         double alpha_z_alpha,
                                         double alpha_z_beta,
                                         // break model param
		                                     double mcut,
		                                     double alphap_m_rho0,
                                         double alphap_m_alpha,
                                         double alphap_m_beta,
		                                     double alpha_c_rho0,
                                         double alpha_c_alpha,
                                         double alpha_c_beta,
                                         // end break model param
                                         double gamma,
                                         double xc,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz);


double get_second_order_bias_at_z_and_nu(double z,
                                         double nu,
                                         struct class_sz_structure * pclass_sz,
                                         struct background * pba);

double get_first_order_bias_at_z_and_nu(double z,
                                         double nu,
                                         struct class_sz_structure * pclass_sz);

double get_ng_bias_contribution_at_z_and_k(double z,
                                           double k,
                                           double bh,
                                           struct background * pba,
                                           struct perturbs * ppt,
                                           struct class_sz_structure * pclass_sz);

double get_scale_dependent_bias_at_z_and_k(double z,
                                           double k,
                                           double bh,
                                           struct class_sz_structure *pclass_sz);


double get_vrms2_at_z(double z,
                      struct class_sz_structure * pclass_sz);

double get_dlnsigma_dlnR_at_z_and_m(double z,
                                    double m,
                                    struct class_sz_structure * pclass_sz,
                                    struct background * pba);
double get_sigma_at_z_and_m(double z,
                            double m,
                            struct class_sz_structure * pclass_sz,
                            struct background * pba);
double get_sigma8_at_z(double z,
                      struct class_sz_structure * pclass_sz,
                      struct background * pba);
double get_nu_at_z_and_m(double z,
                         double m,
                         struct class_sz_structure * pclass_sz,
                         struct background * pba);

// this is r_200c*P_200c
double get_1e6xdy_from_battaglia_pressure_at_x_z_and_m200c(double z,
                                                           double m,
                                                           double x,
                                                           struct background * pba,
                                                           struct class_sz_structure * pclass_sz);

double get_upp_from_gnfw_pressure_at_x_z_and_m500c(double z,
                                                      double m,
                                                      double x,
                                                      double d,
                                                      struct background * pba,
                                                      struct class_sz_structure * pclass_sz);


double get_1e6xdy_from_gnfw_pressure_at_x_z_and_m500c(double z,
                                                      double m,
                                                      double x,
                                                      double d,
                                                      struct background * pba,
                                                      struct class_sz_structure * pclass_sz);

double get_pressure_P_over_P_delta_at_x_M_z_b12_200c(double x_asked,
                                                     double m_asked,
                                                     double z_asked,
                                                     double c_asked,
                                                     double A_P0,
                                                     double A_xc,
                                                     double A_beta,
                                                     double alpha_m_P0,
                                                     double alpha_m_xc,
                                                     double alpha_m_beta,
                                                     double alpha_z_P0,
                                                     double alpha_z_xc,
                                                     double alpha_z_beta,
                                  							     double mcut,
                                  							     double alphap_m_P0,
                                  							     double alphap_m_xc,
                                  							     double alphap_m_beta,
                                  							     double alpha_c_P0,
                                  							     double alpha_c_xc,
                                  							     double alpha_c_beta,
                                                     double alpha,
                                                     double gamma,
                                                     struct background * pba,
                                                     struct class_sz_structure * tsz);

double get_pressure_P_over_P_delta_at_x_gnfw_500c(double x_asked,
                                                      double P0GNFW,
                                                      double alphaGNFW,
                                                      double betaGNFW,
                                                      double gammaGNFW,
                                                      double c500,
                                                      struct background * pba,
                                                      struct class_sz_structure * tsz);


struct Parameters_for_integrand_gas_density_profile_2h{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  struct perturbs * ppt;
  double * pvecback;
  double * pvectsz;
  double z;
  double k;
};


struct Parameters_for_integrand_gas_pressure_profile_2h{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  struct perturbs * ppt;
  double * pvecback;
  double * pvectsz;
  double z;
  double k;
};


struct Parameters_for_integrand_kSZ2_X_at_theta{
struct nonlinear * pnl;
struct primordial * ppm;
struct class_sz_structure * pclass_sz;
struct background * pba;
double * Pvecback;
double * Pvectsz;
double theta;
int index_ell_3;
double * b_l1_l2_l_1d;
double * ln_ell;
};



struct Parameters_for_integrand_kSZ2_X{
struct nonlinear * pnl;
struct primordial * ppm;
struct class_sz_structure * pclass_sz;
struct background * pba;
double * Pvecback;
double * Pvectsz;
int index_ell_3;
double * b_l1_l2_l_1d;
double * ln_ell;
};


struct Parameters_for_integrand_kSZ2_X_lensing_term_at_theta{
struct nonlinear * pnl;
struct primordial * ppm;
struct class_sz_structure * pclass_sz;
struct background * pba;
// double * Pvecback;
// double * Pvectsz;
double theta;
int index_ell;
double * integrand_l_lprime_phi;
double * ln_ellprime;
};



struct Parameters_for_integrand_kSZ2_X_lensing_term{
struct nonlinear * pnl;
struct primordial * ppm;
struct class_sz_structure * pclass_sz;
struct background * pba;
// double * Pvecback;
// double * Pvectsz;
int index_ell;
double * integrand_l_lprime_phi;
double * ln_ellprime;
};





#ifdef __cplusplus
}
#endif

#endif
