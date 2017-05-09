#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ecdm_abcount_old(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_abcounts(SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_bijectionvector(SEXP);
extern SEXP ecdm_ClassbyQmat(SEXP);
extern SEXP ecdm_cond_threshold(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_dina_Gibbs_Q(SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_ETAmat(SEXP, SEXP, SEXP);
extern SEXP ecdm_ETAmat_nok(SEXP);
extern SEXP ecdm_ETAmat_nok_one_m_ac(SEXP);
extern SEXP ecdm_identify_check(SEXP);
extern SEXP ecdm_inv_bijectionvector(SEXP, SEXP);
extern SEXP ecdm_llj(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_lnlik_dina(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_lnlik_dina_condclass(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_OddsRatio(SEXP, SEXP, SEXP);
extern SEXP ecdm_parm_update_nomiss(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_pYit(SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_pYjeq1(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_random_Q(SEXP, SEXP);
extern SEXP ecdm_sim_Y_dina(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_updateQ_DINA(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ecdm_updateQ_DINA_new(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ecdm_abcount_old",          (DL_FUNC) &ecdm_abcount_old,           5},
    {"ecdm_abcounts",             (DL_FUNC) &ecdm_abcounts,              4},
    {"ecdm_bijectionvector",      (DL_FUNC) &ecdm_bijectionvector,       1},
    {"ecdm_ClassbyQmat",          (DL_FUNC) &ecdm_ClassbyQmat,           1},
    {"ecdm_cond_threshold",       (DL_FUNC) &ecdm_cond_threshold,       11},
    {"ecdm_dina_Gibbs_Q",         (DL_FUNC) &ecdm_dina_Gibbs_Q,          4},
    {"ecdm_ETAmat",               (DL_FUNC) &ecdm_ETAmat,                3},
    {"ecdm_ETAmat_nok",           (DL_FUNC) &ecdm_ETAmat_nok,            1},
    {"ecdm_ETAmat_nok_one_m_ac",  (DL_FUNC) &ecdm_ETAmat_nok_one_m_ac,   1},
    {"ecdm_identify_check",       (DL_FUNC) &ecdm_identify_check,        1},
    {"ecdm_inv_bijectionvector",  (DL_FUNC) &ecdm_inv_bijectionvector,   2},
    {"ecdm_llj",                  (DL_FUNC) &ecdm_llj,                   6},
    {"ecdm_lnlik_dina",           (DL_FUNC) &ecdm_lnlik_dina,            8},
    {"ecdm_lnlik_dina_condclass", (DL_FUNC) &ecdm_lnlik_dina_condclass,  8},
    {"ecdm_OddsRatio",            (DL_FUNC) &ecdm_OddsRatio,             3},
    {"ecdm_parm_update_nomiss",   (DL_FUNC) &ecdm_parm_update_nomiss,   10},
    {"ecdm_pYit",                 (DL_FUNC) &ecdm_pYit,                  4},
    {"ecdm_pYjeq1",               (DL_FUNC) &ecdm_pYjeq1,                5},
    {"ecdm_random_Q",             (DL_FUNC) &ecdm_random_Q,              2},
    {"ecdm_sim_Y_dina",           (DL_FUNC) &ecdm_sim_Y_dina,            6},
    {"ecdm_updateQ_DINA",         (DL_FUNC) &ecdm_updateQ_DINA,          5},
    {"ecdm_updateQ_DINA_new",     (DL_FUNC) &ecdm_updateQ_DINA_new,     12},
    {NULL, NULL, 0}
};

void R_init_ecdm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
