/* Thin C shim over libxc: hides xc_func_type allocation/layout from Rust. */
#include <stdlib.h>
#include <xc.h>

void *qc_xc_init(int id) {
    xc_func_type *f = (xc_func_type *)malloc(sizeof(xc_func_type));
    if (!f) return NULL;
    if (xc_func_init(f, id, XC_UNPOLARIZED) != 0) {
        free(f);
        return NULL;
    }
    return (void *)f;
}

/* LDA: energy density (exc, per particle) and potential (vxc). */
void qc_xc_lda_exc_vxc(void *f, int n, const double *rho, double *exc, double *vxc) {
    xc_lda_exc_vxc((const xc_func_type *)f, (size_t)n, rho, exc, vxc);
}

/* GGA: exc, ∂E/∂ρ (vrho), ∂E/∂σ (vsigma) where σ = |∇ρ|². */
void qc_xc_gga_exc_vxc(void *f, int n, const double *rho, const double *sigma,
                       double *exc, double *vrho, double *vsigma) {
    xc_gga_exc_vxc((const xc_func_type *)f, (size_t)n, rho, sigma, exc, vrho, vsigma);
}

/* Exact-exchange mixing coefficient of a hybrid functional (0 for pure DFT). */
double qc_xc_hyb_exx(void *f) {
    return xc_hyb_exx_coef((const xc_func_type *)f);
}

/* Range-separated (CAM) parameters: exact exchange = α/r + β·erf(ωr)/r. */
void qc_xc_cam_coef(void *f, double *omega, double *alpha, double *beta) {
    xc_hyb_cam_coef((const xc_func_type *)f, omega, alpha, beta);
}

/* Spin-polarized init (needed for the spin-resolved kernel). */
void *qc_xc_init_pol(int id) {
    xc_func_type *f = (xc_func_type *)malloc(sizeof(xc_func_type));
    if (!f) return NULL;
    if (xc_func_init(f, id, XC_POLARIZED) != 0) {
        free(f);
        return NULL;
    }
    return (void *)f;
}

/* LDA second derivative: v2rho2 = [f_uu, f_ud, f_dd] per point (polarized). */
void qc_xc_lda_fxc(void *f, int n, const double *rho, double *v2rho2) {
    xc_lda_fxc((const xc_func_type *)f, (size_t)n, rho, v2rho2);
}

void qc_xc_free(void *f) {
    if (f) {
        xc_func_end((xc_func_type *)f);
        free(f);
    }
}
