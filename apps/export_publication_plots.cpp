// patched section only
#pragma omp parallel for
for (int i = 0; i < args.n_mt; ++i) {
    const double mt = mt_raw[i];
    prim[i] = dN_dmt_primordial(mt, 0.0, p.mu, p.mass, p.g);
    dirac[i] = p.dirac_mt(mt, 0.0);
    bw[i] = p.bw_mt(mt, 0.0);
    ps[i] = p.ps_mt(mt, 0.0);
    total_dirac[i] = prim[i] + dirac[i];
    total_bw[i] = prim[i] + bw[i];
    total_ps[i] = prim[i] + ps[i];
}
