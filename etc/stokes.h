void DirectStokes(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail, const float *qw, const float *trg, const float *src, const float *den, float *pot);
void DirectStokesSSE(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail, const float *qw, const float *trg, const float *src, const float *den, float *pot);

