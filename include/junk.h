int DmSH(T* vecIn,int du,int dv, T* vecOut);
/// Spherical harmonic differentiation


int shAna(T* vecIn, T* vecOut);
/// Spherical harmonics transform

int shSyn(T* vecIn, T* vecOut);
/// Spherical harmonics inverse transform

int interpsh(T* vecIn, T* vecOut);
/// Interpolation function via spherical harmonics

int filtersh(T* vecIn, T* vecOut);
/// Filtering via spherical harmonics
