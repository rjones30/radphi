#ifndef DEFINE_SQR_ON_STANDARD_TYPES
#define DEFINE_SQR_ON_STANDARD_TYPES

inline unsigned int sqr(unsigned int x) { return x*x; }
inline Int_t sqr(Int_t x) { return x*x; }
inline Float_t sqr(Float_t x) { return x*x; }
inline Double_t sqr(Double_t x) { return x*x; }
//inline LDouble_t sqr(LDouble_t x) { return x*x; }
//inline Complex_t sqr(Complex_t x) { return x*x; }

#endif
