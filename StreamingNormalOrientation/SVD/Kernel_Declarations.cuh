//#####################################################################
// Copyright (c) 2010-2011, Eftychios Sifakis.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or
//     other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
// BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//#####################################################################

//###########################################################
// Local variable declarations
//###########################################################

#ifdef PRINT_DEBUGGING_OUTPUT

#ifdef USE_SSE_IMPLEMENTATION
float buf[4];
float A11, A21, A31, A12, A22, A32, A13, A23, A33;
float S11, S21, S31, S22, S32, S33;
#ifdef COMPUTE_V_AS_QUATERNION
float QVS, QVVX, QVVY, QVVZ;
#endif
#ifdef COMPUTE_V_AS_MATRIX
float V11, V21, V31, V12, V22, V32, V13, V23, V33;
#endif
#ifdef COMPUTE_U_AS_QUATERNION
float QUS, QUVX, QUVY, QUVZ;
#endif
#ifdef COMPUTE_U_AS_MATRIX
float U11, U21, U31, U12, U22, U32, U13, U23, U33;
#endif
#endif

#ifdef USE_AVX_IMPLEMENTATION
float buf[8];
float A11, A21, A31, A12, A22, A32, A13, A23, A33;
float S11, S21, S31, S22, S32, S33;
#ifdef COMPUTE_V_AS_QUATERNION
float QVS, QVVX, QVVY, QVVZ;
#endif
#ifdef COMPUTE_V_AS_MATRIX
float V11, V21, V31, V12, V22, V32, V13, V23, V33;
#endif
#ifdef COMPUTE_U_AS_QUATERNION
float QUS, QUVX, QUVY, QUVZ;
#endif
#ifdef COMPUTE_U_AS_MATRIX
float U11, U21, U31, U12, U22, U32, U13, U23, U33;
#endif
#endif

#endif

const float Four_Gamma_Squared = sqrt(8.) + 3.;
const float Sine_Pi_Over_Eight = .5*sqrt(2. - sqrt(2.));
const float Cosine_Pi_Over_Eight = .5*sqrt(2. + sqrt(2.));

union { float f; unsigned int ui; } Sfour_gamma_squared;
union { float f; unsigned int ui; } Ssine_pi_over_eight;
union { float f; unsigned int ui; } Scosine_pi_over_eight;
union { float f; unsigned int ui; } Sone_half;
union { float f; unsigned int ui; } Sone;
union { float f; unsigned int ui; } Stiny_number;
union { float f; unsigned int ui; } Ssmall_number;

Sfour_gamma_squared.f = Four_Gamma_Squared;
Ssine_pi_over_eight.f = Sine_Pi_Over_Eight;
Scosine_pi_over_eight.f = Cosine_Pi_Over_Eight;
Sone_half.f = .5;
Sone.f = 1.;
Stiny_number.f = 1.e-20;
Ssmall_number.f = 1.e-12;

union { float f; unsigned int ui; } Sa11;
union { float f; unsigned int ui; } Sa21;
union { float f; unsigned int ui; } Sa31;
union { float f; unsigned int ui; } Sa12;
union { float f; unsigned int ui; } Sa22;
union { float f; unsigned int ui; } Sa32;
union { float f; unsigned int ui; } Sa13;
union { float f; unsigned int ui; } Sa23;
union { float f; unsigned int ui; } Sa33;

#ifdef COMPUTE_V_AS_MATRIX
union { float f; unsigned int ui; } Sv11;
union { float f; unsigned int ui; } Sv21;
union { float f; unsigned int ui; } Sv31;
union { float f; unsigned int ui; } Sv12;
union { float f; unsigned int ui; } Sv22;
union { float f; unsigned int ui; } Sv32;
union { float f; unsigned int ui; } Sv13;
union { float f; unsigned int ui; } Sv23;
union { float f; unsigned int ui; } Sv33;
#endif

#ifdef COMPUTE_V_AS_QUATERNION
union { float f; unsigned int ui; } Sqvs;
union { float f; unsigned int ui; } Sqvvx;
union { float f; unsigned int ui; } Sqvvy;
union { float f; unsigned int ui; } Sqvvz;
#endif

#ifdef COMPUTE_U_AS_MATRIX
union { float f; unsigned int ui; } Su11;
union { float f; unsigned int ui; } Su21;
union { float f; unsigned int ui; } Su31;
union { float f; unsigned int ui; } Su12;
union { float f; unsigned int ui; } Su22;
union { float f; unsigned int ui; } Su32;
union { float f; unsigned int ui; } Su13;
union { float f; unsigned int ui; } Su23;
union { float f; unsigned int ui; } Su33;
#endif

#ifdef COMPUTE_U_AS_QUATERNION
union { float f; unsigned int ui; } Squs;
union { float f; unsigned int ui; } Squvx;
union { float f; unsigned int ui; } Squvy;
union { float f; unsigned int ui; } Squvz;
#endif

union { float f; unsigned int ui; } Sc;
union { float f; unsigned int ui; } Ss;
union { float f; unsigned int ui; } Sch;
union { float f; unsigned int ui; } Ssh;
union { float f; unsigned int ui; } Stmp1;
union { float f; unsigned int ui; } Stmp2;
union { float f; unsigned int ui; } Stmp3;
union { float f; unsigned int ui; } Stmp4;
union { float f; unsigned int ui; } Stmp5;
