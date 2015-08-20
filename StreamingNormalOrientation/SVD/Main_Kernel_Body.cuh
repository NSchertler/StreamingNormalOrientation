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

// #define USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION
// #define PERFORM_STRICT_QUATERNION_RENORMALIZATION

{ // Begin block : Scope of qV (if not maintained)

#ifndef COMPUTE_V_AS_QUATERNION
	union { float f; unsigned int ui; } Sqvs;
	union { float f; unsigned int ui; } Sqvvx;
	union { float f; unsigned int ui; } Sqvvy;
	union { float f; unsigned int ui; } Sqvvz;
#endif

	{ // Begin block : Symmetric eigenanalysis

		union { float f; unsigned int ui; } Ss11;
		union { float f; unsigned int ui; } Ss21;
		union { float f; unsigned int ui; } Ss31;
		union { float f; unsigned int ui; } Ss22;
		union { float f; unsigned int ui; } Ss32;
		union { float f; unsigned int ui; } Ss33;

		Sqvs.f = 1.;
		Sqvvx.f = 0.;
		Sqvvy.f = 0.;
		Sqvvz.f = 0.;

		//######################################################
		// Compute normal equations matrix
		//######################################################

		Ss11.f = Sa11.f*Sa11.f;
		Stmp1.f = Sa21.f*Sa21.f;
		Ss11.f = Stmp1.f + Ss11.f;
		Stmp1.f = Sa31.f*Sa31.f;
		Ss11.f = Stmp1.f + Ss11.f;

		Ss21.f = Sa12.f*Sa11.f;
		Stmp1.f = Sa22.f*Sa21.f;
		Ss21.f = Stmp1.f + Ss21.f;
		Stmp1.f = Sa32.f*Sa31.f;
		Ss21.f = Stmp1.f + Ss21.f;

		Ss31.f = Sa13.f*Sa11.f;
		Stmp1.f = Sa23.f*Sa21.f;
		Ss31.f = Stmp1.f + Ss31.f;
		Stmp1.f = Sa33.f*Sa31.f;
		Ss31.f = Stmp1.f + Ss31.f;

		Ss22.f = Sa12.f*Sa12.f;
		Stmp1.f = Sa22.f*Sa22.f;
		Ss22.f = Stmp1.f + Ss22.f;
		Stmp1.f = Sa32.f*Sa32.f;
		Ss22.f = Stmp1.f + Ss22.f;

		Ss32.f = Sa13.f*Sa12.f;
		Stmp1.f = Sa23.f*Sa22.f;
		Ss32.f = Stmp1.f + Ss32.f;
		Stmp1.f = Sa33.f*Sa32.f;
		Ss32.f = Stmp1.f + Ss32.f;

		Ss33.f = Sa13.f*Sa13.f;
		Stmp1.f = Sa23.f*Sa23.f;
		Ss33.f = Stmp1.f + Ss33.f;
		Stmp1.f = Sa33.f*Sa33.f;
		Ss33.f = Stmp1.f + Ss33.f;

		//######################################################
		// Solve symmetric eigenproblem using Jacobi iteration
		//######################################################

		for (int sweep = 1; sweep <= 4; sweep++){

			// First Jacobi conjugation

#define SS11 Ss11
#define SS21 Ss21
#define SS31 Ss31
#define SS22 Ss22
#define SS32 Ss32
#define SS33 Ss33
#define SQVVX Sqvvx
#define SQVVY Sqvvy
#define SQVVZ Sqvvz
#define STMP1 Stmp1
#define STMP2 Stmp2
#define STMP3 Stmp3

#define VS11 Vs11
#define VS21 Vs21
#define VS31 Vs31
#define VS22 Vs22
#define VS32 Vs32
#define VS33 Vs33
#define VQVVX Vqvvx
#define VQVVY Vqvvy
#define VQVVZ Vqvvz
#define VTMP1 Vtmp1
#define VTMP2 Vtmp2
#define VTMP3 Vtmp3

#include "Jacobi_Conjugation_Kernel.cuh"

#undef SS11
#undef SS21
#undef SS31
#undef SS22
#undef SS32
#undef SS33
#undef SQVVX
#undef SQVVY
#undef SQVVZ
#undef STMP1
#undef STMP2
#undef STMP3

#undef VS11
#undef VS21
#undef VS31
#undef VS22
#undef VS32
#undef VS33
#undef VQVVX
#undef VQVVY
#undef VQVVZ
#undef VTMP1
#undef VTMP2
#undef VTMP3

			// Second Jacobi conjugation

#define SS11 Ss22
#define SS21 Ss32
#define SS31 Ss21
#define SS22 Ss33
#define SS32 Ss31
#define SS33 Ss11
#define SQVVX Sqvvy
#define SQVVY Sqvvz
#define SQVVZ Sqvvx
#define STMP1 Stmp2
#define STMP2 Stmp3
#define STMP3 Stmp1

#define VS11 Vs22
#define VS21 Vs32
#define VS31 Vs21
#define VS22 Vs33
#define VS32 Vs31
#define VS33 Vs11
#define VQVVX Vqvvy
#define VQVVY Vqvvz
#define VQVVZ Vqvvx
#define VTMP1 Vtmp2
#define VTMP2 Vtmp3
#define VTMP3 Vtmp1

#include "Jacobi_Conjugation_Kernel.cuh"

#undef SS11
#undef SS21
#undef SS31
#undef SS22
#undef SS32
#undef SS33
#undef SQVVX
#undef SQVVY
#undef SQVVZ
#undef STMP1
#undef STMP2
#undef STMP3

#undef VS11
#undef VS21
#undef VS31
#undef VS22
#undef VS32
#undef VS33
#undef VQVVX
#undef VQVVY
#undef VQVVZ
#undef VTMP1
#undef VTMP2
#undef VTMP3

			// Third Jacobi conjugation

#define SS11 Ss33
#define SS21 Ss31
#define SS31 Ss32
#define SS22 Ss11
#define SS32 Ss21
#define SS33 Ss22
#define SQVVX Sqvvz
#define SQVVY Sqvvx
#define SQVVZ Sqvvy
#define STMP1 Stmp3
#define STMP2 Stmp1
#define STMP3 Stmp2

#define VS11 Vs33
#define VS21 Vs31
#define VS31 Vs32
#define VS22 Vs11
#define VS32 Vs21
#define VS33 Vs22
#define VQVVX Vqvvz
#define VQVVY Vqvvx
#define VQVVZ Vqvvy
#define VTMP1 Vtmp3
#define VTMP2 Vtmp1
#define VTMP3 Vtmp2

#include "Jacobi_Conjugation_Kernel.cuh"

#undef SS11
#undef SS21
#undef SS31
#undef SS22
#undef SS32
#undef SS33
#undef SQVVX
#undef SQVVY
#undef SQVVZ
#undef STMP1
#undef STMP2
#undef STMP3

#undef VS11
#undef VS21
#undef VS31
#undef VS22
#undef VS32
#undef VS33
#undef VQVVX
#undef VQVVY
#undef VQVVZ
#undef VTMP1
#undef VTMP2
#undef VTMP3
		}


	} // End block : Symmetric eigenanalysis

	//###########################################################
	// Normalize quaternion for matrix V
	//###########################################################

#if !defined(USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION) || defined(PERFORM_STRICT_QUATERNION_RENORMALIZATION)
	Stmp2.f = Sqvs.f*Sqvs.f;
	Stmp1.f = Sqvvx.f*Sqvvx.f;
	Stmp2.f = Stmp1.f + Stmp2.f;
	Stmp1.f = Sqvvy.f*Sqvvy.f;
	Stmp2.f = Stmp1.f + Stmp2.f;
	Stmp1.f = Sqvvz.f*Sqvvz.f;
	Stmp2.f = Stmp1.f + Stmp2.f;

	Stmp1.f = rsqrt(Stmp2.f);
	Stmp4.f = Stmp1.f*Sone_half.f;
	Stmp3.f = Stmp1.f*Stmp4.f;
	Stmp3.f = Stmp1.f*Stmp3.f;
	Stmp3.f = Stmp2.f*Stmp3.f;
	Stmp1.f = Stmp1.f + Stmp4.f;
	Stmp1.f = Stmp1.f - Stmp3.f;

	Sqvs.f = Sqvs.f*Stmp1.f;
	Sqvvx.f = Sqvvx.f*Stmp1.f;
	Sqvvy.f = Sqvvy.f*Stmp1.f;
	Sqvvz.f = Sqvvz.f*Stmp1.f;
#endif

	{ // Begin block : Conjugation with V

		//###########################################################
		// Transform quaternion to matrix V
		//###########################################################

#ifndef COMPUTE_V_AS_MATRIX
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

		Stmp1.f = Sqvvx.f*Sqvvx.f;
		Stmp2.f = Sqvvy.f*Sqvvy.f;
		Stmp3.f = Sqvvz.f*Sqvvz.f;
		Sv11.f = Sqvs.f*Sqvs.f;
		Sv22.f = Sv11.f - Stmp1.f;
		Sv33.f = Sv22.f - Stmp2.f;
		Sv33.f = Sv33.f + Stmp3.f;
		Sv22.f = Sv22.f + Stmp2.f;
		Sv22.f = Sv22.f - Stmp3.f;
		Sv11.f = Sv11.f + Stmp1.f;
		Sv11.f = Sv11.f - Stmp2.f;
		Sv11.f = Sv11.f - Stmp3.f;
		Stmp1.f = Sqvvx.f + Sqvvx.f;
		Stmp2.f = Sqvvy.f + Sqvvy.f;
		Stmp3.f = Sqvvz.f + Sqvvz.f;
		Sv32.f = Sqvs.f*Stmp1.f;
		Sv13.f = Sqvs.f*Stmp2.f;
		Sv21.f = Sqvs.f*Stmp3.f;
		Stmp1.f = Sqvvy.f*Stmp1.f;
		Stmp2.f = Sqvvz.f*Stmp2.f;
		Stmp3.f = Sqvvx.f*Stmp3.f;
		Sv12.f = Stmp1.f - Sv21.f;
		Sv23.f = Stmp2.f - Sv32.f;
		Sv31.f = Stmp3.f - Sv13.f;
		Sv21.f = Stmp1.f + Sv21.f;
		Sv32.f = Stmp2.f + Sv32.f;
		Sv13.f = Stmp3.f + Sv13.f;


		//###########################################################
		// Multiply (from the right) with V
		//###########################################################

		Stmp2.f = Sa12.f;
		Stmp3.f = Sa13.f;
		Sa12.f = Sv12.f*Sa11.f;
		Sa13.f = Sv13.f*Sa11.f;
		Sa11.f = Sv11.f*Sa11.f;
		Stmp1.f = Sv21.f*Stmp2.f;
		Sa11.f = Sa11.f + Stmp1.f;
		Stmp1.f = Sv31.f*Stmp3.f;
		Sa11.f = Sa11.f + Stmp1.f;
		Stmp1.f = Sv22.f*Stmp2.f;
		Sa12.f = Sa12.f + Stmp1.f;
		Stmp1.f = Sv32.f*Stmp3.f;
		Sa12.f = Sa12.f + Stmp1.f;
		Stmp1.f = Sv23.f*Stmp2.f;
		Sa13.f = Sa13.f + Stmp1.f;
		Stmp1.f = Sv33.f*Stmp3.f;
		Sa13.f = Sa13.f + Stmp1.f;

		Stmp2.f = Sa22.f;
		Stmp3.f = Sa23.f;
		Sa22.f = Sv12.f*Sa21.f;
		Sa23.f = Sv13.f*Sa21.f;
		Sa21.f = Sv11.f*Sa21.f;
		Stmp1.f = Sv21.f*Stmp2.f;
		Sa21.f = Sa21.f + Stmp1.f;
		Stmp1.f = Sv31.f*Stmp3.f;
		Sa21.f = Sa21.f + Stmp1.f;
		Stmp1.f = Sv22.f*Stmp2.f;
		Sa22.f = Sa22.f + Stmp1.f;
		Stmp1.f = Sv32.f*Stmp3.f;
		Sa22.f = Sa22.f + Stmp1.f;
		Stmp1.f = Sv23.f*Stmp2.f;
		Sa23.f = Sa23.f + Stmp1.f;
		Stmp1.f = Sv33.f*Stmp3.f;
		Sa23.f = Sa23.f + Stmp1.f;

		Stmp2.f = Sa32.f;
		Stmp3.f = Sa33.f;
		Sa32.f = Sv12.f*Sa31.f;
		Sa33.f = Sv13.f*Sa31.f;
		Sa31.f = Sv11.f*Sa31.f;
		Stmp1.f = Sv21.f*Stmp2.f;
		Sa31.f = Sa31.f + Stmp1.f;
		Stmp1.f = Sv31.f*Stmp3.f;
		Sa31.f = Sa31.f + Stmp1.f;
		Stmp1.f = Sv22.f*Stmp2.f;
		Sa32.f = Sa32.f + Stmp1.f;
		Stmp1.f = Sv32.f*Stmp3.f;
		Sa32.f = Sa32.f + Stmp1.f;
		Stmp1.f = Sv23.f*Stmp2.f;
		Sa33.f = Sa33.f + Stmp1.f;
		Stmp1.f = Sv33.f*Stmp3.f;
		Sa33.f = Sa33.f + Stmp1.f;

	} // End block : Conjugation with V

} // End block : Scope of qV (if not maintained)

//###########################################################
// Permute columns such that the singular values are sorted
//###########################################################

Stmp1.f = Sa11.f*Sa11.f;
Stmp4.f = Sa21.f*Sa21.f;
Stmp1.f = Stmp1.f + Stmp4.f;
Stmp4.f = Sa31.f*Sa31.f;
Stmp1.f = Stmp1.f + Stmp4.f;

Stmp2.f = Sa12.f*Sa12.f;
Stmp4.f = Sa22.f*Sa22.f;
Stmp2.f = Stmp2.f + Stmp4.f;
Stmp4.f = Sa32.f*Sa32.f;
Stmp2.f = Stmp2.f + Stmp4.f;

Stmp3.f = Sa13.f*Sa13.f;
Stmp4.f = Sa23.f*Sa23.f;
Stmp3.f = Stmp3.f + Stmp4.f;
Stmp4.f = Sa33.f*Sa33.f;
Stmp3.f = Stmp3.f + Stmp4.f;

// Swap columns 1-2 if necessary

Stmp4.ui = (Stmp1.f < Stmp2.f ? 0xffffffff : 0);
Stmp5.ui = Sa11.ui^Sa12.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sa11.ui = Sa11.ui^Stmp5.ui;
Sa12.ui = Sa12.ui^Stmp5.ui;

Stmp5.ui = Sa21.ui^Sa22.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sa21.ui = Sa21.ui^Stmp5.ui;
Sa22.ui = Sa22.ui^Stmp5.ui;

Stmp5.ui = Sa31.ui^Sa32.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sa31.ui = Sa31.ui^Stmp5.ui;
Sa32.ui = Sa32.ui^Stmp5.ui;

#ifdef COMPUTE_V_AS_MATRIX
Stmp5.ui = Sv11.ui^Sv12.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sv11.ui = Sv11.ui^Stmp5.ui;
Sv12.ui = Sv12.ui^Stmp5.ui;

Stmp5.ui = Sv21.ui^Sv22.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sv21.ui = Sv21.ui^Stmp5.ui;
Sv22.ui = Sv22.ui^Stmp5.ui;

Stmp5.ui = Sv31.ui^Sv32.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sv31.ui = Sv31.ui^Stmp5.ui;
Sv32.ui = Sv32.ui^Stmp5.ui;
#endif

Stmp5.ui = Stmp1.ui^Stmp2.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Stmp1.ui = Stmp1.ui^Stmp5.ui;
Stmp2.ui = Stmp2.ui^Stmp5.ui;

// If columns 1-2 have been swapped, negate 2nd column of A and V so that V is still a rotation

Stmp5.f = -2.;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Stmp4.f = 1.;
Stmp4.f = Stmp4.f + Stmp5.f;

Sa12.f = Sa12.f*Stmp4.f;
Sa22.f = Sa22.f*Stmp4.f;
Sa32.f = Sa32.f*Stmp4.f;

#ifdef COMPUTE_V_AS_MATRIX
Sv12.f = Sv12.f*Stmp4.f;
Sv22.f = Sv22.f*Stmp4.f;
Sv32.f = Sv32.f*Stmp4.f;
#endif

// If columns 1-2 have been swapped, also update quaternion representation of V (the quaternion may become un-normalized after this)

#ifdef COMPUTE_V_AS_QUATERNION
Stmp4.f = Stmp4.f*Sone_half.f;
Stmp4.f = Stmp4.f - Sone_half.f;

Stmp5.f = Stmp4.f*Sqvvz.f;
Stmp5.f = Stmp5.f + Sqvs.f;
Sqvs.f = Sqvs.f*Stmp4.f;
Sqvvz.f = Sqvvz.f - Sqvs.f;
Sqvs.f = Stmp5.f;

Stmp5.f = Stmp4.f*Sqvvx.f;
Stmp5.f = Stmp5.f + Sqvvy.f;
Sqvvy.f = Sqvvy.f*Stmp4.f;
Sqvvx.f = Sqvvx.f - Sqvvy.f;
Sqvvy.f = Stmp5.f;
#endif

// Swap columns 1-3 if necessary

Stmp4.ui = (Stmp1.f < Stmp3.f ? 0xffffffff : 0);
Stmp5.ui = Sa11.ui^Sa13.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sa11.ui = Sa11.ui^Stmp5.ui;
Sa13.ui = Sa13.ui^Stmp5.ui;

Stmp5.ui = Sa21.ui^Sa23.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sa21.ui = Sa21.ui^Stmp5.ui;
Sa23.ui = Sa23.ui^Stmp5.ui;

Stmp5.ui = Sa31.ui^Sa33.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sa31.ui = Sa31.ui^Stmp5.ui;
Sa33.ui = Sa33.ui^Stmp5.ui;

#ifdef COMPUTE_V_AS_MATRIX
Stmp5.ui = Sv11.ui^Sv13.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sv11.ui = Sv11.ui^Stmp5.ui;
Sv13.ui = Sv13.ui^Stmp5.ui;

Stmp5.ui = Sv21.ui^Sv23.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sv21.ui = Sv21.ui^Stmp5.ui;
Sv23.ui = Sv23.ui^Stmp5.ui;

Stmp5.ui = Sv31.ui^Sv33.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sv31.ui = Sv31.ui^Stmp5.ui;
Sv33.ui = Sv33.ui^Stmp5.ui;
#endif

Stmp5.ui = Stmp1.ui^Stmp3.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Stmp1.ui = Stmp1.ui^Stmp5.ui;
Stmp3.ui = Stmp3.ui^Stmp5.ui;

// If columns 1-3 have been swapped, negate 1st column of A and V so that V is still a rotation

Stmp5.f = -2.;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Stmp4.f = 1.;
Stmp4.f = Stmp4.f + Stmp5.f;

Sa11.f = Sa11.f*Stmp4.f;
Sa21.f = Sa21.f*Stmp4.f;
Sa31.f = Sa31.f*Stmp4.f;

#ifdef COMPUTE_V_AS_MATRIX
Sv11.f = Sv11.f*Stmp4.f;
Sv21.f = Sv21.f*Stmp4.f;
Sv31.f = Sv31.f*Stmp4.f;
#endif

// If columns 1-3 have been swapped, also update quaternion representation of V (the quaternion may become un-normalized after this)

#ifdef COMPUTE_V_AS_QUATERNION
Stmp4.f = Stmp4.f*Sone_half.f;
Stmp4.f = Stmp4.f - Sone_half.f;

Stmp5.f = Stmp4.f*Sqvvy.f;
Stmp5.f = Stmp5.f + Sqvs.f;
Sqvs.f = Sqvs.f*Stmp4.f;
Sqvvy.f = Sqvvy.f - Sqvs.f;
Sqvs.f = Stmp5.f;

Stmp5.f = Stmp4.f*Sqvvz.f;
Stmp5.f = Stmp5.f + Sqvvx.f;
Sqvvx.f = Sqvvx.f*Stmp4.f;
Sqvvz.f = Sqvvz.f - Sqvvx.f;
Sqvvx.f = Stmp5.f;
#endif

// Swap columns 2-3 if necessary

Stmp4.ui = (Stmp2.f < Stmp3.f ? 0xffffffff : 0);
Stmp5.ui = Sa12.ui^Sa13.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sa12.ui = Sa12.ui^Stmp5.ui;
Sa13.ui = Sa13.ui^Stmp5.ui;

Stmp5.ui = Sa22.ui^Sa23.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sa22.ui = Sa22.ui^Stmp5.ui;
Sa23.ui = Sa23.ui^Stmp5.ui;

Stmp5.ui = Sa32.ui^Sa33.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sa32.ui = Sa32.ui^Stmp5.ui;
Sa33.ui = Sa33.ui^Stmp5.ui;

#ifdef COMPUTE_V_AS_MATRIX
Stmp5.ui = Sv12.ui^Sv13.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sv12.ui = Sv12.ui^Stmp5.ui;
Sv13.ui = Sv13.ui^Stmp5.ui;

Stmp5.ui = Sv22.ui^Sv23.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sv22.ui = Sv22.ui^Stmp5.ui;
Sv23.ui = Sv23.ui^Stmp5.ui;

Stmp5.ui = Sv32.ui^Sv33.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Sv32.ui = Sv32.ui^Stmp5.ui;
Sv33.ui = Sv33.ui^Stmp5.ui;
#endif

Stmp5.ui = Stmp2.ui^Stmp3.ui;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Stmp2.ui = Stmp2.ui^Stmp5.ui;
Stmp3.ui = Stmp3.ui^Stmp5.ui;

// If columns 2-3 have been swapped, negate 3rd column of A and V so that V is still a rotation

Stmp5.f = -2.;
Stmp5.ui = Stmp5.ui&Stmp4.ui;
Stmp4.f = 1.;
Stmp4.f = Stmp4.f + Stmp5.f;

Sa13.f = Sa13.f*Stmp4.f;
Sa23.f = Sa23.f*Stmp4.f;
Sa33.f = Sa33.f*Stmp4.f;

#ifdef COMPUTE_V_AS_MATRIX
Sv13.f = Sv13.f*Stmp4.f;
Sv23.f = Sv23.f*Stmp4.f;
Sv33.f = Sv33.f*Stmp4.f;
#endif

// If columns 2-3 have been swapped, also update quaternion representation of V (the quaternion may become un-normalized after this)

#ifdef COMPUTE_V_AS_QUATERNION
Stmp4.f = Stmp4.f*Sone_half.f;
Stmp4.f = Stmp4.f - Sone_half.f;

Stmp5.f = Stmp4.f*Sqvvx.f;
Stmp5.f = Stmp5.f + Sqvs.f;
Sqvs.f = Sqvs.f*Stmp4.f;
Sqvvx.f = Sqvvx.f - Sqvs.f;
Sqvs.f = Stmp5.f;

Stmp5.f = Stmp4.f*Sqvvy.f;
Stmp5.f = Stmp5.f + Sqvvz.f;
Sqvvz.f = Sqvvz.f*Stmp4.f;
Sqvvy.f = Sqvvy.f - Sqvvz.f;
Sqvvz.f = Stmp5.f;
#endif


//###########################################################
// Re-normalize quaternion for matrix V
//###########################################################

#ifdef COMPUTE_V_AS_QUATERNION
Stmp2.f = Sqvs.f*Sqvs.f;
Stmp1.f = Sqvvx.f*Sqvvx.f;
Stmp2.f = Stmp1.f + Stmp2.f;
Stmp1.f = Sqvvy.f*Sqvvy.f;
Stmp2.f = Stmp1.f + Stmp2.f;
Stmp1.f = Sqvvz.f*Sqvvz.f;
Stmp2.f = Stmp1.f + Stmp2.f;
Stmp1.f = rsqrt(Stmp2.f;)

#ifdef PERFORM_STRICT_QUATERNION_RENORMALIZATION
Stmp4.f = Stmp1.f*Sone_half.f;
Stmp3.f = Stmp1.f*Stmp4.f;
Stmp3.f = Stmp1.f*Stmp3.f;
Stmp3.f = Stmp2.f*Stmp3.f;
Stmp1.f = Stmp1.f + Stmp4.f;
Stmp1.f = Stmp1.f - Stmp3.f;
#endif

Sqvs.f = Sqvs.f*Stmp1.f;
Sqvvx.f = Sqvvx.f*Stmp1.f;
Sqvvy.f = Sqvvy.f*Stmp1.f;
Sqvvz.f = Sqvvz.f*Stmp1.f;
#endif

//###########################################################
// Construct QR factorization of A*V (=U*D) using Givens rotations
//###########################################################

#ifdef COMPUTE_U_AS_MATRIX
Su11.f = 1.;
Su21.f = 0.;
Su31.f = 0.;
Su12.f = 0.;
Su22.f = 1.;
Su32.f = 0.;
Su13.f = 0.;
Su23.f = 0.;
Su33.f = 1.;
#endif

#ifdef COMPUTE_U_AS_QUATERNION
Squs.f = 1.;
Squvx.f = 0.;
Squvy.f = 0.;
Squvz.f = 0.;
#endif

// First Givens rotation

#define SAPIVOT Sa11
#define SANPIVOT Sa21
#define SA11 Sa11
#define SA21 Sa21
#define SA12 Sa12
#define SA22 Sa22
#define SA13 Sa13
#define SA23 Sa23
#define SU11 Su11
#define SU12 Su12
#define SU21 Su21
#define SU22 Su22
#define SU31 Su31
#define SU32 Su32

#define VAPIVOT Va11
#define VANPIVOT Va21
#define VA11 Va11
#define VA21 Va21
#define VA12 Va12
#define VA22 Va22
#define VA13 Va13
#define VA23 Va23
#define VU11 Vu11
#define VU12 Vu12
#define VU21 Vu21
#define VU22 Vu22
#define VU31 Vu31
#define VU32 Vu32

#include "Givens_QR_Factorization_Kernel.cuh"

#undef SAPIVOT
#undef SANPIVOT
#undef SA11
#undef SA21
#undef SA12
#undef SA22
#undef SA13
#undef SA23
#undef SU11
#undef SU12
#undef SU21
#undef SU22
#undef SU31
#undef SU32

#undef VAPIVOT
#undef VANPIVOT
#undef VA11
#undef VA21
#undef VA12
#undef VA22
#undef VA13
#undef VA23
#undef VU11
#undef VU12
#undef VU21
#undef VU22
#undef VU31
#undef VU32

// Update quaternion representation of U

#ifdef COMPUTE_U_AS_QUATERNION
Squs.f = Sch.f;
Squvx.f = 0.;
Squvy.f = 0.;
Squvz.f = Ssh.f;
#endif

// Second Givens rotation

#define SAPIVOT Sa11
#define SANPIVOT Sa31
#define SA11 Sa11
#define SA21 Sa31
#define SA12 Sa12
#define SA22 Sa32
#define SA13 Sa13
#define SA23 Sa33
#define SU11 Su11
#define SU12 Su13
#define SU21 Su21
#define SU22 Su23
#define SU31 Su31
#define SU32 Su33

#define VAPIVOT Va11
#define VANPIVOT Va31
#define VA11 Va11
#define VA21 Va31
#define VA12 Va12
#define VA22 Va32
#define VA13 Va13
#define VA23 Va33
#define VU11 Vu11
#define VU12 Vu13
#define VU21 Vu21
#define VU22 Vu23
#define VU31 Vu31
#define VU32 Vu33

#include "Givens_QR_Factorization_Kernel.cuh"

#undef SAPIVOT
#undef SANPIVOT
#undef SA11
#undef SA21
#undef SA12
#undef SA22
#undef SA13
#undef SA23
#undef SU11
#undef SU12
#undef SU21
#undef SU22
#undef SU31
#undef SU32

#undef VAPIVOT
#undef VANPIVOT
#undef VA11
#undef VA21
#undef VA12
#undef VA22
#undef VA13
#undef VA23
#undef VU11
#undef VU12
#undef VU21
#undef VU22
#undef VU31
#undef VU32

// Update quaternion representation of U

#ifdef COMPUTE_U_AS_QUATERNION
Squvx.f = Ssh.f*Squvz.f;
Ssh.f = Ssh.f*Squs.f;
Squvy.f = Squvy.f - Ssh.f;
Squs.f = Sch.f*Squs.f;
Squvz.f = Sch.f*Squvz.f;
#endif

// Third Givens rotation

#define SAPIVOT Sa22
#define SANPIVOT Sa32
#define SA11 Sa21
#define SA21 Sa31
#define SA12 Sa22
#define SA22 Sa32
#define SA13 Sa23
#define SA23 Sa33
#define SU11 Su12
#define SU12 Su13
#define SU21 Su22
#define SU22 Su23
#define SU31 Su32
#define SU32 Su33

#define VAPIVOT Va22
#define VANPIVOT Va32
#define VA11 Va21
#define VA21 Va31
#define VA12 Va22
#define VA22 Va32
#define VA13 Va23
#define VA23 Va33
#define VU11 Vu12
#define VU12 Vu13
#define VU21 Vu22
#define VU22 Vu23
#define VU31 Vu32
#define VU32 Vu33

#include "Givens_QR_Factorization_Kernel.cuh"

#undef SAPIVOT
#undef SANPIVOT
#undef SA11
#undef SA21
#undef SA12
#undef SA22
#undef SA13
#undef SA23
#undef SU11
#undef SU12
#undef SU21
#undef SU22
#undef SU31
#undef SU32

#undef VAPIVOT
#undef VANPIVOT
#undef VA11
#undef VA21
#undef VA12
#undef VA22
#undef VA13
#undef VA23
#undef VU11
#undef VU12
#undef VU21
#undef VU22
#undef VU31
#undef VU32

// Update quaternion representation of U

#ifdef COMPUTE_U_AS_QUATERNION
Stmp1.f = Ssh.f*Squvx.f;
Stmp2.f = Ssh.f*Squvy.f;
Stmp3.f = Ssh.f*Squvz.f;
Ssh.f = Ssh.f*Squs.f;
Squs.f = Sch.f*Squs.f;
Squvx.f = Sch.f*Squvx.f;
Squvy.f = Sch.f*Squvy.f;
Squvz.f = Sch.f*Squvz.f;
Squvx.f = Squvx.f + Ssh.f;
Squs.f = Squs.f - Stmp1.f;
Squvy.f = Squvy.f + Stmp3.f;
Squvz.f = Squvz.f - Stmp2.f;
#endif