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
// Compute the Givens angle (and half-angle) 
//###########################################################

Ssh.f = SS21.f*Sone_half.f;
Stmp5.f = SS11.f - SS22.f;

Stmp2.f = Ssh.f*Ssh.f;
Stmp1.ui = (Stmp2.f >= Stiny_number.f ? 0xffffffff : 0);
Ssh.ui = Stmp1.ui&Ssh.ui;
Sch.ui = Stmp1.ui&Stmp5.ui;
Stmp2.ui = ~Stmp1.ui&Sone.ui;
Sch.ui = Sch.ui | Stmp2.ui;

Stmp1.f = Ssh.f*Ssh.f;
Stmp2.f = Sch.f*Sch.f;
Stmp3.f = Stmp1.f + Stmp2.f;
Stmp4.f = rsqrt(Stmp3.f);

#ifdef USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION
Ss.f = Stmp4.f*Sone_half.f;
Sc.f = Stmp4.f*Ss.f;
Sc.f = Stmp4.f*Sc.f;
Sc.f = Stmp3.f*Sc.f;
Stmp4.f = Stmp4.f + Ss.f;
Stmp4.f = Stmp4.f - Sc.f;
#endif

Ssh.f = Stmp4.f*Ssh.f;
Sch.f = Stmp4.f*Sch.f;

Stmp1.f = Sfour_gamma_squared.f*Stmp1.f;
Stmp1.ui = (Stmp2.f <= Stmp1.f ? 0xffffffff : 0);

Stmp2.ui = Ssine_pi_over_eight.ui&Stmp1.ui;
Ssh.ui = ~Stmp1.ui&Ssh.ui;
Ssh.ui = Ssh.ui | Stmp2.ui;
Stmp2.ui = Scosine_pi_over_eight.ui&Stmp1.ui;
Sch.ui = ~Stmp1.ui&Sch.ui;
Sch.ui = Sch.ui | Stmp2.ui;

Stmp1.f = Ssh.f*Ssh.f;
Stmp2.f = Sch.f*Sch.f;
Sc.f = Stmp2.f - Stmp1.f;
Ss.f = Sch.f*Ssh.f;
Ss.f = Ss.f + Ss.f;

//##########################################################
// Perform the actual Givens conjugation
//##########################################################

#ifndef USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION
Stmp3.f = Stmp1.f + Stmp2.f;
SS33.f = SS33.f*Stmp3.f;
SS31.f = SS31.f*Stmp3.f;
SS32.f = SS32.f*Stmp3.f;
SS33.f = SS33.f*Stmp3.f;
#endif

Stmp1.f = Ss.f*SS31.f;
Stmp2.f = Ss.f*SS32.f;
SS31.f = Sc.f*SS31.f;
SS32.f = Sc.f*SS32.f;
SS31.f = Stmp2.f + SS31.f;
SS32.f = SS32.f - Stmp1.f;

Stmp2.f = Ss.f*Ss.f;
Stmp1.f = SS22.f*Stmp2.f;
Stmp3.f = SS11.f*Stmp2.f;
Stmp4.f = Sc.f*Sc.f;
SS11.f = SS11.f*Stmp4.f;
SS22.f = SS22.f*Stmp4.f;
SS11.f = SS11.f + Stmp1.f;
SS22.f = SS22.f + Stmp3.f;
Stmp4.f = Stmp4.f - Stmp2.f;
Stmp2.f = SS21.f + SS21.f;
SS21.f = SS21.f*Stmp4.f;
Stmp4.f = Sc.f*Ss.f;
Stmp2.f = Stmp2.f*Stmp4.f;
Stmp5.f = Stmp5.f*Stmp4.f;
SS11.f = SS11.f + Stmp2.f;
SS21.f = SS21.f - Stmp5.f;
SS22.f = SS22.f - Stmp2.f;

//##########################################################
// Compute the cumulative rotation, in quaternion form
//##########################################################

Stmp1.f = Ssh.f*Sqvvx.f;
Stmp2.f = Ssh.f*Sqvvy.f;
Stmp3.f = Ssh.f*Sqvvz.f;
Ssh.f = Ssh.f*Sqvs.f;

Sqvs.f = Sch.f*Sqvs.f;
Sqvvx.f = Sch.f*Sqvvx.f;
Sqvvy.f = Sch.f*Sqvvy.f;
Sqvvz.f = Sch.f*Sqvvz.f;

SQVVZ.f = SQVVZ.f + Ssh.f;
Sqvs.f = Sqvs.f - STMP3.f;
SQVVX.f = SQVVX.f + STMP2.f;
SQVVY.f = SQVVY.f - STMP1.f;
