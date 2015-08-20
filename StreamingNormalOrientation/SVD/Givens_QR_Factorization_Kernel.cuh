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
// Compute the Givens half-angle, construct the Givens quaternion and the rotation sine/cosine (for the full angle)
//###########################################################

Ssh.f = SANPIVOT.f*SANPIVOT.f;
Ssh.ui = (Ssh.f >= Ssmall_number.f ? 0xffffffff : 0);
Ssh.ui = Ssh.ui&SANPIVOT.ui;

Stmp5.f = 0.;
Sch.f = Stmp5.f - SAPIVOT.f;
Sch.f = fmaxf(Sch.f, SAPIVOT.f);
Sch.f = fmaxf(Sch.f, Ssmall_number.f);
Stmp5.ui = (SAPIVOT.f >= Stmp5.f ? 0xffffffff : 0);

Stmp1.f = Sch.f*Sch.f;
Stmp2.f = Ssh.f*Ssh.f;
Stmp2.f = Stmp1.f + Stmp2.f;
Stmp1.f = rsqrt(Stmp2.f);

Stmp4.f = Stmp1.f*Sone_half.f;
Stmp3.f = Stmp1.f*Stmp4.f;
Stmp3.f = Stmp1.f*Stmp3.f;
Stmp3.f = Stmp2.f*Stmp3.f;
Stmp1.f = Stmp1.f + Stmp4.f;
Stmp1.f = Stmp1.f - Stmp3.f;
Stmp1.f = Stmp1.f*Stmp2.f;

Sch.f = Sch.f + Stmp1.f;

Stmp1.ui = ~Stmp5.ui&Ssh.ui;
Stmp2.ui = ~Stmp5.ui&Sch.ui;
Sch.ui = Stmp5.ui&Sch.ui;
Ssh.ui = Stmp5.ui&Ssh.ui;
Sch.ui = Sch.ui | Stmp1.ui;
Ssh.ui = Ssh.ui | Stmp2.ui;

Stmp1.f = Sch.f*Sch.f;
Stmp2.f = Ssh.f*Ssh.f;
Stmp2.f = Stmp1.f + Stmp2.f;
Stmp1.f = rsqrt(Stmp2.f);

Stmp4.f = Stmp1.f*Sone_half.f;
Stmp3.f = Stmp1.f*Stmp4.f;
Stmp3.f = Stmp1.f*Stmp3.f;
Stmp3.f = Stmp2.f*Stmp3.f;
Stmp1.f = Stmp1.f + Stmp4.f;
Stmp1.f = Stmp1.f - Stmp3.f;

Sch.f = Sch.f*Stmp1.f;
Ssh.f = Ssh.f*Stmp1.f;

Sc.f = Sch.f*Sch.f;
Ss.f = Ssh.f*Ssh.f;
Sc.f = Sc.f - Ss.f;
Ss.f = Ssh.f*Sch.f;
Ss.f = Ss.f + Ss.f;

//##########################################################
// Rotate matrix A
//##########################################################

Stmp1.f = Ss.f*SA11.f;
Stmp2.f = Ss.f*SA21.f;
SA11.f = Sc.f*SA11.f;
SA21.f = Sc.f*SA21.f;
SA11.f = SA11.f + Stmp2.f;
SA21.f = SA21.f - Stmp1.f;

Stmp1.f = Ss.f*SA12.f;
Stmp2.f = Ss.f*SA22.f;
SA12.f = Sc.f*SA12.f;
SA22.f = Sc.f*SA22.f;
SA12.f = SA12.f + Stmp2.f;
SA22.f = SA22.f - Stmp1.f;

Stmp1.f = Ss.f*SA13.f;
Stmp2.f = Ss.f*SA23.f;
SA13.f = Sc.f*SA13.f;
SA23.f = Sc.f*SA23.f;
SA13.f = SA13.f + Stmp2.f;
SA23.f = SA23.f - Stmp1.f;

//##########################################################
// Update matrix U
//##########################################################

#ifdef COMPUTE_U_AS_MATRIX
Stmp1.f = Ss.f*SU11.f;
Stmp2.f = Ss.f*SU12.f;
SU11.f = Sc.f*SU11.f;
SU12.f = Sc.f*SU12.f;
SU11.f = SU11.f + Stmp2.f;
SU12.f = SU12.f - Stmp1.f;

Stmp1.f = Ss.f*SU21.f;
Stmp2.f = Ss.f*SU22.f;
SU21.f = Sc.f*SU21.f;
SU22.f = Sc.f*SU22.f;
SU21.f = SU21.f + Stmp2.f;
SU22.f = SU22.f - Stmp1.f;

Stmp1.f = Ss.f*SU31.f;
Stmp2.f = Ss.f*SU32.f;
SU31.f = Sc.f*SU31.f;
SU32.f = Sc.f*SU32.f;
SU31.f = SU31.f + Stmp2.f;
SU32.f = SU32.f - Stmp1.f;
#endif
