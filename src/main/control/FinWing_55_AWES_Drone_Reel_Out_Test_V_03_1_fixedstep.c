/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep.c
 *
 * Code generated for Simulink model 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep'.
 *
 * Model version                  : 1.155
 * Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
 * C/C++ source code generated on : Tue May  5 13:04:57 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep.h"
#include "FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_private.h"
#include "common/time.h"
#include "flight/imu.h"
/* Block signals (default storage) */
B_FinWing_55_AWES_Drone_Reel__T FinWing_55_AWES_Drone_Reel_Ou_B;

/* Continuous states */
X_FinWing_55_AWES_Drone_Reel__T FinWing_55_AWES_Drone_Reel_Ou_X;

/* Periodic continuous states */
PeriodicIndX_FinWing_55_AWES__T FinWing_55_AWES_Dr_PeriodicIndX;
PeriodicRngX_FinWing_55_AWES__T FinWing_55_AWES_Dr_PeriodicRngX;

/* Block states (default storage) */
DW_FinWing_55_AWES_Drone_Reel_T FinWing_55_AWES_Drone_Reel_O_DW;

/* External inputs (root inport signals with default storage) */
ExtU_FinWing_55_AWES_Drone_Re_T FinWing_55_AWES_Drone_Reel_Ou_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_FinWing_55_AWES_Drone_Re_T FinWing_55_AWES_Drone_Reel_Ou_Y;

/* Real-time model */
RT_MODEL_FinWing_55_AWES_Dron_T FinWing_55_AWES_Drone_Reel_O_M_;
RT_MODEL_FinWing_55_AWES_Dron_T *const FinWing_55_AWES_Drone_Reel_O_M =
  &FinWing_55_AWES_Drone_Reel_O_M_;
uint32_T plook_binx(real_T u, const real_T bp[], uint32_T maxIndex, real_T
                    *fraction)
{
  uint32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0U;
    *fraction = (u - bp[0U]) / (bp[1U] - bp[0U]);
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_u32d(u, bp, maxIndex >> 1U, maxIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex - 1U;
    *fraction = (u - bp[maxIndex - 1U]) / (bp[maxIndex] - bp[maxIndex - 1U]);
  }

  return bpIndex;
}

real_T intrp1d_l_pw(uint32_T bpIndex, real_T frac, const real_T table[])
{
  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  return (table[bpIndex + 1U] - table[bpIndex]) * frac + table[bpIndex];
}

int32_T plook_s32dd_binx(real_T u, const real_T bp[], uint32_T maxIndex, real_T *
  fraction)
{
  int32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0;
    *fraction = (u - bp[0U]) / (bp[1U] - bp[0U]);
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_s32d(u, bp, maxIndex >> 1U, maxIndex);
    *fraction = (u - bp[(uint32_T)bpIndex]) / (bp[bpIndex + 1U] - bp[(uint32_T)
      bpIndex]);
  } else {
    bpIndex = (int32_T)(maxIndex - 1U);
    *fraction = (u - bp[maxIndex - 1U]) / (bp[maxIndex] - bp[maxIndex - 1U]);
  }

  return bpIndex;
}

real_T intrp2d_s32dl_pw(const int32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride)
{
  real_T yL_1d;
  uint32_T offset_1d;

  /* Column-major Interpolation 2-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  offset_1d = (uint32_T)(bpIndex[1U] * (int32_T)stride) + bpIndex[0U];
  yL_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
    table[offset_1d];
  offset_1d += stride;
  return (((table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
           table[offset_1d]) - yL_1d) * frac[1U] + yL_1d;
}

real_T intrp1d_s32dl_pw(int32_T bpIndex, real_T frac, const real_T table[])
{
  uint32_T offset_0d;

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  offset_0d = (uint32_T)bpIndex;
  return (table[offset_0d + 1U] - table[offset_0d]) * frac + table[offset_0d];
}

uint32_T binsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T bpIndex;
  uint32_T iRght;
  uint32_T bpIdx;

  /* Binary Search */
  bpIdx = startIndex;
  bpIndex = 0U;
  iRght = maxIndex;
  while (iRght - bpIndex > 1U) {
    if (u < bp[bpIdx]) {
      iRght = bpIdx;
    } else {
      bpIndex = bpIdx;
    }

    bpIdx = (iRght + bpIndex) >> 1U;
  }

  return bpIndex;
}

int32_T binsearch_s32d(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T iRght;
  uint32_T iLeft;
  uint32_T bpIdx;

  /* Binary Search */
  bpIdx = startIndex;
  iLeft = 0U;
  iRght = maxIndex;
  while (iRght - iLeft > 1U) {
    if (u < bp[bpIdx]) {
      iRght = bpIdx;
    } else {
      iLeft = bpIdx;
    }

    bpIdx = (iRght + iLeft) >> 1U;
  }

  return (int32_T)iLeft;
}

/* State reduction function */
void local_stateReduction(real_T* x, int_T* p, int_T n, real_T* r)
{
  int_T i, j;
  for (i = 0, j = 0; i < n; ++i, ++j) {
    int_T k = p[i];
    real_T lb = r[j++];
    real_T xk = x[k]-lb;
    real_T rk = r[j]-lb;
    int_T q = (int_T) floor(xk/rk);
    if (q) {
      x[k] = xk-q*rk+lb;
    }
  }
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 6;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_step(0);
  FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_step(0);
  FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  local_stateReduction(x, rtsiGetPeriodicContStateIndices(si), 1,
                       rtsiGetPeriodicContStateRanges(si));
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T u0_0;
  int32_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = atan2(u0_0, u1_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

/* Model step function */
void FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_step(timeUs_t currentTimeUs)
{
  real_T rtb_Va;
  real_T rtb_Atan;
  real_T rtb_u2rhoV2;
  uint32_T rtb_Prelookup_o1;
  real_T rtb_CM;
  real_T rtb_Gain_f;
  real_T rtb_UnaryMinus;
  int32_T rtb_Prelookup_o1_g;
  real_T rtb_CM_el;
  real_T rtb_gain2_f;
  real_T rtb_Product1_o;
  real_T rtb_Transpose[9];
  real_T rtb_sincos_o1_f[3];
  real_T rtb_sincos_o2[3];
  real_T rtb_Product_l[3];
  real_T rtb_sincos_o2_a;
  real_T frac[2];
  int32_T bpIndex[2];
  real_T rtb_Gain1_b_tmp;
  real_T rtb_Transpose_0[9];
  real_T rtb_Transpose_tmp;
  real_T rtb_sincos_o2_tmp;
  if (rtmIsMajorTimeStep(FinWing_55_AWES_Drone_Reel_O_M)) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
                          ((FinWing_55_AWES_Drone_Reel_O_M->Timing.clockTick0+1)*
      FinWing_55_AWES_Drone_Reel_O_M->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(FinWing_55_AWES_Drone_Reel_O_M)) {
    FinWing_55_AWES_Drone_Reel_O_M->Timing.t[0] = rtsiGetT
      (&FinWing_55_AWES_Drone_Reel_O_M->solverInfo);
  }

  /* Trigonometry: '<S21>/sincos' incorporates:
   *  Integrator: '<S13>/Theta'
   *  Trigonometry: '<S17>/sincos'
   */
  rtb_sincos_o2_tmp = cos(FinWing_55_AWES_Drone_Reel_Ou_X.theta);
  rtb_u2rhoV2 = sin(FinWing_55_AWES_Drone_Reel_Ou_X.theta);

  /* Fcn: '<S21>/Fcn11' incorporates:
   *  Trigonometry: '<S21>/sincos'
   */
  rtb_Transpose[0] = rtb_sincos_o2_tmp;

  /* Fcn: '<S21>/Fcn21' */
  rtb_Transpose[1] = 0.0 * rtb_u2rhoV2;

  /* Fcn: '<S21>/Fcn31' */
  rtb_Transpose[2] = rtb_u2rhoV2;

  /* Fcn: '<S21>/Fcn12' incorporates:
   *  Trigonometry: '<S21>/sincos'
   */
  rtb_Transpose[3] = rtb_sincos_o2_tmp * 0.0;

  /* Fcn: '<S21>/Fcn22' */
  rtb_Transpose[4] = 0.0 * rtb_u2rhoV2 * 0.0 + 1.0;

  /* Fcn: '<S21>/Fcn32' */
  rtb_Transpose[5] = rtb_u2rhoV2 * 0.0;

  /* Fcn: '<S21>/Fcn13' */
  rtb_Transpose[6] = -rtb_u2rhoV2;

  /* Fcn: '<S21>/Fcn23' incorporates:
   *  Trigonometry: '<S21>/sincos'
   */
  rtb_Transpose[7] = 0.0 * rtb_sincos_o2_tmp;

  /* Fcn: '<S21>/Fcn33' incorporates:
   *  Trigonometry: '<S21>/sincos'
   */
  rtb_Transpose[8] = rtb_sincos_o2_tmp;

  /* Math: '<S18>/Transpose' */
  for (rtb_Prelookup_o1_g = 0; rtb_Prelookup_o1_g < 3; rtb_Prelookup_o1_g++) {
    rtb_Transpose_0[3 * rtb_Prelookup_o1_g] = rtb_Transpose[rtb_Prelookup_o1_g];
    rtb_Transpose_0[3 * rtb_Prelookup_o1_g + 1] =
      rtb_Transpose[rtb_Prelookup_o1_g + 3];
    rtb_Transpose_0[3 * rtb_Prelookup_o1_g + 2] =
      rtb_Transpose[rtb_Prelookup_o1_g + 6];
  }

  memcpy(&rtb_Transpose[0], &rtb_Transpose_0[0], 9U * sizeof(real_T));

  /* End of Math: '<S18>/Transpose' */

  /* Integrator: '<S13>/U,w' incorporates:
   *  Constant: '<S13>/U0'
   *  Constant: '<S13>/w0'
   *  SignalConversion generated from: '<S13>/U,w'
   */
  if (FinWing_55_AWES_Drone_Reel_O_DW.Uw_IWORK != 0) {
    FinWing_55_AWES_Drone_Reel_Ou_X.u[0] = 19.997078706479069;
    FinWing_55_AWES_Drone_Reel_Ou_X.u[1] = 0.34182335625497184;
  }

  /* SignalConversion generated from: '<S14>/sincos' incorporates:
   *  Integrator: '<S13>/Theta'
   */
  rtb_sincos_o1_f[0] = 0.0;
  rtb_sincos_o1_f[1] = FinWing_55_AWES_Drone_Reel_Ou_X.theta;
  rtb_sincos_o1_f[2] = 0.0;
  for (rtb_Prelookup_o1_g = 0; rtb_Prelookup_o1_g < 3; rtb_Prelookup_o1_g++) {
    /* Trigonometry: '<S14>/sincos' */
    rtb_Product_l[rtb_Prelookup_o1_g] = cos(rtb_sincos_o1_f[rtb_Prelookup_o1_g]);

    /* Product: '<S18>/Product' incorporates:
     *  Integrator: '<S13>/U,w'
     *  SignalConversion generated from: '<S18>/Product'
     */
    rtb_sincos_o2[rtb_Prelookup_o1_g] = rtb_Transpose[rtb_Prelookup_o1_g + 6] *
      FinWing_55_AWES_Drone_Reel_Ou_X.u[1] + (rtb_Transpose[rtb_Prelookup_o1_g +
      3] * 0.0 + rtb_Transpose[rtb_Prelookup_o1_g] *
      FinWing_55_AWES_Drone_Reel_Ou_X.u[0]);

    /* Trigonometry: '<S14>/sincos' */
    rtb_sincos_o1_f[rtb_Prelookup_o1_g] = sin(rtb_sincos_o1_f[rtb_Prelookup_o1_g]);
  }

  /* Fcn: '<S14>/Fcn11' */
  rtb_Transpose[0] = rtb_Product_l[1] * rtb_Product_l[0];

  /* Fcn: '<S14>/Fcn21' incorporates:
   *  Fcn: '<S14>/Fcn22'
   */
  rtb_u2rhoV2 = rtb_sincos_o1_f[2] * rtb_sincos_o1_f[1];
  rtb_Transpose[1] = rtb_u2rhoV2 * rtb_Product_l[0] - rtb_Product_l[2] *
    rtb_sincos_o1_f[0];

  /* Fcn: '<S14>/Fcn31' incorporates:
   *  Fcn: '<S14>/Fcn32'
   */
  rtb_Transpose_tmp = rtb_Product_l[2] * rtb_sincos_o1_f[1];
  rtb_Transpose[2] = rtb_Transpose_tmp * rtb_Product_l[0] + rtb_sincos_o1_f[2] *
    rtb_sincos_o1_f[0];

  /* Fcn: '<S14>/Fcn12' */
  rtb_Transpose[3] = rtb_Product_l[1] * rtb_sincos_o1_f[0];

  /* Fcn: '<S14>/Fcn22' */
  rtb_Transpose[4] = rtb_u2rhoV2 * rtb_sincos_o1_f[0] + rtb_Product_l[2] *
    rtb_Product_l[0];

  /* Fcn: '<S14>/Fcn32' */
  rtb_Transpose[5] = rtb_Transpose_tmp * rtb_sincos_o1_f[0] - rtb_sincos_o1_f[2]
    * rtb_Product_l[0];

  /* Fcn: '<S14>/Fcn13' */
  rtb_Transpose[6] = -rtb_sincos_o1_f[1];

  /* Fcn: '<S14>/Fcn23' */
  rtb_Transpose[7] = rtb_sincos_o1_f[2] * rtb_Product_l[1];

  /* Fcn: '<S14>/Fcn33' */
  rtb_Transpose[8] = rtb_Product_l[2] * rtb_Product_l[1];

  /* Product: '<S1>/Product' */
  for (rtb_Prelookup_o1_g = 0; rtb_Prelookup_o1_g < 3; rtb_Prelookup_o1_g++) {
    rtb_Product_l[rtb_Prelookup_o1_g] = rtb_Transpose[rtb_Prelookup_o1_g + 6] *
      FinWing_55_AWES_Drone_Re_ConstB.Gain[2] +
      (rtb_Transpose[rtb_Prelookup_o1_g + 3] *
       FinWing_55_AWES_Drone_Re_ConstB.Gain[1] +
       rtb_Transpose[rtb_Prelookup_o1_g] * FinWing_55_AWES_Drone_Re_ConstB.Gain
       [0]);
  }

  /* End of Product: '<S1>/Product' */

  /* Gain: '<S26>/1//2rhoV^2' incorporates:
   *  Constant: '<S36>/Constant3'
   *  Integrator: '<S13>/U,w'
   *  Product: '<S26>/Product2'
   *  Product: '<S31>/Product'
   *  Product: '<S31>/Product1'
   *  Product: '<S31>/Product2'
   *  Sum: '<S2>/Sum'
   *  Sum: '<S31>/Sum'
   */
  rtb_u2rhoV2 = (((FinWing_55_AWES_Drone_Reel_Ou_X.u[0] - -10.0) *
                  (FinWing_55_AWES_Drone_Reel_Ou_X.u[0] - -10.0) + 0.0 * 0.0) +
                 FinWing_55_AWES_Drone_Reel_Ou_X.u[1] *
                 FinWing_55_AWES_Drone_Reel_Ou_X.u[1]) * 1.00649 * 0.5;

  /* Trigonometry: '<S28>/Incidence' incorporates:
   *  Integrator: '<S13>/U,w'
   */
  rtb_Atan = rt_atan2d_snf(FinWing_55_AWES_Drone_Reel_Ou_X.u[1],
    FinWing_55_AWES_Drone_Reel_Ou_X.u[0]);

  /* Trigonometry: '<S29>/Trigonometric Function' incorporates:
   *  Trigonometry: '<S30>/Trigonometric Function'
   *  Trigonometry: '<S8>/Trigonometric Function'
   */
  rtb_Transpose_tmp = sin(rtb_Atan);
  rtb_Gain1_b_tmp = cos(rtb_Atan);

  /* UnitConversion: '<S33>/Unit Conversion' */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  rtb_CM = 57.295779513082323 * rtb_Atan;

  /* PreLookup: '<S29>/Prelookup' */
  rtb_Prelookup_o1 = plook_binx(rtb_CM, FinWing_55_AWES_Drone_Re_ConstP.pooled3,
    9U, &rtb_CM);

  /* Interpolation_n-D: '<S29>/CD' */
  rtb_Gain_f = intrp1d_l_pw(rtb_Prelookup_o1, rtb_CM,
    FinWing_55_AWES_Drone_Re_ConstP.CD_Table);

  /* Interpolation_n-D: '<S29>/CL' */
  rtb_Product1_o = intrp1d_l_pw(rtb_Prelookup_o1, rtb_CM,
    FinWing_55_AWES_Drone_Re_ConstP.CL_Table);

  /* UnitConversion: '<S34>/Unit Conversion' */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  rtb_UnaryMinus = 57.295779513082323 * rtb_Atan;

  /* Interpolation_n-D: '<S30>/CD_el' incorporates:
   *  PreLookup: '<S30>/Prelookup1'
   */
  bpIndex[0] = plook_s32dd_binx(rtb_UnaryMinus,
    FinWing_55_AWES_Drone_Re_ConstP.pooled3, 9U, &rtb_UnaryMinus);

  /* Gain: '<S5>/Gain1' incorporates:
   *  Inport: '<Root>/ElevonPitch'
   */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  FinWing_55_AWES_Drone_Reel_Ou_U.ElevonPitch = (double) attitude.values.pitch ;
  
  rtb_CM_el = 0.017453292519943295 * (FinWing_55_AWES_Drone_Reel_Ou_U.ElevonPitch / 10.0);

  /* Saturate: '<Root>/Saturation1' */
  if (rtb_CM_el > 0.35) {
    rtb_CM_el = 0.35;
  } else {
    if (rtb_CM_el < -0.35) {
      rtb_CM_el = -0.35;
    }
  }

  /* End of Saturate: '<Root>/Saturation1' */

  /* UnitConversion: '<S35>/Unit Conversion' */
  rtb_CM_el *= 57.295779513082323;

  /* PreLookup: '<S30>/Prelookup' */
  rtb_Prelookup_o1_g = plook_s32dd_binx(rtb_CM_el,
    FinWing_55_AWES_Drone_Re_ConstP.Prelookup_BreakpointsData, 4U, &rtb_CM_el);

  /* Interpolation_n-D: '<S30>/CD_el' */
  frac[0] = rtb_UnaryMinus;
  frac[1] = rtb_CM_el;
  bpIndex[1] = rtb_Prelookup_o1_g;
  rtb_UnaryMinus = intrp2d_s32dl_pw(bpIndex, frac,
    FinWing_55_AWES_Drone_Re_ConstP.CD_el_Table, 10U);

  /* Gain: '<S30>/coeffAdjust1' */
  rtb_UnaryMinus = -rtb_UnaryMinus;

  /* Interpolation_n-D: '<S30>/CL_el' */
  rtb_gain2_f = intrp1d_s32dl_pw(rtb_Prelookup_o1_g, rtb_CM_el,
    FinWing_55_AWES_Drone_Re_ConstP.CL_el_Table);

  /* Sum: '<S2>/Fx_tot' incorporates:
   *  Gain: '<S24>/gain1'
   *  Gain: '<S27>/gain1'
   *  Gain: '<S29>/coeffAdjust1'
   *  Gain: '<S29>/coeffAdjust2'
   *  Gain: '<S30>/coeffAdjust'
   *  Product: '<S24>/Product2'
   *  Product: '<S27>/Product2'
   *  Product: '<S29>/Product2'
   *  Product: '<S29>/Product3'
   *  Product: '<S30>/Product2'
   *  Product: '<S30>/Product3'
   *  Sum: '<S29>/Sum1'
   *  Sum: '<S30>/Sum1'
   *  Trigonometry: '<S29>/Trigonometric Function'
   */
  rtb_sincos_o2_a = -(rtb_Gain1_b_tmp * rtb_Gain_f + rtb_Transpose_tmp *
                      -rtb_Product1_o) * rtb_u2rhoV2 * 0.0046 + (rtb_Gain1_b_tmp
    * rtb_UnaryMinus + rtb_Transpose_tmp * -rtb_gain2_f) * rtb_u2rhoV2 * 0.0046;

  /* Gain: '<S27>/gain2 ' incorporates:
   *  Product: '<S27>/Product3'
   *  Product: '<S30>/Product'
   *  Product: '<S30>/Product1'
   *  Sum: '<S30>/Sum'
   */
  rtb_gain2_f = (rtb_Transpose_tmp * rtb_UnaryMinus + rtb_Gain1_b_tmp *
                 rtb_gain2_f) * rtb_u2rhoV2 * 0.0046;

  /* Sum: '<S2>/Fz_tot' incorporates:
   *  Gain: '<S24>/gain2 '
   *  Gain: '<S29>/coeffAdjust'
   *  Product: '<S24>/Product3'
   *  Product: '<S29>/Product'
   *  Product: '<S29>/Product1'
   *  Sum: '<S29>/Sum'
   *  Trigonometry: '<S29>/Trigonometric Function'
   */
  rtb_UnaryMinus = -(rtb_Transpose_tmp * rtb_Gain_f + rtb_Gain1_b_tmp *
                     rtb_Product1_o) * rtb_u2rhoV2 * 0.0046 + rtb_gain2_f;

  /* Trigonometry: '<S1>/Atan' incorporates:
   *  Gain: '<S1>/Gain1'
   *  Integrator: '<S13>/Position'
   *  Product: '<S1>/Divide'
   */
  rtb_Atan = atan(-FinWing_55_AWES_Drone_Reel_Ou_X.Xe[1] /
                  FinWing_55_AWES_Drone_Reel_Ou_X.Xe[0]);

  /* Sum: '<S2>/Va' incorporates:
   *  Constant: '<Root>/Vreel '
   *  Constant: '<S36>/Constant3'
   *  Gain: '<S25>/Gain1'
   *  Gain: '<S2>/Gain'
   *  Product: '<S2>/Divide'
   *  Product: '<S2>/Product'
   *  Sqrt: '<S2>/Sqrt'
   *  Trigonometry: '<S2>/Cos'
   */
  rtb_Va = sqrt(2.0 * rtb_u2rhoV2 / 1.00649) * cos(0.017453292519943295 *
    rtb_Atan) - 10.0;

  /* Product: '<S2>/Divide1' incorporates:
   *  Constant: '<S2>/Add_Cable_Drag'
   *  Product: '<S2>/Product2'
   */
  rtb_Gain_f = rtb_Product1_o / (20.0 * rtb_Gain_f);

  /* Gain: '<S24>/gain2' incorporates:
   *  Math: '<S2>/Square'
   *  Math: '<S2>/Square1'
   *  Product: '<S2>/Product1'
   */
  rtb_Product1_o = rtb_Va * rtb_Va * FinWing_55_AWES_Drone_Re_ConstB.Gain1 *
    rtb_Product1_o * (rtb_Gain_f * rtb_Gain_f) * 0.0046;

  /* Product: '<S11>/Product' incorporates:
   *  Trigonometry: '<S11>/Cos'
   */
  rtb_Gain_f = rtb_Product1_o * cos(rtb_Atan);

  /* Product: '<S11>/Product1' incorporates:
   *  Trigonometry: '<S11>/Sin'
   */
  rtb_Product1_o *= sin(rtb_Atan);

  /* Integrator: '<S13>/q' */
  FinWing_55_AWES_Drone_Reel_Ou_B.q = FinWing_55_AWES_Drone_Reel_Ou_X.q;

  /* Sum: '<S13>/Sum1' incorporates:
   *  Constant: '<S13>/gravity'
   *  Constant: '<S16>/Constant'
   *  Gain: '<S11>/Gain'
   *  Gain: '<S13>/Matrix Gain'
   *  Integrator: '<S13>/Theta'
   *  Integrator: '<S13>/U,w'
   *  Product: '<S13>/Product'
   *  Product: '<S13>/Product1'
   *  Product: '<S17>/Product'
   *  Product: '<S8>/Product'
   *  Product: '<S8>/Product1'
   *  Sum: '<Root>/Sum2'
   *  Sum: '<Root>/Sum3'
   *  Sum: '<S13>/Sum'
   *  Sum: '<S1>/Sum2'
   *  Trigonometry: '<S17>/sincos'
   *  UnaryMinus: '<S17>/Unary Minus'
   */
  FinWing_55_AWES_Drone_Reel_Ou_B.Sum1[0] = ((((rtb_Gain1_b_tmp *
    FinWing_55_AWES_Drone_Re_ConstB.ThrustX + -rtb_Gain_f) + rtb_sincos_o2_a) +
    rtb_Product_l[0]) / 1.5 + -sin(FinWing_55_AWES_Drone_Reel_Ou_X.theta) * 0.0)
    + (0.0 * FinWing_55_AWES_Drone_Reel_Ou_X.u[0] +
       -FinWing_55_AWES_Drone_Reel_Ou_X.u[1]) *
    FinWing_55_AWES_Drone_Reel_Ou_B.q;
  FinWing_55_AWES_Drone_Reel_Ou_B.Sum1[1] = ((((rtb_Transpose_tmp *
    FinWing_55_AWES_Drone_Re_ConstB.ThrustX + rtb_Product1_o) + rtb_UnaryMinus)
    + rtb_Product_l[2]) / 1.5 + rtb_sincos_o2_tmp * 0.0) + (0.0 *
    FinWing_55_AWES_Drone_Reel_Ou_X.u[1] + FinWing_55_AWES_Drone_Reel_Ou_X.u[0])
    * FinWing_55_AWES_Drone_Reel_Ou_B.q;

  /* SignalConversion generated from: '<S13>/Position' */
  FinWing_55_AWES_Drone_Reel_Ou_B.TmpSignalConversionAtPositionIn[0] =
    rtb_sincos_o2[0];
  FinWing_55_AWES_Drone_Reel_Ou_B.TmpSignalConversionAtPositionIn[1] =
    rtb_sincos_o2[2];

  /* Interpolation_n-D: '<S29>/CM' */
  rtb_CM = intrp1d_l_pw(rtb_Prelookup_o1, rtb_CM,
                        FinWing_55_AWES_Drone_Re_ConstP.CM_Table);

  /* Interpolation_n-D: '<S30>/CM_el' */
  rtb_CM_el = intrp1d_s32dl_pw(rtb_Prelookup_o1_g, rtb_CM_el,
    FinWing_55_AWES_Drone_Re_ConstP.CM_el_Table);
  if (rtmIsMajorTimeStep(FinWing_55_AWES_Drone_Reel_O_M)) {
    /* Sum: '<Root>/Sum4' incorporates:
     *  Constant: '<S11>/TetherForce_M'
     */
    FinWing_55_AWES_Drone_Reel_Ou_B.Sum4 = 0.0;
  }

  /* Product: '<S15>/Product2' incorporates:
   *  Constant: '<S16>/Constant1'
   *  Constant: '<S16>/Constant2'
   *  Constant: '<S27>/Constant1'
   *  Gain: '<S24>/gain3'
   *  Gain: '<S27>/gain3'
   *  Product: '<S15>/Product3'
   *  Product: '<S24>/Product4'
   *  Product: '<S27>/Product4'
   *  Product: '<S27>/Product8'
   *  Sum: '<Root>/Sum1'
   *  Sum: '<S15>/Sum1'
   *  Sum: '<S27>/Sum'
   *  Sum: '<S2>/M_tot'
   */
  FinWing_55_AWES_Drone_Reel_Ou_B.Product2 = ((((rtb_u2rhoV2 * rtb_CM_el *
    8.51E-5 + rtb_gain2_f * 0.0) + rtb_u2rhoV2 * rtb_CM * 8.51E-5) +
    FinWing_55_AWES_Drone_Reel_Ou_B.Sum4) - 0.0 *
    FinWing_55_AWES_Drone_Reel_Ou_B.q) / 8.0;

  /* Outport: '<Root>/Altitude' incorporates:
   *  Gain: '<S1>/Gain1'
   *  Integrator: '<S13>/Position'
   */
  FinWing_55_AWES_Drone_Reel_Ou_Y.Altitude =
    -FinWing_55_AWES_Drone_Reel_Ou_X.Xe[1];

  /* S-Function (sfun_tstart): '<S9>/startTime' */
  /* S-Function Block (sfun_tstart): <S9>/startTime */
  FinWing_55_AWES_Drone_Reel_Ou_B.startTime = (0.0);

  /* S-Function (sfun_tstart): '<S10>/startTime' */
  /* S-Function Block (sfun_tstart): <S10>/startTime */
  FinWing_55_AWES_Drone_Reel_Ou_B.startTime_a = (0.0);
  if (rtmIsMajorTimeStep(FinWing_55_AWES_Drone_Reel_O_M)) {
    /* Update for Integrator: '<S13>/U,w' */
    FinWing_55_AWES_Drone_Reel_O_DW.Uw_IWORK = 0;
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(FinWing_55_AWES_Drone_Reel_O_M)) {
    rt_ertODEUpdateContinuousStates(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++FinWing_55_AWES_Drone_Reel_O_M->Timing.clockTick0;
    FinWing_55_AWES_Drone_Reel_O_M->Timing.t[0] = rtsiGetSolverStopTime
      (&FinWing_55_AWES_Drone_Reel_O_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      FinWing_55_AWES_Drone_Reel_O_M->Timing.clockTick1++;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_derivatives(void)
{
  XDot_FinWing_55_AWES_Drone_Re_T *_rtXdot;
  _rtXdot = ((XDot_FinWing_55_AWES_Drone_Re_T *)
             FinWing_55_AWES_Drone_Reel_O_M->derivs);

  /* Derivatives for Integrator: '<S13>/Theta' */
  _rtXdot->theta = FinWing_55_AWES_Drone_Reel_Ou_B.q;

  /* Derivatives for Integrator: '<S13>/U,w' */
  _rtXdot->u[0] = FinWing_55_AWES_Drone_Reel_Ou_B.Sum1[0];

  /* Derivatives for Integrator: '<S13>/Position' */
  _rtXdot->Xe[0] =
    FinWing_55_AWES_Drone_Reel_Ou_B.TmpSignalConversionAtPositionIn[0];

  /* Derivatives for Integrator: '<S13>/U,w' */
  _rtXdot->u[1] = FinWing_55_AWES_Drone_Reel_Ou_B.Sum1[1];

  /* Derivatives for Integrator: '<S13>/Position' */
  _rtXdot->Xe[1] =
    FinWing_55_AWES_Drone_Reel_Ou_B.TmpSignalConversionAtPositionIn[1];

  /* Derivatives for Integrator: '<S13>/q' */
  _rtXdot->q = FinWing_55_AWES_Drone_Reel_Ou_B.Product2;
}

/* Model initialize function */
void FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
                          &FinWing_55_AWES_Drone_Reel_O_M->Timing.simTimeStep);
    rtsiSetTPtr(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo, &rtmGetTPtr
                (FinWing_55_AWES_Drone_Reel_O_M));
    rtsiSetStepSizePtr(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
                       &FinWing_55_AWES_Drone_Reel_O_M->Timing.stepSize0);
    rtsiSetdXPtr(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
                 &FinWing_55_AWES_Drone_Reel_O_M->derivs);
    rtsiSetContStatesPtr(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo, (real_T **)
                         &FinWing_55_AWES_Drone_Reel_O_M->contStates);
    rtsiSetNumContStatesPtr(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
      &FinWing_55_AWES_Drone_Reel_O_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
      &FinWing_55_AWES_Drone_Reel_O_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr
      (&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
       &FinWing_55_AWES_Drone_Reel_O_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr
      (&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
       &FinWing_55_AWES_Drone_Reel_O_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
                          (&rtmGetErrorStatus(FinWing_55_AWES_Drone_Reel_O_M)));
    rtsiSetRTModelPtr(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
                      FinWing_55_AWES_Drone_Reel_O_M);
  }

  rtsiSetSimTimeStep(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,
                     MAJOR_TIME_STEP);
  FinWing_55_AWES_Drone_Reel_O_M->intgData.y =
    FinWing_55_AWES_Drone_Reel_O_M->odeY;
  FinWing_55_AWES_Drone_Reel_O_M->intgData.f[0] =
    FinWing_55_AWES_Drone_Reel_O_M->odeF[0];
  FinWing_55_AWES_Drone_Reel_O_M->intgData.f[1] =
    FinWing_55_AWES_Drone_Reel_O_M->odeF[1];
  FinWing_55_AWES_Drone_Reel_O_M->intgData.f[2] =
    FinWing_55_AWES_Drone_Reel_O_M->odeF[2];
  FinWing_55_AWES_Drone_Reel_O_M->contStates = ((X_FinWing_55_AWES_Drone_Reel__T
    *) &FinWing_55_AWES_Drone_Reel_Ou_X);
  FinWing_55_AWES_Drone_Reel_O_M->periodicContStateIndices = ((int_T*)
    FinWing_55_AWES_Dr_PeriodicIndX);
  FinWing_55_AWES_Drone_Reel_O_M->periodicContStateRanges = ((real_T*)
    FinWing_55_AWES_Dr_PeriodicRngX);
  rtsiSetSolverData(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo, (void *)
                    &FinWing_55_AWES_Drone_Reel_O_M->intgData);
  rtsiSetSolverName(&FinWing_55_AWES_Drone_Reel_O_M->solverInfo,"ode3");
  rtmSetTPtr(FinWing_55_AWES_Drone_Reel_O_M,
             &FinWing_55_AWES_Drone_Reel_O_M->Timing.tArray[0]);
  FinWing_55_AWES_Drone_Reel_O_M->Timing.stepSize0 = 0.001;
  rtmSetFirstInitCond(FinWing_55_AWES_Drone_Reel_O_M, 1);

  /* InitializeConditions for Integrator: '<S13>/Theta' */
  FinWing_55_AWES_Drone_Reel_Ou_X.theta = 0.017092;

  /* InitializeConditions for Integrator: '<S13>/U,w' */
  if (rtmIsFirstInitCond(FinWing_55_AWES_Drone_Reel_O_M)) {
    FinWing_55_AWES_Drone_Reel_Ou_X.u[0] = 19.997078706479069;
    FinWing_55_AWES_Drone_Reel_Ou_X.u[1] = 0.34182335625497184;
  }

  FinWing_55_AWES_Drone_Reel_O_DW.Uw_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S13>/U,w' */

  /* InitializeConditions for Integrator: '<S13>/Position' */
  FinWing_55_AWES_Drone_Reel_Ou_X.Xe[0] = 0.0;
  FinWing_55_AWES_Drone_Reel_Ou_X.Xe[1] = -500.0;

  /* InitializeConditions for Integrator: '<S13>/q' */
  FinWing_55_AWES_Drone_Reel_Ou_X.q = 0.0;

  /* InitializeConditions for root-level periodic continuous states */
  {
    int_T rootPeriodicContStateIndices[1] = { 0 };

    real_T rootPeriodicContStateRanges[2] = { -3.1415926535897931,
      3.1415926535897931 };

    (void) memcpy((void*)FinWing_55_AWES_Dr_PeriodicIndX,
                  rootPeriodicContStateIndices,
                  1*sizeof(int_T));
    (void) memcpy((void*)FinWing_55_AWES_Dr_PeriodicRngX,
                  rootPeriodicContStateRanges,
                  2*sizeof(real_T));
  }

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(FinWing_55_AWES_Drone_Reel_O_M)) {
    rtmSetFirstInitCond(FinWing_55_AWES_Drone_Reel_O_M, 0);
  }
}

/* Model terminate function */
void FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
