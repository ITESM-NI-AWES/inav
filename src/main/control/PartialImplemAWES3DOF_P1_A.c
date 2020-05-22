/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: PartialImplemAWES3DOF_P1_A.c
 *
 * Code generated for Simulink model 'PartialImplemAWES3DOF_P1_A'.
 *
 * Model version                  : 1.154
 * Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
 * C/C++ source code generated on : Fri May 22 17:37:58 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#include "PartialImplemAWES3DOF_P1_A.h"
#include "PartialImplemAWES3DOF_P1_A_private.h"

/* Block signals (default storage) */
B_PartialImplemAWES3DOF_P1_A_T PartialImplemAWES3DOF_P1_A_B;

/* Continuous states */
X_PartialImplemAWES3DOF_P1_A_T PartialImplemAWES3DOF_P1_A_X;

/* Periodic continuous states */
PeriodicIndX_PartialImplemAWE_T PartialImplemAWES3_PeriodicIndX;
PeriodicRngX_PartialImplemAWE_T PartialImplemAWES3_PeriodicRngX;

/* Block states (default storage) */
DW_PartialImplemAWES3DOF_P1_A_T PartialImplemAWES3DOF_P1_A_DW;

/* External outputs (root outports fed by signals with default storage) */
ExtY_PartialImplemAWES3DOF_P1_T PartialImplemAWES3DOF_P1_A_Y;

/* Real-time model */
RT_MODEL_PartialImplemAWES3DO_T PartialImplemAWES3DOF_P1_A_M_;
RT_MODEL_PartialImplemAWES3DO_T *const PartialImplemAWES3DOF_P1_A_M =
  &PartialImplemAWES3DOF_P1_A_M_;
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
  PartialImplemAWES3DOF_P1_A_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  PartialImplemAWES3DOF_P1_A_step();
  PartialImplemAWES3DOF_P1_A_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  PartialImplemAWES3DOF_P1_A_step();
  PartialImplemAWES3DOF_P1_A_derivatives();

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
void PartialImplemAWES3DOF_P1_A_step(void)
{
  uint32_T rtb_Prelookup_o1;
  real_T rtb_Product_m;
  real_T rtb_UnaryMinus;
  int32_T bpIndex[2];
  real_T rtb_Product1_c_tmp;
  int32_T iU;
  real_T rtb_Transpose_tmp;
  real_T rtb_sincos_o2_tmp;
  if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_A_M)) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
                          ((PartialImplemAWES3DOF_P1_A_M->Timing.clockTick0+1)*
      PartialImplemAWES3DOF_P1_A_M->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(PartialImplemAWES3DOF_P1_A_M)) {
    PartialImplemAWES3DOF_P1_A_M->Timing.t[0] = rtsiGetT
      (&PartialImplemAWES3DOF_P1_A_M->solverInfo);
  }

  /* Trigonometry: '<S15>/sincos' incorporates:
   *  Integrator: '<S7>/Theta'
   *  Trigonometry: '<S11>/sincos'
   */
  rtb_sincos_o2_tmp = cos(PartialImplemAWES3DOF_P1_A_X.theta);
  PartialImplemAWES3DOF_P1_A_B.u2rhoV2 = sin(PartialImplemAWES3DOF_P1_A_X.theta);

  /* Fcn: '<S15>/Fcn11' incorporates:
   *  Trigonometry: '<S15>/sincos'
   */
  PartialImplemAWES3DOF_P1_A_B.Transpose[0] = rtb_sincos_o2_tmp;

  /* Fcn: '<S15>/Fcn21' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[1] = 0.0 *
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2;

  /* Fcn: '<S15>/Fcn31' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[2] =
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2;

  /* Fcn: '<S15>/Fcn12' incorporates:
   *  Trigonometry: '<S15>/sincos'
   */
  PartialImplemAWES3DOF_P1_A_B.Transpose[3] = rtb_sincos_o2_tmp * 0.0;

  /* Fcn: '<S15>/Fcn22' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[4] = 0.0 *
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2 * 0.0 + 1.0;

  /* Fcn: '<S15>/Fcn32' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[5] =
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2 * 0.0;

  /* Fcn: '<S15>/Fcn13' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[6] =
    -PartialImplemAWES3DOF_P1_A_B.u2rhoV2;

  /* Fcn: '<S15>/Fcn23' incorporates:
   *  Trigonometry: '<S15>/sincos'
   */
  PartialImplemAWES3DOF_P1_A_B.Transpose[7] = 0.0 * rtb_sincos_o2_tmp;

  /* Fcn: '<S15>/Fcn33' incorporates:
   *  Trigonometry: '<S15>/sincos'
   */
  PartialImplemAWES3DOF_P1_A_B.Transpose[8] = rtb_sincos_o2_tmp;

  /* Math: '<S12>/Transpose' */
  for (iU = 0; iU < 3; iU++) {
    PartialImplemAWES3DOF_P1_A_B.rtb_Transpose_m[3 * iU] =
      PartialImplemAWES3DOF_P1_A_B.Transpose[iU];
    PartialImplemAWES3DOF_P1_A_B.rtb_Transpose_m[3 * iU + 1] =
      PartialImplemAWES3DOF_P1_A_B.Transpose[iU + 3];
    PartialImplemAWES3DOF_P1_A_B.rtb_Transpose_m[3 * iU + 2] =
      PartialImplemAWES3DOF_P1_A_B.Transpose[iU + 6];
  }

  memcpy(&PartialImplemAWES3DOF_P1_A_B.Transpose[0],
         &PartialImplemAWES3DOF_P1_A_B.rtb_Transpose_m[0], 9U * sizeof(real_T));

  /* End of Math: '<S12>/Transpose' */

  /* Integrator: '<S7>/U,w' incorporates:
   *  Constant: '<S7>/U0'
   *  Constant: '<S7>/w0'
   *  SignalConversion generated from: '<S7>/U,w'
   */
  if (PartialImplemAWES3DOF_P1_A_DW.Uw_IWORK != 0) {
    PartialImplemAWES3DOF_P1_A_X.u[0] = 19.997078706479069;
    PartialImplemAWES3DOF_P1_A_X.u[1] = 0.34182335625497184;
  }

  /* SignalConversion generated from: '<S8>/sincos' incorporates:
   *  Integrator: '<S7>/Theta'
   */
  PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[0] = 0.0;
  PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[1] =
    PartialImplemAWES3DOF_P1_A_X.theta;
  PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[2] = 0.0;
  for (iU = 0; iU < 3; iU++) {
    /* Trigonometry: '<S8>/sincos' */
    PartialImplemAWES3DOF_P1_A_B.Product_b[iU] = cos
      (PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[iU]);

    /* Product: '<S12>/Product' incorporates:
     *  Integrator: '<S7>/U,w'
     *  SignalConversion generated from: '<S12>/Product'
     */
    PartialImplemAWES3DOF_P1_A_B.sincos_o2[iU] =
      PartialImplemAWES3DOF_P1_A_B.Transpose[iU + 6] *
      PartialImplemAWES3DOF_P1_A_X.u[1] +
      (PartialImplemAWES3DOF_P1_A_B.Transpose[iU + 3] * 0.0 +
       PartialImplemAWES3DOF_P1_A_B.Transpose[iU] *
       PartialImplemAWES3DOF_P1_A_X.u[0]);

    /* Trigonometry: '<S8>/sincos' */
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[iU] = sin
      (PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[iU]);
  }

  /* Fcn: '<S8>/Fcn11' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[0] =
    PartialImplemAWES3DOF_P1_A_B.Product_b[1] *
    PartialImplemAWES3DOF_P1_A_B.Product_b[0];

  /* Fcn: '<S8>/Fcn21' incorporates:
   *  Fcn: '<S8>/Fcn22'
   */
  PartialImplemAWES3DOF_P1_A_B.u2rhoV2 =
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[2] *
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[1];
  PartialImplemAWES3DOF_P1_A_B.Transpose[1] =
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2 *
    PartialImplemAWES3DOF_P1_A_B.Product_b[0] -
    PartialImplemAWES3DOF_P1_A_B.Product_b[2] *
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[0];

  /* Fcn: '<S8>/Fcn31' incorporates:
   *  Fcn: '<S8>/Fcn32'
   */
  rtb_Transpose_tmp = PartialImplemAWES3DOF_P1_A_B.Product_b[2] *
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[1];
  PartialImplemAWES3DOF_P1_A_B.Transpose[2] = rtb_Transpose_tmp *
    PartialImplemAWES3DOF_P1_A_B.Product_b[0] +
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[2] *
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[0];

  /* Fcn: '<S8>/Fcn12' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[3] =
    PartialImplemAWES3DOF_P1_A_B.Product_b[1] *
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[0];

  /* Fcn: '<S8>/Fcn22' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[4] =
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2 *
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[0] +
    PartialImplemAWES3DOF_P1_A_B.Product_b[2] *
    PartialImplemAWES3DOF_P1_A_B.Product_b[0];

  /* Fcn: '<S8>/Fcn32' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[5] = rtb_Transpose_tmp *
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[0] -
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[2] *
    PartialImplemAWES3DOF_P1_A_B.Product_b[0];

  /* Fcn: '<S8>/Fcn13' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[6] =
    -PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[1];

  /* Fcn: '<S8>/Fcn23' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[7] =
    PartialImplemAWES3DOF_P1_A_B.sincos_o1_j[2] *
    PartialImplemAWES3DOF_P1_A_B.Product_b[1];

  /* Fcn: '<S8>/Fcn33' */
  PartialImplemAWES3DOF_P1_A_B.Transpose[8] =
    PartialImplemAWES3DOF_P1_A_B.Product_b[2] *
    PartialImplemAWES3DOF_P1_A_B.Product_b[1];

  /* Product: '<S1>/Product' */
  for (iU = 0; iU < 3; iU++) {
    PartialImplemAWES3DOF_P1_A_B.Product_b[iU] =
      PartialImplemAWES3DOF_P1_A_B.Transpose[iU + 6] *
      PartialImplemAWES3DOF_P1_ConstB.Gain[2] +
      (PartialImplemAWES3DOF_P1_A_B.Transpose[iU + 3] *
       PartialImplemAWES3DOF_P1_ConstB.Gain[1] +
       PartialImplemAWES3DOF_P1_A_B.Transpose[iU] *
       PartialImplemAWES3DOF_P1_ConstB.Gain[0]);
  }

  /* End of Product: '<S1>/Product' */

  /* Gain: '<S19>/1//2rhoV^2' incorporates:
   *  Constant: '<S29>/Constant3'
   *  Integrator: '<S7>/U,w'
   *  Product: '<S19>/Product2'
   *  Product: '<S24>/Product'
   *  Product: '<S24>/Product1'
   *  Product: '<S24>/Product2'
   *  Sum: '<S24>/Sum'
   *  Sum: '<S2>/Sum'
   */
  PartialImplemAWES3DOF_P1_A_B.u2rhoV2 = (((PartialImplemAWES3DOF_P1_A_X.u[0] -
    -10.0) * (PartialImplemAWES3DOF_P1_A_X.u[0] - -10.0) + 0.0 * 0.0) +
    PartialImplemAWES3DOF_P1_A_X.u[1] * PartialImplemAWES3DOF_P1_A_X.u[1]) *
    1.00649 * 0.5;

  /* Trigonometry: '<S21>/Incidence' incorporates:
   *  Integrator: '<S7>/U,w'
   */
  PartialImplemAWES3DOF_P1_A_B.Gain1 = rt_atan2d_snf
    (PartialImplemAWES3DOF_P1_A_X.u[1], PartialImplemAWES3DOF_P1_A_X.u[0]);

  /* Trigonometry: '<S22>/Trigonometric Function' incorporates:
   *  Trigonometry: '<S23>/Trigonometric Function'
   *  Trigonometry: '<S5>/Trigonometric Function'
   */
  rtb_Transpose_tmp = sin(PartialImplemAWES3DOF_P1_A_B.Gain1);
  rtb_Product1_c_tmp = cos(PartialImplemAWES3DOF_P1_A_B.Gain1);

  /* UnitConversion: '<S26>/Unit Conversion' */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  PartialImplemAWES3DOF_P1_A_B.CM = 57.295779513082323 *
    PartialImplemAWES3DOF_P1_A_B.Gain1;

  /* PreLookup: '<S22>/Prelookup' */
  rtb_Prelookup_o1 = plook_binx(PartialImplemAWES3DOF_P1_A_B.CM,
    PartialImplemAWES3DOF_P1_ConstP.pooled3, 9U,
    &PartialImplemAWES3DOF_P1_A_B.CM);

  /* Interpolation_n-D: '<S22>/CD' */
  PartialImplemAWES3DOF_P1_A_B.CD = intrp1d_l_pw(rtb_Prelookup_o1,
    PartialImplemAWES3DOF_P1_A_B.CM, PartialImplemAWES3DOF_P1_ConstP.CD_Table);

  /* Interpolation_n-D: '<S22>/CL' */
  PartialImplemAWES3DOF_P1_A_B.Gain_j = intrp1d_l_pw(rtb_Prelookup_o1,
    PartialImplemAWES3DOF_P1_A_B.CM, PartialImplemAWES3DOF_P1_ConstP.CL_Table);

  /* UnitConversion: '<S27>/Unit Conversion' */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  PartialImplemAWES3DOF_P1_A_B.Gain1 *= 57.295779513082323;

  /* Interpolation_n-D: '<S23>/CD_el' incorporates:
   *  PreLookup: '<S23>/Prelookup1'
   */
  bpIndex[0] = plook_s32dd_binx(PartialImplemAWES3DOF_P1_A_B.Gain1,
    PartialImplemAWES3DOF_P1_ConstP.pooled3, 9U,
    &PartialImplemAWES3DOF_P1_A_B.Gain1);
  PartialImplemAWES3DOF_P1_A_B.frac[0] = PartialImplemAWES3DOF_P1_A_B.Gain1;
  PartialImplemAWES3DOF_P1_A_B.frac[1] =
    PartialImplemAWES3DOF_P1_ConstB.Prelookup_o2;
  bpIndex[1] = PartialImplemAWES3DOF_P1_ConstB.Prelookup_o1;
  PartialImplemAWES3DOF_P1_A_B.Gain1 = intrp2d_s32dl_pw(bpIndex,
    PartialImplemAWES3DOF_P1_A_B.frac,
    PartialImplemAWES3DOF_P1_ConstP.CD_el_Table, 10U);

  /* Sum: '<S2>/Fx_tot' incorporates:
   *  Gain: '<S18>/gain1'
   *  Gain: '<S20>/gain1'
   *  Gain: '<S22>/coeffAdjust2'
   *  Gain: '<S23>/coeffAdjust1'
   *  Product: '<S18>/Product2'
   *  Product: '<S20>/Product2'
   *  Product: '<S22>/Product2'
   *  Product: '<S22>/Product3'
   *  Product: '<S23>/Product2'
   *  Product: '<S23>/Product3'
   *  Sum: '<S22>/Sum1'
   *  Sum: '<S23>/Sum1'
   *  Trigonometry: '<S22>/Trigonometric Function'
   */
  rtb_UnaryMinus = (-(rtb_Product1_c_tmp * PartialImplemAWES3DOF_P1_A_B.CD) +
                    rtb_Transpose_tmp * PartialImplemAWES3DOF_P1_A_B.Gain_j) *
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2 * 0.0046 + (-(rtb_Product1_c_tmp *
    PartialImplemAWES3DOF_P1_A_B.Gain1) + rtb_Transpose_tmp *
    PartialImplemAWES3DOF_P1_ConstB.coeffAdjust) *
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2 * 0.0046;

  /* Gain: '<S20>/gain2 ' incorporates:
   *  Gain: '<S23>/coeffAdjust2'
   *  Product: '<S20>/Product3'
   *  Product: '<S23>/Product'
   *  Product: '<S23>/Product1'
   *  Sum: '<S23>/Sum'
   */
  PartialImplemAWES3DOF_P1_A_B.Gain1 = (-(rtb_Transpose_tmp *
    PartialImplemAWES3DOF_P1_A_B.Gain1) + rtb_Product1_c_tmp *
    PartialImplemAWES3DOF_P1_ConstB.CL_el) *
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2 * 0.0046;

  /* Trigonometry: '<S1>/Trigonometric Function' incorporates:
   *  Gain: '<S1>/Gain1'
   *  Integrator: '<S7>/Position'
   *  Product: '<S1>/Divide'
   */
  PartialImplemAWES3DOF_P1_A_B.TrigonometricFunction = atan
    (-PartialImplemAWES3DOF_P1_A_X.Xe[1] / PartialImplemAWES3DOF_P1_A_X.Xe[0]);

  /* Trigonometry: '<S2>/Cos' incorporates:
   *  Trigonometry: '<S6>/Cos'
   */
  rtb_Product_m = cos(PartialImplemAWES3DOF_P1_A_B.TrigonometricFunction);

  /* Sum: '<S2>/Va' incorporates:
   *  Constant: '<Root>/Vreel '
   *  Constant: '<S29>/Constant3'
   *  Gain: '<S2>/Gain'
   *  Product: '<S2>/Divide'
   *  Product: '<S2>/Product'
   *  Sqrt: '<S2>/Sqrt'
   *  Trigonometry: '<S2>/Cos'
   */
  PartialImplemAWES3DOF_P1_A_B.sincos_o2_b = sqrt(2.0 *
    PartialImplemAWES3DOF_P1_A_B.u2rhoV2 / 1.00649) * rtb_Product_m - 10.0;

  /* Math: '<S2>/Square1' */
  PartialImplemAWES3DOF_P1_A_B.Square1 =
    PartialImplemAWES3DOF_P1_A_B.sincos_o2_b *
    PartialImplemAWES3DOF_P1_A_B.sincos_o2_b;

  /* Product: '<S2>/Divide1' incorporates:
   *  Constant: '<S2>/Add_Cable_Drag'
   *  Product: '<S2>/Product2'
   */
  PartialImplemAWES3DOF_P1_A_B.sincos_o2_b = PartialImplemAWES3DOF_P1_A_B.Gain_j
    / (20.0 * PartialImplemAWES3DOF_P1_A_B.CD);

  /* Gain: '<S18>/gain2' incorporates:
   *  Math: '<S2>/Square'
   *  Product: '<S2>/Product1'
   */
  PartialImplemAWES3DOF_P1_A_B.sincos_o2_b =
    PartialImplemAWES3DOF_P1_A_B.Square1 * PartialImplemAWES3DOF_P1_ConstB.Gain1
    * PartialImplemAWES3DOF_P1_A_B.Gain_j *
    (PartialImplemAWES3DOF_P1_A_B.sincos_o2_b *
     PartialImplemAWES3DOF_P1_A_B.sincos_o2_b) * 0.0046;

  /* Product: '<S6>/Product' */
  rtb_Product_m *= PartialImplemAWES3DOF_P1_A_B.sincos_o2_b;

  /* Product: '<S6>/Product1' incorporates:
   *  Trigonometry: '<S6>/Sin'
   */
  PartialImplemAWES3DOF_P1_A_B.sincos_o2_b *= sin
    (PartialImplemAWES3DOF_P1_A_B.TrigonometricFunction);

  /* Integrator: '<S7>/q' */
  PartialImplemAWES3DOF_P1_A_B.q = PartialImplemAWES3DOF_P1_A_X.q;

  /* Sum: '<S7>/Sum1' incorporates:
   *  Constant: '<S10>/Constant'
   *  Constant: '<S7>/gravity'
   *  Gain: '<S18>/gain2 '
   *  Gain: '<S22>/coeffAdjust1'
   *  Gain: '<S22>/coeffAdjust3'
   *  Gain: '<S6>/Gain'
   *  Gain: '<S7>/Matrix Gain'
   *  Integrator: '<S7>/Theta'
   *  Integrator: '<S7>/U,w'
   *  Product: '<S11>/Product'
   *  Product: '<S18>/Product3'
   *  Product: '<S22>/Product'
   *  Product: '<S22>/Product1'
   *  Product: '<S5>/Product1'
   *  Product: '<S5>/Product3'
   *  Product: '<S7>/Product'
   *  Product: '<S7>/Product1'
   *  Sum: '<Root>/Sum2'
   *  Sum: '<Root>/Sum3'
   *  Sum: '<S1>/Sum2'
   *  Sum: '<S22>/Sum'
   *  Sum: '<S2>/Fz_tot'
   *  Sum: '<S7>/Sum'
   *  Trigonometry: '<S11>/sincos'
   *  Trigonometry: '<S22>/Trigonometric Function'
   *  UnaryMinus: '<S11>/Unary Minus'
   */
  PartialImplemAWES3DOF_P1_A_B.Sum1[0] = ((((rtb_Product1_c_tmp *
    PartialImplemAWES3DOF_P1_ConstB.Thrust + -rtb_Product_m) + rtb_UnaryMinus) +
    PartialImplemAWES3DOF_P1_A_B.Product_b[0]) / 1.5 + -sin
    (PartialImplemAWES3DOF_P1_A_X.theta) * 0.0) + (0.0 *
    PartialImplemAWES3DOF_P1_A_X.u[0] + -PartialImplemAWES3DOF_P1_A_X.u[1]) *
    PartialImplemAWES3DOF_P1_A_B.q;
  PartialImplemAWES3DOF_P1_A_B.Sum1[1] = (((((-(rtb_Transpose_tmp *
    PartialImplemAWES3DOF_P1_A_B.CD) + -(rtb_Product1_c_tmp *
    PartialImplemAWES3DOF_P1_A_B.Gain_j)) * PartialImplemAWES3DOF_P1_A_B.u2rhoV2
    * 0.0046 + PartialImplemAWES3DOF_P1_A_B.Gain1) + (rtb_Transpose_tmp *
    PartialImplemAWES3DOF_P1_ConstB.Thrust +
    PartialImplemAWES3DOF_P1_A_B.sincos_o2_b)) +
    PartialImplemAWES3DOF_P1_A_B.Product_b[2]) / 1.5 + rtb_sincos_o2_tmp * 0.0)
    + (0.0 * PartialImplemAWES3DOF_P1_A_X.u[1] + PartialImplemAWES3DOF_P1_A_X.u
       [0]) * PartialImplemAWES3DOF_P1_A_B.q;

  /* SignalConversion generated from: '<S7>/Position' */
  PartialImplemAWES3DOF_P1_A_B.TmpSignalConversionAtPositionIn[0] =
    PartialImplemAWES3DOF_P1_A_B.sincos_o2[0];
  PartialImplemAWES3DOF_P1_A_B.TmpSignalConversionAtPositionIn[1] =
    PartialImplemAWES3DOF_P1_A_B.sincos_o2[2];

  /* Interpolation_n-D: '<S22>/CM' */
  PartialImplemAWES3DOF_P1_A_B.CM = intrp1d_l_pw(rtb_Prelookup_o1,
    PartialImplemAWES3DOF_P1_A_B.CM, PartialImplemAWES3DOF_P1_ConstP.CM_Table);
  if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_A_M)) {
    /* Sum: '<Root>/Sum4' incorporates:
     *  Constant: '<S6>/TetherForce_M'
     */
    PartialImplemAWES3DOF_P1_A_B.Sum4 = 0.0;
  }

  /* Product: '<S9>/Product2' incorporates:
   *  Constant: '<S10>/Constant1'
   *  Constant: '<S10>/Constant2'
   *  Constant: '<S20>/Constant1'
   *  Gain: '<S18>/gain3'
   *  Gain: '<S20>/gain3'
   *  Product: '<S18>/Product4'
   *  Product: '<S20>/Product4'
   *  Product: '<S20>/Product8'
   *  Product: '<S9>/Product3'
   *  Sum: '<Root>/Sum1'
   *  Sum: '<S20>/Sum'
   *  Sum: '<S2>/M_tot'
   *  Sum: '<S9>/Sum1'
   */
  PartialImplemAWES3DOF_P1_A_B.Product2 =
    ((((PartialImplemAWES3DOF_P1_A_B.u2rhoV2 *
        PartialImplemAWES3DOF_P1_ConstB.CM_el * 8.51E-5 +
        PartialImplemAWES3DOF_P1_A_B.Gain1 * 0.0) +
       PartialImplemAWES3DOF_P1_A_B.u2rhoV2 * PartialImplemAWES3DOF_P1_A_B.CM *
       8.51E-5) + PartialImplemAWES3DOF_P1_A_B.Sum4) - 0.0 *
     PartialImplemAWES3DOF_P1_A_B.q) / 8.0;

  /* Outport: '<Root>/Altitude' incorporates:
   *  Gain: '<S1>/Gain1'
   *  Integrator: '<S7>/Position'
   */
  PartialImplemAWES3DOF_P1_A_Y.Altitude = -PartialImplemAWES3DOF_P1_A_X.Xe[1];
  if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_A_M)) {
    /* Update for Integrator: '<S7>/U,w' */
    PartialImplemAWES3DOF_P1_A_DW.Uw_IWORK = 0;
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_A_M)) {
    rt_ertODEUpdateContinuousStates(&PartialImplemAWES3DOF_P1_A_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++PartialImplemAWES3DOF_P1_A_M->Timing.clockTick0;
    PartialImplemAWES3DOF_P1_A_M->Timing.t[0] = rtsiGetSolverStopTime
      (&PartialImplemAWES3DOF_P1_A_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      PartialImplemAWES3DOF_P1_A_M->Timing.clockTick1++;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void PartialImplemAWES3DOF_P1_A_derivatives(void)
{
  XDot_PartialImplemAWES3DOF_P1_T *_rtXdot;
  _rtXdot = ((XDot_PartialImplemAWES3DOF_P1_T *)
             PartialImplemAWES3DOF_P1_A_M->derivs);

  /* Derivatives for Integrator: '<S7>/Theta' */
  _rtXdot->theta = PartialImplemAWES3DOF_P1_A_B.q;

  /* Derivatives for Integrator: '<S7>/U,w' */
  _rtXdot->u[0] = PartialImplemAWES3DOF_P1_A_B.Sum1[0];

  /* Derivatives for Integrator: '<S7>/Position' */
  _rtXdot->Xe[0] = PartialImplemAWES3DOF_P1_A_B.TmpSignalConversionAtPositionIn
    [0];

  /* Derivatives for Integrator: '<S7>/U,w' */
  _rtXdot->u[1] = PartialImplemAWES3DOF_P1_A_B.Sum1[1];

  /* Derivatives for Integrator: '<S7>/Position' */
  _rtXdot->Xe[1] = PartialImplemAWES3DOF_P1_A_B.TmpSignalConversionAtPositionIn
    [1];

  /* Derivatives for Integrator: '<S7>/q' */
  _rtXdot->q = PartialImplemAWES3DOF_P1_A_B.Product2;
}

/* Model initialize function */
void PartialImplemAWES3DOF_P1_A_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
                          &PartialImplemAWES3DOF_P1_A_M->Timing.simTimeStep);
    rtsiSetTPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo, &rtmGetTPtr
                (PartialImplemAWES3DOF_P1_A_M));
    rtsiSetStepSizePtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
                       &PartialImplemAWES3DOF_P1_A_M->Timing.stepSize0);
    rtsiSetdXPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
                 &PartialImplemAWES3DOF_P1_A_M->derivs);
    rtsiSetContStatesPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo, (real_T **)
                         &PartialImplemAWES3DOF_P1_A_M->contStates);
    rtsiSetNumContStatesPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
      &PartialImplemAWES3DOF_P1_A_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
      &PartialImplemAWES3DOF_P1_A_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
      &PartialImplemAWES3DOF_P1_A_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
      &PartialImplemAWES3DOF_P1_A_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
                          (&rtmGetErrorStatus(PartialImplemAWES3DOF_P1_A_M)));
    rtsiSetRTModelPtr(&PartialImplemAWES3DOF_P1_A_M->solverInfo,
                      PartialImplemAWES3DOF_P1_A_M);
  }

  rtsiSetSimTimeStep(&PartialImplemAWES3DOF_P1_A_M->solverInfo, MAJOR_TIME_STEP);
  PartialImplemAWES3DOF_P1_A_M->intgData.y = PartialImplemAWES3DOF_P1_A_M->odeY;
  PartialImplemAWES3DOF_P1_A_M->intgData.f[0] =
    PartialImplemAWES3DOF_P1_A_M->odeF[0];
  PartialImplemAWES3DOF_P1_A_M->intgData.f[1] =
    PartialImplemAWES3DOF_P1_A_M->odeF[1];
  PartialImplemAWES3DOF_P1_A_M->intgData.f[2] =
    PartialImplemAWES3DOF_P1_A_M->odeF[2];
  PartialImplemAWES3DOF_P1_A_M->contStates = ((X_PartialImplemAWES3DOF_P1_A_T *)
    &PartialImplemAWES3DOF_P1_A_X);
  PartialImplemAWES3DOF_P1_A_M->periodicContStateIndices = ((int_T*)
    PartialImplemAWES3_PeriodicIndX);
  PartialImplemAWES3DOF_P1_A_M->periodicContStateRanges = ((real_T*)
    PartialImplemAWES3_PeriodicRngX);
  rtsiSetSolverData(&PartialImplemAWES3DOF_P1_A_M->solverInfo, (void *)
                    &PartialImplemAWES3DOF_P1_A_M->intgData);
  rtsiSetSolverName(&PartialImplemAWES3DOF_P1_A_M->solverInfo,"ode3");
  rtmSetTPtr(PartialImplemAWES3DOF_P1_A_M,
             &PartialImplemAWES3DOF_P1_A_M->Timing.tArray[0]);
  PartialImplemAWES3DOF_P1_A_M->Timing.stepSize0 = 0.001;
  rtmSetFirstInitCond(PartialImplemAWES3DOF_P1_A_M, 1);

  /* InitializeConditions for Integrator: '<S7>/Theta' */
  PartialImplemAWES3DOF_P1_A_X.theta = 0.017092;

  /* InitializeConditions for Integrator: '<S7>/U,w' */
  if (rtmIsFirstInitCond(PartialImplemAWES3DOF_P1_A_M)) {
    PartialImplemAWES3DOF_P1_A_X.u[0] = 19.997078706479069;
    PartialImplemAWES3DOF_P1_A_X.u[1] = 0.34182335625497184;
  }

  PartialImplemAWES3DOF_P1_A_DW.Uw_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S7>/U,w' */

  /* InitializeConditions for Integrator: '<S7>/Position' */
  PartialImplemAWES3DOF_P1_A_X.Xe[0] = -10.0;
  PartialImplemAWES3DOF_P1_A_X.Xe[1] = -500.0;

  /* InitializeConditions for Integrator: '<S7>/q' */
  PartialImplemAWES3DOF_P1_A_X.q = 0.0;

  /* ConstCode for Outport: '<Root>/Elevator Command' incorporates:
   *  Constant: '<Root>/elevator'
   */
  PartialImplemAWES3DOF_P1_A_Y.ElevatorCommand = -20.0;

  /* InitializeConditions for root-level periodic continuous states */
  {
    int_T rootPeriodicContStateIndices[1] = { 0 };

    real_T rootPeriodicContStateRanges[2] = { -3.1415926535897931,
      3.1415926535897931 };

    (void) memcpy((void*)PartialImplemAWES3_PeriodicIndX,
                  rootPeriodicContStateIndices,
                  1*sizeof(int_T));
    (void) memcpy((void*)PartialImplemAWES3_PeriodicRngX,
                  rootPeriodicContStateRanges,
                  2*sizeof(real_T));
  }

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(PartialImplemAWES3DOF_P1_A_M)) {
    rtmSetFirstInitCond(PartialImplemAWES3DOF_P1_A_M, 0);
  }
}

/* Model terminate function */
void PartialImplemAWES3DOF_P1_A_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
