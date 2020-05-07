/*
 * Airframe_6DOF_v001.c
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "Airframe_6DOF_v001".
 *
 * Model version              : 1.0
 * Simulink Coder version : 9.3 (R2020a) 18-Nov-2019
 * C source code generated on : Wed Apr 29 23:00:24 2020
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#include "Airframe_6DOF_v001.h"
#include "Airframe_6DOF_v001_private.h"
#include "common/time.h"
#include "flight/imu.h"
/* Block signals (default storage) */
B_Airframe_6DOF_v001_T Airframe_6DOF_v001_B;

/* Continuous states */
X_Airframe_6DOF_v001_T Airframe_6DOF_v001_X;

/* Periodic continuous states */
PeriodicIndX_Airframe_6DOF_v0_T Airframe_6DOF_v001_PeriodicIndX;
PeriodicRngX_Airframe_6DOF_v0_T Airframe_6DOF_v001_PeriodicRngX;

/* Block states (default storage) */
DW_Airframe_6DOF_v001_T Airframe_6DOF_v001_DW;

/* External inputs (root inport signals with default storage) */
ExtU_Airframe_6DOF_v001_T Airframe_6DOF_v001_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_Airframe_6DOF_v001_T Airframe_6DOF_v001_Y;

/* Real-time model */
RT_MODEL_Airframe_6DOF_v001_T Airframe_6DOF_v001_M_;
RT_MODEL_Airframe_6DOF_v001_T *const Airframe_6DOF_v001_M =
  &Airframe_6DOF_v001_M_;
static void rate_scheduler(void);

/*     Initialize pressure and temperature tables. */
void InitCalcAtmosCOESA(real_T *temperature76, real_T *pressureRatio76)
{
  if (temperature76[0] != TEMPERATURE0 ) {
    int_T k;
    temperature76[0] = TEMPERATURE0;
    pressureRatio76[0] = 1.0;

    /* set up the data at the 1976 altitude breakpoints */
    for (k=0; k<(NUM1976PTS-1); k++) {
      if (tempGradient76[k] != 0.0) {
        temperature76[k+1] = temperature76[k] +
          tempGradient76[k]*(altitude76[k+1] - altitude76[k]);
        pressureRatio76[k+1] = pressureRatio76[k] *
          exp(log(temperature76[k]/temperature76[k+1]) * GMR/tempGradient76[k]);
      } else {
        temperature76[k+1] = temperature76[k];
        pressureRatio76[k+1] = pressureRatio76[k] *
          exp((-GMR)*(altitude76[k+1] - altitude76[k])/temperature76[k]);
      }
    }
  }
}

/*
 *     Using cached pressure and temperature tables, find the
 *     working interval and perform logarithmic interpolation.
 */
void CalcAtmosCOESA(const real_T *altitude, real_T *temp, real_T *pressure,
                    real_T *density, real_T *speedofsound, real_T *temperature76,
                    real_T *pressureRatio76, int_T numPoints)
{
  int_T i;
  for (i=0; i < numPoints; i++) {
    int_T bottom = 0;
    int_T top = NUM1976PTS-1;
    int_T idx;

    /* Find altitude interval using binary search
     *
     * Deal with the extreme cases first:
     *   if altitude <= altitude76[bottom] then return idx = bottom
     *   if altitude >= altitude76[top]    then return idx = top
     */
    if (altitude[i] <= altitude76[bottom]) {
      idx = bottom;
    } else if (altitude[i] >= altitude76[top]) {
      idx = NUM1976PTS-2;
    } else {
      for (;;) {
        idx = (bottom + top)/2;
        if (altitude[i] < altitude76[idx]) {
          top = idx - 1;
        } else if (altitude[i] >= altitude76[idx+1]) {
          bottom = idx + 1;
        } else {
          /* we have altitude76[idx] <= altitude[i] < altitude76[idx+1],
           * so break and just use idx
           */
          break;
        }
      }
    }

    /* Interval has been obtained, now do linear temperature
     * interpolation and log pressure interpolation.
     */
    if (tempGradient76[idx] != 0.0 ) {
      temp[i] = temperature76[idx] +
        tempGradient76[idx] * (altitude[i] - altitude76[idx]);
      pressure[i] = PRESSURE0 * pressureRatio76[idx] *
        (rt_powd_snf(temperature76[idx]/temp[i], GMR/tempGradient76[idx]));
    } else {
      temp[i] = temperature76[idx];
      pressure[i] = PRESSURE0 * pressureRatio76[idx] *
        exp((-GMR)*(altitude[i] - altitude76[idx]) / temperature76[idx]);
    }

    density[i] = pressure[i] / ((R_HAT/MOL_WT)*temp[i]);
    speedofsound[i] = sqrt(GAMMA*temp[i]*(R_HAT/MOL_WT));
  }
}

uint32_T plook_bincpa(real_T u, const real_T bp[], uint32_T maxIndex, real_T
                      *fraction, uint32_T *prevIndex)
{
  uint32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'on'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0U;
    *fraction = 0.0;
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_u32d_prevIdx(u, bp, *prevIndex, maxIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex;
    *fraction = 0.0;
  }

  *prevIndex = bpIndex;
  return bpIndex;
}

real_T intrp2d_la(const uint32_T bpIndex[], const real_T frac[], const real_T
                  table[], const uint32_T stride, const uint32_T maxIndex[])
{
  real_T y;
  real_T yR_1d;
  uint32_T offset_1d;

  /* Column-major Interpolation 2-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'on'
     Overflow mode: 'wrapping'
   */
  offset_1d = bpIndex[1U] * stride + bpIndex[0U];
  if (bpIndex[0U] == maxIndex[0U]) {
    y = table[offset_1d];
  } else {
    y = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] + table[offset_1d];
  }

  if (bpIndex[1U] == maxIndex[1U]) {
  } else {
    offset_1d += stride;
    if (bpIndex[0U] == maxIndex[0U]) {
      yR_1d = table[offset_1d];
    } else {
      yR_1d = (table[offset_1d + 1U] - table[offset_1d]) * frac[0U] +
        table[offset_1d];
    }

    y += (yR_1d - y) * frac[1U];
  }

  return y;
}

uint32_T binsearch_u32d_prevIdx(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T bpIndex;
  uint32_T iRght;
  uint32_T iLeft;
  uint32_T found;

  /* Binary Search using Previous Index */
  bpIndex = startIndex;
  iLeft = 0U;
  iRght = maxIndex;
  found = 0U;
  while (found == 0U) {
    if (u < bp[bpIndex]) {
      iRght = bpIndex - 1U;
      bpIndex = (iRght + iLeft) >> 1U;
    } else if (u < bp[bpIndex + 1U]) {
      found = 1U;
    } else {
      iLeft = bpIndex + 1U;
      bpIndex = (iRght + iLeft) >> 1U;
    }
  }

  return bpIndex;
}

/*
 *   This function updates active task flag for each subrate.
 * The function is called at model base rate, hence the
 * generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (Airframe_6DOF_v001_M->Timing.TaskCounters.TID[2])++;
  if ((Airframe_6DOF_v001_M->Timing.TaskCounters.TID[2]) > 99) {/* Sample time: [0.1s, 0.0s] */
    Airframe_6DOF_v001_M->Timing.TaskCounters.TID[2] = 0;
  }
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
  int_T nXc = 28;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  Airframe_6DOF_v001_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  Airframe_6DOF_v001_step(0);
  Airframe_6DOF_v001_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  Airframe_6DOF_v001_step(0);
  Airframe_6DOF_v001_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  local_stateReduction(x, rtsiGetPeriodicContStateIndices(si), 3,
                       rtsiGetPeriodicContStateRanges(si));
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  uint32_T lo;
  uint32_T hi;

  /* Uniform random number generator (random number between 0 and 1)

     #define IA      16807                      magic multiplier = 7^5
     #define IM      2147483647                 modulus = 2^31-1
     #define IQ      127773                     IM div IA
     #define IR      2836                       IM modulo IA
     #define S       4.656612875245797e-10      reciprocal of 2^31-1
     test = IA * (seed % IQ) - IR * (seed/IQ)
     seed = test < 0 ? (test + IM) : test
     return (seed*S)
   */
  lo = *u % 127773U * 16807U;
  hi = *u / 127773U * 2836U;
  if (lo < hi) {
    *u = 2147483647U - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return (real_T)*u * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  real_T y;
  real_T sr;
  real_T si;

  /* Normal (Gaussian) random number generator */
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = sqrt(-2.0 * log(si) / si) * sr;
  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
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

void rt_mrdivide_U1d1x3_U2d_9vOrDY9Z(const real_T u0[3], const real_T u1[9],
  real_T y[3])
{
  real_T A[9];
  int32_T r1;
  int32_T r2;
  int32_T r3;
  real_T maxval;
  real_T a21;
  int32_T rtemp;
  memcpy(&A[0], &u1[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabs(u1[0]);
  a21 = fabs(u1[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (fabs(u1[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  A[r2] = u1[r2] / u1[r1];
  A[r3] /= A[r1];
  A[r2 + 3] -= A[r1 + 3] * A[r2];
  A[r3 + 3] -= A[r1 + 3] * A[r3];
  A[r2 + 6] -= A[r1 + 6] * A[r2];
  A[r3 + 6] -= A[r1 + 6] * A[r3];
  if (fabs(A[r3 + 3]) > fabs(A[r2 + 3])) {
    rtemp = r2 + 1;
    r2 = r3;
    r3 = rtemp - 1;
  }

  A[r3 + 3] /= A[r2 + 3];
  A[r3 + 6] -= A[r3 + 3] * A[r2 + 6];
  y[r1] = u0[0] / A[r1];
  y[r2] = u0[1] - A[r1 + 3] * y[r1];
  y[r3] = u0[2] - A[r1 + 6] * y[r1];
  y[r2] /= A[r2 + 3];
  y[r3] -= A[r2 + 6] * y[r2];
  y[r3] /= A[r3 + 6];
  y[r2] -= A[r3 + 3] * y[r3];
  y[r1] -= y[r3] * A[r3];
  y[r1] -= y[r2] * A[r2];
}

real_T rt_modd_snf(real_T u0, real_T u1)
{
  real_T y;
  boolean_T yEq;
  real_T q;
  y = u0;
  if (u1 == 0.0) {
    if (u0 == 0.0) {
      y = u1;
    }
  } else if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (u0 == 0.0) {
    y = 0.0 / u1;
  } else if (rtIsInf(u1)) {
    if ((u1 < 0.0) != (u0 < 0.0)) {
      y = u1;
    }
  } else {
    y = fmod(u0, u1);
    yEq = (y == 0.0);
    if ((!yEq) && (u1 > floor(u1))) {
      q = fabs(u0 / u1);
      yEq = !(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q);
    }

    if (yEq) {
      y = u1 * 0.0;
    } else {
      if ((u0 < 0.0) != (u1 < 0.0)) {
        y += u1;
      }
    }
  }

  return y;
}

/* Model step function */
void Airframe_6DOF_v001_step(timeUs_t currentTimeUs)
{
  UNUSED(currentTimeUs);
  real_T qo;
  boolean_T rEQ0;
  real_T rtb_Sum_ji;
  real_T rtb_Airspeed;
  real_T rtb_sigma_ugsigma_vg;
  real_T rtb_Product_ob[3];
  real_T rtb_ixj;
  real_T rtb_Sum2_j[6];
  real_T rtb_Product_m[3];
  real_T rtb_Product3_h[3];
  real_T rtb_UnitConversion;
  real_T rtb_LowAltitudeScaleLength;
  real_T rtb_VectorConcatenate[9];
  int8_T rtAction;
  real_T rtb_Transpose[9];
  real_T rtb_Sum_br;
  real_T rtb_sincos_o1[3];
  real_T rtb_UnitConversion_p;
  real_T frac[2];
  uint32_T bpIndex[2];
  real_T rtb_MediumHighAltitudeIntensity;
  int32_T i;
  real_T tmp[3];
  real_T tmp_0[3];
  real_T rtb_Product_ig[3];
  real_T rtb_Product_ig_0[3];
  real_T rtb_Sum2_i[3];
  real_T rtb_Product_ig_1[3];
  real_T rtb_rgw_p_idx_1;
  real_T tmp_1;
  real_T rtb_LowAltitudeScaleLength_tmp;


  if (attitude.values.pitch > 200) attitude.values.pitch = 200;
  if (attitude.values.pitch < -200) attitude.values.pitch = -200;
  if (attitude.values.roll > 200) attitude.values.roll = 200;
  if (attitude.values.roll < -200) attitude.values.roll = -200;
  if (attitude.values.yaw > 200) attitude.values.yaw = 200;
  if (attitude.values.yaw < -200) attitude.values.yaw = -200;
  
  Airframe_6DOF_v001_U.Aileron = (double)(attitude.values.pitch * 0.001745329);
  Airframe_6DOF_v001_U.Elevator = (double)(attitude.values.roll * 0.001745329);
  Airframe_6DOF_v001_U.Rudder = (double)(attitude.values.yaw * 0.001745329);

  if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M)) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&Airframe_6DOF_v001_M->solverInfo,
                          ((Airframe_6DOF_v001_M->Timing.clockTick0+1)*
      Airframe_6DOF_v001_M->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(Airframe_6DOF_v001_M)) {
    Airframe_6DOF_v001_M->Timing.t[0] = rtsiGetT
      (&Airframe_6DOF_v001_M->solverInfo);
  }

  /* Integrator: '<S7>/phi theta psi' */
  rtb_sincos_o1[0] = Airframe_6DOF_v001_X.phi[0];
  rtb_sincos_o1[1] = Airframe_6DOF_v001_X.phi[1];
  rtb_sincos_o1[2] = Airframe_6DOF_v001_X.phi[2];

  /* Sum: '<S52>/Sum1' incorporates:
   *  Integrator: '<S6>/xe,ye,ze'
   *  UnaryMinus: '<S52>/Ze2height'
   */
  Airframe_6DOF_v001_B.Sum1 = -Airframe_6DOF_v001_X.Xe[2];

  /* UnitConversion: '<S63>/Unit Conversion' */
  /* Unit Conversion - from: m to: ft
     Expression: output = (3.28084*input) + (0) */
  rtb_UnitConversion = 3.280839895013123 * Airframe_6DOF_v001_B.Sum1;

  /* UnitConversion: '<S69>/Unit Conversion' incorporates:
   *  Integrator: '<S6>/ub,vb,wb'
   *  Product: '<S121>/Product'
   *  Product: '<S121>/Product1'
   *  Sqrt: '<S53>/Airspeed'
   *  Sum: '<S121>/Sum'
   */
  /* Unit Conversion - from: m/s to: ft/s
     Expression: output = (3.28084*input) + (0) */
  rtb_UnitConversion_p = sqrt(Airframe_6DOF_v001_X.Ubody[0] *
    Airframe_6DOF_v001_X.Ubody[0] + Airframe_6DOF_v001_X.Ubody[2] *
    Airframe_6DOF_v001_X.Ubody[2]) * 3.280839895013123;

  /* Saturate: '<S96>/Limit Function 10ft to 1000ft' incorporates:
   *  Saturate: '<S79>/Limit Height h<1000ft'
   */
  if (rtb_UnitConversion > 1000.0) {
    rtb_Airspeed = 1000.0;
    rtb_Sum_ji = 1000.0;
  } else {
    if (rtb_UnitConversion < 10.0) {
      rtb_Airspeed = 10.0;
    } else {
      rtb_Airspeed = rtb_UnitConversion;
    }

    if (rtb_UnitConversion < 0.0) {
      rtb_Sum_ji = 0.0;
    } else {
      rtb_Sum_ji = rtb_UnitConversion;
    }
  }

  /* End of Saturate: '<S96>/Limit Function 10ft to 1000ft' */

  /* Fcn: '<S96>/Low Altitude Scale Length' */
  rtb_LowAltitudeScaleLength = rtb_Airspeed / rt_powd_snf(0.000823 *
    rtb_Airspeed + 0.177, 1.2);

  /* Product: '<S79>/sigma_ug, sigma_vg' incorporates:
   *  Fcn: '<S79>/Low Altitude Intensity'
   */
  rtb_sigma_ugsigma_vg = 1.0 / rt_powd_snf(0.000823 * rtb_Sum_ji + 0.177, 0.4) *
    1.6404199475065615;

  /* Interpolation_n-D: '<S78>/Medium//High Altitude Intensity' incorporates:
   *  PreLookup: '<S78>/PreLook-Up Index Search  (altitude)'
   */
  bpIndex[0] = plook_bincpa(rtb_UnitConversion,
    Airframe_6DOF_v001_ConstP.PreLookUpIndexSearchaltitude_Br, 11U, &rtb_Sum_ji,
    &Airframe_6DOF_v001_DW.PreLookUpIndexSearchaltitude_DW);
  frac[0] = rtb_Sum_ji;
  frac[1] = 0.0;
  bpIndex[1] = 2U;
  rtb_MediumHighAltitudeIntensity = intrp2d_la(bpIndex, frac,
    Airframe_6DOF_v001_ConstP.MediumHighAltitudeIntensity_Tab, 12U,
    Airframe_6DOF_v001_ConstP.MediumHighAltitudeIntensity_max);
  if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
      Airframe_6DOF_v001_M->Timing.TaskCounters.TID[2] == 0) {
    /* Product: '<S71>/Product' incorporates:
     *  RandomNumber: '<S71>/White Noise'
     */
    Airframe_6DOF_v001_B.Product[0] = 5.6049912163979281 *
      Airframe_6DOF_v001_DW.NextOutput[0];
    Airframe_6DOF_v001_B.Product[1] = 5.6049912163979281 *
      Airframe_6DOF_v001_DW.NextOutput[1];
    Airframe_6DOF_v001_B.Product[2] = 5.6049912163979281 *
      Airframe_6DOF_v001_DW.NextOutput[2];
    Airframe_6DOF_v001_B.Product[3] = 5.6049912163979281 *
      Airframe_6DOF_v001_DW.NextOutput[3];
  }

  /* Outputs for Enabled SubSystem: '<S62>/Hvgw(s)' incorporates:
   *  EnablePort: '<S76>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S62>/Hugw(s)' incorporates:
   *  EnablePort: '<S75>/Enable'
   */
  if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
      Airframe_6DOF_v001_M->Timing.TaskCounters.TID[1] == 0) {
    if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
        Airframe_6DOF_v001_DW.Hugws_MODE) {
      /* Disable for Outport: '<S75>/ugw' */
      Airframe_6DOF_v001_B.w1_a[0] = 0.0;
      Airframe_6DOF_v001_B.w1_a[1] = 0.0;
      Airframe_6DOF_v001_DW.Hugws_MODE = false;
    }

    if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
        Airframe_6DOF_v001_DW.Hvgws_MODE) {
      /* Disable for Outport: '<S76>/vgw' */
      Airframe_6DOF_v001_B.w1[0] = 0.0;
      Airframe_6DOF_v001_B.w1[1] = 0.0;
      Airframe_6DOF_v001_DW.Hvgws_MODE = false;
    }
  }

  /* End of Outputs for SubSystem: '<S62>/Hvgw(s)' */
  if (Airframe_6DOF_v001_DW.Hugws_MODE) {
    /* Product: '<S75>/Lug//V' */
    rtb_Sum_ji = rtb_LowAltitudeScaleLength / rtb_UnitConversion_p;
    rtb_rgw_p_idx_1 = Airframe_6DOF_v001_ConstB.UnitConversion_a /
      rtb_UnitConversion_p;

    /* Product: '<S75>/w' incorporates:
     *  Gain: '<S75>/(2//pi)'
     *  Integrator: '<S75>/ug_p'
     *  Product: '<S75>/Lug//V1'
     *  Sqrt: '<S75>/sqrt'
     *  Sum: '<S75>/Sum'
     */
    Airframe_6DOF_v001_B.w_a[0] = (sqrt(0.63661977236758138 * rtb_Sum_ji) *
      Airframe_6DOF_v001_B.Product[0] - Airframe_6DOF_v001_X.ug_p_CSTATE[0]) /
      rtb_Sum_ji;
    Airframe_6DOF_v001_B.w_a[1] = (sqrt(0.63661977236758138 * rtb_rgw_p_idx_1) *
      Airframe_6DOF_v001_B.Product[0] - Airframe_6DOF_v001_X.ug_p_CSTATE[1]) /
      rtb_rgw_p_idx_1;

    /* Product: '<S75>/w1' incorporates:
     *  Integrator: '<S75>/ug_p'
     */
    Airframe_6DOF_v001_B.w1_a[0] = Airframe_6DOF_v001_X.ug_p_CSTATE[0] *
      rtb_sigma_ugsigma_vg;
    Airframe_6DOF_v001_B.w1_a[1] = Airframe_6DOF_v001_X.ug_p_CSTATE[1] *
      rtb_MediumHighAltitudeIntensity;
  }

  /* End of Outputs for SubSystem: '<S62>/Hugw(s)' */

  /* Outputs for Enabled SubSystem: '<S62>/Hvgw(s)' incorporates:
   *  EnablePort: '<S76>/Enable'
   */
  if (Airframe_6DOF_v001_DW.Hvgws_MODE) {
    /* Product: '<S76>/Lvg//V' incorporates:
     *  Gain: '<S68>/Lv'
     */
    rtb_Sum_ji = rtb_LowAltitudeScaleLength / rtb_UnitConversion_p;
    rtb_rgw_p_idx_1 = Airframe_6DOF_v001_ConstB.UnitConversion_a /
      rtb_UnitConversion_p;

    /* Product: '<S76>/w' incorporates:
     *  Gain: '<S76>/(1//pi)'
     *  Integrator: '<S76>/vg_p1'
     *  Product: '<S76>/Lug//V1'
     *  Sqrt: '<S76>/sqrt'
     *  Sum: '<S76>/Sum'
     */
    Airframe_6DOF_v001_B.w_n[0] = (sqrt(0.31830988618379069 * rtb_Sum_ji) *
      Airframe_6DOF_v001_B.Product[1] - Airframe_6DOF_v001_X.vg_p1_CSTATE[0]) /
      rtb_Sum_ji;

    /* Product: '<S76>/w ' incorporates:
     *  Gain: '<S76>/(1//pi)'
     *  Gain: '<S76>/sqrt(3)'
     *  Integrator: '<S76>/vg_p1'
     *  Integrator: '<S76>/vgw_p2'
     *  Product: '<S76>/Lvg//V '
     *  Sum: '<S76>/Sum1'
     */
    Airframe_6DOF_v001_B.w_d[0] = (Airframe_6DOF_v001_B.w_n[0] * rtb_Sum_ji *
      1.7320508075688772 + (Airframe_6DOF_v001_X.vg_p1_CSTATE[0] -
      Airframe_6DOF_v001_X.vgw_p2_CSTATE[0])) / rtb_Sum_ji;

    /* Product: '<S76>/w' incorporates:
     *  Gain: '<S76>/(1//pi)'
     *  Integrator: '<S76>/vg_p1'
     *  Product: '<S76>/Lug//V1'
     *  Sqrt: '<S76>/sqrt'
     *  Sum: '<S76>/Sum'
     */
    Airframe_6DOF_v001_B.w_n[1] = (sqrt(0.31830988618379069 * rtb_rgw_p_idx_1) *
      Airframe_6DOF_v001_B.Product[1] - Airframe_6DOF_v001_X.vg_p1_CSTATE[1]) /
      rtb_rgw_p_idx_1;

    /* Product: '<S76>/w ' incorporates:
     *  Gain: '<S76>/(1//pi)'
     *  Gain: '<S76>/sqrt(3)'
     *  Integrator: '<S76>/vg_p1'
     *  Integrator: '<S76>/vgw_p2'
     *  Product: '<S76>/Lvg//V '
     *  Sum: '<S76>/Sum1'
     */
    Airframe_6DOF_v001_B.w_d[1] = (Airframe_6DOF_v001_B.w_n[1] * rtb_rgw_p_idx_1
      * 1.7320508075688772 + (Airframe_6DOF_v001_X.vg_p1_CSTATE[1] -
      Airframe_6DOF_v001_X.vgw_p2_CSTATE[1])) / rtb_rgw_p_idx_1;

    /* Product: '<S76>/w 1' incorporates:
     *  Integrator: '<S76>/vgw_p2'
     */
    Airframe_6DOF_v001_B.w1[0] = rtb_sigma_ugsigma_vg *
      Airframe_6DOF_v001_X.vgw_p2_CSTATE[0];
    Airframe_6DOF_v001_B.w1[1] = rtb_MediumHighAltitudeIntensity *
      Airframe_6DOF_v001_X.vgw_p2_CSTATE[1];
  }

  /* End of Outputs for SubSystem: '<S62>/Hvgw(s)' */

  /* Gain: '<S68>/Lw' */
  frac[0] = rtb_Airspeed;

  /* Outputs for Enabled SubSystem: '<S62>/Hwgw(s)' incorporates:
   *  EnablePort: '<S77>/Enable'
   */
  if ((rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
       Airframe_6DOF_v001_M->Timing.TaskCounters.TID[1] == 0) &&
      rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
      Airframe_6DOF_v001_DW.Hwgws_MODE) {
    /* Disable for Outport: '<S77>/wgw' */
    Airframe_6DOF_v001_B.LwgV1[0] = 0.0;
    Airframe_6DOF_v001_B.LwgV1[1] = 0.0;
    Airframe_6DOF_v001_DW.Hwgws_MODE = false;
  }

  if (Airframe_6DOF_v001_DW.Hwgws_MODE) {
    /* Product: '<S77>/Lwg//V 1' incorporates:
     *  Integrator: '<S77>/wg_p2'
     */
    Airframe_6DOF_v001_B.LwgV1[0] = 1.6404199475065615 *
      Airframe_6DOF_v001_X.wg_p2_CSTATE[0];
    Airframe_6DOF_v001_B.LwgV1[1] = rtb_MediumHighAltitudeIntensity *
      Airframe_6DOF_v001_X.wg_p2_CSTATE[1];

    /* Product: '<S77>/Lwg//V' incorporates:
     *  Gain: '<S68>/Lw'
     */
    rtb_sigma_ugsigma_vg = rtb_Airspeed / rtb_UnitConversion_p;

    /* Product: '<S77>/w' incorporates:
     *  Gain: '<S77>/1//pi'
     *  Integrator: '<S77>/wg_p1'
     *  Product: '<S77>/Lug//V1'
     *  Sqrt: '<S77>/sqrt1'
     *  Sum: '<S77>/Sum'
     */
    Airframe_6DOF_v001_B.w[0] = (sqrt(0.31830988618379069 * rtb_sigma_ugsigma_vg)
      * Airframe_6DOF_v001_B.Product[2] - Airframe_6DOF_v001_X.wg_p1_CSTATE[0]) /
      rtb_sigma_ugsigma_vg;

    /* Product: '<S77>/w ' incorporates:
     *  Integrator: '<S77>/wg_p1'
     *  Integrator: '<S77>/wg_p2'
     *  Product: '<S77>/Lwg//V '
     *  Sum: '<S77>/Sum1'
     */
    Airframe_6DOF_v001_B.w_c[0] = (Airframe_6DOF_v001_B.w[0] *
      1.7320508075688772 * rtb_sigma_ugsigma_vg +
      (Airframe_6DOF_v001_X.wg_p1_CSTATE[0] - Airframe_6DOF_v001_X.wg_p2_CSTATE
       [0])) / rtb_sigma_ugsigma_vg;

    /* Product: '<S77>/Lwg//V' incorporates:
     *  Gain: '<S68>/Lw'
     */
    rtb_sigma_ugsigma_vg = Airframe_6DOF_v001_ConstB.UnitConversion_a /
      rtb_UnitConversion_p;

    /* Product: '<S77>/w' incorporates:
     *  Gain: '<S77>/1//pi'
     *  Integrator: '<S77>/wg_p1'
     *  Product: '<S77>/Lug//V1'
     *  Sqrt: '<S77>/sqrt1'
     *  Sum: '<S77>/Sum'
     */
    Airframe_6DOF_v001_B.w[1] = (sqrt(0.31830988618379069 * rtb_sigma_ugsigma_vg)
      * Airframe_6DOF_v001_B.Product[2] - Airframe_6DOF_v001_X.wg_p1_CSTATE[1]) /
      rtb_sigma_ugsigma_vg;

    /* Product: '<S77>/w ' incorporates:
     *  Integrator: '<S77>/wg_p1'
     *  Integrator: '<S77>/wg_p2'
     *  Product: '<S77>/Lwg//V '
     *  Sum: '<S77>/Sum1'
     */
    Airframe_6DOF_v001_B.w_c[1] = (Airframe_6DOF_v001_B.w[1] *
      1.7320508075688772 * rtb_sigma_ugsigma_vg +
      (Airframe_6DOF_v001_X.wg_p1_CSTATE[1] - Airframe_6DOF_v001_X.wg_p2_CSTATE
       [1])) / rtb_sigma_ugsigma_vg;
  }

  /* End of Outputs for SubSystem: '<S62>/Hwgw(s)' */

  /* Trigonometry: '<S54>/sincos' incorporates:
   *  Integrator: '<S7>/phi theta psi'
   */
  rtb_Product3_h[0] = sin(Airframe_6DOF_v001_X.phi[0]);
  rtb_Product_ob[0] = cos(Airframe_6DOF_v001_X.phi[0]);
  rtb_Product3_h[1] = sin(Airframe_6DOF_v001_X.phi[1]);
  rtb_Product_ob[1] = cos(Airframe_6DOF_v001_X.phi[1]);
  rtb_Product3_h[2] = sin(Airframe_6DOF_v001_X.phi[2]);
  rtb_Product_ob[2] = cos(Airframe_6DOF_v001_X.phi[2]);

  /* Fcn: '<S54>/Fcn11' */
  rtb_VectorConcatenate[0] = rtb_Product_ob[1] * rtb_Product_ob[0];

  /* Fcn: '<S54>/Fcn21' incorporates:
   *  Fcn: '<S54>/Fcn22'
   */
  rtb_sigma_ugsigma_vg = rtb_Product3_h[2] * rtb_Product3_h[1];
  rtb_VectorConcatenate[1] = rtb_sigma_ugsigma_vg * rtb_Product_ob[0] -
    rtb_Product_ob[2] * rtb_Product3_h[0];

  /* Fcn: '<S54>/Fcn31' incorporates:
   *  Fcn: '<S54>/Fcn32'
   */
  rtb_Airspeed = rtb_Product_ob[2] * rtb_Product3_h[1];
  rtb_VectorConcatenate[2] = rtb_Airspeed * rtb_Product_ob[0] + rtb_Product3_h[2]
    * rtb_Product3_h[0];

  /* Fcn: '<S54>/Fcn12' */
  rtb_VectorConcatenate[3] = rtb_Product_ob[1] * rtb_Product3_h[0];

  /* Fcn: '<S54>/Fcn22' */
  rtb_VectorConcatenate[4] = rtb_sigma_ugsigma_vg * rtb_Product3_h[0] +
    rtb_Product_ob[2] * rtb_Product_ob[0];

  /* Fcn: '<S54>/Fcn32' */
  rtb_VectorConcatenate[5] = rtb_Airspeed * rtb_Product3_h[0] - rtb_Product3_h[2]
    * rtb_Product_ob[0];

  /* Fcn: '<S54>/Fcn13' */
  rtb_VectorConcatenate[6] = -rtb_Product3_h[1];

  /* Fcn: '<S54>/Fcn23' */
  rtb_VectorConcatenate[7] = rtb_Product3_h[2] * rtb_Product_ob[1];

  /* Fcn: '<S54>/Fcn33' */
  rtb_VectorConcatenate[8] = rtb_Product_ob[2] * rtb_Product_ob[1];

  /* If: '<S67>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M)) {
    if (rtb_UnitConversion <= 1000.0) {
      rtAction = 0;
    } else if (rtb_UnitConversion >= 2000.0) {
      rtAction = 1;
    } else {
      rtAction = 2;
    }

    Airframe_6DOF_v001_DW.ifHeightMaxlowaltitudeelseifHei = rtAction;
  } else {
    rtAction = Airframe_6DOF_v001_DW.ifHeightMaxlowaltitudeelseifHei;
  }

  switch (rtAction) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S67>/Low altitude  velocities' incorporates:
     *  ActionPort: '<S89>/Action Port'
     */
    /* Sum: '<S95>/Sum' incorporates:
     *  Product: '<S95>/Product1'
     *  Product: '<S95>/Product2'
     */
    rtb_Product3_h[0] = Airframe_6DOF_v001_B.w1_a[0] - 0.0 *
      Airframe_6DOF_v001_B.w1[0];

    /* Sum: '<S95>/Sum1' incorporates:
     *  Product: '<S95>/Product1'
     *  Product: '<S95>/Product2'
     */
    rtb_Product3_h[1] = 0.0 * Airframe_6DOF_v001_B.w1_a[0] +
      Airframe_6DOF_v001_B.w1[0];

    /* Reshape: '<S94>/Reshape1' incorporates:
     *  Product: '<S94>/Product'
     *  SignalConversion generated from: '<S94>/Vector Concatenate'
     */
    for (i = 0; i < 3; i++) {
      rtb_Product_ob[i] = rtb_VectorConcatenate[i + 6] *
        Airframe_6DOF_v001_B.LwgV1[0] + (rtb_VectorConcatenate[i + 3] *
        rtb_Product3_h[1] + rtb_VectorConcatenate[i] * rtb_Product3_h[0]);
    }

    /* End of Reshape: '<S94>/Reshape1' */
    /* End of Outputs for SubSystem: '<S67>/Low altitude  velocities' */
    break;

   case 1:
    /* Outputs for IfAction SubSystem: '<S67>/Medium//High  altitude velocities' incorporates:
     *  ActionPort: '<S90>/Action Port'
     */
    /* Gain: '<S90>/Gain' */
    rtb_Product_ob[0] = Airframe_6DOF_v001_B.w1_a[1];
    rtb_Product_ob[1] = Airframe_6DOF_v001_B.w1[1];
    rtb_Product_ob[2] = Airframe_6DOF_v001_B.LwgV1[1];

    /* End of Outputs for SubSystem: '<S67>/Medium//High  altitude velocities' */
    break;

   case 2:
    /* Outputs for IfAction SubSystem: '<S67>/Interpolate  velocities' incorporates:
     *  ActionPort: '<S88>/Action Port'
     */
    /* Sum: '<S93>/Sum' incorporates:
     *  Product: '<S93>/Product1'
     *  Product: '<S93>/Product2'
     */
    rtb_Product_ob[0] = Airframe_6DOF_v001_B.w1_a[0] - 0.0 *
      Airframe_6DOF_v001_B.w1[0];

    /* Sum: '<S93>/Sum1' incorporates:
     *  Product: '<S93>/Product1'
     *  Product: '<S93>/Product2'
     */
    rtb_Product_ob[1] = 0.0 * Airframe_6DOF_v001_B.w1_a[0] +
      Airframe_6DOF_v001_B.w1[0];

    /* Product: '<S92>/Product' incorporates:
     *  SignalConversion generated from: '<S92>/Vector Concatenate'
     */
    for (i = 0; i < 3; i++) {
      rtb_Product3_h[i] = rtb_VectorConcatenate[i + 6] *
        Airframe_6DOF_v001_B.LwgV1[0] + (rtb_VectorConcatenate[i + 3] *
        rtb_Product_ob[1] + rtb_VectorConcatenate[i] * rtb_Product_ob[0]);
    }

    /* End of Product: '<S92>/Product' */

    /* Sum: '<S88>/Sum3' incorporates:
     *  Constant: '<S88>/max_height_low'
     *  Product: '<S88>/Product1'
     *  Sum: '<S88>/Sum1'
     *  Sum: '<S88>/Sum2'
     */
    rtb_Product_ob[0] = (Airframe_6DOF_v001_B.w1_a[1] - rtb_Product3_h[0]) *
      (rtb_UnitConversion - 1000.0) / 1000.0 + rtb_Product3_h[0];
    rtb_Product_ob[1] = (Airframe_6DOF_v001_B.w1[1] - rtb_Product3_h[1]) *
      (rtb_UnitConversion - 1000.0) / 1000.0 + rtb_Product3_h[1];
    rtb_Product_ob[2] = (Airframe_6DOF_v001_B.LwgV1[1] - rtb_Product3_h[2]) *
      (rtb_UnitConversion - 1000.0) / 1000.0 + rtb_Product3_h[2];

    /* End of Outputs for SubSystem: '<S67>/Interpolate  velocities' */
    break;
  }

  /* End of If: '<S67>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */

  /* UnitConversion: '<S51>/Unit Conversion' incorporates:
   *  Integrator: '<S6>/ub,vb,wb'
   *  Sum: '<S2>/Sum1'
   */
  /* Unit Conversion - from: ft/s to: m/s
     Expression: output = (0.3048*input) + (0) */
  rtb_Product_ob[0] = Airframe_6DOF_v001_X.Ubody[0] - 0.3048 * rtb_Product_ob[0];
  rtb_Product_ob[1] = Airframe_6DOF_v001_X.Ubody[1] - 0.3048 * rtb_Product_ob[1];

  /* Sum: '<S2>/Sum1' incorporates:
   *  Integrator: '<S6>/ub,vb,wb'
   *  UnitConversion: '<S51>/Unit Conversion'
   */
  rtb_LowAltitudeScaleLength = Airframe_6DOF_v001_X.Ubody[2] - 0.3048 *
    rtb_Product_ob[2];

  /* S-Function (saeroatmos): '<S50>/S-Function' */
  {
    /* S-Function Block: <S50>/S-Function */
    real_T *temp_table = (real_T *) &Airframe_6DOF_v001_DW.SFunction_temp_table
      [0];
    real_T *pres_table = (real_T *) &Airframe_6DOF_v001_DW.SFunction_pres_table
      [0];

    /* COESA */
    CalcAtmosCOESA( &Airframe_6DOF_v001_B.Sum1,
                   &Airframe_6DOF_v001_B.SFunction_o1,
                   &Airframe_6DOF_v001_B.SFunction_o3,
                   &Airframe_6DOF_v001_B.SFunction_o4,
                   &Airframe_6DOF_v001_B.SFunction_o2, temp_table, pres_table, 1);
  }

  /* Sum: '<S46>/Sum' incorporates:
   *  Product: '<S46>/Product'
   *  Product: '<S46>/Product1'
   *  Product: '<S46>/Product2'
   *  Sum: '<S49>/Sum'
   */
  rtb_Airspeed = (rtb_Product_ob[0] * rtb_Product_ob[0] + rtb_Product_ob[1] *
                  rtb_Product_ob[1]) + rtb_LowAltitudeScaleLength *
    rtb_LowAltitudeScaleLength;

  /* Gain: '<S27>/reference area' incorporates:
   *  Gain: '<S28>/1//2rhoV^2'
   *  Product: '<S28>/Product2'
   *  Sum: '<S46>/Sum'
   */
  rtb_sigma_ugsigma_vg = rtb_Airspeed * Airframe_6DOF_v001_B.SFunction_o4 * 0.5 *
    23.23;

  /* Sqrt: '<S29>/Airspeed' */
  rtb_Airspeed = sqrt(rtb_Airspeed);

  /* Product: '<S29>/Product' */
  rtb_Sum_ji = rtb_Product_ob[1] / rtb_Airspeed;

  /* Trigonometry: '<S29>/Sideslip' */
  if (rtb_Sum_ji > 1.0) {
    rtb_Sum_ji = 1.0;
  } else {
    if (rtb_Sum_ji < -1.0) {
      rtb_Sum_ji = -1.0;
    }
  }

  rtb_Sum_ji = asin(rtb_Sum_ji);

  /* End of Trigonometry: '<S29>/Sideslip' */

  /* Trigonometry: '<S29>/Incidence' */
  rtb_rgw_p_idx_1 = rt_atan2d_snf(rtb_LowAltitudeScaleLength, rtb_Product_ob[0]);

  /* Product: '<S36>/Product18' */
  rtb_LowAltitudeScaleLength = rtb_rgw_p_idx_1 * rtb_rgw_p_idx_1;

  /* Product: '<S36>/Product19' */
  rtb_Sum_br = rtb_LowAltitudeScaleLength * rtb_rgw_p_idx_1;

  /* Sum: '<S36>/Sum1' incorporates:
   *  Constant: '<S36>/Constant10'
   *  Constant: '<S36>/Constant11'
   *  Constant: '<S36>/Constant9'
   *  Product: '<S36>/Product10'
   *  Product: '<S36>/Product6'
   *  Product: '<S36>/Product7'
   */
  rtb_ixj = (0.00292 * rtb_rgw_p_idx_1 + 5.459 * rtb_LowAltitudeScaleLength) +
    -5.162 * rtb_Sum_br;

  /* Sum: '<S36>/Sum2' incorporates:
   *  Constant: '<S36>/Constant12'
   *  Constant: '<S36>/Constant14'
   *  Product: '<S36>/Product11'
   *  Product: '<S36>/Product14'
   */
  rtb_Sum_br = rtb_Sum_br * 3.442 + rtb_rgw_p_idx_1 * -5.578;

  /* Product: '<S36>/Product9' incorporates:
   *  Product: '<S32>/Product4'
   */
  rtb_LowAltitudeScaleLength_tmp = rtb_Sum_ji * rtb_Sum_ji;

  /* Outputs for Enabled SubSystem: '<S61>/Hqgw' incorporates:
   *  EnablePort: '<S73>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S61>/Hpgw' incorporates:
   *  EnablePort: '<S72>/Enable'
   */
  if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
      Airframe_6DOF_v001_M->Timing.TaskCounters.TID[1] == 0) {
    if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
        Airframe_6DOF_v001_DW.Hpgw_MODE) {
      /* Disable for Outport: '<S72>/pgw' */
      Airframe_6DOF_v001_B.sigma_w[0] = 0.0;
      Airframe_6DOF_v001_B.sigma_w[1] = 0.0;
      Airframe_6DOF_v001_DW.Hpgw_MODE = false;
    }

    if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
        Airframe_6DOF_v001_DW.Hqgw_MODE) {
      /* Disable for Outport: '<S73>/qgw' */
      Airframe_6DOF_v001_B.w_cz[0] = 0.0;
      Airframe_6DOF_v001_B.w_cz[1] = 0.0;
      Airframe_6DOF_v001_DW.Hqgw_MODE = false;
    }
  }

  /* End of Outputs for SubSystem: '<S61>/Hqgw' */
  if (Airframe_6DOF_v001_DW.Hpgw_MODE) {
    /* Fcn: '<S72>/sqrt(0.8//V)' */
    qo = 0.8 / rtb_UnitConversion_p;
    if (qo < 0.0) {
      qo = -sqrt(-qo);
    } else {
      qo = sqrt(qo);
    }

    /* Product: '<S72>/sigma_w' incorporates:
     *  Integrator: '<S72>/pgw_p'
     */
    Airframe_6DOF_v001_B.sigma_w[0] = 1.6404199475065615 *
      Airframe_6DOF_v001_X.pgw_p_CSTATE[0];
    Airframe_6DOF_v001_B.sigma_w[1] = rtb_MediumHighAltitudeIntensity *
      Airframe_6DOF_v001_X.pgw_p_CSTATE[1];

    /* Product: '<S72>/w3' */
    rtb_MediumHighAltitudeIntensity = rtb_UnitConversion_p *
      0.059847340050885565;

    /* Product: '<S72>/w' incorporates:
     *  Fcn: '<S72>/sqrt(0.8//V)'
     *  Integrator: '<S72>/pgw_p'
     *  Math: '<S72>/L^1//3'
     *  Product: '<S72>/Lug//V1'
     *  Product: '<S72>/w1'
     *  Product: '<S72>/w2'
     *  Sum: '<S72>/Sum'
     */
    Airframe_6DOF_v001_B.w_d5[0] = (qo / rt_powd_snf(frac[0],
      0.33333333333333331) * 0.62542342293925646 * Airframe_6DOF_v001_B.Product
      [3] - Airframe_6DOF_v001_X.pgw_p_CSTATE[0]) *
      rtb_MediumHighAltitudeIntensity;

    /* Math: '<S72>/L^1//3' incorporates:
     *  Gain: '<S68>/Lw'
     */
    if (Airframe_6DOF_v001_ConstB.UnitConversion_a < 0.0) {
      tmp_1 = -rt_powd_snf(-Airframe_6DOF_v001_ConstB.UnitConversion_a,
                           0.33333333333333331);
    } else {
      tmp_1 = rt_powd_snf(Airframe_6DOF_v001_ConstB.UnitConversion_a,
                          0.33333333333333331);
    }

    /* Product: '<S72>/w' incorporates:
     *  Fcn: '<S72>/sqrt(0.8//V)'
     *  Integrator: '<S72>/pgw_p'
     *  Product: '<S72>/Lug//V1'
     *  Product: '<S72>/w1'
     *  Product: '<S72>/w2'
     *  Sum: '<S72>/Sum'
     */
    Airframe_6DOF_v001_B.w_d5[1] = (qo / tmp_1 * 0.62542342293925646 *
      Airframe_6DOF_v001_B.Product[3] - Airframe_6DOF_v001_X.pgw_p_CSTATE[1]) *
      rtb_MediumHighAltitudeIntensity;
  }

  /* End of Outputs for SubSystem: '<S61>/Hpgw' */

  /* Outputs for Enabled SubSystem: '<S61>/Hqgw' incorporates:
   *  EnablePort: '<S73>/Enable'
   */
  if (Airframe_6DOF_v001_DW.Hqgw_MODE) {
    /* Gain: '<S73>/pi//4' */
    rtb_MediumHighAltitudeIntensity = 0.78539816339744828 * rtb_UnitConversion_p;

    /* Product: '<S73>/w' incorporates:
     *  Integrator: '<S73>/qgw_p'
     *  Product: '<S73>/wg//V'
     *  Sum: '<S73>/Sum'
     */
    Airframe_6DOF_v001_B.w_cz[0] = (Airframe_6DOF_v001_B.LwgV1[0] /
      rtb_UnitConversion_p - Airframe_6DOF_v001_X.qgw_p_CSTATE[0]) *
      (rtb_MediumHighAltitudeIntensity / 13.123359580052492);
    Airframe_6DOF_v001_B.w_cz[1] = (Airframe_6DOF_v001_B.LwgV1[1] /
      rtb_UnitConversion_p - Airframe_6DOF_v001_X.qgw_p_CSTATE[1]) *
      (rtb_MediumHighAltitudeIntensity / 13.123359580052492);
  }

  /* End of Outputs for SubSystem: '<S61>/Hqgw' */

  /* Outputs for Enabled SubSystem: '<S61>/Hrgw' incorporates:
   *  EnablePort: '<S74>/Enable'
   */
  if ((rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
       Airframe_6DOF_v001_M->Timing.TaskCounters.TID[1] == 0) &&
      rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
      Airframe_6DOF_v001_DW.Hrgw_MODE) {
    /* Disable for Outport: '<S74>/rgw' */
    Airframe_6DOF_v001_B.UnaryMinus[0] = 0.0;
    Airframe_6DOF_v001_B.UnaryMinus[1] = 0.0;
    Airframe_6DOF_v001_DW.Hrgw_MODE = false;
  }

  if (Airframe_6DOF_v001_DW.Hrgw_MODE) {
    /* Gain: '<S74>/pi//3' */
    rtb_MediumHighAltitudeIntensity = 1.0471975511965976 * rtb_UnitConversion_p;

    /* Product: '<S74>/w' incorporates:
     *  Integrator: '<S74>/rgw_p'
     *  Product: '<S74>/vg//V'
     *  Sum: '<S74>/Sum'
     */
    Airframe_6DOF_v001_B.w_l[0] = (Airframe_6DOF_v001_B.w1[0] /
      rtb_UnitConversion_p - Airframe_6DOF_v001_X.rgw_p_CSTATE[0]) *
      (rtb_MediumHighAltitudeIntensity / 13.123359580052492);

    /* UnaryMinus: '<S74>/Unary Minus' */
    Airframe_6DOF_v001_B.UnaryMinus[0] = -Airframe_6DOF_v001_B.w_l[0];

    /* Product: '<S74>/w' incorporates:
     *  Integrator: '<S74>/rgw_p'
     *  Product: '<S74>/vg//V'
     *  Sum: '<S74>/Sum'
     */
    Airframe_6DOF_v001_B.w_l[1] = (Airframe_6DOF_v001_B.w1[1] /
      rtb_UnitConversion_p - Airframe_6DOF_v001_X.rgw_p_CSTATE[1]) *
      (rtb_MediumHighAltitudeIntensity / 13.123359580052492);

    /* UnaryMinus: '<S74>/Unary Minus' */
    Airframe_6DOF_v001_B.UnaryMinus[1] = -Airframe_6DOF_v001_B.w_l[1];
  }

  /* End of Outputs for SubSystem: '<S61>/Hrgw' */

  /* If: '<S66>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M)) {
    if (rtb_UnitConversion <= 1000.0) {
      rtAction = 0;
    } else if (rtb_UnitConversion >= 2000.0) {
      rtAction = 1;
    } else {
      rtAction = 2;
    }

    Airframe_6DOF_v001_DW.ifHeightMaxlowaltitudeelseifH_k = rtAction;
  } else {
    rtAction = Airframe_6DOF_v001_DW.ifHeightMaxlowaltitudeelseifH_k;
  }

  switch (rtAction) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S66>/Low altitude  rates' incorporates:
     *  ActionPort: '<S81>/Action Port'
     */
    /* Sum: '<S87>/Sum' incorporates:
     *  Product: '<S87>/Product1'
     *  Product: '<S87>/Product2'
     */
    rtb_Product_ob[0] = Airframe_6DOF_v001_B.sigma_w[0] - 0.0 *
      Airframe_6DOF_v001_B.w_cz[0];

    /* Sum: '<S87>/Sum1' incorporates:
     *  Product: '<S87>/Product1'
     *  Product: '<S87>/Product2'
     */
    rtb_Product_ob[1] = 0.0 * Airframe_6DOF_v001_B.sigma_w[0] +
      Airframe_6DOF_v001_B.w_cz[0];

    /* Reshape: '<S86>/Reshape1' incorporates:
     *  Product: '<S86>/Product'
     *  SignalConversion generated from: '<S86>/Vector Concatenate'
     */
    for (i = 0; i < 3; i++) {
      rtb_Product_m[i] = rtb_VectorConcatenate[i + 6] *
        Airframe_6DOF_v001_B.UnaryMinus[0] + (rtb_VectorConcatenate[i + 3] *
        rtb_Product_ob[1] + rtb_VectorConcatenate[i] * rtb_Product_ob[0]);
    }

    /* End of Reshape: '<S86>/Reshape1' */
    /* End of Outputs for SubSystem: '<S66>/Low altitude  rates' */
    break;

   case 1:
    /* Outputs for IfAction SubSystem: '<S66>/Medium//High  altitude rates' incorporates:
     *  ActionPort: '<S82>/Action Port'
     */
    /* Gain: '<S82>/Gain' */
    rtb_Product_m[0] = Airframe_6DOF_v001_B.sigma_w[1];
    rtb_Product_m[1] = Airframe_6DOF_v001_B.w_cz[1];
    rtb_Product_m[2] = Airframe_6DOF_v001_B.UnaryMinus[1];

    /* End of Outputs for SubSystem: '<S66>/Medium//High  altitude rates' */
    break;

   case 2:
    /* Outputs for IfAction SubSystem: '<S66>/Interpolate  rates' incorporates:
     *  ActionPort: '<S80>/Action Port'
     */
    /* Sum: '<S85>/Sum' incorporates:
     *  Product: '<S85>/Product1'
     *  Product: '<S85>/Product2'
     */
    rtb_Product_m[0] = Airframe_6DOF_v001_B.sigma_w[0] - 0.0 *
      Airframe_6DOF_v001_B.w_cz[0];

    /* Sum: '<S85>/Sum1' incorporates:
     *  Product: '<S85>/Product1'
     *  Product: '<S85>/Product2'
     */
    rtb_Product_m[1] = 0.0 * Airframe_6DOF_v001_B.sigma_w[0] +
      Airframe_6DOF_v001_B.w_cz[0];

    /* Product: '<S84>/Product' incorporates:
     *  SignalConversion generated from: '<S84>/Vector Concatenate'
     */
    for (i = 0; i < 3; i++) {
      rtb_Product_ob[i] = rtb_VectorConcatenate[i + 6] *
        Airframe_6DOF_v001_B.UnaryMinus[0] + (rtb_VectorConcatenate[i + 3] *
        rtb_Product_m[1] + rtb_VectorConcatenate[i] * rtb_Product_m[0]);
    }

    /* End of Product: '<S84>/Product' */

    /* Sum: '<S80>/Sum1' incorporates:
     *  Constant: '<S80>/max_height_low'
     */
    rtb_UnitConversion -= 1000.0;

    /* Sum: '<S80>/Sum3' incorporates:
     *  Product: '<S80>/Product1'
     *  Sum: '<S80>/Sum2'
     */
    rtb_Product_m[0] = (Airframe_6DOF_v001_B.sigma_w[1] - rtb_Product_ob[0]) *
      rtb_UnitConversion / 1000.0 + rtb_Product_ob[0];
    rtb_Product_m[1] = (Airframe_6DOF_v001_B.w_cz[1] - rtb_Product_ob[1]) *
      rtb_UnitConversion / 1000.0 + rtb_Product_ob[1];
    rtb_Product_m[2] = (Airframe_6DOF_v001_B.UnaryMinus[1] - rtb_Product_ob[2]) *
      rtb_UnitConversion / 1000.0 + rtb_Product_ob[2];

    /* End of Outputs for SubSystem: '<S66>/Interpolate  rates' */
    break;
  }

  /* End of If: '<S66>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */

  /* Sum: '<S2>/Sum3' incorporates:
   *  Integrator: '<S6>/p,q,r '
   */
  rtb_Product_m[0] += Airframe_6DOF_v001_X.p[0];
  rtb_Product_m[1] += Airframe_6DOF_v001_X.p[1];
  rtb_UnitConversion = rtb_Product_m[2] + Airframe_6DOF_v001_X.p[2];

  /* Gain: '<S35>/Reference Span' incorporates:
   *  Constant: '<S37>/Constant1'
   *  Product: '<S35>/Product'
   */
  rtb_UnitConversion_p = 0.0 * rtb_Product_m[0] * 7.315;

  /* Gain: '<S35>/Reference Span ' incorporates:
   *  Constant: '<S39>/Constant1'
   *  Product: '<S35>/Product2'
   */
  rtb_MediumHighAltitudeIntensity = 0.0 * rtb_UnitConversion * 7.315;

  /* Sum: '<S2>/Sum2' incorporates:
   *  Constant: '<S25>/Constant'
   *  Constant: '<S33>/Constant5'
   *  Constant: '<S34>/Constant1'
   *  Constant: '<S38>/Constant8'
   *  Gain: '<S35>/Reference Length'
   *  Gain: '<S35>/Reference Span'
   *  Gain: '<S35>/Reference Span '
   *  Inport: '<Root>/RudderCmd'
   *  Product: '<S33>/Product1'
   *  Product: '<S33>/Product4'
   *  Product: '<S34>/Product5'
   *  Product: '<S35>/Product1'
   *  Product: '<S35>/Product3'
   *  Sum: '<S25>/Sum'
   *  Sum: '<S26>/Sum'
   *  Sum: '<S35>/Sum1'
   *  Sum: '<S36>/Sum'
   */
  rtb_Sum2_j[0] = (((rtb_UnitConversion_p + -0.6748 * rtb_Product_m[1] * 1.5875)
                    + rtb_MediumHighAltitudeIntensity) / rtb_Airspeed + (rtb_ixj
    + -0.03554)) + (0.0 * rtb_rgw_p_idx_1 * 1.106 + Airframe_6DOF_v001_U.Rudder *
                    0.03412);

  /* Gain: '<S35>/Reference Length' incorporates:
   *  Constant: '<S38>/Constant1'
   *  Product: '<S35>/Product1'
   */
  rtb_ixj = 0.0 * rtb_Product_m[1] * 1.5875;

  /* Sum: '<S2>/Sum2' incorporates:
   *  Constant: '<S25>/Constant'
   *  Constant: '<S31>/Constant4'
   *  Constant: '<S31>/Constant5'
   *  Constant: '<S31>/Constant6'
   *  Constant: '<S31>/Constant7'
   *  Constant: '<S32>/Constant5'
   *  Constant: '<S32>/Constant6'
   *  Constant: '<S32>/Constant7'
   *  Constant: '<S33>/Constant1'
   *  Constant: '<S34>/Constant4'
   *  Constant: '<S34>/Constant5'
   *  Constant: '<S34>/Constant6'
   *  Constant: '<S34>/Constant7'
   *  Constant: '<S36>/Constant13'
   *  Constant: '<S36>/Constant15'
   *  Constant: '<S36>/Constant16'
   *  Constant: '<S36>/Constant17'
   *  Constant: '<S36>/Constant18'
   *  Constant: '<S36>/Constant19'
   *  Constant: '<S36>/Constant8'
   *  Constant: '<S37>/Constant4'
   *  Constant: '<S37>/Constant5'
   *  Constant: '<S37>/Constant8'
   *  Constant: '<S38>/Constant3'
   *  Constant: '<S38>/Constant4'
   *  Constant: '<S38>/Constant5'
   *  Constant: '<S39>/Constant3'
   *  Constant: '<S39>/Constant4'
   *  Constant: '<S39>/Constant5'
   *  Constant: '<S39>/Constant6'
   *  Gain: '<S35>/Reference Length'
   *  Gain: '<S35>/Reference Span'
   *  Gain: '<S35>/Reference Span '
   *  Inport: '<Root>/AileronCmd'
   *  Inport: '<Root>/ElevatorCmd'
   *  Inport: '<Root>/RudderCmd'
   *  Product: '<S31>/Product1'
   *  Product: '<S31>/Product2'
   *  Product: '<S31>/Product3'
   *  Product: '<S31>/Product4'
   *  Product: '<S31>/Product5'
   *  Product: '<S32>/Product1'
   *  Product: '<S32>/Product2'
   *  Product: '<S32>/Product3'
   *  Product: '<S32>/Product4'
   *  Product: '<S33>/Product2'
   *  Product: '<S33>/Product5'
   *  Product: '<S34>/Product1'
   *  Product: '<S34>/Product2'
   *  Product: '<S34>/Product3'
   *  Product: '<S34>/Product4'
   *  Product: '<S34>/Product6'
   *  Product: '<S35>/Product'
   *  Product: '<S35>/Product1'
   *  Product: '<S35>/Product2'
   *  Product: '<S35>/Product3'
   *  Product: '<S36>/Product12'
   *  Product: '<S36>/Product13'
   *  Product: '<S36>/Product15'
   *  Product: '<S36>/Product16'
   *  Product: '<S36>/Product17'
   *  Product: '<S36>/Product20'
   *  Product: '<S36>/Product5'
   *  Product: '<S36>/Product8'
   *  Product: '<S36>/Product9'
   *  Sum: '<S25>/Sum'
   *  Sum: '<S26>/Sum'
   *  Sum: '<S31>/Sum1'
   *  Sum: '<S32>/Sum'
   *  Sum: '<S34>/Sum1'
   *  Sum: '<S35>/Sum1'
   *  Sum: '<S36>/Sum'
   *  Sum: '<S36>/Sum3'
   *  Sum: '<S36>/Sum4'
   */
  
  

  rtb_Sum2_j[1] = (((-0.124 * rtb_Product_m[0] * 7.315 + rtb_ixj) + 0.3666 *
                    rtb_UnitConversion * 7.315) / rtb_Airspeed + (rtb_Sum_ji *
    -0.7678 + -0.002226)) + ((rtb_rgw_p_idx_1 * Airframe_6DOF_v001_U.Rudder *
    0.5238 + 0.1158 * Airframe_6DOF_v001_U.Rudder) +
    Airframe_6DOF_v001_U.Aileron * -0.02956);
  rtb_Sum2_j[2] = (((rtb_UnitConversion_p + -2.988 * rtb_Product_m[1] * 1.5875)
                    + rtb_MediumHighAltitudeIntensity) / rtb_Airspeed +
                   (rtb_Sum_br + -0.05504)) + ((rtb_LowAltitudeScaleLength_tmp *
    Airframe_6DOF_v001_U.Elevator * -15.93 + -0.398 *
    Airframe_6DOF_v001_U.Elevator) + rtb_rgw_p_idx_1 * 0.0 * -1.261);
  rtb_Sum2_j[3] = (((-0.5045 * rtb_Product_m[0] * 7.315 + rtb_ixj) + 0.1695 *
                    rtb_UnitConversion * 7.315) / rtb_Airspeed + (rtb_Sum_ji *
    -0.0618 + 0.000591)) + ((rtb_rgw_p_idx_1 * Airframe_6DOF_v001_U.Aileron *
    -0.08269 + -0.09917 * Airframe_6DOF_v001_U.Aileron) +
    Airframe_6DOF_v001_U.Rudder * 0.006934);
  rtb_Sum2_j[4] = (((rtb_UnitConversion_p + -15.56 * rtb_Product_m[1] * 1.5875)
                    + -0.3118 * rtb_UnitConversion * 7.315) / rtb_Airspeed +
                   (((rtb_rgw_p_idx_1 * -0.6028 + rtb_LowAltitudeScaleLength *
                      -2.14) + rtb_LowAltitudeScaleLength_tmp * 0.6921) +
                    0.09448)) + (Airframe_6DOF_v001_U.Elevator * -1.921 +
    Airframe_6DOF_v001_ConstB.Product3);
  rtb_Sum2_j[5] = (((-0.1585 * rtb_Product_m[0] * 7.315 + 0.1595 *
                     rtb_Product_m[1] * 1.5875) + -0.1112 * rtb_UnitConversion *
                    7.315) / rtb_Airspeed + ((rtb_LowAltitudeScaleLength_tmp *
    rtb_Sum_ji * 0.1373 + rtb_Sum_ji * 0.006719) + -0.003117)) +
    (Airframe_6DOF_v001_U.Aileron * -0.003872 + Airframe_6DOF_v001_U.Rudder *
     -0.08265);
  for (i = 0; i < 3; i++) {
    /* Product: '<S27>/Product' */
    rtb_Product_m[i] = rtb_sigma_ugsigma_vg * rtb_Sum2_j[i];

    /* Product: '<S19>/Product' incorporates:
     *  Integrator: '<S6>/p,q,r '
     */
    rtb_Product_ob[i] = Airframe_6DOF_v001_ConstB.Selector[i + 6] *
      Airframe_6DOF_v001_X.p[2] + (Airframe_6DOF_v001_ConstB.Selector[i + 3] *
      Airframe_6DOF_v001_X.p[1] + Airframe_6DOF_v001_ConstB.Selector[i] *
      Airframe_6DOF_v001_X.p[0]);
  }

  /* Sum: '<S18>/Sum' incorporates:
   *  Integrator: '<S6>/p,q,r '
   *  Product: '<S21>/i x j'
   *  Product: '<S21>/j x k'
   *  Product: '<S21>/k x i'
   *  Product: '<S22>/i x k'
   *  Product: '<S22>/j x i'
   *  Product: '<S22>/k x j'
   */
  tmp[0] = Airframe_6DOF_v001_X.p[1] * rtb_Product_ob[2];
  tmp[1] = Airframe_6DOF_v001_X.p[2] * rtb_Product_ob[0];
  tmp[2] = Airframe_6DOF_v001_X.p[0] * rtb_Product_ob[1];
  tmp_0[0] = Airframe_6DOF_v001_X.p[2] * rtb_Product_ob[1];
  tmp_0[1] = Airframe_6DOF_v001_X.p[0] * rtb_Product_ob[2];
  tmp_0[2] = Airframe_6DOF_v001_X.p[1] * rtb_Product_ob[0];

  /* Sum: '<S40>/Sum' incorporates:
   *  Product: '<S44>/i x j'
   *  Product: '<S44>/j x k'
   *  Product: '<S44>/k x i'
   *  Product: '<S45>/i x k'
   *  Product: '<S45>/j x i'
   *  Product: '<S45>/k x j'
   */
  rtb_Product_ig[0] = rtb_Product_m[1] * Airframe_6DOF_v001_ConstB.Sum[2];
  rtb_Product_ig[1] = rtb_Product_m[2] * Airframe_6DOF_v001_ConstB.Sum[0];
  rtb_Product_ig[2] = rtb_Product_m[0] * Airframe_6DOF_v001_ConstB.Sum[1];
  rtb_Product_ig_0[0] = rtb_Product_m[2] * Airframe_6DOF_v001_ConstB.Sum[1];
  rtb_Product_ig_0[1] = rtb_Product_m[0] * Airframe_6DOF_v001_ConstB.Sum[2];
  rtb_Product_ig_0[2] = rtb_Product_m[1] * Airframe_6DOF_v001_ConstB.Sum[0];

  /* Product: '<S27>/Product1' incorporates:
   *  Constant: '<S27>/Constant'
   *  Constant: '<S27>/Constant1'
   *  Product: '<S27>/Product3'
   */
  rtb_Sum2_i[0] = 14.63 * rtb_sigma_ugsigma_vg * rtb_Sum2_j[3];
  rtb_Sum2_i[1] = 1.5875 * rtb_sigma_ugsigma_vg * rtb_Sum2_j[4];
  rtb_Sum2_i[2] = 14.63 * rtb_sigma_ugsigma_vg * rtb_Sum2_j[5];

  /* SignalConversion generated from: '<S15>/sincos' incorporates:
   *  Integrator: '<S7>/phi theta psi'
   */
  rtb_Product3_h[0] = Airframe_6DOF_v001_X.phi[2];
  rtb_Product3_h[1] = Airframe_6DOF_v001_X.phi[1];
  rtb_Product3_h[2] = Airframe_6DOF_v001_X.phi[0];
  for (i = 0; i < 3; i++) {
    /* Sum: '<S8>/Sum2' incorporates:
     *  Integrator: '<S6>/p,q,r '
     *  Product: '<S20>/Product'
     *  Sum: '<S18>/Sum'
     *  Sum: '<S27>/Sum1'
     *  Sum: '<S40>/Sum'
     */
    rtb_Product_ig_1[i] = (((rtb_Product_ig[i] - rtb_Product_ig_0[i]) +
      rtb_Sum2_i[i]) - (0.0 * Airframe_6DOF_v001_X.p[2] + (0.0 *
      Airframe_6DOF_v001_X.p[1] + 0.0 * Airframe_6DOF_v001_X.p[0]))) - (tmp[i] -
      tmp_0[i]);

    /* Sum: '<S18>/Sum' incorporates:
     *  Trigonometry: '<S15>/sincos'
     */
    rtb_Product_ob[i] = cos(rtb_Product3_h[i]);

    /* Trigonometry: '<S15>/sincos' */
    rtb_Product3_h[i] = sin(rtb_Product3_h[i]);
  }

  /* Product: '<S8>/Product2' */
  rt_mrdivide_U1d1x3_U2d_9vOrDY9Z(rtb_Product_ig_1,
    Airframe_6DOF_v001_ConstB.Selector2, Airframe_6DOF_v001_B.Product2);

  /* Fcn: '<S15>/Fcn11' */
  rtb_Transpose[0] = rtb_Product_ob[1] * rtb_Product_ob[0];

  /* Fcn: '<S15>/Fcn21' incorporates:
   *  Fcn: '<S15>/Fcn22'
   */
  rtb_sigma_ugsigma_vg = rtb_Product3_h[2] * rtb_Product3_h[1];
  rtb_Transpose[1] = rtb_sigma_ugsigma_vg * rtb_Product_ob[0] - rtb_Product_ob[2]
    * rtb_Product3_h[0];

  /* Fcn: '<S15>/Fcn31' incorporates:
   *  Fcn: '<S15>/Fcn32'
   */
  rtb_UnitConversion = rtb_Product_ob[2] * rtb_Product3_h[1];
  rtb_Transpose[2] = rtb_UnitConversion * rtb_Product_ob[0] + rtb_Product3_h[2] *
    rtb_Product3_h[0];

  /* Fcn: '<S15>/Fcn12' */
  rtb_Transpose[3] = rtb_Product_ob[1] * rtb_Product3_h[0];

  /* Fcn: '<S15>/Fcn22' */
  rtb_Transpose[4] = rtb_sigma_ugsigma_vg * rtb_Product3_h[0] + rtb_Product_ob[2]
    * rtb_Product_ob[0];

  /* Fcn: '<S15>/Fcn32' */
  rtb_Transpose[5] = rtb_UnitConversion * rtb_Product3_h[0] - rtb_Product3_h[2] *
    rtb_Product_ob[0];

  /* Fcn: '<S15>/Fcn13' */
  rtb_Transpose[6] = -rtb_Product3_h[1];

  /* Fcn: '<S15>/Fcn23' */
  rtb_Transpose[7] = rtb_Product3_h[2] * rtb_Product_ob[1];

  /* Fcn: '<S15>/Fcn33' */
  rtb_Transpose[8] = rtb_Product_ob[2] * rtb_Product_ob[1];

  /* UnitConversion: '<S102>/Unit Conversion' incorporates:
   *  Integrator: '<S6>/xe,ye,ze'
   *  Product: '<S101>/rad lat'
   *  Product: '<S101>/x*cos'
   *  Product: '<S101>/y*sin'
   *  Sum: '<S101>/Sum'
   */
  /* Unit Conversion - from: rad to: deg
     Expression: output = (57.2958*input) + (0) */
  frac[0] = (Airframe_6DOF_v001_X.Xe[0] - Airframe_6DOF_v001_X.Xe[1] * 0.0) *
    1.5784225029068334E-7 * 57.295779513082323;

  /* Switch: '<S106>/Switch' incorporates:
   *  Abs: '<S106>/Abs'
   *  Bias: '<S106>/Bias'
   *  Bias: '<S106>/Bias1'
   *  Constant: '<S106>/Constant2'
   *  Constant: '<S107>/Constant'
   *  Math: '<S106>/Math Function1'
   *  RelationalOperator: '<S107>/Compare'
   *  Sum: '<S52>/Sum'
   */
  if (fabs(frac[0]) > 180.0) {
    rtb_LowAltitudeScaleLength = rt_modd_snf(frac[0] + 180.0, 360.0) + -180.0;
  } else {
    rtb_LowAltitudeScaleLength = frac[0];
  }

  /* End of Switch: '<S106>/Switch' */

  /* Abs: '<S103>/Abs1' */
  rtb_UnitConversion = fabs(rtb_LowAltitudeScaleLength);

  /* Switch: '<S103>/Switch' incorporates:
   *  Bias: '<S103>/Bias'
   *  Bias: '<S103>/Bias1'
   *  Constant: '<S105>/Constant'
   *  Gain: '<S103>/Gain'
   *  Product: '<S103>/Divide1'
   *  RelationalOperator: '<S105>/Compare'
   *  Switch: '<S99>/Switch1'
   */
  if (rtb_UnitConversion > 90.0) {
    /* Signum: '<S103>/Sign1' */
    if (rtb_LowAltitudeScaleLength < 0.0) {
      rtb_LowAltitudeScaleLength = -1.0;
    } else if (rtb_LowAltitudeScaleLength > 0.0) {
      rtb_LowAltitudeScaleLength = 1.0;
    } else if (rtb_LowAltitudeScaleLength == 0.0) {
      rtb_LowAltitudeScaleLength = 0.0;
    } else {
      rtb_LowAltitudeScaleLength = (rtNaN);
    }

    /* End of Signum: '<S103>/Sign1' */
    rtb_LowAltitudeScaleLength *= -(rtb_UnitConversion + -90.0) + 90.0;
  }

  /* End of Switch: '<S103>/Switch' */

  /* GravityWGS84: '<S3>/WGS84 Gravity Model  ' */
  rtb_Sum_ji = rtb_LowAltitudeScaleLength * 0.017453292519943295;
  rtb_UnitConversion_p = fabs(rtb_Sum_ji);
  rtb_UnitConversion = 1.0;
  if (rtb_UnitConversion_p > 3.1415926535897931) {
    if (rtb_Sum_ji < -3.1415926535897931) {
      rtb_UnitConversion = -1.0;
    }

    if (rtIsInf(rtb_UnitConversion_p + 3.1415926535897931)) {
      rtb_LowAltitudeScaleLength = (rtNaN);
    } else {
      rtb_LowAltitudeScaleLength = fmod(rtb_UnitConversion_p +
        3.1415926535897931, 6.2831853071795862);
      rEQ0 = (rtb_LowAltitudeScaleLength == 0.0);
      if (!rEQ0) {
        rtb_sigma_ugsigma_vg = (rtb_UnitConversion_p + 3.1415926535897931) /
          6.2831853071795862;
        rEQ0 = !(fabs(rtb_sigma_ugsigma_vg - floor(rtb_sigma_ugsigma_vg + 0.5)) >
                 2.2204460492503131E-16 * rtb_sigma_ugsigma_vg);
      }

      if (rEQ0) {
        rtb_LowAltitudeScaleLength = 0.0;
      }
    }

    rtb_Sum_ji = (rtb_LowAltitudeScaleLength - 3.1415926535897931) *
      rtb_UnitConversion;
    rtb_UnitConversion_p = fabs(rtb_Sum_ji);
  }

  if (rtb_UnitConversion_p > 1.5707963267948966) {
    if (rtb_Sum_ji > 1.5707963267948966) {
      rtb_Sum_ji = 1.5707963267948966 - (rtb_UnitConversion_p -
        1.5707963267948966);
    }

    if (rtb_Sum_ji < -1.5707963267948966) {
      rtb_Sum_ji = -(1.5707963267948966 - (rtb_UnitConversion_p -
        1.5707963267948966));
    }
  }

  rtb_sigma_ugsigma_vg = sin(rtb_Sum_ji);
  rtb_UnitConversion = rtb_sigma_ugsigma_vg * rtb_sigma_ugsigma_vg;

  /* Gain: '<S1>/Gain' incorporates:
   *  GravityWGS84: '<S3>/WGS84 Gravity Model  '
   */
  rtb_sigma_ugsigma_vg = ((1.0 - (1.006802597171564 - 2.0 * rtb_UnitConversion /
    298.257223563) * 2.0 * Airframe_6DOF_v001_B.Sum1 / 6.378137E+6) + 3.0 *
    Airframe_6DOF_v001_B.Sum1 * Airframe_6DOF_v001_B.Sum1 / 4.0680631590769E+13)
    * ((0.00193185265241 * rtb_UnitConversion + 1.0) * 9.7803253359 / sqrt(1.0 -
        0.00669437999014 * rtb_UnitConversion)) * 2288.231;

  /* Sum: '<Root>/Sum2' */
  tmp[0] = rtb_Product_m[0] + Airframe_6DOF_v001_ConstB.ThrustX;
  tmp[1] = rtb_Product_m[1];
  tmp[2] = rtb_Product_m[2];

  /* Outport: '<Root>/StatesOut' incorporates:
   *  Integrator: '<S6>/p,q,r '
   *  Integrator: '<S6>/ub,vb,wb'
   *  Integrator: '<S6>/xe,ye,ze'
   *  Integrator: '<S7>/phi theta psi'
   */
  Airframe_6DOF_v001_Y.StatesOut[0] = Airframe_6DOF_v001_X.Xe[0];
  Airframe_6DOF_v001_Y.StatesOut[1] = Airframe_6DOF_v001_X.Xe[2];
  Airframe_6DOF_v001_Y.StatesOut[2] = Airframe_6DOF_v001_X.phi[1];
  Airframe_6DOF_v001_Y.StatesOut[3] = Airframe_6DOF_v001_X.Ubody[0];
  Airframe_6DOF_v001_Y.StatesOut[4] = Airframe_6DOF_v001_X.Ubody[2];
  Airframe_6DOF_v001_Y.StatesOut[5] = Airframe_6DOF_v001_X.p[1];
  Airframe_6DOF_v001_Y.StatesOut[6] = Airframe_6DOF_v001_B.Product2[1];
  Airframe_6DOF_v001_Y.StatesOut[9] = Airframe_6DOF_v001_X.Xe[1];
  Airframe_6DOF_v001_Y.StatesOut[10] = Airframe_6DOF_v001_X.phi[0];
  Airframe_6DOF_v001_Y.StatesOut[11] = Airframe_6DOF_v001_X.phi[2];
  Airframe_6DOF_v001_Y.StatesOut[12] = Airframe_6DOF_v001_X.Ubody[1];
  Airframe_6DOF_v001_Y.StatesOut[13] = Airframe_6DOF_v001_X.p[0];
  Airframe_6DOF_v001_Y.StatesOut[14] = Airframe_6DOF_v001_X.p[2];
  Airframe_6DOF_v001_Y.StatesOut[15] = Airframe_6DOF_v001_B.Product2[0];
  Airframe_6DOF_v001_Y.StatesOut[16] = Airframe_6DOF_v001_B.Product2[2];
  for (i = 0; i < 3; i++) {
    /* Trigonometry: '<S16>/sincos' */
    rtb_Product3_h[i] = cos(rtb_sincos_o1[i]);

    /* Product: '<S6>/Product' incorporates:
     *  Math: '<S6>/Transpose'
     *  Product: '<S1>/Product'
     */
    rtb_UnitConversion = rtb_Transpose[i + 3];
    rtb_UnitConversion_p = rtb_Transpose[i + 6];

    /* Sum: '<Root>/Sum2' incorporates:
     *  Constant: '<S9>/Constant'
     *  Product: '<S1>/Product'
     *  Product: '<S6>/Product'
     *  Sum: '<S1>/Sum2'
     */
    rtb_Product_m[i] = (((rtb_UnitConversion * 0.0 + rtb_Transpose[i] * 0.0) +
                         rtb_UnitConversion_p * rtb_sigma_ugsigma_vg) + tmp[i]) /
      2288.231;

    /* Math: '<S6>/Transpose' */
    rtb_VectorConcatenate[3 * i] = rtb_Transpose[i];
    rtb_VectorConcatenate[3 * i + 1] = rtb_UnitConversion;
    rtb_VectorConcatenate[3 * i + 2] = rtb_UnitConversion_p;

    /* Trigonometry: '<S16>/sincos' */
    rtb_sincos_o1[i] = sin(rtb_sincos_o1[i]);
  }

  /* Outport: '<Root>/StatesOut' */
  Airframe_6DOF_v001_Y.StatesOut[7] = rtb_Product_m[0];
  Airframe_6DOF_v001_Y.StatesOut[8] = rtb_Product_m[2];
  Airframe_6DOF_v001_Y.StatesOut[17] = rtb_Product_m[1];

  /* Fcn: '<S16>/phidot' incorporates:
   *  Fcn: '<S16>/psidot'
   *  Integrator: '<S6>/p,q,r '
   */
  rtb_sigma_ugsigma_vg = Airframe_6DOF_v001_X.p[1] * rtb_sincos_o1[0] +
    Airframe_6DOF_v001_X.p[2] * rtb_Product3_h[0];

  /* SignalConversion generated from: '<S7>/phi theta psi' incorporates:
   *  Fcn: '<S16>/phidot'
   *  Fcn: '<S16>/psidot'
   *  Fcn: '<S16>/thetadot'
   *  Integrator: '<S6>/p,q,r '
   */
  Airframe_6DOF_v001_B.TmpSignalConversionAtphithetaps[0] = rtb_sigma_ugsigma_vg
    * (rtb_sincos_o1[1] / rtb_Product3_h[1]) + Airframe_6DOF_v001_X.p[0];
  Airframe_6DOF_v001_B.TmpSignalConversionAtphithetaps[1] =
    Airframe_6DOF_v001_X.p[1] * rtb_Product3_h[0] - Airframe_6DOF_v001_X.p[2] *
    rtb_sincos_o1[0];
  Airframe_6DOF_v001_B.TmpSignalConversionAtphithetaps[2] = rtb_sigma_ugsigma_vg
    / rtb_Product3_h[1];

  /* Sum: '<S6>/Sum' incorporates:
   *  Integrator: '<S6>/p,q,r '
   *  Integrator: '<S6>/ub,vb,wb'
   *  Product: '<S23>/i x j'
   *  Product: '<S23>/j x k'
   *  Product: '<S23>/k x i'
   *  Product: '<S24>/i x k'
   *  Product: '<S24>/j x i'
   *  Product: '<S24>/k x j'
   *  Sum: '<S10>/Sum'
   */
  Airframe_6DOF_v001_B.Sum[0] = (Airframe_6DOF_v001_X.Ubody[1] *
    Airframe_6DOF_v001_X.p[2] - Airframe_6DOF_v001_X.Ubody[2] *
    Airframe_6DOF_v001_X.p[1]) + rtb_Product_m[0];
  Airframe_6DOF_v001_B.Sum[1] = (Airframe_6DOF_v001_X.Ubody[2] *
    Airframe_6DOF_v001_X.p[0] - Airframe_6DOF_v001_X.Ubody[0] *
    Airframe_6DOF_v001_X.p[2]) + rtb_Product_m[1];
  Airframe_6DOF_v001_B.Sum[2] = (Airframe_6DOF_v001_X.Ubody[0] *
    Airframe_6DOF_v001_X.p[1] - Airframe_6DOF_v001_X.Ubody[1] *
    Airframe_6DOF_v001_X.p[0]) + rtb_Product_m[2];

  /* Math: '<S6>/Transpose' */
  memcpy(&rtb_Transpose[0], &rtb_VectorConcatenate[0], 9U * sizeof(real_T));

  /* Product: '<S14>/Product' incorporates:
   *  Integrator: '<S6>/ub,vb,wb'
   */
  for (i = 0; i < 3; i++) {
    Airframe_6DOF_v001_B.Product_b[i] = 0.0;
    Airframe_6DOF_v001_B.Product_b[i] += rtb_Transpose[i] *
      Airframe_6DOF_v001_X.Ubody[0];
    Airframe_6DOF_v001_B.Product_b[i] += rtb_Transpose[i + 3] *
      Airframe_6DOF_v001_X.Ubody[1];
    Airframe_6DOF_v001_B.Product_b[i] += rtb_Transpose[i + 6] *
      Airframe_6DOF_v001_X.Ubody[2];
  }

  /* End of Product: '<S14>/Product' */
  if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M)) {
    if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M) &&
        Airframe_6DOF_v001_M->Timing.TaskCounters.TID[2] == 0) {
      /* Update for RandomNumber: '<S71>/White Noise' */
      Airframe_6DOF_v001_DW.NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf
        (&Airframe_6DOF_v001_DW.RandSeed[0]);
      Airframe_6DOF_v001_DW.NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf
        (&Airframe_6DOF_v001_DW.RandSeed[1]);
      Airframe_6DOF_v001_DW.NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf
        (&Airframe_6DOF_v001_DW.RandSeed[2]);
      Airframe_6DOF_v001_DW.NextOutput[3] = rt_nrand_Upu32_Yd_f_pw_snf
        (&Airframe_6DOF_v001_DW.RandSeed[3]);
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(Airframe_6DOF_v001_M)) {
    rt_ertODEUpdateContinuousStates(&Airframe_6DOF_v001_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++Airframe_6DOF_v001_M->Timing.clockTick0;
    Airframe_6DOF_v001_M->Timing.t[0] = rtsiGetSolverStopTime
      (&Airframe_6DOF_v001_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      Airframe_6DOF_v001_M->Timing.clockTick1++;
    }

    rate_scheduler();
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void Airframe_6DOF_v001_derivatives(void)
{
  XDot_Airframe_6DOF_v001_T *_rtXdot;
  _rtXdot = ((XDot_Airframe_6DOF_v001_T *) Airframe_6DOF_v001_M->derivs);

  /* Derivatives for Integrator: '<S6>/xe,ye,ze' */
  _rtXdot->Xe[0] = Airframe_6DOF_v001_B.Product_b[0];

  /* Derivatives for Integrator: '<S7>/phi theta psi' */
  _rtXdot->phi[0] = Airframe_6DOF_v001_B.TmpSignalConversionAtphithetaps[0];

  /* Derivatives for Integrator: '<S6>/ub,vb,wb' */
  _rtXdot->Ubody[0] = Airframe_6DOF_v001_B.Sum[0];

  /* Derivatives for Integrator: '<S6>/p,q,r ' */
  _rtXdot->p[0] = Airframe_6DOF_v001_B.Product2[0];

  /* Derivatives for Integrator: '<S6>/xe,ye,ze' */
  _rtXdot->Xe[1] = Airframe_6DOF_v001_B.Product_b[1];

  /* Derivatives for Integrator: '<S7>/phi theta psi' */
  _rtXdot->phi[1] = Airframe_6DOF_v001_B.TmpSignalConversionAtphithetaps[1];

  /* Derivatives for Integrator: '<S6>/ub,vb,wb' */
  _rtXdot->Ubody[1] = Airframe_6DOF_v001_B.Sum[1];

  /* Derivatives for Integrator: '<S6>/p,q,r ' */
  _rtXdot->p[1] = Airframe_6DOF_v001_B.Product2[1];

  /* Derivatives for Integrator: '<S6>/xe,ye,ze' */
  _rtXdot->Xe[2] = Airframe_6DOF_v001_B.Product_b[2];

  /* Derivatives for Integrator: '<S7>/phi theta psi' */
  _rtXdot->phi[2] = Airframe_6DOF_v001_B.TmpSignalConversionAtphithetaps[2];

  /* Derivatives for Integrator: '<S6>/ub,vb,wb' */
  _rtXdot->Ubody[2] = Airframe_6DOF_v001_B.Sum[2];

  /* Derivatives for Integrator: '<S6>/p,q,r ' */
  _rtXdot->p[2] = Airframe_6DOF_v001_B.Product2[2];

  /* Derivatives for Enabled SubSystem: '<S62>/Hugw(s)' */
  if (Airframe_6DOF_v001_DW.Hugws_MODE) {
    /* Derivatives for Integrator: '<S75>/ug_p' */
    _rtXdot->ug_p_CSTATE[0] = Airframe_6DOF_v001_B.w_a[0];
    _rtXdot->ug_p_CSTATE[1] = Airframe_6DOF_v001_B.w_a[1];
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_Airframe_6DOF_v001_T *) Airframe_6DOF_v001_M->derivs)
             ->ug_p_CSTATE[0]);
      for (i=0; i < 2; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S62>/Hugw(s)' */

  /* Derivatives for Enabled SubSystem: '<S62>/Hvgw(s)' */
  if (Airframe_6DOF_v001_DW.Hvgws_MODE) {
    /* Derivatives for Integrator: '<S76>/vg_p1' */
    _rtXdot->vg_p1_CSTATE[0] = Airframe_6DOF_v001_B.w_n[0];

    /* Derivatives for Integrator: '<S76>/vgw_p2' */
    _rtXdot->vgw_p2_CSTATE[0] = Airframe_6DOF_v001_B.w_d[0];

    /* Derivatives for Integrator: '<S76>/vg_p1' */
    _rtXdot->vg_p1_CSTATE[1] = Airframe_6DOF_v001_B.w_n[1];

    /* Derivatives for Integrator: '<S76>/vgw_p2' */
    _rtXdot->vgw_p2_CSTATE[1] = Airframe_6DOF_v001_B.w_d[1];
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_Airframe_6DOF_v001_T *) Airframe_6DOF_v001_M->derivs)
             ->vg_p1_CSTATE[0]);
      for (i=0; i < 4; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S62>/Hvgw(s)' */

  /* Derivatives for Enabled SubSystem: '<S62>/Hwgw(s)' */
  if (Airframe_6DOF_v001_DW.Hwgws_MODE) {
    /* Derivatives for Integrator: '<S77>/wg_p1' */
    _rtXdot->wg_p1_CSTATE[0] = Airframe_6DOF_v001_B.w[0];

    /* Derivatives for Integrator: '<S77>/wg_p2' */
    _rtXdot->wg_p2_CSTATE[0] = Airframe_6DOF_v001_B.w_c[0];

    /* Derivatives for Integrator: '<S77>/wg_p1' */
    _rtXdot->wg_p1_CSTATE[1] = Airframe_6DOF_v001_B.w[1];

    /* Derivatives for Integrator: '<S77>/wg_p2' */
    _rtXdot->wg_p2_CSTATE[1] = Airframe_6DOF_v001_B.w_c[1];
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_Airframe_6DOF_v001_T *) Airframe_6DOF_v001_M->derivs)
             ->wg_p1_CSTATE[0]);
      for (i=0; i < 4; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S62>/Hwgw(s)' */

  /* Derivatives for Enabled SubSystem: '<S61>/Hpgw' */
  if (Airframe_6DOF_v001_DW.Hpgw_MODE) {
    /* Derivatives for Integrator: '<S72>/pgw_p' */
    _rtXdot->pgw_p_CSTATE[0] = Airframe_6DOF_v001_B.w_d5[0];
    _rtXdot->pgw_p_CSTATE[1] = Airframe_6DOF_v001_B.w_d5[1];
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_Airframe_6DOF_v001_T *) Airframe_6DOF_v001_M->derivs)
             ->pgw_p_CSTATE[0]);
      for (i=0; i < 2; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S61>/Hpgw' */

  /* Derivatives for Enabled SubSystem: '<S61>/Hqgw' */
  if (Airframe_6DOF_v001_DW.Hqgw_MODE) {
    /* Derivatives for Integrator: '<S73>/qgw_p' */
    _rtXdot->qgw_p_CSTATE[0] = Airframe_6DOF_v001_B.w_cz[0];
    _rtXdot->qgw_p_CSTATE[1] = Airframe_6DOF_v001_B.w_cz[1];
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_Airframe_6DOF_v001_T *) Airframe_6DOF_v001_M->derivs)
             ->qgw_p_CSTATE[0]);
      for (i=0; i < 2; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S61>/Hqgw' */

  /* Derivatives for Enabled SubSystem: '<S61>/Hrgw' */
  if (Airframe_6DOF_v001_DW.Hrgw_MODE) {
    /* Derivatives for Integrator: '<S74>/rgw_p' */
    _rtXdot->rgw_p_CSTATE[0] = Airframe_6DOF_v001_B.w_l[0];
    _rtXdot->rgw_p_CSTATE[1] = Airframe_6DOF_v001_B.w_l[1];
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_Airframe_6DOF_v001_T *) Airframe_6DOF_v001_M->derivs)
             ->rgw_p_CSTATE[0]);
      for (i=0; i < 2; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S61>/Hrgw' */
}

/* Model initialize function */
void Airframe_6DOF_v001_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)Airframe_6DOF_v001_M, 0,
                sizeof(RT_MODEL_Airframe_6DOF_v001_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Airframe_6DOF_v001_M->solverInfo,
                          &Airframe_6DOF_v001_M->Timing.simTimeStep);
    rtsiSetTPtr(&Airframe_6DOF_v001_M->solverInfo, &rtmGetTPtr
                (Airframe_6DOF_v001_M));
    rtsiSetStepSizePtr(&Airframe_6DOF_v001_M->solverInfo,
                       &Airframe_6DOF_v001_M->Timing.stepSize0);
    rtsiSetdXPtr(&Airframe_6DOF_v001_M->solverInfo,
                 &Airframe_6DOF_v001_M->derivs);
    rtsiSetContStatesPtr(&Airframe_6DOF_v001_M->solverInfo, (real_T **)
                         &Airframe_6DOF_v001_M->contStates);
    rtsiSetNumContStatesPtr(&Airframe_6DOF_v001_M->solverInfo,
      &Airframe_6DOF_v001_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&Airframe_6DOF_v001_M->solverInfo,
      &Airframe_6DOF_v001_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&Airframe_6DOF_v001_M->solverInfo,
      &Airframe_6DOF_v001_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&Airframe_6DOF_v001_M->solverInfo,
      &Airframe_6DOF_v001_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&Airframe_6DOF_v001_M->solverInfo, (&rtmGetErrorStatus
      (Airframe_6DOF_v001_M)));
    rtsiSetRTModelPtr(&Airframe_6DOF_v001_M->solverInfo, Airframe_6DOF_v001_M);
  }

  rtsiSetSimTimeStep(&Airframe_6DOF_v001_M->solverInfo, MAJOR_TIME_STEP);
  Airframe_6DOF_v001_M->intgData.y = Airframe_6DOF_v001_M->odeY;
  Airframe_6DOF_v001_M->intgData.f[0] = Airframe_6DOF_v001_M->odeF[0];
  Airframe_6DOF_v001_M->intgData.f[1] = Airframe_6DOF_v001_M->odeF[1];
  Airframe_6DOF_v001_M->intgData.f[2] = Airframe_6DOF_v001_M->odeF[2];
  Airframe_6DOF_v001_M->contStates = ((X_Airframe_6DOF_v001_T *)
    &Airframe_6DOF_v001_X);
  Airframe_6DOF_v001_M->periodicContStateIndices = ((int_T*)
    Airframe_6DOF_v001_PeriodicIndX);
  Airframe_6DOF_v001_M->periodicContStateRanges = ((real_T*)
    Airframe_6DOF_v001_PeriodicRngX);
  rtsiSetSolverData(&Airframe_6DOF_v001_M->solverInfo, (void *)
                    &Airframe_6DOF_v001_M->intgData);
  rtsiSetSolverName(&Airframe_6DOF_v001_M->solverInfo,"ode3");
  rtmSetTPtr(Airframe_6DOF_v001_M, &Airframe_6DOF_v001_M->Timing.tArray[0]);
  Airframe_6DOF_v001_M->Timing.stepSize0 = 0.001;

  /* block I/O */
  (void) memset(((void *) &Airframe_6DOF_v001_B), 0,
                sizeof(B_Airframe_6DOF_v001_T));

  /* states (continuous) */
  {
    (void) memset((void *)&Airframe_6DOF_v001_X, 0,
                  sizeof(X_Airframe_6DOF_v001_T));
  }

  /* Periodic continuous states */
  {
    (void) memset((void*) Airframe_6DOF_v001_PeriodicIndX, 0,
                  3*sizeof(int_T));
    (void) memset((void*) Airframe_6DOF_v001_PeriodicRngX, 0,
                  6*sizeof(real_T));
  }

  /* states (dwork) */
  (void) memset((void *)&Airframe_6DOF_v001_DW, 0,
                sizeof(DW_Airframe_6DOF_v001_T));

  /* external inputs */
  (void)memset(&Airframe_6DOF_v001_U, 0, sizeof(ExtU_Airframe_6DOF_v001_T));

  /* external outputs */
  (void) memset(&Airframe_6DOF_v001_Y.StatesOut[0], 0,
                18U*sizeof(real_T));

  /* Start for If: '<S67>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  Airframe_6DOF_v001_DW.ifHeightMaxlowaltitudeelseifHei = -1;

  /* Start for S-Function (saeroatmos): '<S50>/S-Function' */
  {
    real_T *temp_table = (real_T *) &Airframe_6DOF_v001_DW.SFunction_temp_table
      [0];
    real_T *pres_table = (real_T *) &Airframe_6DOF_v001_DW.SFunction_pres_table
      [0];

    /* COESA */
    /*
     * Initialize COESA pressure and temperature tables.
     */
    InitCalcAtmosCOESA( temp_table, pres_table );
  }

  /* Start for If: '<S66>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  Airframe_6DOF_v001_DW.ifHeightMaxlowaltitudeelseifH_k = -1;

  /* InitializeConditions for Integrator: '<S6>/xe,ye,ze' */
  Airframe_6DOF_v001_X.Xe[0] = 0.0;

  /* InitializeConditions for Integrator: '<S7>/phi theta psi' */
  Airframe_6DOF_v001_X.phi[0] = 0.0;

  /* InitializeConditions for Integrator: '<S6>/ub,vb,wb' */
  Airframe_6DOF_v001_X.Ubody[0] = 44.54;

  /* InitializeConditions for Integrator: '<S6>/p,q,r ' */
  Airframe_6DOF_v001_X.p[0] = 0.0;

  /* InitializeConditions for Integrator: '<S6>/xe,ye,ze' */
  Airframe_6DOF_v001_X.Xe[1] = 0.0;

  /* InitializeConditions for Integrator: '<S7>/phi theta psi' */
  Airframe_6DOF_v001_X.phi[1] = 0.1308996938995747;

  /* InitializeConditions for Integrator: '<S6>/ub,vb,wb' */
  Airframe_6DOF_v001_X.Ubody[1] = 2.714;

  /* InitializeConditions for Integrator: '<S6>/p,q,r ' */
  Airframe_6DOF_v001_X.p[1] = 0.0;

  /* InitializeConditions for Integrator: '<S6>/xe,ye,ze' */
  Airframe_6DOF_v001_X.Xe[2] = -2202.0;

  /* InitializeConditions for Integrator: '<S7>/phi theta psi' */
  Airframe_6DOF_v001_X.phi[2] = 0.0;

  /* InitializeConditions for Integrator: '<S6>/ub,vb,wb' */
  Airframe_6DOF_v001_X.Ubody[2] = 5.836;

  /* InitializeConditions for Integrator: '<S6>/p,q,r ' */
  Airframe_6DOF_v001_X.p[2] = 0.0;

  /* InitializeConditions for RandomNumber: '<S71>/White Noise' */
  Airframe_6DOF_v001_DW.RandSeed[0] = 1529675776U;
  Airframe_6DOF_v001_DW.NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf
    (&Airframe_6DOF_v001_DW.RandSeed[0]);
  Airframe_6DOF_v001_DW.RandSeed[1] = 1529741312U;
  Airframe_6DOF_v001_DW.NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf
    (&Airframe_6DOF_v001_DW.RandSeed[1]);
  Airframe_6DOF_v001_DW.RandSeed[2] = 1529806848U;
  Airframe_6DOF_v001_DW.NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf
    (&Airframe_6DOF_v001_DW.RandSeed[2]);
  Airframe_6DOF_v001_DW.RandSeed[3] = 1529872384U;
  Airframe_6DOF_v001_DW.NextOutput[3] = rt_nrand_Upu32_Yd_f_pw_snf
    (&Airframe_6DOF_v001_DW.RandSeed[3]);

  /* SystemInitialize for Enabled SubSystem: '<S62>/Hugw(s)' */
  /* InitializeConditions for Integrator: '<S75>/ug_p' */
  Airframe_6DOF_v001_X.ug_p_CSTATE[0] = 0.0;

  /* SystemInitialize for Outport: '<S75>/ugw' */
  Airframe_6DOF_v001_B.w1_a[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S62>/Hugw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S62>/Hvgw(s)' */
  /* InitializeConditions for Integrator: '<S76>/vg_p1' */
  Airframe_6DOF_v001_X.vg_p1_CSTATE[0] = 0.0;

  /* InitializeConditions for Integrator: '<S76>/vgw_p2' */
  Airframe_6DOF_v001_X.vgw_p2_CSTATE[0] = 0.0;

  /* SystemInitialize for Outport: '<S76>/vgw' */
  Airframe_6DOF_v001_B.w1[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S62>/Hvgw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S62>/Hwgw(s)' */
  /* InitializeConditions for Integrator: '<S77>/wg_p1' */
  Airframe_6DOF_v001_X.wg_p1_CSTATE[0] = 0.0;

  /* InitializeConditions for Integrator: '<S77>/wg_p2' */
  Airframe_6DOF_v001_X.wg_p2_CSTATE[0] = 0.0;

  /* SystemInitialize for Outport: '<S77>/wgw' */
  Airframe_6DOF_v001_B.LwgV1[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S62>/Hwgw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S61>/Hpgw' */
  /* InitializeConditions for Integrator: '<S72>/pgw_p' */
  Airframe_6DOF_v001_X.pgw_p_CSTATE[0] = 0.0;

  /* SystemInitialize for Outport: '<S72>/pgw' */
  Airframe_6DOF_v001_B.sigma_w[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S61>/Hpgw' */

  /* SystemInitialize for Enabled SubSystem: '<S61>/Hqgw' */
  /* InitializeConditions for Integrator: '<S73>/qgw_p' */
  Airframe_6DOF_v001_X.qgw_p_CSTATE[0] = 0.0;

  /* SystemInitialize for Outport: '<S73>/qgw' */
  Airframe_6DOF_v001_B.w_cz[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S61>/Hqgw' */

  /* SystemInitialize for Enabled SubSystem: '<S61>/Hrgw' */
  /* InitializeConditions for Integrator: '<S74>/rgw_p' */
  Airframe_6DOF_v001_X.rgw_p_CSTATE[0] = 0.0;

  /* SystemInitialize for Outport: '<S74>/rgw' */
  Airframe_6DOF_v001_B.UnaryMinus[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S61>/Hrgw' */

  /* SystemInitialize for Enabled SubSystem: '<S62>/Hugw(s)' */
  /* InitializeConditions for Integrator: '<S75>/ug_p' */
  Airframe_6DOF_v001_X.ug_p_CSTATE[1] = 0.0;

  /* SystemInitialize for Outport: '<S75>/ugw' */
  Airframe_6DOF_v001_B.w1_a[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S62>/Hugw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S62>/Hvgw(s)' */
  /* InitializeConditions for Integrator: '<S76>/vg_p1' */
  Airframe_6DOF_v001_X.vg_p1_CSTATE[1] = 0.0;

  /* InitializeConditions for Integrator: '<S76>/vgw_p2' */
  Airframe_6DOF_v001_X.vgw_p2_CSTATE[1] = 0.0;

  /* SystemInitialize for Outport: '<S76>/vgw' */
  Airframe_6DOF_v001_B.w1[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S62>/Hvgw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S62>/Hwgw(s)' */
  /* InitializeConditions for Integrator: '<S77>/wg_p1' */
  Airframe_6DOF_v001_X.wg_p1_CSTATE[1] = 0.0;

  /* InitializeConditions for Integrator: '<S77>/wg_p2' */
  Airframe_6DOF_v001_X.wg_p2_CSTATE[1] = 0.0;

  /* SystemInitialize for Outport: '<S77>/wgw' */
  Airframe_6DOF_v001_B.LwgV1[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S62>/Hwgw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S61>/Hpgw' */
  /* InitializeConditions for Integrator: '<S72>/pgw_p' */
  Airframe_6DOF_v001_X.pgw_p_CSTATE[1] = 0.0;

  /* SystemInitialize for Outport: '<S72>/pgw' */
  Airframe_6DOF_v001_B.sigma_w[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S61>/Hpgw' */

  /* SystemInitialize for Enabled SubSystem: '<S61>/Hqgw' */
  /* InitializeConditions for Integrator: '<S73>/qgw_p' */
  Airframe_6DOF_v001_X.qgw_p_CSTATE[1] = 0.0;

  /* SystemInitialize for Outport: '<S73>/qgw' */
  Airframe_6DOF_v001_B.w_cz[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S61>/Hqgw' */

  /* SystemInitialize for Enabled SubSystem: '<S61>/Hrgw' */
  /* InitializeConditions for Integrator: '<S74>/rgw_p' */
  Airframe_6DOF_v001_X.rgw_p_CSTATE[1] = 0.0;

  /* SystemInitialize for Outport: '<S74>/rgw' */
  Airframe_6DOF_v001_B.UnaryMinus[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S61>/Hrgw' */

  /* InitializeConditions for root-level periodic continuous states */
  {
    int_T rootPeriodicContStateIndices[3] = { 3, 4, 5 };

    real_T rootPeriodicContStateRanges[6] = { -3.1415926535897931,
      3.1415926535897931, -3.1415926535897931, 3.1415926535897931,
      -3.1415926535897931, 3.1415926535897931 };

    (void) memcpy((void*)Airframe_6DOF_v001_PeriodicIndX,
                  rootPeriodicContStateIndices,
                  3*sizeof(int_T));
    (void) memcpy((void*)Airframe_6DOF_v001_PeriodicRngX,
                  rootPeriodicContStateRanges,
                  6*sizeof(real_T));
  }
}

/* Model terminate function */
void Airframe_6DOF_v001_terminate(void)
{
  /* (no terminate code required) */
}
