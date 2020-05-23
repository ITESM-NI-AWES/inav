/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: PartialImplemAWES3DOF_P1_B.c
 *
 * Code generated for Simulink model 'PartialImplemAWES3DOF_P1_B'.
 *
 * Model version                  : 1.21
 * Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
 * C/C++ source code generated on : Fri May 22 19:31:49 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "PartialImplemAWES3DOF_P1_B.h"
#include "PartialImplemAWES3DOF_P1_B_private.h"

/* Block signals (default storage) */
B_PartialImplemAWES3DOF_P1_B_T PartialImplemAWES3DOF_P1_B_B;

/* Continuous states */
X_PartialImplemAWES3DOF_P1_B_T PartialImplemAWES3DOF_P1_B_X;

/* Block states (default storage) */
DW_PartialImplemAWES3DOF_P1_B_T PartialImplemAWES3DOF_P1_B_DW;

/* External outputs (root outports fed by signals with default storage) */
ExtY_PartialImplemAWES3DOF_P1_T PartialImplemAWES3DOF_P1_B_Y;

/* Real-time model */
RT_MODEL_PartialImplemAWES3DO_T PartialImplemAWES3DOF_P1_B_M_;
RT_MODEL_PartialImplemAWES3DO_T *const PartialImplemAWES3DOF_P1_B_M =
  &PartialImplemAWES3DOF_P1_B_M_;

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
  int_T nXc = 5;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  PartialImplemAWES3DOF_P1_B_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  PartialImplemAWES3DOF_P1_B_step();
  PartialImplemAWES3DOF_P1_B_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  PartialImplemAWES3DOF_P1_B_step();
  PartialImplemAWES3DOF_P1_B_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model step function */
void PartialImplemAWES3DOF_P1_B_step(void)
{
  real_T rtb_Plant;
  real_T rtb_Subtract1;
  real_T u0;
  int32_T tmp;
  if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_B_M)) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
                          ((PartialImplemAWES3DOF_P1_B_M->Timing.clockTick0+1)*
      PartialImplemAWES3DOF_P1_B_M->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(PartialImplemAWES3DOF_P1_B_M)) {
    PartialImplemAWES3DOF_P1_B_M->Timing.t[0] = rtsiGetT
      (&PartialImplemAWES3DOF_P1_B_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_B_M)) {
    /* Step: '<Root>/Positive_Altitude_Step' */
    if (((PartialImplemAWES3DOF_P1_B_M->Timing.clockTick1) * 0.001) < 0.01) {
      tmp = 0;
    } else {
      tmp = 50;
    }

    /* End of Step: '<Root>/Positive_Altitude_Step' */

    /* Sum: '<Root>/Subtract' */
    PartialImplemAWES3DOF_P1_B_B.Subtract =
      PartialImplemAWES3DOF_P1_ConstB.Altitude_Selector2 + (real_T)tmp;
  }

  /* TransferFcn: '<Root>/Plant' */
  rtb_Plant = 0.0 * PartialImplemAWES3DOF_P1_B_X.Plant_CSTATE[0] + -400000.0 *
    PartialImplemAWES3DOF_P1_B_X.Plant_CSTATE[1];

  /* Sum: '<Root>/Sum' */
  PartialImplemAWES3DOF_P1_B_B.Sum = PartialImplemAWES3DOF_P1_B_B.Subtract +
    rtb_Plant;

  /* Sum: '<S1>/Subtract1' incorporates:
   *  ZeroPole: '<S1>/Controller'
   */
  rtb_Subtract1 = (-0.63110383296 *
                   PartialImplemAWES3DOF_P1_B_X.Controller_CSTATE + 0.033728 *
                   PartialImplemAWES3DOF_P1_B_B.Sum) -
    PartialImplemAWES3DOF_P1_ConstB.Gain1;

  /* Gain: '<S39>/Filter Coefficient' incorporates:
   *  Gain: '<S30>/Derivative Gain'
   *  Integrator: '<S31>/Filter'
   *  Sum: '<S31>/SumD'
   */
  PartialImplemAWES3DOF_P1_B_B.FilterCoefficient = (2.0 * rtb_Subtract1 -
    PartialImplemAWES3DOF_P1_B_X.Filter_CSTATE) * 100.0;

  /* Sum: '<S45>/Sum' incorporates:
   *  Gain: '<S41>/Proportional Gain'
   *  Integrator: '<S36>/Integrator'
   */
  PartialImplemAWES3DOF_P1_B_B.Sum_i = (-300.0 * rtb_Subtract1 +
    PartialImplemAWES3DOF_P1_B_X.Integrator_CSTATE) +
    PartialImplemAWES3DOF_P1_B_B.FilterCoefficient;
  if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_B_M)) {
    /* Sum: '<S3>/Diff' incorporates:
     *  UnitDelay: '<S3>/UD'
     *
     * Block description for '<S3>/Diff':
     *
     *  Add in CPU
     *
     * Block description for '<S3>/UD':
     *
     *  Store in Global RAM
     */
    u0 = PartialImplemAWES3DOF_P1_B_B.Sum_i -
      PartialImplemAWES3DOF_P1_B_DW.UD_DSTATE;

    /* Saturate: '<S1>/Saturation' */
    if (u0 > 0.35) {
      /* Outport: '<Root>/Altitude Command' */
      PartialImplemAWES3DOF_P1_B_Y.AltitudeCommand = 0.35;
    } else if (u0 < -0.35) {
      /* Outport: '<Root>/Altitude Command' */
      PartialImplemAWES3DOF_P1_B_Y.AltitudeCommand = -0.35;
    } else {
      /* Outport: '<Root>/Altitude Command' */
      PartialImplemAWES3DOF_P1_B_Y.AltitudeCommand = u0;
    }

    /* End of Saturate: '<S1>/Saturation' */
  }

  /* Gain: '<S33>/Integral Gain' */
  PartialImplemAWES3DOF_P1_B_B.IntegralGain = -15000.0 * rtb_Subtract1;

  /* Outport: '<Root>/Altitude Error' */
  PartialImplemAWES3DOF_P1_B_Y.AltitudeError = PartialImplemAWES3DOF_P1_B_B.Sum;

  /* Outport: '<Root>/Altitude' */
  PartialImplemAWES3DOF_P1_B_Y.Altitude = rtb_Plant;
  if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_B_M)) {
    if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_B_M)) {
      /* Update for UnitDelay: '<S3>/UD'
       *
       * Block description for '<S3>/UD':
       *
       *  Store in Global RAM
       */
      PartialImplemAWES3DOF_P1_B_DW.UD_DSTATE =
        PartialImplemAWES3DOF_P1_B_B.Sum_i;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(PartialImplemAWES3DOF_P1_B_M)) {
    rt_ertODEUpdateContinuousStates(&PartialImplemAWES3DOF_P1_B_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++PartialImplemAWES3DOF_P1_B_M->Timing.clockTick0;
    PartialImplemAWES3DOF_P1_B_M->Timing.t[0] = rtsiGetSolverStopTime
      (&PartialImplemAWES3DOF_P1_B_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      PartialImplemAWES3DOF_P1_B_M->Timing.clockTick1++;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void PartialImplemAWES3DOF_P1_B_derivatives(void)
{
  XDot_PartialImplemAWES3DOF_P1_T *_rtXdot;
  _rtXdot = ((XDot_PartialImplemAWES3DOF_P1_T *)
             PartialImplemAWES3DOF_P1_B_M->derivs);

  /* Derivatives for TransferFcn: '<Root>/Plant' incorporates:
   *  Outport: '<Root>/Altitude Command'
   */
  _rtXdot->Plant_CSTATE[0] = 0.0;
  _rtXdot->Plant_CSTATE[0] += -1200.0 *
    PartialImplemAWES3DOF_P1_B_X.Plant_CSTATE[0];
  _rtXdot->Plant_CSTATE[1] = 0.0;
  _rtXdot->Plant_CSTATE[0] += -2400.0 *
    PartialImplemAWES3DOF_P1_B_X.Plant_CSTATE[1];
  _rtXdot->Plant_CSTATE[1] += PartialImplemAWES3DOF_P1_B_X.Plant_CSTATE[0];
  _rtXdot->Plant_CSTATE[0] += PartialImplemAWES3DOF_P1_B_Y.AltitudeCommand;

  /* Derivatives for ZeroPole: '<S1>/Controller' */
  _rtXdot->Controller_CSTATE = 0.0;
  _rtXdot->Controller_CSTATE += -25.7546 *
    PartialImplemAWES3DOF_P1_B_X.Controller_CSTATE;
  _rtXdot->Controller_CSTATE += PartialImplemAWES3DOF_P1_B_B.Sum;

  /* Derivatives for Integrator: '<S36>/Integrator' */
  _rtXdot->Integrator_CSTATE = PartialImplemAWES3DOF_P1_B_B.IntegralGain;

  /* Derivatives for Integrator: '<S31>/Filter' */
  _rtXdot->Filter_CSTATE = PartialImplemAWES3DOF_P1_B_B.FilterCoefficient;
}

/* Model initialize function */
void PartialImplemAWES3DOF_P1_B_initialize(void)
{
  /* Registration code */
  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
                          &PartialImplemAWES3DOF_P1_B_M->Timing.simTimeStep);
    rtsiSetTPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo, &rtmGetTPtr
                (PartialImplemAWES3DOF_P1_B_M));
    rtsiSetStepSizePtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
                       &PartialImplemAWES3DOF_P1_B_M->Timing.stepSize0);
    rtsiSetdXPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
                 &PartialImplemAWES3DOF_P1_B_M->derivs);
    rtsiSetContStatesPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo, (real_T **)
                         &PartialImplemAWES3DOF_P1_B_M->contStates);
    rtsiSetNumContStatesPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
      &PartialImplemAWES3DOF_P1_B_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
      &PartialImplemAWES3DOF_P1_B_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
      &PartialImplemAWES3DOF_P1_B_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
      &PartialImplemAWES3DOF_P1_B_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
                          (&rtmGetErrorStatus(PartialImplemAWES3DOF_P1_B_M)));
    rtsiSetRTModelPtr(&PartialImplemAWES3DOF_P1_B_M->solverInfo,
                      PartialImplemAWES3DOF_P1_B_M);
  }

  rtsiSetSimTimeStep(&PartialImplemAWES3DOF_P1_B_M->solverInfo, MAJOR_TIME_STEP);
  PartialImplemAWES3DOF_P1_B_M->intgData.y = PartialImplemAWES3DOF_P1_B_M->odeY;
  PartialImplemAWES3DOF_P1_B_M->intgData.f[0] =
    PartialImplemAWES3DOF_P1_B_M->odeF[0];
  PartialImplemAWES3DOF_P1_B_M->intgData.f[1] =
    PartialImplemAWES3DOF_P1_B_M->odeF[1];
  PartialImplemAWES3DOF_P1_B_M->intgData.f[2] =
    PartialImplemAWES3DOF_P1_B_M->odeF[2];
  PartialImplemAWES3DOF_P1_B_M->contStates = ((X_PartialImplemAWES3DOF_P1_B_T *)
    &PartialImplemAWES3DOF_P1_B_X);
  rtsiSetSolverData(&PartialImplemAWES3DOF_P1_B_M->solverInfo, (void *)
                    &PartialImplemAWES3DOF_P1_B_M->intgData);
  rtsiSetSolverName(&PartialImplemAWES3DOF_P1_B_M->solverInfo,"ode3");
  rtmSetTPtr(PartialImplemAWES3DOF_P1_B_M,
             &PartialImplemAWES3DOF_P1_B_M->Timing.tArray[0]);
  PartialImplemAWES3DOF_P1_B_M->Timing.stepSize0 = 0.001;

  /* InitializeConditions for TransferFcn: '<Root>/Plant' */
  PartialImplemAWES3DOF_P1_B_X.Plant_CSTATE[0] = 0.0;
  PartialImplemAWES3DOF_P1_B_X.Plant_CSTATE[1] = 0.0;

  /* InitializeConditions for ZeroPole: '<S1>/Controller' */
  PartialImplemAWES3DOF_P1_B_X.Controller_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S36>/Integrator' */
  PartialImplemAWES3DOF_P1_B_X.Integrator_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S31>/Filter' */
  PartialImplemAWES3DOF_P1_B_X.Filter_CSTATE = 0.0;
}

/* Model terminate function */
void PartialImplemAWES3DOF_P1_B_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
