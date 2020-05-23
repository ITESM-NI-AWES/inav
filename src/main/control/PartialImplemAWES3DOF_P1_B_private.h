/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: PartialImplemAWES3DOF_P1_B_private.h
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

#ifndef RTW_HEADER_PartialImplemAWES3DOF_P1_B_private_h_
#define RTW_HEADER_PartialImplemAWES3DOF_P1_B_private_h_
#include "rtwtypes.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

/* private model entry point functions */
extern void PartialImplemAWES3DOF_P1_B_derivatives(void);

#endif                    /* RTW_HEADER_PartialImplemAWES3DOF_P1_B_private_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
