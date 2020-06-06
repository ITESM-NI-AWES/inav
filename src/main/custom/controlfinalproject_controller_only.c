/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: controlfinalproject_controller_only.c
 *
 * Code generated for Simulink model 'controlfinalproject_controller_only'.
 *
 * Model version                  : 1.10
 * Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
 * C/C++ source code generated on : Fri Jun  5 22:51:29 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "controlfinalproject_controller_only.h"
#include "controlfinalproject_controller_only_private.h"

/* Block states (default storage) */
DW_controlfinalproject_contro_T controlfinalproject_controll_DW;

/* External inputs (root inport signals with default storage) */
ExtU_controlfinalproject_cont_T controlfinalproject_controlle_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_controlfinalproject_cont_T controlfinalproject_controlle_Y;

/* Real-time model */
RT_MODEL_controlfinalproject__T controlfinalproject_controll_M_;
RT_MODEL_controlfinalproject__T *const controlfinalproject_controll_M =
  &controlfinalproject_controll_M_;

/* Model step function */
void controlfinalproject_controller_only_step(void)
{
  real_T DiscreteAutopilotForwardDiffe_h;

  /* DiscreteTransferFcn: '<Root>/Discrete Autopilot Forward Differencing' incorporates:
   *  Inport: '<Root>/Error Signal'
   */
  DiscreteAutopilotForwardDiffe_h = controlfinalproject_controlle_U.ErrorSignal
    - -0.9 * controlfinalproject_controll_DW.DiscreteAutopilotForwardDiffere;

  /* Outport: '<Root>/Manipulated Variable' incorporates:
   *  DiscreteTransferFcn: '<Root>/Discrete Autopilot Forward Differencing'
   */
  controlfinalproject_controlle_Y.ManipulatedVariable = 2.8 *
    DiscreteAutopilotForwardDiffe_h + -2.7159999999999997 *
    controlfinalproject_controll_DW.DiscreteAutopilotForwardDiffere;

  /* Update for DiscreteTransferFcn: '<Root>/Discrete Autopilot Forward Differencing' */
  controlfinalproject_controll_DW.DiscreteAutopilotForwardDiffere =
    DiscreteAutopilotForwardDiffe_h;
}

/* Model initialize function */
void controlfinalproject_controller_only_initialize(void)
{
  /* (no initialization code required) */
}

/* Model terminate function */
void controlfinalproject_controller_only_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
