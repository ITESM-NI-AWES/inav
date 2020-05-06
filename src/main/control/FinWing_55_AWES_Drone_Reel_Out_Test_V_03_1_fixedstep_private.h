/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_private.h
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

#ifndef RTW_HEADER_FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_private_h_
#define RTW_HEADER_FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_private_h_
#include "rtwtypes.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetFirstInitCond
# define rtmSetFirstInitCond(rtm, val) ((rtm)->Timing.firstInitCondFlag = (val))
#endif

#ifndef rtmIsFirstInitCond
# define rtmIsFirstInitCond(rtm)       ((rtm)->Timing.firstInitCondFlag)
#endif

#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern uint32_T plook_binx(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction);
extern real_T intrp1d_l_pw(uint32_T bpIndex, real_T frac, const real_T table[]);
extern int32_T plook_s32dd_binx(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction);
extern real_T intrp2d_s32dl_pw(const int32_T bpIndex[], const real_T frac[],
  const real_T table[], const uint32_T stride);
extern real_T intrp1d_s32dl_pw(int32_T bpIndex, real_T frac, const real_T table[]);
extern uint32_T binsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex);
extern int32_T binsearch_s32d(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex);

/* private model entry point functions */
extern void FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_derivatives
  (void);

#endif
/* RTW_HEADER_FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_private_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
