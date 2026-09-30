---
name: machine-kinematics
description: "Current Apex printer is R C Z B only — plate angle C is imposed by the point position, B tilts one side 0-185°; an extra linear axis comes with the post-Maker-Faire upgrade"
metadata:
  node_type: memory
  type: project
  originSessionId: df56fa55-73de-4689-8817-bf5ea9359353
  modified: 2026-09-28T16:36:45.042Z
---

Current machine (Maker Faire Rome, Oct 2026): axes R (radius), C (plate rotation), Z, B (head tilt) only — no lateral linear axis. To reach a point, C = polar angle of the point (imposed); B is the ONLY orientation DOF, tilting in the radial plane, one side only, B ∈ [0°, 185°] (B > 0 tilts the nozzle toward the centre, layer normal / overhang toward +R). Planned upgrade (several hundred euros, full machine refactor, NOT before Maker Faire) adds a linear axis, making C free for orientation. The slicer must support BOTH kinematics. The machine currently has a cardboard bed; only one elbow-tube print has succeeded so far (2026-09-28).

**Why:** until 2026-09-28 Apex step8 and the research P2 assumed C free (full sphere of directions), which the current machine cannot execute; measured: head within 5° of the surface normal on only 7-11 % of the smoothing area with C imposed. The claim that Bird already fits this machine is FALSE in part (measured 2026-09-29, Benchy body normals: 16 % lean toward the centre, 2.3 % beyond 45°).

**How to apply:** any head-orientation code (step8, step9/10, P2 accessibility, P5, G-code) must take the kinematic variant from the machine profile; never assume a free C. User tolerance for head vs surface normal on smoothing layers: parameter, provisional 45°. Future machine types announced by the user (2026-09-28): (1) + a "D axis" = second head rotation (perpendicular to the Y axis rotation): C still position-coupled, B+D give orientation; (2) + a plate sliding along X (linear axis): C free; (3) both. Head (B) limits given by the user: max 750°/s, max acceleration 300°/s² (stored in `main._default_machine_profile`). Modelled by `MachineProfile.rotation_axes` + `position_coupled_axes` ([0] when C is imposed, [] when a linear axis frees it). Related: [[maker-faire-rome-deadline]].

D axis (user, 2026-09-29): UNLIMITED rotation (homing by roller endstop). Measured on Benchy body normals: current 2.3 % beyond 45°; + D unlimited 0 % (full coverage, like + Y). Y could be sold as an option. The Bird tilt constraint to reachable directions must be a SELECTABLE parameter (not forced).

Plate: 300 mm diameter (radius 150). Part placement: the UI assembles all models on the plate into one, origin of the final model at the plate centre, then feeds the pipeline; an explicit XY offset parameter is still wanted for tests and the experimental Apex (research: `bird_loop.py --offset=dx,dy`, `<prefix>_shift.npy`).
