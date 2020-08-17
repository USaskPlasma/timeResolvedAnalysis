# Time-Resolved Langmuir analysis
MATLAB code for time-resolved Langmuir analysis <br />
Written for use with SIGLENT SDS 1104X-E scope. <br />
 
### How to take time-resolved measurements
* Connect oscilloscope channels as follows:
  * Ch 1 ---> HV pulser voltage
  * Ch 2 ---> HV pulser current (rogowski coil)
  * Ch 3 ---> Probe current
  * Ch 4 ---> Probe voltage
* Vary the probe voltage from ~-40 to ~+40 V  in steps of 0.5 V (depending on conditions of discharge).
* Record data at each voltage step as .csv using scope.
* Load dataset into code.
* Plot I-V curves in 3D and 2D to identify areas of bad fitting.
* Make changes to excluded points in fit until good fit is achieved for most of the timesteps.

<br />
Output examples: I-V curve, 3d plot and plasma parameters <br />
<p align="center">
<img src="https://github.com/USaskPlasma/timeResolvedAnalysis/blob/master/I-Vexample.png" width="400"> 
<img src="https://github.com/USaskPlasma/timeResolvedAnalysis/blob/master/3dplotexample.png" width="400">
<img src="https://github.com/USaskPlasma/timeResolvedAnalysis/blob/master/Outputexample.png" width="400"> <br />
</p>

