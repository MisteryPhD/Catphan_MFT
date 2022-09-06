# Catphan_MFT
Modulation Transfer Function evaluation based on Catphan phantom CBCT images

**Given data**: *Catphan phantom  CBCT images*.

**The GOAL**: is to compute *MTF* (Modulation Transfer Function) for the apparatus used to obtain the given set of CBCT images of Catphan phantom.

Catphan phantom contains several types of gauges - one of them (that is imaged on 66, 67 and 68 CBCT slices) are High Resolution test gauges constructed from equidistantly separated metal plates forming several groups with different "spatial frequency" (67 slice is presented on the figure below - fig.1)

![Alt text](localization_trap_all.png?raw=true "Figure 1 Catphan phantom - CBCT slice #67 (high resolution gauge marked by a red line trapezoid)")

$MTF(f) = \frac{\pi\sqrt{2}}{2} \frac{M(f)}{M_0}$, $M(f) = \sqrt{SD_{gauge}^2-SD^2}$, SD^2 = \frac{SD_{metal}^2+SD_{background}^2}{2}
```
section
```
