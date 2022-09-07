# Catphan_MFT
Modulation Transfer Function evaluation based on Catphan phantom CBCT images

**Given data**: *Catphan phantom  CBCT images*.

**The GOAL**: is to compute *MTF* (Modulation Transfer Function) for the apparatus used to obtain the given set of CBCT images of Catphan phantom.

Catphan phantom contains several types of gauges - one of them (that is imaged on 66, 67 and 68 CBCT slices) are High Resolution test gauges constructed from equidistantly separated metal plates forming several groups with different "spatial frequency" (67 slice is presented on the figure below - fig.1)

![Alt text](localization_trap_all.png?raw=true "Figure 1 Catphan phantom - CBCT slice #67 (high resolution gauge marked by a red line trapezoid)")

$$MTF(f) = \frac{\pi\sqrt{2}}{2} \frac{M(f)}{M_0},    M(f) = \sqrt{SD_{gauge}^2-SD^2}$$

$$SD^2 = \frac{SD_{metal}^2+SD_{background}^2}{2},    M_0 = \frac{|M_{metal}-M_{background}|}{2}$$


Where “SD” means the standard deviation, “M” means the mean value, “gauge” index means the rectangular region on the image that contains the gauge, “metal” index means the region on the image that contains only the metal (here the largest plate on the image is used), “background” index means the region on the image where are no gauge located (here the 20x20 pixels region in the center of the image is taken).  
The spatial frequency for each gauge, “f”, is computed simply as the amount of plates in the gauge divided by its length (the length is the distance between the first and the last plate in the gauge).

$$f = N_{plates}/L$$

Assuming that amount of plates in each gauge is a known data the task is transforming into the “gauge 
localization” task:

$$N_{plates} = [5 5 5 5 5 5 5 5 5 5 4 4 4 3 2]$$ (in respect to the fig.1)
