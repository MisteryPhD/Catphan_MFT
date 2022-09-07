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

$N_{plates}$ = [5 5 5 5 5 5 5 5 5 5 4 4 4 3 2] (in respect to the fig.1)

To simplify the coding a bit, let's assume a “polar rectangular” instead of straight rectangular regions for each gauge – “polar rectangular” is a circular section bounded by two angles and two radius values. In the given case each gauge has a small angular size so, this polar rectangular are very close to the straight rectangular (the figure bounded by two x-axis values and two y-axis values), so the using of the “Droege-Morin method” interpretation given above is reasonable. 

To localize each gauge region, it was decided to move a ray starting from the bottom (near the gauge #1 - -90 degrees) in the counterclockwise direction (from the gauge #1 to gauge #15) making a full circle. Moving the ray with the smallest reasonable angle step (angular size of one pixel), the program analyzes values of the pixels that lie on the ray path (fig.2 – fig.3): if some pixels have a value greater than the value of the “free-of-objects” region mean value (20x20 pixels in the center) by 10 standard deviation of the same “free-of-objects” region or more it is classified as a “gauge”. The program memorizes each such “hit” – the first point of the ray with the “large” value, the last point of the ray with the “large” value (“radius begin” and “radius end”), the mean value between the first and the last point and the angular position of the ray. Later, analyzing these hits the program will identify “polar rectangles” for each gauge: “radius begin”, “radius end”, “angle begin”, “angle end” (fig.4). 

Knowing the boundaries (two angles and two radiuses) for each gauge – it is easy to identify the set of 
image points (pixels) that lies inside these boundaries and so to identify the set of values for each 
specific gauge (fig.5). Using these values, it is easy to compute its standard deviation or the mean value 
that is used for MTF computing ( $SD_{gauge}$ ).

To estimate the $M_{metal}$ and $SD_{metal}$ values (that is defined as Mean and Standard Deviation of the 
pixels values at the largest plate) it is required to identify the set of image pixels that belong only to the 
largest plate on the image – to do so, the set of pixels that have been found for the largest (#15) gauge is 
split to two classes using a k-means clustering analysis, and points with the higher values is assumed to 
be the “plates” points.

![Alt text](ray_on_free_region.png?raw=true "Figure 2 Gauge localization: Ray out of the gauge")

![Alt text](ray_on_gauge.png?raw=true "Figure 3 Gauge localization: Ray hit the gauge")

![Alt text](ray_hits.png?raw=true "Figure 4 Ray hits: radius end AND the mean value")

![Alt text](localization_trap.png?raw=true "Figure 5 Gauge localization: find the gauge pixels")
