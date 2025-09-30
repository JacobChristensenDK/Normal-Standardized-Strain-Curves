This code interpolates and creates EDS and LDS values.

strain_DATA must be a struct organized as in the following way: strain_DATA -->
patient ID --> 3 fields called x2CH (2 chamber), x4CH (four chamber) and
APLAX (3 chamber). --> up to 9 dobules with Time, 1-6 segment names,left_marker and right_marker.

There must be valve times including AVO, MVO, and AVC.
