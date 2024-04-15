# Reverb models

Reverberation modelling using image-source modelling and plate modelling.

### image_source method

 Calculates the impulse response of a room using the image-source method using three different algorithms
 
 * basic: contant absorption for all surfaces
 * advanced: variable absorption across surfaces and accounts for stereo separation based on ear distance
 * vectorised: advanced algorithm but vectorised using a 3D array

 ### plate_reverb 

 Calculates the impulse response of a metal plate. Material options are steel or aluminium.
 IR can be convolved with an input signal.