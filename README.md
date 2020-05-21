# standard-PIV-image-generator
Matlab based image generator for standard PIV (Image Processing Toolbox required).

Each Matlab file in this repository uses a different velocity field to generate PIV particle images:
- a second order velocity field,
- a Rankine vortex,
- a Rankine vortex with shear effect.

You can easily change image size, number of particle images, particle image density, amount of synthetic noise and a number of other parameters. 

All the scripts use the same process of image generation. It allows particles to enter and leave the camera field of view between first and second frame.

See "details.pdf" for details and figures about the structure of velocity fields and the generation process.


