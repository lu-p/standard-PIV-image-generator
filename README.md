# standard-PIV-image-generator
Matlab based image generator for standard PIV (Image Processing Toolbox required).

Each Matlab file in this repository uses a different velocity field to generate PIV particle images:
- a second order velocity field,
- a Rankine vortex,
- a Rankine vortex with shear effect.

You can easily change the image size, the number of particle images, the particle image density, the amount of synthetic noise and a number of other parameters. Particles are allowed to enter and leave the camera field of view.

All the scripts use the same process of image generation.

See the pdf file for details about structure of velocity fields and generation process.
