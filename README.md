# ITK patch morphology tools: denoising, segmentation, grading

## implements methods from:
  * __denoising__: Pierrick Coupé, Pierre Yger, Christian Barillot 
                     "Fast Non Local Means Denoising for 3D MR Images" 
                     http://dx.doi.org/10.1007/11866763_5
  * __adaptative denoising__: José V. Manjón PhD,Pierrick Coupé PhD, Luis Martí-Bonmatí PhD, D. Louis Collins PhD and Montserrat Robles PhD 
                     "Adaptive non-local means denoising of MR images with spatially varying noise levels" 
                     http://dx.doi.org/10.1002/jmri.22003
  * __segmentation__: Pierrick Coupé, , José V. Manjón, Vladimir Fonov, Jens Pruessner, Montserrat Robles, D. Louis Collins 
                     "Patch-based segmentation using expert priors: Application to hippocampus and ventricle segmentation" 
                     http://dx.doi.org/10.1016/j.neuroimage.2010.09.018
                     
                     "MICCAI 2012 Workshop on Multi-Atlas Labeling" ISBN-10: 1479126187 entry BIC-IPL and BIC-IPL-HR,
                     https://masi.vuse.vanderbilt.edu/workshop2012/images/c/c8/MICCAI_2012_Workshop_v2.pdf
                     
                     Katrin Weier, Vladimir Fonov, Karyne Lavoie, Julien Doyon and D. Louis Collins
                     "Rapid automatic segmentation of the human cerebellum and its lobules (RASCAL) 
                     Implementation and application of the patch-based label-fusion technique with a template library to segment the human cerebellum"
                     http://dx.doi.org/10.1002/hbm.22529
                     
  * __grading__: Pierrick Coupé, Simon F Eskildsen, José V Manjón, Vladimir S Fonov, D Louis Collins, Alzheimer's disease Neuroimaging Initiative 
                     "Simultaneous segmentation and grading of anatomical structures for patient's classification: application to Alzheimer's disease" 
                     http://dx.doi.org/10.1016/j.neuroimage.2011.10.080
                     
## Building:
* Building require ITK version 4.XX and cmake version 2.6 or later

## Usage (concise):
General notes, these programs were originally designed to be used with MINC files, but should work with any file format supported by ITK version 4
Most tools have an optional paramter called search radius - which specifies how non-local the search should be (in voxels)
and patch radius - the radius of the local patch used to extract features.

### Denoising
Denoising requires specifying noise level (```$sigma```), and optionall search radius ```$search_radius``` and patch radius ```$patch_radius``` 
```sh
itk_minc_nonlocal_filter input.mnc output.mnc --noise $sigma --search $search_radius --patch $patch_radius
```

### Adaptative denoising
Adaptative denoising  have optional parameters: search radius ```$search_radius``` and patch radius ```$patch_radius``` 
```sh
itk_minc_nonlocal_filter input.mnc output.mnc --search $search_radius --patch $patch_radius --anlm
```

### Segmentation
Segmetnation tool require library of presegmented examples ```$train```, number of classes including background ```$classes``` and optionally search radius and patch radius
Training examples are referenced in a comma separated file in the format: ```<image.mnc>,<labels.mnc>```

```
itk_patch_morphology input.mnc output_labels.mnc --discrete $classes --search $search_radius --patch $patch_radius --train $train
```

Several high level scripts are included in ```scripts``` directory:
 * ```ventricles_segmentation_pipeline.pl``` - segmentation script for latera ventricle segmentation, uses volume_patches program from legacy directory
 * ```hcag_segmentation_pipeline.pl``` - Hippocampus and Amygdala segmentation script
 * ```miccai2012_segmentation_minipipe.pl``` - whole head segmentation pipeline, used in MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling
 * ```patch_segmentation_pipeline.pl``` - generic segmetnation script used in RASCAL paper


### Grading
Similar to segmetnation tool, grading require library of presegmented examples ```$train``` and optionally search radius and patch radius
Training examples are referenced in a comma separated file in the format: ```<image.mnc>,<labels.mnc>,<grading>```
Two training libraries can be provided , which are both loaded (and optionally each is used independently for pre-selection)
```
itk_patch_morphology input.mnc --grading output_grading.mnc  --search $search_radius --patch $patch_radius --train $train --train2 $train2

```

High level script for simultaneous grading and segmentation: ```scripts/snipe_grading_pipeline.pl```
