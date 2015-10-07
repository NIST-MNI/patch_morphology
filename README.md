# ITK patch morphology tools: denoising, segmentation, grading

## implements methods from:
  * ** denoising **: Pierrick Coupé, Pierre Yger, Christian Barillot 
                     "Fast Non Local Means Denoising for 3D MR Images" 
                     http://dx.doi.org/10.1007/11866763_5
  * ** adaptative denoising **: José V. Manjón PhD,Pierrick Coupé PhD, Luis Martí-Bonmatí PhD, D. Louis Collins PhD and Montserrat Robles PhD 
                     "Adaptive non-local means denoising of MR images with spatially varying noise levels" 
                     http://dx.doi.org/10.1002/jmri.22003
  * ** segmentation **: Pierrick Coupé, , José V. Manjón, Vladimir Fonov, Jens Pruessner, Montserrat Robles, D. Louis Collins 
                     "Patch-based segmentation using expert priors: Application to hippocampus and ventricle segmentation" 
                     http://dx.doi.org/10.1016/j.neuroimage.2010.09.018
  * ** grading **: Pierrick Coupé, Simon F Eskildsen, José V Manjón, Vladimir S Fonov, D Louis Collins, Alzheimer's disease Neuroimaging Initiative 
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

### Grading
Similar to segmetnation tool, grading require library of presegmented examples ```$train``` and optionally search radius and patch radius
Training examples are referenced in a comma separated file in the format: ```<image.mnc>,<labels.mnc>,<grading>```
Two training libraries can be provided , which are both loaded (and optionally each is used independently for pre-selection)
```
itk_patch_morphology input.mnc --grading output_grading.mnc  --search $search_radius --patch $patch_radius --train $train --train2 $train2

```