# Sensor Selection with Composite Features (SS-CF)
 Matlab implementation of the sensor selection algorithm in the paper:
 * [T. Yu, A. Mohammadi, Y. Tan, P. Choong and D. Oetomo, "Sensor Selection With Composite Features in Identifying User-Intended Poses for Human-Prosthetic Interfaces," in IEEE Transactions on Neural Systems and Rehabilitation Engineering, vol. 31, pp. 1732-1742, 2023, doi: 10.1109/TNSRE.2023.3258225.](https://ieeexplore.ieee.org/document/10073539)  

## Quick start
 - Run *Main_SSCF.m* for an example.
 - Sensor selection function: *runSSCF.m*. 

## About Dataset
 * dataset.mat contains one non-disabled subject performing forward reaching to spatial target points with their upper limb. 
 * The goal is to select $q$ number of sensors that can best differentiate the three elbow poses required by the spatial target points. 
 * Upper limb and upper body joint kinematics and sEMG signals are recorded with features extracted.
 * Full dataset of 10 subjects available at: [https://doi.org/10.26188/23294693](https://doi.org/10.26188/23294693).

## Acknowledgement
We would appreciate your acknowledgement by citing our paper:

``` 
@article{Yu2023,
author={Yu, Tianshi and Mohammadi, Alireza and Tan, Ying and Choong, Peter and Oetomo, Denny},journal={IEEE Transactions on Neural Systems and Rehabilitation Engineering}, 
title={Sensor Selection With Composite Features in Identifying User-Intended Poses for Human-Prosthetic Interfaces}, 
year={2023},
volume={31},
number={},
pages={1732-1742},
doi={10.1109/TNSRE.2023.3258225}}
``` 

## License
[![License: CC BY 4.0](https://licensebuttons.net/l/by/4.0/80x15.png)](https://creativecommons.org/licenses/by/4.0/) CC BY 4.0 
