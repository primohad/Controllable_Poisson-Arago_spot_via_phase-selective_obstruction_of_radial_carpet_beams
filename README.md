# Controllable Poisson–Arago Spot via Phase-Selective Obstruction of Radial Carpet Beams

This repository contains the MATLAB code used for simulations, data analysis, and figure generation in the paper:

> **"Controllable Poisson–Arago spot via phase-selective obstruction of radial carpet beams"**
> Mohaddeseh Mohammadi Masouleh, Somaye Fathollazade, Saifollah Rasouli
> (2026)

## Repository Structure

```
code/
├── 01_initial_idea/                                    # Fig. 1
├── 02_generation_initial_idea_with_normalization/       # Fig. 2
├── 02_generation_initial_idea_without_normalization/    # Fig. 2 & Visualizations 1 & 2
├── 03_slm_characterization/                             # Supplement Document Sec. 2
├── 04_intensity_detection_plane/                        # Fig. 4
├── 05_phase_patterns_second_grating/                    # Fig. 5
├── 06_power_in_phase_regions/                            # Supplement Document Sec. 3
├── 07_calibration_theta_effect/                         # Supplement Document Sec. 4
└── 08_power_of_spots/                                   # Fig. 6 & Visualization 3
```

## Requirements

- MATLAB R2020a or later
- Toolboxes: Image Processing Toolbox, Signal Processing Toolbox

## Usage

Navigate to the desired subfolder and run `main.m`:

```matlab
cd code/01_initial_idea
main
```

Each subfolder is self-contained and regenerates the corresponding figure(s) or supplementary result listed above.

## Citation

If you use this code, please cite:

```bibtex
@article{masouleh2026poissonarago,
  title   = {Controllable Poisson--Arago spot via phase-selective obstruction of radial carpet beams},
  author  = {Mohammadi Masouleh, Mohaddeseh and Fathollazade, Somaye and Rasouli, Saifollah},
  year    = {2026}
}
```

## License

MIT License

Copyright (c) 2026 Mohaddeseh Mohammadi Masouleh, Somaye Fathollazade, Saifollah Rasouli

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
