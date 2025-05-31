# Matryoshka.jl
_Global trajectory optimization using Matryoshka mapping. [1]_

This is a Julia version of previously published Python's [ng_trajectory](https://github.com/jara001/ng_trajectory). Configuration files and logs should be compatible between both programs. Whereas ng was used mostly for testing, Julia version implements only the more promissing approaches.

## Implemented features
_Note: Bold features are default options._

- [ ] Optimizers
  - [X] **matryoshka**
  - [ ] ~~braghin~~
- [ ] Criterions
  - [ ] curvature
  - [X] jazar_model
  - [ ] length
  - [X] **profile**
  - [X] manual
- [X] Interpolators
  - [X] none (partial)
  - [X] **cubic_spline**
- [ ] Segmentators
  - [ ] euclidean
  - [X] **flood_fill**
- [ ] Selectors
  - [ ] curvature
  - [ ] uniform
  - [ ] curvature_sample
  - [X] uniform_distance
  - [ ] fixed
  - [X] **uniform_time**
  - [ ] curvature2
- [ ] Penalizers
  - [ ] segment
  - [X] **curvature**
  - [ ] centerline
  - [X] count
  - [X] none
  - [ ] borderlines
- [X] Cascade evaluation
- [ ] Dynamic plotting
- [ ] Plotting of Matryoshka segments
- [ ] Plotting of flood_fill segmentator
- [X] Resume from logs
- [ ] Selector rotation
- [X] Variating parameter
