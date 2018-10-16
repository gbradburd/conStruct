# conStruct 1.0.1

## Major changes
 + following move to C++14 by Stan
 + `structure2conStruct` now works for multiple STRUCTURE file formats

## Bug fixes
 + users can now specify their own custom plotting colors in `make.all.the.plots` 

# conStruct 1.0.0

## Major changes
 + stan model blocks are now compiled at package installation instead of at a call to `conStruct` or `x.validation`.
 + `x.validation` is now parallelizable
 + new `model-comparison` vignette (see `vignette("model-comparison",package="conStruct")`)
 + alphaD parameter is now rescaled to reflect non-normalized geographic distances
 + compare.two.runs function added

## Bug fixes
 + Removed large files in git history on repo

# conStruct 0.0.0.9000

## Beta release