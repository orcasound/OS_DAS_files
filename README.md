## These Python scripts create modeled DAS fibers and modeled arrival time differences for chosen sources.

### DAS_defineFiber.py   -----> outputs a binary (pickle) summary of the fiber

## Multiple signal sources can be used to improve the localization of the DAS hydrophones.

### DAS_locateArray.py   -----> outputs plots and adds optimized hydrophone locations to the binary summary of the fiber


#### Here is an output for a "L" shaped fiber. The fiber descends at 30 degrees to a depth of 100 m and then turns "North"
#### Red is the 'actual' (i.e. assumed fiber) and Blue is the located array after localization using a number of locating signals.
![3d PLOT OF FIBER](LshapedFiber.png)

## A localized fiber is then used to calculate source locations.

### DAS_locateSource.py

#### An assumed source is used to determine time differences recorded on the fiber.
#### These time differences and the localized fiber are used to localize the 'unknown' source
#### Calculated Source location is  [[501.04589546 -99.55370506  -1.29196534]]
#### 'Actual' source location is  [500, -100, 0]

## DAS_classes.py has a number of class and helper functions.
