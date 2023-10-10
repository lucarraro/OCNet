# OCNet 1.1.0.9000

## Minor changes

* `create_general_contour_OCN`: documentation updated.
* `OCN_to_AEM`: error message updated.
* `OCN_to_AEM`: `moranI` option added.

# OCNet 1.1.0

## Major changes

* `paths_OCN`: `whichNodes` argument and relative example added.
* `OCN_to_AEM`: function added.

## Bugs fixed

* `paths_OCN`: bug fixed for rivers with multiple outlets.

# OCNet 1.0.1

## Minor changes

* `landscape_OCN`: faster algorithm
* `aggregate_OCN`: faster algorithm (via `Rcpp` implementation) 
* `paths_OCN`: faster algorithm (via `Rcpp` implementation) 

## Bugs fixed

* `aggregate_OCN`: bug fixed for `thrA = 0`.

# OCNet 1.0.0

## Major changes

* OCNs are objects of the `river` class. 
* `plot`, `show`, `names`, `"$<-"`, `"$"`, `"[["`, `"[[<-"` methods for `river` class defined and examples added. 
* `aggregate_OCN`: bug in use of option `maxReachLength` corrected and example added. Note that location and number of AG nodes in long reaches may be altered with respect to previous package versions. 
* `draw_thematic_OCN`: default value for argument `chooseAggregation` changed.

# OCNet 0.7.1

## Minor changes

* `draw_elev2D_OCN`, `draw_elev3D_OCN`: `args_imagePlot` argument added.
* `draw_subcatchments_OCN`: `add` and `args_imagePlot` arguments added.
* `draw_thematic_OCN`: `add`, `args_imagePlot`, `args_legend` arguments and relative examples added. 

# OCNet 0.7.0

## Major changes

* `create_OCN`: initialization algorithm improved.
* `landscape_OCN`: algorithm improved.
* `aggregate_OCN`: `displayUpdates` argument added. Improved algorithm.
* `paths_OCN`: `level` argument added.
* `draw_thematic_OCN`: algorithm improved.
* `draw_subcatchments_OCN`: `theme` argument and relative example added.
* `min_lwd` and `max_lwd` added to all drawing functions.

## Bugs fixed

* `paths_OCN`: Fixed bug when number of nodes is very large

# OCNet 0.6.0

## Major changes

* `create_general_contour_OCN`: function added. 
* `draw_thematic_OCN`: `...` argument and relative example added.
* `aggregate_OCN`: `breakpoints` argument and relative example added.
* `OCN_to_SSN`: `obsSites`, `predSites`, `randomAllocation` arguments added.

## Minor changes

* All functions have been updated for consistence/compatibility with `create_general_contour_OCN`.
* `draw_thematic_OCN`: default option for `theme` added; inverted order of input variables (with backward compatibility maintained); functioning of `colLevels` improved.

## Bugs fixed

* `paths_OCN`: corrected bug for `includeDownstreamNode = TRUE`.

# OCNet 0.5.1

## Bugs fixed

* `paths_OCN`: corrected bug in export of `downstreamPath`. Option `displayUpdates` added. Example in `continue_OCN` modified.

# OCNet 0.5.0

## Major changes

* `paths_OCN`: improved algorithm. Option `includeUnconnectedPaths` added. Option `includePaths` replaces `pathsRN`.

# OCNet 0.4.0

## Major changes

* `continue_OCN`: function added

## Bugs fixed

* `paths_OCN`: Fixed bug when number of nodes is very large

## Others

* `create_OCN`: `OCN$energyInit`, `OCN$FD$perm` are added to output

# OCNet 0.3.2

* * `CITATION` updated

# OCNet 0.3.1

## Bugs fixed

* `draw_contour_OCN`: corrected issue in example 2

# OCNet 0.3.0

## Major changes

* `aggregate_OCN`, `find_area_theshold_OCN`: functioning of `maxReachLength` option has changed 

## Bugs fixed

* `aggregate_OCN`: fixed issue in calculation of `RN$toCM` and `AG$toCM`
* `aggregate_OCN`, `find_area_theshold_OCN`: `maxReachLength` and `streamOrderType` are returned 
* `draw_elev2D_OCN`: added option `addLegend`
* `draw_thematic_OCN`: added option `nodeType`, documentation updated

# OCNet 0.2.0

## Major changes

* `create_OCN`: changed `coolingRate` definition 
* `create_OCN`: perimeteral pixels next to outlets cannot be rewired
* `create_OCN`: default values for `nIter`, `coolingRate` changed
* `create_OCN`: changed initial network state when `nOutlet > 1`
* all uploaded data have been re-created in accordance with the new settings

## Bugs fixed

* `create_OCN`: allowed `nIter < 2`
* `draw_thematic_OCN`: fixed `backgroundColor` issue when `chooseCM = TRUE`
* `landscape_OCN`: fixed `displayUpdates` issue
* `paths_OCN`: added option `includeDownstreamNode` 
* `find_area_theshold_OCN`: allowed `nNodesRN < 2`
* `create_OCN`: added option `easyDraw`
* `draw_thematic_OCN`: added option `chooseAggregation`

## Others

* CITATION added
* vignette: metapopulation model updated
