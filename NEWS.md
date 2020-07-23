# OCNet 0.3.1.9999

* `CITATION` added

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
