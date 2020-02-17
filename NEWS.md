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
* vignette: metapopulation model updated
* `find_area_theshold_OCN`: allowed `nNodesRN < 2`
* `create_OCN`: added option `easyDraw`
* `draw_thematic_OCN`: added option `chooseAggregation`


