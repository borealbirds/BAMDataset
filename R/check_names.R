# Adds columns to a data.frame if necessary
#
# dat: data.frame
# nms: character vector of proposed names to add; if they already exist in "dat", nothing happens
# fill: value to assign to each column
#
# Returns a data.frame
check_names = function(dat,
                       nms,
                       fill = NA) {
  
  for (nm in nms) {
    if (nm %notin% names(dat)) dat[[nm]] = fill
  }
  
  dat
  
}