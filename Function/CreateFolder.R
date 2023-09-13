#______________________________________________________________________________#
#     Function to check directory and create if don't exist                    #
#______________________________________________________________________________#

CreateFolder <- function(SubDir){
  if (! dir.exists(SubDir)){        ## Check exist directory
    dir.create(SubDir)            ## Create directory
  }
}
