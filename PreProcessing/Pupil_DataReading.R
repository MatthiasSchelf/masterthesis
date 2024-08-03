library(data.table)
library(R.matlab)

# Load data
sub_091 <- fread("D:/UGent_gerelateerd/Masterproef/Data/Pupilprepro/sub-091_task-memory_pupil.tsv")

# Extracting variables 2, 7 from DT
variables_to_export <- sub_091[,c(1,2,7)]

# Naming the variables using the names retrieved earlier
variable_names <- names(sub_091)[c(1,2,7)]

# Assigning names to the columns of variables_to_export
setnames(variables_to_export, variable_names)

# Exporting the variables as a MATLAB .mat file
writeMat("sub-091_selected.mat", variables_to_export = variables_to_export)


