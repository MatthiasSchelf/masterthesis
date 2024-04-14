library(data.table)
library(R.matlab)

# Load data
sub_098 <- fread("D:/UGent_gerelateerd/Masterproef/Data/Pupilprepro/sub-098_task-memory_pupil.tsv")

# Extracting variables 2, 7 from DT
variables_to_export <- sub_098[, c(2, 7)]

# Naming the variables using the names retrieved earlier
variable_names <- names(sub_098)[c(2, 7)]

# Assigning names to the columns of variables_to_export
setnames(variables_to_export, variable_names)

# Exporting the variables as a MATLAB .mat file
writeMat("sub_098.mat", variables_to_export = variables_to_export)


