
#====================================================#
#=== step 1:  specify output folder =================#
#====================================================#

today = paste0("10152024")
output_folder = paste0("output_",today)
if (!file.exists(output_folder)){dir.create(output_folder)}

## SLURM_ARRAY_TASK_ID on bluehive
### this is a unique id for each task when you submit it; if your #SBATCH -a 1-500, this number would be 1 to 500;
### I usually use this number as seed for data generation
SLURM_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

source("fun_USMART_V1.R")

seed = 10000+SLURM_ID
set.seed(seed)

sample_size = 1000
patient_id = sample(1:pop_size, sample_size, replace = FALSE)
sim1 = population.data[patient_id,]

results = tryCatch(simul( sim1, seed), error = function(e) {print(e); NA})

save(results, sample_size, seed, 
     file = paste0(output_folder,"/USMART_sim",seed,"setting_1.rda"))





