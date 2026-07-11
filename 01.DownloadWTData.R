# ---
# title: BAM dataset - download WildTrax data
# author: Elly Knight
# created: March 2, 2026
# ---

#NOTES################################

#PURPOSE: This script downloads WildTrax data for the BAM dataset. This secript and the subsequent sets of script should be run biannually to update the BAM dataset for modelling and use by contributing scientists.

#We use a for loop because it allows for debugging and evaluation of project download status

#There are a handful of projects that are not downloading properly via wildRtrax due to special characters. An issue is open on this. These projects are listed in the error.log object. This output should be reviewed and any additional datasets desired should be downloaded manually and incorporated into the dataset in script 02. Note that many of these datasets are single species call rate datasets that will be of less interest for BAM's purposes.

#Projects that are removed for exclusion are projects for which we have access to the entire organization (e.g., CWS), but where the organization managers would prefer we exclude the projects from the dataset. This is instead of providing access to individual projects within the organization, which can be unwieldy for large organization. Hopefully, an "exclude" button will be available within WildTrax sharing in the future.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(wildrtrax) #to download data from wildtrax
library(readxl) #to read excel files

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"

#3. Login to WildTrax----
source("WTlogin.R")
wt_auth()

#INVENTORY#################

#See printed emails in "Data Assessment/CWS engagement" on the shared drive for further details on the preferences of each CWS region for how to handle data use.

#1. Get project list -----
proj.aru <- wt_get_projects("ARU")
proj.pc <- wt_get_projects("PC")

#2. Get list for exclusion for CWS-PRA ----
ex <- read_excel(file.path(root, "Dataset Assessment", "Exclusion", "CWS-PRA_projects_Apr2026.xlsx")) |> 
  dplyr::filter(Sharing=="Exclude")

#3. Filter project list ----
#remove certain project statuses for CWS-ONT as per preferences
proj <- rbind(proj.aru, proj.pc) |> 
  dplyr::filter(project_status!="Test Only",
                tasks_completed > 0,
                !project_id %in% ex$project_id,
                !(organization_name=="CWS-ONT" & project_status %in% c("Active", "Published = Private")))

#DOWNLOAD ###############

#1. Set up loop ----
n_proj = nrow(proj)
aru.list <- vector(mode = "list", length = n_proj)
pc.list <- vector(mode = "list", length = n_proj)
error.log <- data.frame()
for(i in seq_len(n_proj)){
  
  #authenticate each time because this loop takes forever
  wt_auth()
  
  #2. Download ----
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  if(proj$project_sensor[i]=="ARU"){
    
    dat.try <- try(wt_download_report(project_id = proj$project_id[i], sensor_id = "ARU", report = "main"))
    
    if("data.frame" %in% class(dat.try)){
      aru.list[[i]] <- dat.try
    }
    
  }
  
  if(proj$project_sensor[i]=="PC"){
    
    dat.try <- try(wt_download_report(project_id = proj$project_id[i], "PC", report="main"))
    
    if("data.frame" %in% class(dat.try)){
      pc.list[[i]] <- dat.try
    }
    
  }
  
  #Log projects that error
  if(!"data.frame" %in% class(dat.try)){
    error.log <- rbind(error.log, 
                       proj[i,])
    
  }
  
  print(paste0("Finished dataset ", proj$project[i], " : ", i, " of ", nrow(proj), " projects"))
  
}

for (i in seq_len(n_proj)) {
  
  if (!is.data.frame(pc.list[[i]]) && !is.data.frame(aru.list[[i]])) {
    message("retrying i = ", i)
    
    if(proj$project_sensor[i]=="ARU"){
      
      dat.try <- try(wt_download_report(project_id = proj$project_id[i], sensor_id = "ARU", report = "main"))
      
      if("data.frame" %in% class(dat.try)){
        aru.list[[i]] <- dat.try
      }
      
    }
    
    if(proj$project_sensor[i]=="PC"){
      
      dat.try <- try(wt_download_report(project_id = proj$project_id[i], "PC", report="main"))
      
      if("data.frame" %in% class(dat.try)){
        pc.list[[i]] <- dat.try
      }
      
    }
    
    if ("data.frame" %in% class(dat.try)) error.log = error.log %>% dplyr::filter(project_id != proj$project_id[i])
    
  }
  
}

#PUT TOGETHER##############

#1. Name the lists ----
names(aru.list) <- proj$project_id[1:length(aru.list)]
names(pc.list) <- proj$project_id[1:length(pc.list)]

#2. Take out the empty objects ----
aru.wt <- aru.list[!sapply(aru.list, is.null)]
pc.wt <- pc.list[!sapply(pc.list, is.null)]

#3. Get the list of errored projects ----
#Note most are published-map only projects that we don't have organization-level access. We don't filter these out earlier because we do have access to some published-map only projects for which we have organization access.

error <- tryCatch({
  data.frame(file = list.files(file.path(root, "WildTrax", "website downloads (error projects)"))) %>% 
    separate(file, into=c("sensor", "project_id", "report"), remove=FALSE)},
  error = function(e) data.frame())

#4. Read in errored files ----
for(i in seq_len(nrow(error))){
  
  if(error$sensor[i]=="ARU"){
    aru.wt[[length(aru.wt)+1]] <- read.csv(file.path(root, "WildTrax", "website downloads (error projects)", error$file[i]))
    names(aru.wt)[[length(aru.wt)]] <- error$project_id[i]
  }
  
  if(error$sensor[i]=="PC"){
    pc.wt[[length(pc.wt)+1]] <- read.csv(file.path(root, "WildTrax", "website downloads (error projects)", error$file[i]))
    names(pc.wt)[[length(pc.wt)]] <- error$project_id[i]
  }
  
  print(paste0("Finished dataset ", i, " of ", nrow(error), " projects"))
  
}

#5. Save date stamped data & project list----
out_dir = file.path(root, "WildTrax", Sys.Date())
if (!dir.exists(out_dir)) dir.create(out_dir)
save(aru.wt, pc.wt, proj, error.log, file=file.path(out_dir, paste0("01_wildtrax_raw_", Sys.Date(), ".Rdata")))
