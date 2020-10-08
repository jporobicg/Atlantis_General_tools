## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~               Transformation from csv to NC            ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
library(ncdf4)
## diccionary
dic <- read.csv('dic_var.csv')
## reading csv files
folder <- paste0('Folders/',dir('Folders/'))
for ( fol in folder){
    files   <- dir(fol)
    cat(paste0('the numer of files in ',fol, 'is ',length(files),'\n'))
    for(f in files){
        ## creating the name for the new file
        t.name   <- unlist(strsplit(gsub('.csv', '.nc', f), '_'))
        if(any(t.name == "ipsl")){ ## for specific naming
            t.name[which(t.name == "ipsl")] <- 'ipsl-cm6a-lr'
        }
        t.name[which(t.name %in% c('region','regional'))]   <- t.name[1]
        out.name <- paste0('Files4upload/',paste(t.name[-1], collapse = '_'))
        ## reading csv files
        tmp.f    <- read.csv(paste0(fol,'/',f))
        ## -- Crating Netcdf files -- ##
        ##----------------------------##
        ## getting the dimensions
        Month_vector <- tmp.f$timestep
        ## Define dimensions
        monthly      <- ncdim_def('Monthly', 'month', Month_vector)
        ## defining the data
        name.var <- colnames(tmp.f)[2]
        longname <- as.character(dic[which(dic[, 2] == name.var), 1])
        variable <- ncvar_def(name.var, 'g m-2', monthly, missval = 1.e+20, longname)
        nc.out   <- nc_create(out.name, list(variable))
        ncvar_put(nc.out, variable, tmp.f[, 2])
        nc_close(nc.out)
    }
}
