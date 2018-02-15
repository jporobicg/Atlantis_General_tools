bords <- function(dat = NULL, lon = NULL, lat = NULL, K = FALSE, distance = NULL, quiet = TRUE, ...){
    ################################################################
    # This fuction can obtain the limits point a certain distance  #
    # Creator   :   Javier Porobic                                 #
    # date      :   10-12-2014                                     #
    # Place     :   CSIRO                                          #
    #
    # -  -  -  -  -  -  - Inputs -  -  -  -  -  -  -  -  -  -  -  -#
    # dat  : data object, with lon and lat
    # lon  : lon coordinates
    # lat  : lat coordinates
    # K    : if the distance is in kilometer
    # the distance by defect is in miles but you can put in kilometers
    # but you need to set up K as TRUE
    R        <- 6371                        #  Average radius of the earth in
    bearing  <- seq(1, 360, 10) * pi / 180   #  Bearing in radians
    if(!is.null(dat)){
        lon <- dat[, 1] * pi / 180          #  In radians
        lat <- dat[, 2] * pi / 180          #  In radians
    } else {
        if(is.null(lat) || is.null(lon)) stop('\n Is necessary to enter a value for Latitude and Longitude')
        if(length(lat) != length(lon))   stop('lon and lat need to have the same length')
        lat <- lat * pi / 180               #  In radians
        lon <- lon * pi / 180               #  In radians
    }
    if(!isTRUE(K)){
        distance <- distance  * 1.852       #  km
    }
    tag <- 'km'
    lat2 <- lon2 <- matrix(NA, length(lat), length(bearing))
    for(j in 1 : length(bearing)){
        lat2[, j]  <- asin(sin(lat) * cos(distance / R) + cos(lat) * sin(distance / R) * cos(bearing[j]))
        lon2[, j]  <- lon + atan2(sin(bearing[j]) * sin(distance / R ) * cos(lat), cos(distance / R) - sin(lat) * sin(lat2[, j]))
    }
    lat2 <- lat2 * 180 / pi                 #  Radians to Degrees
    lon2 <- lon2 * 180 / pi                 #  Radians to Degrees
    lat  <- lat * 180 / pi                  #  Radians to Degrees
    lon  <- lon * 180 / pi                  #  Radians to Degrees
    lat2 <- as.vector(lat2)                 #  From Matrix to vector
    lon2 <- as.vector(lon2)                 #  From Matrix to vector
    pos  <- c(chull(x = lon2, y = lat2), chull(x = lon2, y = lat2)[1])
    poly <-data.frame(latitude = lat2[pos], longitude = lon2[pos])
    if(!isTRUE(K)){
        distance <- distance / 1.852       #  km to Nm
        tag <- 'nm'
    }

    #  Graphics
    if(!isTRUE(quiet)){
        plot(lon, lat, xlab = 'Longitude', ylab = 'Latitude', bty = 'n', pch = 19,
             ylim = c(min(lat2) - mean(diff(lat2)),  max(lat2)) + mean(diff(lat2)),
             xlim = c(min(lon2) - mean(diff(lon2)),  max(lon2)) + mean(diff(lon2)),
             main = paste('Distance from the boundary line \n', round(distance, 3), tag, sep = ' ') )
        with(poly, lines(longitude, latitude, col = 2, lwd = 1.2))
    }

    return(poly)
}

stage <- function(arr.fg, size.col = 2, age.lim = 0){
    ##---------------------------------------------------------------------------##
    ## This function add the division between adults an juvenils in the selected ##
    ## dataframe.                                                                ##
    #----------------------------------------------------------------------------##
    arr.fg$Stage <- ifelse(arr.fg[,size.col] >= age.lim, 'Adult', 'Juvenil')
    arr.fg       <- split(arr.fg, arr.fg$Stage)
    return(arr.fg)
}

hist.deph <- function(list.l, col.val = 1, dep.lim = NULL, n.file = 'Rfile', save = TRUE, relative = FALSE){
    ## list.l : List witht the depth for ear functional groups
    count = 1
    output <- list()
    if(isTRUE(save)) png(paste(n.file, '.png', sep = ''), width = 1200, height = 1200)
    ifelse((length(list.l) > 2 ), par(mfrow = c(3, 4), cex = 1.4), par(mfrow = c(1, 2), cex = 1.4))
    names.out <- ''
    for( i in 1 : length(list.l)){
        data  <- list.l[[i]]
        min.d <- dep.lim[i, 2]
        max.d <- dep.lim[i, 3]
        if(is.null(dim(data))){
            ## For FG with juvenil and adult stages
            for(j in 1 : 2){
                data2             <- data[[j]]
                depth.dist        <- hist(data2[, col.val], breaks = seq(from = min.d, to = max.d, length = 9), plot = FALSE)
                depth.dist$counts <- with(depth.dist, counts / sum(counts))
                with(depth.dist, plot(counts,  - mids, type = 'b', pch = 19, las = 1, bty = 'n',
                                      ylab = ifelse(isTRUE(relative), 'Layer', 'Depth'), xlab  =  'Density', main = paste(data2$FG[1], data2$Stage[1], sep = '\n'),
                                      xlim = c(0, max(counts) * 1.01), ylim = c( - max.d, 0)))
                with(depth.dist, abline(h =  - breaks, lty = 2, col = 'grey'))
                output[[count]]   <- rev(with(depth.dist, counts))
                names.out         <- c(names.out, paste(data2$FG[1], data2$Stage[1], sep = '-'))
                count             <- count + 1
            }
        } else {
            depth.dist            <- hist(data[, col.val], breaks = seq(from = min.d, to = max.d, length = 9), plot = FALSE)
            depth.dist$counts     <- with(depth.dist, counts / sum(counts))
            with(depth.dist, plot(counts,  - mids, type = 'b', pch = 19, las = 1, bty = 'n',
                                  ylab = ifelse(isTRUE(relative), 'Layer', 'Depth'), xlab  =  'Density', main = data$FG[1],
                                  xlim = c(0, max(counts) * 1.01), ylim = c( - max.d, 0)))
            with(depth.dist, abline(h =  - breaks, lty = 2, col = 'grey'))
            output[[count]]       <- rev(with(depth.dist, counts))
            count                 <- count + 1
            names.out             <- c(names.out, data$FG[1])
        }
    }
    if(isTRUE(save)) dev.off()
    output <- matrix(unlist(output), ncol = (count - 1))
    colnames(output) <- names.out[ - 1]
    return(output)
}

points.img <- function(x.lim = NA, y.lim = NA, filename = NULL, grid = FALSE, convert = FALSE, superpose = FALSE, poly = NULL){
    #######################################################
    # This function is usefull to extract data from plots #
    ## x.vals  = minimum and maximum values of X          #
    ## Y.vals  = minimum and maximum values of Y          #
    ## filename  =  Name of the image file                #
    #######################################################
    library(pixmap)
    ## convertion on linux
    system(paste('convert', filename, 'temporal.ppm', sep = ' '))
    #windows user should fine the way to convert
    ## read the pnm file
    img <- read.pnm("temporal.ppm")         #some warnings related with X axis

    ## locator to identified the min and the max
    if(isTRUE(convert)){
    flip    <- function(matriz) t(matriz)[,nrow(matriz):1]
    red.mat <- matrix(NA, img@size[1], img@size[2])
    red.mat <- img@red

    image(flip(red.mat), bty = 'n' ,axes = FALSE)
    } else {
        plot(img)
    }
    if(grid == TRUE) grid(10, col = 1, lwd = .3, lty = 1)

    ## locator to identified the min and the max
    cat('set the lower value of X and the maximun\n')
    x.axs <- locator(2)
    cat('set the lower value of Y and the maximun\n\n\n')
    y.axs  <- locator(2)


    xseq    <- seq(x.axs$x[1], x.axs$x[2], length.out = 10000)
    yseq    <- seq(y.axs$y[1], y.axs$y[2], length.out = 10000)
    x.u.pos <- seq(x.lim[1], x.lim[2], length.out = 10000)
    y.u.pos <- seq(y.lim[1], y.lim[2], length.out = 10000)

    if(isTRUE(superpose)){
        if(is.null(poly)) stop('\n It is necesary provide the values to plot')
        if(!is.list(poly)){
            x.pol <- apply(as.matrix(poly[, 1], 1), 1, function(x)which.min(abs(x - x.u.pos)))
            y.pol <- apply(as.matrix(poly[, 2], 1), 1, function(x)which.min(abs(x - y.u.pos)))
            lines(xseq[x.pol], yseq[y.pol])
        } else {
            for( list in 1 : nrow(poly[[1]])){
                x.temp <- poly[[1]][list, ]
                y.temp <- poly[[2]][list, ]
                x.temp <- matrix(x.temp[!is.na(x.temp)],1)
                y.temp <- matrix(y.temp[!is.na(y.temp)],1)
                x.pol <- apply(x.temp, 2, function(x)which.min(abs(x - x.u.pos)))
                y.pol <- apply(y.temp, 2, function(x)which.min(abs(x - y.u.pos)))
                lines(xseq[x.pol], yseq[y.pol])
                readline("Press <return to continue")
            }
        }
    }

    cat('Select the point,  when ready right - click\n')
    pos     <- locator()
    lines(pos)

    ## find the values
    value.y <- value.x <- rep(NA, length(pos$x))

    for( i in 1 : length(value.x)){
        value.x[i] <- which.min(abs(pos$x[i] - xseq))
        value.y[i] <- which.min(abs(pos$y[i] - yseq))
    }

    value.x <- x.u.pos[value.x]
    value.y <- y.u.pos[value.y]

    values <- data.frame(Y = value.y, X = value.x)

    return(values)
}


convex<-function(data = NULL, lon = NULL, lat = NULL, positive=FALSE){
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # data is a lon-lat object [dataframe]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(!is.null(data)){
        lon <- data[,1]
        lat <- data[,2]
    } else {
        if(is.null(lat) || is.null(lon)) stop('\n Is necessary to enter a value for Latitude and Longitude')
        if(length(lat) != length(lon))   stop('lon and lat need to have the same length')
    }
    if(isTRUE(positive)){
        lon <- lon *  - 1
        lat <- lat *  - 1
    }
    pos  <- chull(x = lon, y = lat)
    pos  <- c(pos, pos[1])
    poly <-data.frame(latitude = lat[pos], longitude = lon[pos])
    return(poly)
}

FG <- function(data, column=1){
    data$FG <- ifelse(data[, column] %in% c('Pez Escorpion', 'Callanthias platei', 'Pampanito', 'Chromis meridiana', 'Paratrimma nigrimenta', 'Scartichthys variolatus','Suezichthys sp.','Malapterus reticulatus','Pseudolabrus gayi','Kyphosus cinerascens','Graniento', 'graniento', 'Paratrachichthys fernandezianus', 'Paralichthys fernandezianus'), 'SPF',
               ifelse(data[, column] %in% c('Corvina', 'Jurel', 'Jurelillo', 'Sierra', 'Umbrina reedi', 'Pez Volador', 'Atun'), 'LPF',
               ifelse(data[, column] ==  'Vidriola', 'VID',
               ifelse(data[, column] %in% c('Chironemus delfini', 'Chalaco', 'Chironemus bicornis', 'Colorado', 'Col', 'Cabrilla', 'Chancharro', 'Pez mariposa', 'Aseraggodes bahamondei', 'Pez Mariposa','Jerguilla', 'Monocentris reedi', 'Parapercis dockinsi', 'Scorpaena thomsoni'), 'SBF',
               ifelse(data[, column] == 'Alfonsino', 'ALF',
               ifelse(data[, column] %in% c('Bacalao', 'Lenguado', 'Paralabrax sp.', 'Aseraggodes bahamondei', 'Congrio', 'Lotella fernandeziana', 'Salmon', 'Jerguillon'), 'LBF',
               ifelse(data[, column] == 'Tollo', 'CHO',
               ifelse(data[, column] == 'Anguila', 'ANG',
               ifelse(data[, column] %in% c('Caracol Blanco', 'Ostra', 'Loco', 'Chiton', 'Gastropod', 'Slug', 'Snail', 'Bivalve'), 'MOL',
               ifelse(data[, column] %in% c('Crab', 'Ce', 'Centolla', 'Jaiva', 'Jaiba', 'Estomatopodo', 'shrimp'), 'SCR',
               ifelse(data[, column] %in% c('Pulpo', 'Octopus'), 'OCT',
               ifelse(data[, column] == 'Lobo', 'OTA',
               ifelse(data[, column] %in% c('Nemadactylus gayi', 'Breca'), 'BRC',
               ifelse(data[, column] %in% c('Delphinus_delphi', 'Stenella_coeruleoalba', 'Tasmacetus_shepherdii', 'tursiops_truncatus'), 'DOL',
               ifelse(data[, column] %in% c('Whale', 'Balaenoptera_borealis', 'Balaenoptera_inusculus', 'Balaenoptera_physalus_quoyi', 'Ziphius_cavirostris','Globicephala_macrorhynchus', 'Kogia_breviceps', 'Mesoplodon_bahamondi', 'Mesoplodones', 'Orcinus_orca', 'Physeter_catodon1980', 'Physeter_catodon'), 'CET',
               ifelse(data[, column] == 'Orange roughy', 'ORO',
               ifelse(data[, column] == 'Langosta', 'SPL',
               ifelse(data[, column] == 'Cangrejo dorado', 'GCR',
               ifelse(data[, column] %in% c('Actinia', 'Amphiura', 'Anemone', 'Erizo', 'Estrella de mar', 'Echinoide'), 'SUR',
               ifelse(data[, column] %in% c('Hidrozoo', 'Esponja','Bernacle', 'Black hydrozoan', 'Briozoo', 'Cucumber', 'Encrusting bryozoan', 'Hydrozoo', 'Polychaete', 'Worm', 'Sponge'), 'BFF',
               ifelse(data[, column] %in% c('Bare Rock', 'Ruble', 'Cobble', 'Rubble'), 'ROC',
               ifelse(data[, column] == 'Sand', 'SAN',
               ifelse(data[, column] %in% c('Brown Alga', 'Green Alga', 'Red Alga', 'Green Microalga', 'Brown Microalga'), 'MA',
               ifelse(data[, column] %in% c('Coral','Coral negro', 'Coral blanco'), 'COR',
               ifelse(data[, column] == 'Orange Roughy', 'ORO', NA)))))))))))))))))))))))))
    return(data)
}

read.poly <- function(data){
    #############################################
    #  This function reads csv files from Qgis  #
    #############################################
    attrib <- read.csv(data)[, - 1]
    coord  <- read.csv(data)
    ## Reeplace the characters from Qgis
    ver1   <- gsub("POLYGON \\(\\(", "", coord[, 1])
    ver2   <- gsub("\\)\\)", "", ver1)
    ## Read the coordinates
    for( i in 1 : length(ver2)){
        ver3 <- read.table(text = gsub(" ", ",", ver2[i]), sep = ',')
        if ( i == 1){
            ver.f <- ver3
        } else {
            if ( dim(ver3)[2] < dim(ver.f)[2]){
                ver3[, (dim(ver3)[2] + 1) : dim(ver.f)[2]] <- NA
            } else if (dim(ver3)[2] > dim(ver.f)[2]) {
                ver.f[, (dim(ver.f)[2] + 1) : dim(ver3)[2]] <- NA
            }
            ver.f <- rbind(ver.f, ver3)
        }
    }
    return(list(coor = ver.f, attrib = attrib))
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Distribution of the species by polygon
##' @param poly List object witht the polygons
##' @param lista List witht he functional groups
##' @param la ID of the column where to find the latitude
##' @param lo ID of the column where to find the longitude
##' @param mult ID of the column where to find the multiplicative
##' @param sea ID of the column where to find the seasons
##' @param only.s  Indicate if is only one seaon
##' @param is.Crus Logic,  inidcating if the species are crustaceon. this is due to the use of traps in JFR
##' @param quiet Logic,  False :  will gave an outpur of each step
##' @return A table witht the relative abundance by polygon of each functional group
##' @author Demiurgo
by.pol <- function(poly, lista, la = 4, lo = 5, mult = 6, sea = 10, only.s = FALSE, is.Crus=FALSE, quiet = TRUE){
    ## Function to clasified the species by polygon
    library(sp)
    lon <- as.matrix(poly$coor[, c(seq(from = 1, to = dim(poly$coor)[2], by = 2))])
    lat <- as.matrix(poly$coor[, c(seq(from = 2, to = dim(poly$coor)[2], by = 2))])
    output <- list()
    lista2 <- lista
    for(season in 1 : ifelse(isTRUE(only.s), 1, 4)){                             # By Season
        for( i in 1 : nrow(poly$attrib)){        # By Polygon
            for(j in 1 : length(lista)){         # By Sample
                lista       <- lista2
                total.traps <- length(unique(lista[[j]][ , la]))
                if(!isTRUE(only.s)){
                    lista[[j]] <- lista[[j]][which(lista[[j]][, sea]  == season), ] # Only samples from the same Season
                }
                ## Vector of 0 and 1 with te location of the samples inpol=1 if not inpol = 0
                inpol <- ifelse(point.in.polygon( - lista[[j]][, lo],  - lista[[j]][, la], lon[i,], lat[i, ]) > 0, 1, 0)
                scal  <- inpol * lista[[j]][, mult]
                if(!quiet & sum(scal) > 0){
                    cat('\n\n =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = ')
                    cat('\n\nSpecies : ', lista[[j]]$FG[1])
                    cat('\nNumber in polygon', sum(scal))
                    cat('\nSeason',  season, '\tPolygon(i)', poly[[2]][i, 1], '\tSample(j)', j)
                }
                if(!isTRUE(is.Crus)){
                    poly$attrib[i, j + 6] <- sum(scal, na.rm = TRUE) / sum(scal > 0)
                #cat('\nsum :', sum(scal), '\t weight :', sum(scal>0), '\tCPUE', poly$attrib[i, j + 6])
                } else  {
                    ## Devide the total catch by the number of traps
                    n.traps   <- length(levels(as.factor(lista[[j]][which(inpol > 0), la])))
                    abun.Weig <- (sum(scal) / n.traps) * (n.traps / total.traps)
                    poly$attrib[i, j + 6] <- abun.Weig
                    if(!quiet & sum(scal) > 0){
                        cat('\nTraps ', n.traps)
                        cat('\nTotal Traps ', total.traps)
                        cat('\nRelative Abundance',  sum(scal) / n.traps)
                        cat('\nRelative Abundance Weigthed',  abun.Weig)
                        readline("\nPress <return to continue")
                    }
                }
            }
        }
        names(poly$attrib)[7 : ((7 + length(lista)) - 1)] <- names(lista)
        output[[season]] <- poly
    }
    return(output)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Relocation of species from inside of a polygon to another closest
##' @param eddies Location of eddies
##' @param centroids Centroids of the polygons
##' @param polygon IDPolygons to remove
##' @return the new position of the eddies
##' @author Demiurgo
move.point <- function(points, centroid, island){
    library(sp)
    ##########################################################
    ##  This function move the points from inside on an area #
    ##  to another area closer                               #
    ##########################################################
    in.pol <- point.in.polygon(points[, 1], points[, 2], island[, 1], island[, 2])
    if(all(in.pol == 0)) return(points)
    pos           <- which(in.pol == 1)
    in.island     <- points[pos, ]
    in.island.rad <- deg2rad(in.island)
    centroid.rad  <- deg2rad(centroid)
    for(i.ed in 1: length(pos)){
        dist.pos <- which.min(gcd.hf(pos1 = matrix(in.island.rad[i.ed, ], ncol = 2), pos2 = -centroid.rad, mat = TRUE))
        points[pos[i.ed], 1 : 2] <-  unlist(centroid[dist.pos, 1 : 2]  *  - 1, use.names = FALSE)
    }
    return(points)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Degrees to radians
##' @param deg Coordinates in degrees
##' @return Object with the Lat and Lon information in Radians
##' @author Demiurgo
deg2rad <- function(deg){
    rad <- deg * pi / 180
    return(rad)
}

##' @description Calculates the geodesic distance between two points specified by radian latitude/longitude using the Haversine formula (hf)
##' @title Haversine distance formula
##' @param long1 Longitude of the first point in radians
##' @param lat1 Latitude of the second point in radians
##' @param long2 Longitude of the first point in radians
##' @param lat2 Latitude of the second point in radians
##' @param mat TRUE or FALSE,  identified if the information is provided in two diferent matrix, one for the initial position and a second matrix for the final position
##' @param pos1 If mat = TRUE this matrix represent the first position witht the first colum for the Longitude and the second for the Latitude
##' @param pos2 If mat = TRUE this matrix represent the second position witht the first colum for the Longitude and the second for the Latitude
##' @return the distance in m for eatch points
##' @author Demiurgo
gcd.hf <- function(long1 = NULL, lat1 = NULL, long2 = NULL, lat2 = NULL, mat = FALSE, pos1 = NULL, pos2 = NULL) {
    if(isTRUE(mat)){
        long1 <- pos1[, 1]
        lat1  <- pos1[, 2]
        long2 <- pos2[, 1]
        lat2  <- pos2[, 2]
    }
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat / 2)^ 2 + cos(lat1) * cos(lat2) * sin(delta.long / 2)^ 2
    if(isTRUE(mat)){
        sq               <- sqrt(a)
        sq[which(sq >1)] <- 1
        c                <- 2 * asin(sq)
    }else{
        c <- 2 * asin(min(1, sqrt(a)))
    }
  d <- R * c
  d <- d * 1000
  return(d) # Distance in m
}

plot.map <- function(poly, xlim = NULL, ylim = NULL, sc = 4, OnlyPoly = TRUE, leg.color = 'normal', specie=FALSE, leg = NULL, save = FALSE, name = NULL, PAR = TRUE, ...){
    #==================================================================#
                                        #     # xlim     = vector of longitude limits                            #
    # ylim     = vector of latitude limits                   island     # sc       = Scale factor (between lat and lon)                    #
                                        # OnlyPoly = Information for lon and lat is obtained from the     # leg.col  = Color of the legend base on the atributes             #
    # leg      = Legend                                      island     # =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = #

    ## Lybraries
    library(maps)
    library(mapdata)
    library(mapproj)
    library(RColorBrewer)

    ## Separe Lat and Lon
    lon <- as.matrix(poly$coor[, c(seq(from = 1, to = dim(poly$coor)[2], by = 2))])
    lat <- as.matrix(poly$coor[, c(seq(from = 2, to = dim(poly$coor)[2], by = 2))])
    ## Color
    if (leg.color == 'normal') {
        color <- rgb(t(col2rgb("royalblue")), alpha = 30, maxColorValue = 255)
        etiq  <- NULL
    } else if (any(leg.color == names(poly$attrib)) && leg.color != 'box_id'){
        var   <- unlist(poly$attrib[leg.color])
        var[which(is.na(var))] <- 0
        div   <- 1
        while(any(var > 1000)){
            div <- div * 10
            var <- round(var / div, 3)
        }
        if(specie == TRUE){
            #browser()
            max.var <- max(var)
            r.m.var <- max(var / sum(var))
            var     <- var / max.var
            max.var <- max(var)
            quart   <- 1 : length(var)
            color   <- rgb(t(col2rgb("royalblue")), alpha = var * 255, maxColorValue = 255)
            etiq     <- ceiling(seq(from = 0, to = r.m.var * 100, length = 4))
            color.et <- rgb(t(col2rgb("royalblue")), alpha = (etiq / 100) * 255, maxColorValue = 255)
        } else {
            quart <- cut(var, breaks = unique(quantile(var, probs = c(1, 0.9, 0.7, 0.5, 0.3, 0))),
                         include.lowest = TRUE)
            color <- brewer.pal(length(levels(quart)), "BuGn")
            etiq  <- round(quantile(var, probs = c(1, 0.9, 0.7, 0.5, 0.3, 0)), 3)
        }
        end   <- length(etiq)
        for( i in 1 : (end - 1)){
            if (i == length(etiq)){
                etiq[i] <- paste(etiq[i], '0', sep = '-')
            } else {
                etiq[i] <- paste(etiq[i], etiq[i + 1], sep = '-')
            }
        }
        etiq <- etiq[-end]
    } else if (leg.color == 'box_id'){
        color <- 'gray90'
        etiq  <- NULL
        lon.box <- apply(lon, 1, function(x) mean(x, na.rm = TRUE))
        lat.box <- apply(lat, 1, function(x) mean(x, na.rm = TRUE))
        lab.box <- unlist(poly$attrib['box_id'])
    }
    ## Limits
    if (isTRUE(OnlyPoly)){
        ylim <- range(lat, na.rm = TRUE)
        xlim <- range(lon, na.rm = TRUE)
        ylim <- c(ylim[1] * 1.001, ylim[2] * 0.999)
        xlim <- c(xlim[1] * 1.001, xlim[2] * 0.999)
    } else {
        if (is.null(xlim) | is.null(ylim)){
            cat('You need to set xlim and ylim for the map boundaries')
        }
    }
    ## Boundaries frame
    i.ylim <- c(ylim[1] * 1.001, ylim[2] * 0.999)
    i.xlim <- c(xlim[1] * 1.001, xlim[2] * 0.999)
    yrec   <- seq(from = i.ylim[1], to = i.ylim[2], by = sc)
    xrec   <- seq(from = i.xlim[1], to = i.xlim[2], by = sc)
    if (tail(xrec, 1) != i.xlim[2]) xrec <- c(xrec, i.xlim[2])
    if (tail(yrec, 1) != i.ylim[2]) yrec <- c(yrec, i.ylim[2])
    if (isTRUE(save)){
        ratio <- diff(ylim) / diff(xlim)
        if (ratio >=  1) {
            width <- 600
            height <- 600  * ratio
        } else {
            width <- 600 / ratio
            height <- 600
        }
        png(paste(name, '_plot.png', sep = ''), width = width, height = height)
    }
    ## Plot with coastal map
    if (!isTRUE(OnlyPoly)){
        ## Plot with land
        par(mar = c(4, 2, 4, 4) + 0.1, oma = c(3, 3, 3, 3))
        map('worldHires', ylim = i.ylim, xlim = i.xlim, fill = FALSE,
            xaxs = 'i', yaxs = 'i', boundary = F, col = 'white', ...)
        map('worldHires', ylim = ylim, xlim = xlim, col='black',
            fill=TRUE, xaxs='i', yaxs = 'i', add = T, ...)
        ## Polygons
        for( i in 1 : dim(lat)[1]){
            polygon(lon[i, ], lat[i, ], col  = ifelse(length(color) == 1, color, color[as.numeric(quart[i])]),  border = 'gray9')
        }
    } else {
        ## Only polygon
        if(PAR) par(mar=c(4, 4, 4, 4) + 0.1,oma=c(3,3,3,3), ...)
        plot(1, type = "n", axes = F, xlab = "", ylab = "", xlim = xlim, ylim = ylim,)

        ## Plot Polygons
        for( i in 1 : dim(lat)[1]){
            polygon(lon[i, ], lat[i, ], col  = ifelse(length(color) == 1, color, color[as.numeric(quart[i])]),  border = 'gray9')
        }

    }

    ## plot ID
    if(leg.color == 'box_id'){
        text(lon.box, lat.box, label = as.vector(lab.box))
    }

    ## Sum and inf
    for(i in 1 : (length(xrec) - 1)){
        rect(xrec[i], i.ylim[2], xrec[i + 1], i.ylim[2] + (diff(xlim) * 0.02), col = ifelse(i%%2 == 0, 'white', 'black'))
        rect(xrec[i], i.ylim[1], xrec[i + 1], i.ylim[1] - (diff(xlim) * 0.02), col = ifelse(i%%2 == 0, 'black', 'white'))
    }

    ## Right and left

    for(i in 1 : (length(yrec) - 1)){
        rect(i.xlim[2], yrec[i], i.xlim[2]  -  (diff(ylim) * 0.02), yrec[i + 1], col = ifelse(i%%2 == 0, 'white', 'black'))
        rect(i.xlim[1], yrec[i], i.xlim[1]  +  (diff(ylim) * 0.02), yrec[i + 1], col = ifelse(i%%2 == 0, 'black', 'white'))
    }


    xtick   <- xrec[seq(from = 1, to = length(xrec), by = 2)]
    xlabels <- parse(text = paste( -ceiling(xtick),'*degree*',floor(-(xtick - ceiling(xtick)) * 60 ),'*minute*W',  sep = ''))
    ytick   <- yrec[seq(from = 2, to = length(yrec), by = 2)]
    ylabels <- parse(text = paste( -ceiling(ytick),'*degree*',floor(-(ytick - ceiling(ytick)) * 60 ),'*minute*S',  sep = ''))
    axis(1, at = xtick, labels = xlabels, cex.axis = 1.5, lwd = 0)
    axis(2, at = ytick, labels = ylabels, cex.axis = 1.5, lwd = 0)
    mtext('Longitude', side = 1, outer = TRUE, cex = 1.5)
    mtext('Latitude', side = 2, outer = TRUE, line = 1.2, cex = 1.5)
    abline(h = yrec, v = xrec, lty = 3)

    if (!is.null(leg)){
        position <- seq(from = ylim[1], to = ylim[2], length.out = 30)
        arrows (leg, position[27],y1 = position[29] , lwd = 3)
        text(leg, position[26] + 0.005, 'N',  cex = 2, font = 2)
        map.scale(x = xlim[2] - 1.5, y = position[3], ratio = FALSE, relwidth = 0.1)
        if (!is.null(etiq)){
            legend(x = leg + .5, y = position[29],
                   title = ifelse(specie == FALSE, ifelse(div > 1, paste(leg.color, 'x', div), leg.color), paste(leg.color, '(%)')) ,
                   legend = etiq, pch = 15, col = if(specie == FALSE) { rev(color)} else {color.et},
                   bty = 'n', cex = 1.2 )
        }
    }
    if (isTRUE(save))  dev.off()
}

correct <- function(real.pos, x = 2, y = 1){
    ## Correct the location of the points
    ## x is the column for the X-position
    ## y is the column for the Y-position
    cat('\nclick the points that you want to change\n')
    pos <- locator()
    if(is.null(pos)) return(real.pos)
    repl <- vector()
    for(i in 1 : length(pos$x)){
        repl[[i]]   <- which.min(abs(real.pos[,  x] - pos$x[i]))
    }
    cat('\nYou want to change', length(pos$x), 'points, choose the new position\n')
    new.pos <- locator(n = length(pos$x))
    real.pos[repl, x] <- new.pos$x
    real.pos[repl, y] <- new.pos$y
    return(real.pos)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title von Bertalanfy growth equation
##' @param linf Asymptotic length or L infinity (if you work in fisheries)
##' @param k Individual growth rate
##' @param to Initial size of the organism and is defined as age at which the organisms would have had zero size
##' @param age Current age of the organism
##' @return Return the length of the organism at the 'age'
##' @author Demiurgo
vb <- function(linf, k, to, age){
    l <- linf * (1 - exp(-k * (age - to)))
    return(l)
}


dup.row <- function(Dat, Expand){
    ## Repeat the rows in a data frame based on values in a specific column ##
    ## Dat     : Is the data frame and expand in the number of the column
    ## Expand  : The column to expand the rows
    int.dup.row <- function(Dat, Expand) {
      # browser()
        rr <- data.frame(Dat)
        rr <- cbind(rr, n = rep(NA, Expand), row.names = NULL)
        nc <- ncol(rr)
        rr <- rr[ ,  - nc]
        return(rr)
    }
    ## Expande rows (as a list)
    exp.r <- lapply(1 : nrow(Dat), function(x) int.dup.row(Dat[x, ], Dat[x, Expand]))
    ## Bind the List
    n.Dat <- do.call(rbind, exp.r)
    return(n.Dat)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Estimation of structural and reserve weight
##' @param FG Vector with the name of the functional groups
##' @param weight vector witht the weigths
##' @param metric the default metric for the weight is 'mg' (miligrams), the other option are: 'g' grams; 'Kg' kilograms and 'gN' grams of nitrogen
##' @param wet Logical vector for the type weight, the default the is wet weight
##' @return a dataframe witht the fucntional groups and the structural and reserve weight
##' @author Demiurgo
weights <- function(FG, weight, metric = 'mg', wet = TRUE){
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~                      estimation of reserve and structural weight                   ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    if(metric == 'g')  weight <- weight * 1000
    if(metric == 'Kg') weight <- weight * 1000000

    ## Wet to dry weight
    if(wet) weight <- (weight / 20)

    ## From dry weight to nitrogen
    weight  <- weight / 5.7

    ## ratio between structural and reserve weight 1 : 2.65
    factor  <- weight / 3.65
    out     <- data.frame(Functional_Groups = FG, structural = factor, reserve = factor * 2.65)
    return(out)
}


##' \description{Function to estimate clearance based on biological parameters.}
##' @title Clearance
##' @param fg Vector with the name of the functional groups
##' @param speed Vector with the speed in mh-1. If the Speed is nor available ('NA') the speed would be calculated using the length at age class of the FG;
##' @param len Vector witht he length of the at Age class of the functional groups
##' @param height Vector of the ratio height compare to the length at age of the functinal group, if is not provided the stimation would be using 1/5 of the total length
##' @param ratio is the ratio between the height and the with of the individual. In other words, if the ratio is .5, the width is onlye the hal of the height for thar specie. If is not provided the value of width would be the same than height
##' @param time.l Vector proportion of time in a day invested for the specie searching for food.
##' @param max.speed Maximum speed reported for the functional group. That set the high boundary for the calculation of the speed
##' @param by.group Arrange the output for a easy manipulation
##' @return A vector witht the values of clereances for each functional group in mgNm3d-1
##' @author Demiurgo
clearance <- function(fg = NULL, speed, len, height = NA, ratio = NULL, time.l = NULL, max.speed = NULL, alfa = NULL, beta = NULL, by.group = TRUE){
    ## General assumption, If I dont have the speed I will use the mass to calculate the speed
    mass <- alometric((len * 100), alfa, beta) / 1000  ## in kilograms

    ## Assumption based on Sato et al 2007 Swimming speed
    speed <- ifelse(is.na(speed), mass ^ 0.27 * 3600, speed) # speed mh-1
    speed <- ifelse(is.na(max.speed), speed, ifelse(max.speed > speed, speed, max.speed))
    ## I assume the the height is at least 1/5 of the length

    height <- ifelse(is.na(height), len / 5, len * height)
    ## the with can be the same than the eight, but can have a deformation (flat fish)
    if(!is.null(ratio)){
        width  <- height * ratio
    } else {
        width  <- height
    }
    ## The time that the fish is looking for food
    if(is.null(time.l)) time.l <- 0.5
    ## all the time is in hours, the idea is to have it per day and per area (m3)
    out <- speed * height * width * 24 * time.l
    if(!is.null(fg)){
        out <- data.frame(FG = fg, Clearance = out)
        if(isTRUE(by.group)){
            temp.out <- split(out, out$FG)
            m        <- as.vector(table(fg))
            out      <- matrix(NA, nrow = max(m), ncol = length(m))
            for(i in 1 : length(temp.out)){
                out[, i] <- c(temp.out[[i]][, 2], rep(NA, 10 - length(temp.out[[i]][, 1])))
            }
            out           <- t(out)
            rownames(out) <- names(temp.out)
            colnames(out) <- c(1 : (ncol(out)))
        }
    }
    return(out)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Estimation of the allometric relationship
##' @param length Vector of values of length
##' @param alfa Parameter
##' @param beta Parameter
##' @return a vector with the weight
##' @author Demiurgo
alometric <- function(length, alfa, beta){
    weight <- alfa  * length ^ beta
    return(weight)
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Transform from conventional wet weight to gram of nitrogen
##' @param weight Vector of weight
##' @param metric The default metric for the wet weight is 'mg' (miligrams), the other options are:  'g' for grams and 'Kg' for Kilograms
##' @param wet Logical vector for the type weight, the default the is wet weight
##' @return a vector witht the weight in miligrams of nitrogen (mgN)
##' @author Demiurgo
weight2N <- function(weight, metric = 'mg', wet = TRUE, reserve = TRUE){
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~            transform weight to Ng        ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    if(metric == 'g')  weight <- weight * 1000
    if(metric == 'Kg') weight <- weight * 1000000
    ## Wet to dry weight
    if(wet) weight <- (weight / 20)
    ## From dry weight to nitrogen
    weight  <- weight / 5.7
    ## Structural weight
    weight  <- weight / 3.65
    ## If is reserve
    if(isTRUE(reserve)){
        weight  <- weight * 2.65
    }
    return(weight)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Estimation of groth rate 'mum' for the holling type II equation implemented in Atlantis for each functional groups
##' @param length Vecto of length in cm
##' @param weight Vector of weight in gr
##' @param metric The default metric for the wet weight is 'mg' (miligrams), the other options are:  'g' for grams and 'Kg' for Kilograms
##' @param wet Logical vector for the type weight, the default the is wet weight
##' @param spw.rate Spawning rate, vector of values that represent the increase in waeight due to reproduction
##' @param mature Logical, vector of TRUE or FALSE to descrive if the FG at that weight/length its mature
##' @param AgeClass Age classes by cohort
##' @return A vector with the values of mum for each length (mgNd-1)
##' @author Demiurgo
mum.f <- function(length, weight, metric  = 'mg', wet = TRUE, spw.rate, mature, AgeClass){
    we  <- c(0, weight)
    we  <- weight2N(we, metric = metric)
    len <- c(0, length)
    d.w <- vector('numeric')
    for(i in 2 : length(we)){
        wet        <- ifelse(isTRUE(mature[i - 1]), we[i] * spw.rate, we[i])
        d.w[i - 1] <- wet - we[i - 1]
    }
    d.l         <- diff(len)
    growth.rate <- d.w / d.l
    growth.rate <- growth.rate / (AgeClass * 365)
    return(growth.rate)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Weight to number for the data base based on the average weight
##' @param db Data base output
##' @param colW Number of the  column with the weight
##' @param colSp Number of the column with the species
##' @param colInd Identification of data type (weight or number)
##' @return
##' @author Demiurgo
WtoN <- function(db, colW, colSp, colInd){
    ubi  <- which(db[, colInd] %in% c('p','P'))
    temp <- ifelse(db[ubi, colSp] == 'OCT', round(db[ubi, colW] / 723),
            ifelse(db[ubi, colSp] == 'BRC', round(db[ubi, colW] / 746),
            ifelse(db[ubi, colSp] == 'ANG', round(db[ubi, colW] / 1013), NA)))
    nvec <- db[, colW]
    nvec[ubi] <- temp
    db$Count <- nvec
    return(db)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Griffiths model to estimate the temperature correction scalar
##' @param alfa Default  =  0.851. Species-specific coefficient
##' @param beta Defaul  =  1.066. Clobal temperature coefficient
##' @param tempC Default = 1. Global coefficient,  a multiplicative scalar
##' @param temp Enviromental temperature. is ambient water temperature
##' @param t.opt Is a species-specific optimum temperature
##' @param t.exp Defaul = 3. Exponent parameter
##' @param corr Default = 1000. Q10 correction parameter
##' @param vector Default = FALSE. Specified if the data of optimal temperature is a vector or not.
##' @return
##' @author Demiurgo
temp.corr <- function(alfa = 0.851, beta = 1.066, tempC = 1, temp, t.opt, t.exp = 3, corr, vector = FALSE){
    if(vector){
        out[fg, ] <- log(2) * alfa[fg] * (beta ^ temp) * exp( (- tempC * (abs(temp - t.opt[fg])) ^ t.exp) / corr[fg])
    } else {
        out <- array(NA, c(length(alfa), length(temp)))
        for(fg in 1 :  length(alfa)){
            out[fg, ] <- log(2) * alfa[fg] * (beta ^ temp) * exp( (- tempC * (abs(temp - t.opt[fg])) ^ t.exp) / corr[fg])
        }
    }
    return(out)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param x Vector of values
##' @param lag Lag for the average
##' @return a vector with the moving average
##' @author Demiurgo
movavg <- function(x, lag){
    dimen.i <- length(x)
    cero        <- rep(0, floor(lag / 2))
    x           <- c(cero, x, cero)
    #if(lag %% 2 > 1) c(x, 0)
    x[is.na(x)] <- 0
    dimen       <- length(x)
    out         <- vector('numeric', dimen)
    count       <- 1
    for (i in ceiling(lag / 2) : (dimen - (lag - 1) / 2)) {
        serie    <- count : (lag + count - 1)
        out[i]   <- mean(x[serie], na.rm = TRUE)
        count    <- count+1
    }

    out <- out[(floor(lag / 2) + 1) : (floor(lag / 2) + dimen.i)]
    out[is.na(out)] <- 0
    return(out)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Abunance initial condition
##' @param data Abundance observed by Atlantis box
##' @param in.bios Initial Biomass
##' @param boxes Boxes information from BGM files
##' @param groups Functional groups information from the csv file
##' @param lfd Lenght frequency by cohort
##' @param cover.d Data frame witht the proportion of cover by box
##' @return Matrix with number or biomass (N) by Functional group
##' @author Demiurgo
init.number <- function(data, in.bios, boxes, groups, lfd, cover.d, m.weight){
    ## if the Biomass or number is empty the code will use the distribution assuming that
    ## the proportion by box is given or the proportion of cover.
    area  <- boxes[order(boxes$box_id), ]$area  # m2
    depth <- -boxes[order(boxes$box_id), ]$botz
    depth[depth < 0] <- NA
    Volumen <- area * depth
    ## noBio <- which(is.na(in.bios$BioNumber))
    ## idx   <- which(colnames(data) %in% in.bios$FG[noBio])
    ## for(i in 1 : length(noBio)){
    ##     in.bios$BioNumber[noBio[i]] <- sum(data[, idx[i]])
    ## }
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~             BIOMASSS POOLS ESTIMATION          ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ##~ Biomass pools,  stimation of initial abundance by box
    b.pools <- groups$Code[which(groups$NumCohorts == 1)]
    ##~ which functional groups are biomass pool
    i.pools <- which(colnames(data) %in% b.pools)
    ##~ Proportion by box
    dat.tmp <- data[, i.pools]
    prp.box <- dat.tmp * NA
    for( i in 1 : ncol(dat.tmp)){
        prp.box[, i] <- dat.tmp[, i] / sum(dat.tmp[, i])
    }
    ##~ transfor from number to Biomass using the mean weight
    n.pools                    <- which(in.bios$FG %in% b.pools & in.bios$type == 'N')
    in.bios$BioNumber[n.pools] <- in.bios$BioNumber[n.pools] * in.bios$weight[n.pools]
    ##~ Transform the weight in Nitrogen for the Biomass
    only.pools <- which(in.bios$FG %in% b.pools)
    N.pool     <-  with(in.bios, data.frame(FG  = as.character(FG[only.pools]),
                                            N  = weight2N(BioNumber[only.pools], 'g', reserve = FALSE)))
    ##~ match the same order of columns
    len.ord <- sum(colnames(prp.box) %in% N.pool$FG)
    #ord.pool   <- which(colnames(prp.box) %in% N.pool$FG)
    ##~ data in a dataframe
    N          <- array(NA, dim = c(len.ord, nrow(prp.box)))
    namN       <- vector('character')
    for( i in 1 : len.ord){
        ord.prp <- which(colnames(prp.box) == N.pool$FG[i])
        prop   <- prp.box[, ord.prp] / sum(prp.box[, ord.prp])
        N[i, ] <- as.numeric(prop * N.pool$N[i])
        if(N.pool$FG[i] %in% c('OCT', 'SQD')){
            N[i, ] <- N[i, ] / (Volumen)
        } else if(!(N.pool$FG[i] %in% c('LPH', 'SPH','SZO', 'MZO', 'LZO'))){
            N[i, ] <- N[i, ] / area
        }
        N[i, is.na(N[i, ])] <- 0

        pos.n  <- with(groups, which(Code %in% N.pool$FG[i]))
        namN   <- c(namN, paste(groups$Name[pos.n], '_N', sep = ''))
    }
    ## ##~ added the name of each row
    row.names(N) <- namN

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~               Adding the cover information           ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    cov.pos <- which(groups$Code %in% colnames(cover.d))
    cov.nam <- c(paste(groups$Name[cov.pos], '_Cover', sep = ''), 'reef', 'soft')
    colnames(cover.d)[2 : ncol(cover.d)] <- cov.nam
    cover.d$Rugosity <- rowSums(cover.d[, - c(1, 5)])
    cov.dat <- t(cover.d)[ - 1, ]

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~         NUMBER ESTIMATION      ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## I will Use the same abundance for each class. I can provided the age structure but I',  azy to do it
    ## Convert Biomass to number
    numb    <- groups$Code[which(groups$NumCohorts > 1)]
    cohorts <- groups$NumCohorts[which(groups$NumCohorts > 1)]
    ##~ which functional groups are biomass pool
    o.numb  <- which(in.bios$FG %in% numb)
    ##~ change from Biomass to numbers
    in.bios$fg.Nums <- in.bios$BioNumber
    in.bios$RN      <- in.bios$SN <- NA
    ## numbers
    ##~ Transform into number by cohort and area
    nam <- vector('character')
    for(i in o.numb){
        loc         <- which(lfd$FG %in% in.bios$FG[i])
        l.prop      <- which(colnames(data) %in% in.bios$FG[i])
        prop.box    <- data[, l.prop] / sum(data[, l.prop])
        ## By cohort
        N.at.Cohort  <- as.numeric(lfd[loc, 2 :  ncol(lfd)] * in.bios$fg.Nums[i])
        if(in.bios$type[i] == 'B'){ # from Biomass at age to numbers
            p.w  <- which(m.weight$FG %in% in.bios$FG[i])
            N.at.Cohort <-  as.numeric(N.at.Cohort / m.weight[p.w, 2 : ncol(m.weight)])
            N.at.Cohort[is.na(N.at.Cohort)] <- 0
        }
        #RN.at.Cohort <- as.numeric(lfd[loc, 2 :  ncol(lfd)] * in.bios$RN[i])
        #SN.at.Cohort <- as.numeric(lfd[loc, 2 :  ncol(lfd)] * in.bios$SN[i])
        ## By Area
        tmp.N  <- array(NA, dim = c(groups$NumCohorts[which(groups$Code %in% in.bios$FG[i])], length(prop.box)))
        #tmp.SN <- tmp.RN <- tmp.N
        for(j in 1 : sum(N.at.Cohort > 0)){
            #prop        <- prp.box[,  l.prop] / sum(prp.box[,  l.prop])
            tmp.N[j, ]  <- prop.box * N.at.Cohort[j]
            #tmp.RN[j, ] <- prop * RN.at.Cohort[j]
            #tmp.SN[j, ] <- prop * SN.at.Cohort[j]
            pos.n       <- with(groups, which(Code %in% in.bios$FG[i]))
            nam         <- c(nam, paste(groups$Name[pos.n], j, '_Nums', sep = ''))
        }
        if(i == o.numb[1]){
            out    <- tmp.N
        } else {
            out <- rbind(out, tmp.N)
        }
    }
    nrow(out)
    length(nam)
    rownames(out)    <- nam
    output           <- rbind(N, out, cov.dat)
    return(output)
}

read.dat <- function(file){
    options(warn=-1)
    ##-----------------------------------------------------------#
    ##  Function to read data file in a ".dat" files             #
    ##        Creador:  Steve Martell  Mod: Javier Porobic       #
    ##-----------------------------------------------------------#
    ifile <- scan(file, what = "character", flush = TRUE, blank.lines.skip = FALSE, quiet = TRUE)
    idx   <- sapply(as.double(ifile), is.na)
    vnam  <- ifile[idx]	                    #list names
    nv    <- length(vnam)	                    #number of objects
    A     <- list()
    ir    <- 0
    for(i in 1 : nv)
    {
        ir <- match(vnam[i], ifile)
        if(i != nv) irr <- match(vnam[i + 1], ifile) else irr = length(ifile) + 1 #next row
        dum <- NA
        if(irr-ir==2) dum <- as.double(scan(file, skip = ir, nlines = 1, quiet = TRUE, what = ""))
        if(irr-ir>2)  dum <- as.matrix(read.table(file, skip = ir, nrow = irr-ir-1, fill = TRUE))
        if(is.numeric(dum))#Logical test to ensure dealing with numbers
        {
            A[[ vnam[i ]]] <- dum
        }
    }
    return(A)
}

##' @title Parameter file reader
##' @param text Biological parametar file for Atlatnis
##' @param pattern Text that you are looking
##' @param FG Name of the functional groups
##' @param Vector Logic argument, if the data is on vectors or not
##' @return A matrix with the values from the .prm file
##' @author Demiurgo
text2num <- function(text, pattern, FG = NULL, Vector = FALSE){
    if(!isTRUE(Vector)){
        text <- text[grep(pattern = pattern, text)]
        txt  <- gsub(pattern = '[[:space:]]+' ,  '|',  text)
        col1 <- col2 <- vector()
        for( i in 1 : length(txt)){
            tmp     <- unlist(strsplit(txt[i], split = '|', fixed = TRUE))
            tmp2    <- unlist(strsplit(tmp[1], split = '_'))
            if(FG[1] == 'look') {
                col1[i] <- tmp2[1]
            } else {
                id.co   <- which(tmp2 %in% FG )
                col1[i] <- tmp2[id.co]
            }
            col2[i] <- as.numeric(tmp[2])
        }
        if(is.null(FG)) col1 <- rep('FG', length(col2))
        return(data.frame(FG = col1, Value = col2))
    } else {
        l.pat <- grep(pattern = pattern, text)
        nam   <- gsub(pattern = '[ ]+' ,  '|',  text[l.pat])
        fg    <- vector()
        pos   <- 1
        for( i in 1 : length(nam)){
            tmp     <- unlist(strsplit(nam[i], split = '|', fixed = TRUE))
            if(tmp[1] %in% c('#','##', '###')) next  ## check this part!!
            fg[pos] <- tmp[1]
            if(pos == 1) {
                pp.mat <- matrix(as.numeric(unlist(strsplit(text[l.pat[i] + 1], split = ' ', fixed = TRUE))), nrow = 1)
                pos    <- pos + 1
            } else {
                pp.tmp <- matrix(as.numeric(unlist(strsplit(text[l.pat[i] + 1], split = ' ', fixed = TRUE))), nrow = 1)
                if(ncol(pp.mat) != ncol(pp.tmp)) stop('\nError: The pPrey vector for', tmp[1], ' has ', ncol(pp.tmp))
                pp.mat <- rbind(pp.mat, pp.tmp)
                pos    <- pos + 1
            }
        }
        if(all(is.na(pp.mat[, 1]))) pp.mat <- pp.mat[, - 1]
        row.names(pp.mat)                  <- fg
        return(pp.mat)
    }
}


## ~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~       Diet Checker   ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~ ##

diet.chk <- function(group.file, prm.file, nc.file){
    library(reshape)
    library(ggplot2)
    ## Reading files
    groups.csv <- read.csv(grp.file)
    prm        <- readLines(prm.file, warn = FALSE)

    ## Gape size and adult and young age
    KLP                     <- text2num(prm, 'KLP', FG = as.character(groups.csv$Code))
    KUP                     <- text2num(prm, 'KUP',  FG = as.character(groups.csv$Code))
    age                     <- text2num(prm, '_age_mat', FG = as.character(groups.csv$Code))
    Gape                    <- data.frame(FG = KLP$FG, KLP = KLP$Value, KUP = KUP$Value, Age.Adult = NA)
    pos.Age                 <- which(Gape$FG %in% age$FG)
    Gape$Age.Adult[pos.Age] <- age$Value
    Gape$Age.Young          <- Gape$Age.Adult - 1
    ## availability matrix
    Ava.mat            <- text2num(prm, 'pPREY', Vector=TRUE)
    colnames(Ava.mat)  <- c(as.character(groups.csv$Code), 'DLsed', 'DRsed', 'DCsed')
    Is.off             <- which(groups.csv$IsTurnedOn == 0)
    ## getting the RN and SN from the NC file
    nc.out <- nc_open(nc.file)
    FG     <- as.character(groups.csv$Name)
    Biom.N <- array(data = NA, dim = c(length(FG), max(groups.csv$NumCohorts)))
    Struct <- Biom.N
    for(code in 1 : length(FG)){
        if(groups.csv$NumCohorts[code] == 1 && groups.csv$IsTurnedOn[code] == 1){
            N.tot <- ncvar_get(nc.out, paste(FG[code], "_N", sep = ""))
            if(all(is.na(N.tot)) || all(N.tot == 0) || sum(N.tot, na.rm = TRUE) == 0){
                Biom.N[code, 1] <- ncatt_get(nc.out, varid = paste(FG[code], "_N", sep = ""), attname = "_FillValue")$value
            } else {
                if(length(dim(N.tot)) > 3){
                    N.tot <- N.tot[, , 1]
                }
                Biom.N[code, 1] <- sum(N.tot, na.rm = TRUE)
            }
        } else if(groups.csv$NumCohorts[code] > 1 && groups.csv$IsTurnedOn[code] == 1) {
            for(cohort in 1 : groups.csv$NumCohorts[code]){
                StructN <- ncvar_get(nc.out, paste(FG[code], as.character(cohort), "_StructN", sep = ""))
                ReservN <- ncvar_get(nc.out, paste(FG[code], as.character(cohort), "_ResN", sep = ""))
                Numb    <- ncvar_get(nc.out, paste(FG[code], as.character(cohort), "_Nums", sep = ""))
                if(length(dim(ReservN)) > 3){
                    StructN <- StructN[, , 1]
                    ReservN <- ReservN[, , 1]
                }
                if(length(dim(Numb)) > 3){
                    Numb    <- Numb[, , 1]
                }
                Biom.N[code, cohort] <- (max(colSums(StructN,  na.rm = TRUE), na.rm = TRUE)  +
                                         max(colSums(ReservN,  na.rm = TRUE), na.rm = TRUE)) *
                    sum(Numb, na.rm = TRUE)
                Struct[code, cohort] <- (max(colSums(ReservN,  na.rm = TRUE), na.rm = TRUE)) #* sum(Numb, na.rm = TRUE)
            }
        }
    }
    nc_close(nc.out)
    row.names(Biom.N) <- as.character(groups.csv$Code)
    row.names(Struct) <- as.character(groups.csv$Code)
    if(length(Is.off) > 0){
        ## Remove groups that are not On in the model
        Struct            <- Struct[ - Is.off, ]
        Biom.N            <- Biom.N[ - Is.off, ]
    }
    Gape$Age.Young    <- ifelse(Gape$Age.Young == 0,  1, Gape$Age.Young)
    ## Pre-Calculations
    ## Be sure that the FG have the same order
    Biom.N        <- Biom.N[order(row.names(Biom.N)), ]
    Struct        <- Struct[order(row.names(Struct)), ]
    Gape          <- Gape[order(Gape$FG), ]
    G.pos         <- which(row.names(Struct) %in% Gape$FG)
    Gape$juv.Min  <- Struct[G.pos, 1] * Gape$KLP
    for( i in 1 : length(G.pos)){
        Gape$adult.Min[i]  <- Struct[G.pos[i], Gape$Age.Adult[i]] * Gape$KLP[i]
        Gape$adult.Max[i]  <- Struct[G.pos[i], length(sum(!is.na(Struct[G.pos[i], ])))] * Gape$KLP[i]
        Gape$juv.Max[i]    <- Struct[G.pos[i], Gape$Age.Young[i]] * Gape$KUP[i]
        Gape$JminS[i]      <- Struct[G.pos[i], 1]
        Gape$AminS[i]      <- Struct[G.pos[i], Gape$Age.Adult[i]]
        Gape$JmaxS[i]      <- Struct[G.pos[i], Gape$Age.Young[i]]
        Gape$AmaxS[i]      <- Struct[G.pos[i], length(sum(!is.na(Struct[G.pos[i], ])))]
    }
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~        Overlap-Matrix    ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    Over.mat <- Ava.mat * 0
    Prey     <- colnames(Ava.mat)
    Pred     <- row.names(Ava.mat)
    cp <- 2
    for( py in 1: length(Prey)){
        for( pd in 1: length(Pred)){
            c.pred      <- unlist(strsplit(Pred[pd],'pPREY'))[2]
            predator    <- gsub(pattern = "[[:digit:]]+", '\\1', c.pred)
            a.pred.prey <- as.numeric(unlist(strsplit(c.pred, predator)))
            pry.loc     <- which(Gape$FG %in% Prey[py])
            prd.loc     <- which(Gape$FG %in% predator)
            if(length(pry.loc) == 0 || is.na(a.pred.prey)){
                Over.mat [pd, py] <- 1
            } else {
                if(a.pred.prey[1] == 1){
                    ## Young Predator
                    if(a.pred.prey[2] == 1){
                        ## Young Prey
                        Over.mat [pd, py]  <- ifelse(Gape$JminS[pry.loc] >= Gape$juv.Min[prd.loc],
                                              ifelse(Gape$JminS[pry.loc] <= Gape$juv.Max[prd.loc], 1, 0),
                                              ifelse(Gape$JmaxS[pry.loc] >= Gape$juv.Min[prd.loc],
                                              ifelse(Gape$JmaxS[pry.loc] <= Gape$juv.Max[prd.loc], 1, 1), 0))
                        if(is.na(Over.mat[pd, py])) Over.mat[pd, py] <- 1
                    } else {
                        ## Adult Prey
                        Over.mat [pd, py]  <- ifelse(Gape$AminS[pry.loc] >= Gape$juv.Min[prd.loc],
                                              ifelse(Gape$AminS[pry.loc] <= Gape$juv.Max[prd.loc], 1, 0),
                                              ifelse(Gape$AmaxS[pry.loc] >= Gape$juv.Min[prd.loc],
                                              ifelse(Gape$AmaxS[pry.loc] <= Gape$juv.Max[prd.loc], 1, 1), 0))
                        if(is.na(Over.mat[pd, py])) Over.mat[pd, py] <- 1
                    }
                } else {
                    ## Adult Predator
                    if(a.pred.prey[2] == 1){
                        ## Young Prey
                        Over.mat [pd, py]  <- ifelse(Gape$JminS[pry.loc] >= Gape$adult.Min[prd.loc],
                                              ifelse(Gape$JminS[pry.loc] <= Gape$adult.Max[prd.loc], 1, 0),
                                              ifelse(Gape$JmaxS[pry.loc] >= Gape$adult.Min[prd.loc],
                                              ifelse(Gape$JmaxS[pry.loc] <= Gape$adult.Max[prd.loc], 1, 1), 0))
                        if(is.na(Over.mat[pd, py])) Over.mat[pd, py] <- 1
                    } else {
                        ## Adult Prey
                        Over.mat [pd, py]  <- ifelse(Gape$AminS[pry.loc] >= Gape$adult.Min[prd.loc],
                                              ifelse(Gape$AminS[pry.loc] <= Gape$adult.Max[prd.loc], 1, 0),
                                              ifelse(Gape$AmaxS[pry.loc] >= Gape$adult.Min[prd.loc],
                                              ifelse(Gape$AmaxS[pry.loc] <= Gape$adult.Max[prd.loc], 1, 1), 0))
                        if(is.na(Over.mat[pd, py])) Over.mat[pd, py] <- 1
                    }
                }
            }
        }
    }
    ## total biomasss by Juv and Adults
    Biom.N  <- Biom.N[order(row.names(Biom.N)), ]
    fg      <- row.names(Biom.N)
    bio.juv <- bio.adl <- matrix(NA, ncol = 2, nrow = nrow(Biom.N))
    for( i in 1 : nrow(Biom.N)){
        l.age <- which(age$FG == fg[i])
        if(length(l.age) !=  0){
            bio.juv[i, ] <- c(fg[i], sum(Biom.N[i, 1 : (age$Value[l.age] - 1)], na.rm = TRUE))
            bio.adl[i, ] <- c(fg[i], sum(Biom.N[i, age$Value[l.age] : ncol(Biom.N)], na.rm = TRUE))
        } else {
            ## For Biomass pool,  only adult
            bio.juv[i, ] <- c(fg[i], sum(Biom.N[i, 1], na.rm = TRUE))
            bio.adl[i, ] <- c(fg[i], sum(Biom.N[i, 1], na.rm = TRUE))
        }
    }
    ## Sort based on the order of the prey in the Availavility matrix   ##
    or.prey <- match(colnames(Over.mat), bio.juv[, 1])
    bio.juv <- bio.juv[or.prey, ]
    bio.adl <- bio.adl[or.prey, ]
    Biom.N
    ## Real feeding
    is.feeding <- Over.mat * Ava.mat
    real.feed  <- is.feeding * NA
    pred       <- row.names(Over.mat)
    for( pd in 1 : nrow(is.feeding)){
        ## Getting the number of biomass needed by each functional group
        c.pred      <- unlist(strsplit(pred[pd],'pPREY'))[2]
        predator    <- gsub(pattern = "[[:digit:]]+", '\\1', c.pred)
        a.pred.prey <- as.numeric(unlist(strsplit(c.pred, predator)))
        pry.loc     <- which(bio.adl[, 1] %in% predator)
        if(length(a.pred.prey) == 0 || is.na(a.pred.prey)) a.pred.prey[2] <- 2
        ## Young Predator
        if(a.pred.prey[2] == 1){
            ## Young Prey
            real.feed[pd, ] <- (is.feeding[pd, ] * as.numeric(bio.juv[, 2]))
        } else {
            ## Adult Prey
            real.feed[pd, ] <- (is.feeding[pd, ] * as.numeric(bio.adl[, 2]))
        }
    }
    ## Plot output
    real.vec.pprey <- melt(t(log(real.feed)))
    ggplot(data = real.vec.pprey,
           aes(x = X1, y = X2, fill = value)) + geom_tile() +
        scale_fill_gradient(limits=c(0, max(real.vec.pprey$value, na.rm = TRUE)), name = 'Predation value', low="white", high="red", na.value = 'white')  +
        theme(panel.background = element_blank()) + labs(x = 'Prey', y = 'Predator')+scale_x_discrete(position = "top")
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Clean data for spatial distribution of species
##' @param ##' @param nomb name of the colums or functional groups to be clean
##' @return a clean poly structure
##' @author Demiurgo
clean <- function(poly, nomb){
    cir.se <- c(4, 1 : 4, 1)
    for(fg in nomb){
        dbs  <- cbind(poly[[1]]$attrib[, fg], poly[[2]]$attrib[, fg],
                      poly[[3]]$attrib[, fg], poly[[4]]$attrib[, fg])
        rsea <- apply(dbs, 2, function(x) all(is.na(x)))
        if(sum(rsea) > 2){ ## At least 3 seasons need to have data
            mean.s <- apply(dbs, 1, function(x) sum(x, na.rm = TRUE))
            for(sea in 1 : 4){
                poly[[sea]]$attrib[, fg] <- mean.s
            }
        } else if(sum(rsea)  == 1){
            na.s                      <- which(rsea)
            poly[[na.s]]$attrib[, fg] <- rowMeans(dbs[, cir.se[c(na.s, na.s + 2)]], na.rm = TRUE)
        }
        ## reading again the fixed distribution
        dbs  <- cbind(poly[[1]]$attrib[, fg], poly[[2]]$attrib[, fg],
                      poly[[3]]$attrib[, fg], poly[[4]]$attrib[, fg])
        for(box in 1 : nrow(dbs)){
            if(all(is.na(dbs[box, ]))) next()
            na.box <- which(is.na(dbs[box, ]))
            if(length(na.box) == 1){
                dbs[box, na.box] <- mean(dbs[box, cir.se[c(na.box, na.box + 2)]], na.rm = TRUE)
            } else {
                dbs[box, ] <- mean(dbs[box, ], na.rm = TRUE)
            }
        }
        for(sea in 1 : 4){
            poly[[sea]]$attrib[, fg] <- dbs[, sea]
        }
    }
    return(poly)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Selectivity of Atlatnis
##' @param len Length
##' @param selb Selectividad parameter
##' @param lsm 50% average
##' @return the proportion
##' @author Demiurgo
selec <- function(len, selb, lsm){
    psel <- 1 / (1 + exp( - selb * (len - lsm)))
    return(psel)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param var.name Name of the Variable
##' @param data vector with time and
##' @param long.name
##' @param units
##' @param time.unit
##' @return
##' @author Demiurgo
write.ts <- function(var.name, ext.nam = NULL, data, long.name, units, time.unit = 'seconds since 1900-01-01 00:00:0.0', n.var = 1){
    if(n.var > 1){
        sink(paste0(ext.nam, '.ts'))
    } else {
        sink(paste0(var.name[1], ext.nam, '.ts'))
    }
    cat('# Time serie of', long.name[1], 'Created using the write.ts function\n')
    cat('# -------------------------------------------------------------------\n')
    cat('#\n')
    cat(paste0('## COLUMNS ', 1 + n.var, '\n'))
    cat('##\n')
    cat('## COLUMN1.name Time\n')
    cat('## COLUMN1.long_name Time\n')
    cat('## COLUMN1.units', time.unit, '\n')
    cat('## COLUMN1.missing_value 0\n')
    cat('##\n')
    for(i in 1 : n.var){
        val.var = 1 + i
        cat(paste0('## COLUMN', val.var, '.name'), var.name[i],'\n')
        cat(paste0('## COLUMN', val.var, '.long_name'), long.name[i],'\n')
        cat(paste0('## COLUMN', val.var, '.units'), units,'\n')
        cat(paste0('## COLUMN', val.var, '.missing_value'), '0\n')
        cat('##\n')
    }
    for(nr in 1 : nrow(data)){
        cat(data[nr, ],'\n')
    }
    sink()
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Summary of catch and CPU for bycatch
##' @param list list with the information of Catch and trips per season/Island and Year
##' @return data.frame with catch per trip per season per year
##' @author Demiurgo
sumry <- function(list){
    list  <- data.frame(list)
    trips <- data.frame(table(list$Season))
    out   <- NULL
    for(i in trips[, 1]){
        list2 <- list[which(list$Season == i),  ]
        tmp   <- unlist(tapply(list2$Total_unit, list2$FG, sum, na.rm = TRUE))
        nam   <- sort(unique(list2$FG))
        tmp   <- data.frame(Season = i, FG = nam, total_c= tmp, trips = trips[which(trips[, 1] == i), 2])
        out   <- rbind(out, tmp)
    }
    out$cperT <- with(out, total_c /trips)
    out$Year  <- unique(list$Year)
    out$Island <- unique(list$Island)
    return(out)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Average by years
##' @param dat data frame
##' @param FGs Names of the functional groups
##' @param start Year that you want to start to calculate the average
##' @param n.years Number of years
##' @param Time Vector with the dates
##' @param years Number of years that
##' @return A dataframe witht he average for the las 10 years
##' @author Demiurgo
avg.fg <- function(dat, FGs, start = 2000, n.years = 10, Time){
    col.n  <- which(colnames(dat) %in% FGs)
    dat2   <- dat[, col.n]
    out    <- rowsum(dat2, format(Time, '%Y'))
    tframe <- seq(from = which(row.names(out) == start), length = n.years)
    out    <- colMeans(dat2[tframe, ])
    return(out)
}

##' @title Exponential growth model
##' @param t time
##' @param No initial abundance
##' @param r growth rate
##' @return Estiamted abundance
##' @author Demiurgo
mod.exp <- function(t, No, r){
    yr1 <- min(t, na.rm = TRUE)
    yr2 <- max(t, na.rm = TRUE)
    nobs <- length(t)
    nyrs <- yr2 - yr1 + 1
    years <- yr1 : yr2
    yobs <- rep(0, nobs)
    x <- rep(0, nyrs)
    x[yr1] <- No
    for(i in yr1 : (yr2 - 1)){
         x[i + 1] = x[i] + r * x[i]
         }
    for(i in 1:nobs) yobs[i] = x[t[i]]
    return(yobs)
   }

##' @title Densodependent growth model
##' @param t time
##' @param No initial abundance
##' @param r growth rate
##' @param K carrying capacity
##' @return Estimated abundance
##' @author Demiurgo
mod.denso <- function(t, No, r, K){
    yr1 <- min(t, na.rm = TRUE)
    yr2 <- max(t, na.rm = TRUE)
    nobs <- length(t)
    nyrs <- yr2 - yr1 + 1
    years <- yr1 : yr2
    yobs <- rep(0, nobs)
    x <- rep(0, nyrs)
    x[yr1] <- No
    for(i in yr1 : (yr2 - 1)){
         x[i + 1] = x[i] + r * x[i] * (1 - (x[i] / K))
         }
     for(i in 1:nobs) yobs[i]= x[t[i]]
    return(yobs)
   }

##' @title Read Atlantis initial condition file
##' @param nc.file netCDF initial confition file
##' @param just.Sum This fucntion only calculate the Sum by box. If you want the mean this should be false
##' @param file.out Name of the output file
##' @return A csv file with the sum by box of the values in the initial condition
##' @author Demiurgo
transform.ini <- function(nc.file = NULL, just.Sum = TRUE, file.out = 'IniCond'){
    nc    <- nc_open(nc.file)
    var2d <- NULL
    var   <- names(nc[['var']])
    ## Ugly loop
    for(v in 1 : length(var)){
        var.tmp <- ncvar_get(nc, variables[v])
        if(length(dim(var.tmp)) >= 2){ ## Avoiding problem with the variables like "numlayers" that only have one dimention
            if(just.Sum){
                temp2d  <- colSums(var.tmp, na.rm = TRUE)
            } else{
                temp2d  <- colMeans(var.tmp, na.rm = TRUE)
            }
            var2d   <- cbind(var2d,  temp2d)
        } else {
            var2d   <- cbind(var2d,  var.tmp)
        }
    }
    colnames(var2d) <- var
    if(just.Sum){
        file <- paste0(file.out, '_Sums.csv')
    } else {
        file <- paste0(file.out, '_Means.csv')
    }
    write.table( var2d,  file,  row.names = FALSE)
    cat(paste('file', file, 'has been created'))
}




##' @title Mean growth
##' @param linf Maximum size in average
##' @param l.cur Current size
##' @param k.g growth
##' @return The growth at a given size
##' @author Demiurgo
mgrowth <- function(linf, l.cur, k.g){
    growth <- (linf - l.cur)*(1-exp(-k.g))
    return(growth)
}

##' @title Estimation of selectivity at size
##' @param l50 size at which 50% is selected
##' @param l95 size at which 95% is selected
##' @param size Vector of size
##' @return A vector with the selectivity
##' @author Demiurgo
selectivity.f <- function(l50, l95, size){
    vect <- 1 / (1 + exp(( - log(19) * (size - l50)) / (l95 - l50)))
    return(vect)
}




##' @title Calculate the different weigt (estructural reserve) and number of indivduals
##' @param nc.out Netcdf out file. this is the traditional .nc file from atlantis
##' @param grp groups csv file
##' @param FG Functional group (selected)
##' @param By option to aggregate the time series in just one or leave the values by cohort
##' @param box.info Internal function
##' @return List witht he weight (reserve and structural), total biomass (tons) and number
##' @author Demiurgo
nitro.weight <- function(nc.out, grp, FG, By = 'Total', box.info, mg2t, x.cn){
    ## Age classes
    pos.fg <- which(grp$Code == FG)
    Bio <- Num <- SN  <- RN  <- list()
    if(grp[pos.fg, 'NumCohorts'] > 1 & grp[pos.fg, 'GroupType'] != 'PWN'){
        n.coh <- grp[pos.fg, 'NumCohorts']
        for(coh in 1 : n.coh){
            name.fg <- paste0(grp$Name[pos.fg], coh)
            resN    <- ncvar_get(nc.out, paste0(name.fg, '_ResN'))
            strN    <- ncvar_get(nc.out, paste0(name.fg, '_StructN'))
            nums    <- ncvar_get(nc.out, paste0(name.fg, '_Nums'))
            b.coh   <- (resN  + strN)  * nums * mg2t * x.cn
            if(By %in% c('Total', 'Cohort')){
                nums    <- apply(nums, 3, sum, na.rm = TRUE)
                b.coh   <- apply(b.coh, 3, sum, na.rm = TRUE)
                strN    <- apply(strN, 3, mean, na.rm = TRUE)
                resN    <- apply(resN, 3, mean, na.rm = TRUE)
            }
                RN[[coh]]  <- resN
                SN[[coh]]  <- strN
                Bio[[coh]] <- b.coh
                Num[[coh]] <- nums
        }
        if(By ==  'Total'){
            RN  <- rowSums(matrix(unlist(RN), ncol = n.coh))
            SN  <- rowSums(matrix(unlist(SN), ncol = n.coh))
            Bio <- rowSums(matrix(unlist(Bio), ncol = n.coh))
            Num <- rowSums(matrix(unlist(Num), ncol = n.coh))
        }
        type <- 'AgeClass'
    } else if (grp[pos.fg, 'NumCohorts'] == 1){ ## Biomass pool
        name.fg <- paste0(grp$Name[pos.fg], '_N')
        biom    <- ncvar_get(nc.out, name.fg)
        if(length(dim(biom)) == 3){
            if(grp$GroupType[pos.fg] == 'LG_INF'){
                biom <- apply(biom, 3, '*', box.info$VolInf)
            }else{
                biom <- apply(biom, 3, '*', box.info$Vol)
            }
            if(By == 'Total'){
                biom <- apply(biom, 2, sum, na.rm = TRUE)
            }
        } else {
            biom <- apply(biom, 2, function(x) x * box.info$info$Area)
            if(By == 'Total'){
                biom <- apply(biom, 2, sum, na.rm = TRUE)
            }
        }
        Bio <- biom
        type <- 'BioPool'
    } else if(grp[pos.fg, 'NumCohorts'] > 1 & grp[pos.fg, 'GroupType']  == 'PWN'){
        ## Some model use Agestructured biomass pools
        for(coh in 1 : pwn.grp[pwn, 'NumCohorts']){
            name.fg <- paste0(grp$Name[pos.fg],'_N', coh)
            biom    <- ncvar_get(nc.out, name.fg)
            biom    <- apply(biom, 3, '*', box.info$Vol)
            if(By == 'Total'){
                biom <- apply(biom, 2, sum, na.rm = TRUE)
            }
            Bio[[coh]] <- biom
        }
        type <- 'AgeBioPool'
    }
    return <- list(Biomass    = Bio,
                   Numbers    = Num,
                   Structural = SN,
                   Reserve    = RN,
                   Type       = type)
}



##' @title Box information
##' @param bgm.file BGM file,  Atlantis input
##' @param depths Cummulative depths (de max depth of each layer)
##' @return a dataframe with information (box-id, area (m), )
##' @author Demiurgo
boxes.prop <- function(bgm.file, depths){
    bgm       <- readLines(bgm.file, warn = FALSE)
    boxes     <- text2num(bgm, 'nbox', FG = 'look')
    out       <- NULL
    depths    <- depths[ - which(depths == 0)]
    max.nlyrs <- length(depths)              ## maximum number of water layers
    vol       <- array(NA, dim = c(boxes$Value, length(depths)))
    for(b in 1 : boxes$Value){
        area      <- text2num(bgm, paste0('box', b - 1,'.area'), FG = 'look')
        z         <- text2num(bgm, paste0('box', b - 1,'.botz'), FG = 'look')
        box.lyrs  <- sum(depths <  - z$Value)
        box.lyrs  <- pmin(box.lyrs, max.nlyrs) # bound by maximum depth
        out       <- rbind(out, data.frame(Boxid = b - 1, Area  = area$Value, Volumen = area$Value *  -z$Value,
                                           Depth =  -z$Value, Layers = box.lyrs))
        vol[b, 1 : box.lyrs] <- area$Value * depths[1 : box.lyrs]
    }
    vol                <- cbind(out$Area, vol) # to include the sediment layer for later calculations
    vol                <- t(vol[, ncol(vol) : 1])
    vol2               <- vol ## arrange not for infauna
    vol2[nrow(vol2), ] <- 0
    if(any(out$Area < 1)) warning('\nOne (or more) of the boxes areas is less than 1m2,  Check if the right BGM file in xy coordinates')
    out[out$Depth <= 0, 2 : ncol(out)] <- 0
    out <- list(info   = out,
                Vol    = vol2,
                VolInf = vol)
    return(out)
}
