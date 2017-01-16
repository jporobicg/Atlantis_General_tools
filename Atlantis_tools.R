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
    ifelse((length(list.l) > 2 ), par(mfrow = c(3, 3), cex = 1.4), par(mfrow = c(1, 2), cex = 1.4))
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
    data$FG <- ifelse(data[, column] %in% c('Callanthias platei', 'Pampanito', 'Chromis meridiana', 'Scartichthys variolatus','Suezichthys sp.','Malapterus reticulatus','Pseudolabrus gayi','Kyphosus cinerascens','Graniento'), 'Small pelagic fish',
               ifelse(data[, column] %in% c('Corvina', 'Jurel', 'Jurelillo'), 'Large pelagic fish',
               ifelse(data[, column] ==  'Vidriola', 'Vidriola',
               ifelse(data[, column] %in% c('Chironemus bicornis', 'Colorado', 'Cabrilla', 'Chancharro','Aseraggodes bahamondei', 'Pez Mariposa','Jerguilla', 'Monocentris reedi'), 'Small benthic fish',
               ifelse(data[, column] == 'Alfonsino', 'Alfonsino',
               ifelse(data[, column] %in% c('Bacalao', 'Lenguado', 'Paralabrax sp.', 'Aseraggodes bahamondei', 'Congrio', 'Lotella fernandeziana', 'Salmon'), 'Large benthic fish',
               ifelse(data[, column] == 'Tollo', 'CHO',
               ifelse(data[, column] == 'Anguila', 'Anguila',
               ifelse(data[, column] %in% c('Caracol Blanco', 'Ostra', 'Loco', 'Chiton', 'Gastropod', 'Slug', 'Snail'), 'Mollusca',
               ifelse(data[, column] %in% c('Centolla', 'Jaiva', 'Estomatopodo', 'shrimp'), 'Small Crustacean',
               ifelse(data[, column] == 'Pulpo', 'Octupus',
               ifelse(data[, column] == 'Lobo', 'OTA',
               ifelse(data[, column] == 'Breca', 'Breca',
               ifelse(data[, column] == 'Langosta', 'Lobster',
               ifelse(data[, column] == 'Cangrejo dorado', 'Golden Crab',
               ifelse(data[, column] %in% c('Actinia', 'Amphiura', 'Anemone', 'Erizo', 'Estrella de mar', 'Echinoide'), 'Benthic Carnivorous',
               ifelse(data[, column] %in% c('Bernacle', 'Black hydrozoan', 'Briozoo', 'Cucumber', 'Encrusting bryozoan', 'Hydrozoo', 'Polychaete', 'Worm', 'Sponge'), 'Benthic feeder',
               ifelse(data[, column] %in% c('Bare Rock', 'Ruble', 'Sand'), 'Sediment',
               ifelse(data[, column] %in% c('Brown Alga', 'Green Alga', 'Red Alga'), 'Macro Algae',
               ifelse(data[, column] %in% c('Coral'), 'Coral', NA))))))))))))))))))))
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
               lista      <- lista2
                if(!isTRUE(only.s)){
                    lista[[j]] <- lista[[j]][which(lista[[j]][, sea]  == season), ] # Only samples from the same Season
                }
                ## Vector of 0 and 1 with te location of the samples inpol=1 if not inpol = 0
                inpol      <- ifelse(point.in.polygon( - lista[[j]][, lo],  - lista[[j]][, la], lon[i,], lat[i, ]) > 0, 1, 0)
                if(!quiet){
                    cat('\n\n =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = ')
                    cat('\n\nEspecies : ', lista[[j]][1, 9])
                    cat('\nnumero de individuos iniciales ', sum(inpol))
                }
                if(poly$attrib[i, 3] >= 0){
                    ## be Sure you don't have samples on land
                    ## If you have will change the position to the more close polygon
                    change <- which(inpol > 0)
                    for(h in change){
                        pos.new           <- which.min(abs(-lista[[j]][h, la]  - lat[i,]))
                        lista[[j]][h, la] <-  - lat[i, pos.new]
                        lista[[j]][h, lo] <-  - lon[i, pos.new]
                    }
                    inpol <- ifelse(point.in.polygon( - lista[[j]][, lo],  - lista[[j]][, la], lon[i,], lat[i, ]) > 0, 0, 0)
                }
                scal  <- inpol * lista[[j]][, mult]
                if(!quiet){
                    cat('\nFinal number', sum(scal))
                    cat('\nSeason',  season, '\tPolygon(i)', poly[[2]][i, 1], '\tSample(j)', j)
                    readline("\nPress <return to continue")
                }
                if(!isTRUE(is.Crus)){
                    poly$attrib[i, j + 6] <- sum(scal, na.rm = TRUE) / sum(scal > 0)
                #cat('\nsum :', sum(scal), '\t weight :', sum(scal>0), '\tCPUE', poly$attrib[i, j + 6])
                } else  {
                    ## Devide the total catch by the number of traps
                    n.traps <- length(levels(as.factor(lista[[j]][which(inpol > 0), la])))
                    poly$attrib[i, j + 6] <- sum(scal) / n.traps
                }
            }
        }
        names(poly$attrib)[7 : ((7 + length(lista)) - 1)] <- names(lista)
        output[[season]] <- poly
    }
    return(output)
}



plot.map <- function(poly, xlim = NULL, ylim = NULL, sc = 4, OnlyPoly = TRUE, leg.color = 'normal', specie=FALSE, leg = NULL, save = FALSE, name = NULL, ...){
    #==================================================================#
    # poly     = Poly data form the read.poly() fuction                #
    # xlim     = vector of longitude limits                            #
    # ylim     = vector of latitude limits                             #
    # sc       = Scale factor (between lat and lon)                    #
    # OnlyPoly = Information for lon and lat is obtained from the poly #
    # leg.col  = Color of the legend base on the atributes             #
    # leg      = Legend                                                #
    # =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = #

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
        color <- 'gray90'
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
            var     <- var / sum(var)
            quart   <- 1 : length(var)
            color   <- rgb(t(col2rgb("royalblue")), alpha = var * 255, maxColorValue = 255)
            etiq     <- seq(from = 0, to = 100,  by = 25)
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
        par(mar=c(4, 4, 4, 4) + 0.1,oma=c(3,3,3,3), ...)
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
##' @param height Vector of the height at age of the functinal group, if is not provided the stimation would be using 1/5 of the total length
##' @param ratio is the ratio between the height and the with of the individual. In other words, if the ratio is .5, the width is onlye the hal of the height for thar specie. If is not provided the value of width would be the same than height
##' @param time.l Vector proportion of time invested for the specie searching for food.
##' @param by.group  Arrange the output for a easy manipulation
##' @return A vector witht the values of clereances for each functional group in mgNm3d-1
##' @author Demiurgo
clearance <- function(fg = NULL, speed, len, height = NULL, ratio = NULL, time.l = NULL, by.group = TRUE){
    ## Assumption Is length in 1 seconds
    speed <- ifelse(is.na(speed), len * 3600, speed) #speed mh-1
    ## I assume the the height is at least 1/5 of the length
    if(is.null(height)) height <-  len / 5
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

    ## If is reserve
    if(isTRUE(reserve)){
        factor  <- weight / 3.65
        weight  <- factor * 2.65
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
##' @param AgeClass
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
    temp <- ifelse(db[ubi, colSp] == 'Octupus', round(db[ubi, colW] / 723),
            ifelse(db[ubi, colSp] == 'Breca', round(db[ubi, colW] / 746),
            ifelse(db[ubi, colSp] == 'Anguila', round(db[ubi, colW] / 1013), NA)))
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
