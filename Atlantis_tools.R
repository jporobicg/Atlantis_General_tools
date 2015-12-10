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
    bearing  <- seq(1, 360, 1) * pi / 180   #  Bearing in radians
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

points.img <- function(x.lim = NA, y.lim = NA, filename = NULL, grid = FALSE){
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
    flip    <- function(matriz) t(matriz)[,nrow(matriz):1]
    red.mat <- matrix(NA, img@size[1], img@size[2])
    red.mat <- img@red

    image(flip(red.mat), bty = 'n' ,axes = FALSE)
    if(grid == TRUE) grid(100, col = 2, lwd = .3, lty = 1)

    ## locator to identified the min and the max
    cat('set the lower value of X and the maximun\n')
    x.axs <- locator(2)
    cat('set the lower value of Y and the maximun\n\n\n')
    y.axs  <- locator(2)

    xseq    <- seq(x.axs$x[1], x.axs$x[2], length.out = 10000)
    yseq    <- seq(y.axs$y[1], y.axs$y[2], length.out = 10000)
    x.u.pos <- seq(x.lim[1], x.lim[2], length.out = 10000)
    y.u.pos <- seq(y.lim[1], y.lim[2], length.out = 10000)
    cat('Select the point,  when ready right - click')
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
    data$FG <- ifelse(lc[, column] %in% c('Jurelillo', 'Jurel', 'Pampanito', 'Pez Mariposa', 'Jerguilla'), 'Small pelagic fish',
               ifelse(lc[, column] %in% c('Corvina', 'Vidriola'), 'Large pelagic fish',
               ifelse(lc[, column] %in% c('Colorado', 'Graniento', 'Cabrilla', 'Chancharro'), 'Small benthic fish',
               ifelse(lc[, column] == 'Alfonsino', 'Alfonsino',
               ifelse(lc[, column] %in% c('Bacalao', 'Tollo', 'Anguila', 'Lenguado', 'Congrio'), 'Large benthic fish',
               ifelse(lc[, column] %in% c('Caracol Blanco', 'Ostra', 'Loco'), 'Mollusca',
               ifelse(lc[, column] %in% c('Centolla', 'Jaiva'), 'Small Crustacean',
               ifelse(lc[, column] == 'Pulpo', 'Octupus',
               ifelse(lc[, column] == 'Breca', 'Breca',
               ifelse(lc[, column] %in% c('Erizo', 'Estrella de mar'), 'Sea Urchin' , NA))))))))))
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
by.pol <- function(poly, lista, la = 4, lo = 5, mult = 6, sea = 10){
    ## Function to clasified the species by polygon
    library(sp)
    lon <- as.matrix(poly$coor[, c(seq(from = 1, to = dim(poly$coor)[2], by = 2))])
    lat <- as.matrix(poly$coor[, c(seq(from = 2, to = dim(poly$coor)[2], by = 2))])
    season <- list()
    lista2 <- lista
    for(s in 1 : 4){
        for( i in 1 : nrow(poly$attrib)){
            for(j in 1 : length(lista)){
                lista <- lista2
                lista[[j]] <- lista[[j]][ - which(lista[[j]][, sea] != s), ]
                inpol <- ifelse(point.in.polygon( - lista[[j]][, lo],  - lista[[j]][, la], lon[i,], lat[i, ]) > 0, 1, 0)
                if(poly$attrib[i, 3] > 0){
                    change <- which(inpol>0)
                    for(h in change){
                        pos.new <- which.min( (lista[[j]][h, la] *  - 1) - lat[i,])
                        lista[[j]][h, la] <-  - lat[i, pos.new]
                        lista[[j]][h, lo] <-  - lon[i, pos.new]
                    }
                    inpol <- ifelse(point.in.polygon( - lista[[j]][, lo],  - lista[[j]][, la], lon[i,], lat[i, ]) > 0, 0, 0)
                }
                scal  <- inpol * lista[[j]][, mult]
                poly$attrib[i, j + 6] <- sum( scal)
            }
        }
        names(poly$attrib)[7 : ((7 + length(lista)) - 1)] <- names(lista)
      #  season[[s]] <- list(poly=poly)
        season[[s]] <- poly
    }
    return(season)
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
