uv2rho <- function(var, varname=NULL){
    ##----------------------------------------------------------------------------##
    ## This function transfert a field at u or v points to a field at rho points  ##
    ##----------------------------------------------------------------------------##
    MARGIN <- 1
    V_R    <- 2
    if(varname == 'u') {
        V_R    <- 1
        MARGIN <- MARGIN + 1
    }
    incr <- dim(var)[V_R]
    cx   <- apply(var, MARGIN, function(x){
        temp <- 0.5 * (x[1 : (incr - 1)] + x[2 : incr])
        temp <- c(temp[1], temp, temp[incr])
        return(temp)
    }
    )
    if(varname == 'v') cx <- t(cx)
    return(cx)
}

deg2rad <- function(deg){
    rad <- deg*pi/180
    return(rad)
}

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
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

okubo <- function(lat, lon, u, v){
    ## Global variables
    size  <- dim(u)
    distx <- u * NA
    ## transformation to radians
    lat2 <- deg2rad(lat)
    lon2 <- deg2rad(lon)
    ## distance in lon
    for(j in 1 : size[1]){
        for(i in 1 : size[2]){
            distx[j, i] <- gcd.hf(lon2[1,i],lat2[1,i],lon2[j, i],lat2[j, i])
        }
    }

    ## dist in lat
    disty <- v * NA
    for(i in 1 : size[2]){
        for(j in 1 : size[1]){
            disty[j, i] <- gcd.hf(lon2[j, 1], lat2[j, 1],lon2[j, i],lat2[j, i])
        }
    }
    ## Definintion vars
    duy <- dvx <- dvy <- dux <- dx <- dy <- u * NA
    ## dy
    dy <- u * NA
    h1 <- disty[,1 : (size[2] - 2)]
    h2 <- disty[,3 : size[2]]
    dy[,2 : (size[2] - 1)] <- h2 - h1
    ## dx
    h1 <- distx[1 : (size[1] - 2),]
    h2 <- distx[3 : size[1],]
    dx[2 : (size[1] - 1),] <- h2 - h1 # diferencial de u en x
    ## differentials
    ## dudx
    h1 <- u[1 : (size[1] - 2),]
    h2 <- u[3 : size[1],]
    dux[ 2 : (size[1] - 1),] <- h2-h1 # diferencial de u en x
    ##dvdy
    h1 <- v[,1 : (size[2] - 2) ]
    h2 <- v[,3 : size[2] ]
    dvy[,2 : (size[2] - 1)] <- h2 - h1 # diferencial de v en y
    ##dvdx
    h1 <- v[ 1 : (size[1] - 2), ]
    h2 <- v[ 3 : size[1],]
    dvx[2 : (size[1] - 1),] <- h2-h1 # diferencial de v en x
    ##dudy
    h1 <- u[,1 : (size[2] - 2)]
    h2 <- u[,3 : size[2]]
    duy[,2 : (size[2] - 1)] <- h2-h1 # diferencial de u en y
    ## Flux parameters
    sn <- dux/dx - (dvy/dy)  #Componente normal corriente
    ss <- dvx/dx + duy/dy  # Zisalle de la corriente
    w  <- dvx/dx - duy/dy  # vorticidad relativa w
    ow <- sn^2 + ss^2 - w^2  # Parametro de Okubo-Weiss


    out <- list(relat.vort=w, okubo=ow, distx = distx, disty = disty)
    return(out)
}

edd.pos <-function(eddies.pos, distx = NA, disty = NA, min.area = 100){
    #############################################################
    # Patern recognition to look for the position of the eddie  #
    ############################################################
    disty <- uv2rho(apply(disty, 1, diff), varname= 'u')/ 1000
    distx <- uv2rho(apply(distx, 2, diff), varname= 'u') / 1000
    area  <- t(disty) * distx
    pos.list <- area.tot <- n.cells <- list()
    n.eddies <- 1
    while(!is.null(dim(eddies.pos))){
        if(dim(eddies.pos)[1] <=  1) break
        i <- 1
        eddies.pos1 <- matrix(eddies.pos[i, ], ncol = 2)
        while((eddies.pos1[nrow(eddies.pos1), 1] == (eddies.pos[i + 1, 1] - 1)) &
              (eddies.pos1[1, 2] == eddies.pos[i + 1, 2])){
                  eddies.pos1 <- rbind(eddies.pos1, eddies.pos[i + 1, ])
                  eddies.pos <- matrix(eddies.pos[ - c(i, (i + 1)), ], ncol = 2)
                  i <- i + 1
                  if(dim(eddies.pos)[1] <= i | dim(eddies.pos)[1] == 0) break
              }
        newcolloop = TRUE
        while(isTRUE(newcolloop)){
            size.orig  <- dim(eddies.pos1)
            l.col  <- tail(eddies.pos1[, 2], 1)
            samecol = TRUE
            while(isTRUE(samecol) & !is.null(dim(eddies.pos))){
                newval <- eddies.pos[eddies.pos[, 2] == (l.col + 1), ]
                if(is.vector(newval)) newval <- matrix(newval, ncol = 2)
                new.pos <- which(eddies.pos[,2] == (l.col + 1))
                oldval  <- unique(eddies.pos1[which(eddies.pos1[,2] %in% c(l.col, l.col + 1)), 1])
                oldval  <- unique(c(oldval - 1, oldval, oldval + 1))
                mat.val <- which(newval[, 1] %in% oldval)
                if(length(mat.val) > 0){
                    eddies.pos1 <- rbind(eddies.pos1, newval[mat.val, ])
                    eddies.pos <- eddies.pos[ - new.pos[mat.val], ]
                } else{
                    samecol = FALSE
                }
            }
             if(!size.orig[1] < dim(eddies.pos1)[1]) newcolloop = FALSE
        }
        if(dim(eddies.pos1)[1] == 1) eddies.pos <- eddies.pos[ - 1, ] # just in case Have only one point
        temp.area <- sum(area[eddies.pos1])
        if(temp.area >= min.area){
            pos.list[[n.eddies]] <- eddies.pos1
            area.tot[[n.eddies]] <- temp.area
            n.cells[[n.eddies]]  <- dim(eddies.pos1)[1]
            n.eddies             <- n.eddies + 1
        }
    }
    output <- list(position = pos.list, area = area.tot, n.cell = n.cells)
    return(output)
}

eddy.bord <- function(eddy){
    #########################################
    # look for the border of the eddy       #
    #########################################
    levels <- unique(eddy[, 2])
    for( i in 1 : length(levels)){
        if(i==1) {
            bord   <- eddy[which(eddy[, 2] == levels[1]), ]
        } else if(i > 1 & i < length(levels)){
            sub.ed <- eddy[eddy[, 2] == levels[i], ]
            if(!is.matrix(sub.ed)) sub.ed <- matrix(sub.ed, ncol = 2)
            rang   <- unique(range(sub.ed[, 1]))
            bord   <- rbind(bord, sub.ed[sub.ed[, 1] %in% rang,])
        } else{
            bord   <- rbind(bord, eddy[which(eddy[, 2] == tail(levels, 1)), ])
        }
    }
    return(bord)
}

eddies.year <- function(files, folder){
    ################################################################################
    ##  integrated code. Identification of the eddies and the borders by year
    border <- eddy <- list()
    id.tot <- 1
    for(n.f in 1 :  length(files)){
        nc   <- open.nc(paste(folder, files[n.f], sep = '/'))
        if(n.f == 1){
            lat  <- var.get.nc(nc, "lat_rho")
            lon  <- var.get.nc(nc, "lon_rho")
            mask <- var.get.nc(nc, "mask_rho")
        }
        mask[which(mask == 0, arr.ind = TRUE)] <- NA
        u    <- var.get.nc(nc, "u")
        v    <- var.get.nc(nc, "v")
        tid  <- dim(v)[4]
        close.nc(nc)
        for(id in 1 : tid){
            u.new <- u[ , , 38, id]
            u.new <- uv2rho(u.new,'u')
            u.new <- u.new * mask
            v.new <- v[ , , 38, id]
            v.new <- uv2rho(v.new, 'v')
            v.new <- v.new * mask
            ##--------------------------------------------------------------
            ##  Estimation of the Okobo Weis parameter
            ##--------------------------------------------------------------
            ## Okubo parameter estimation
            ok <- okubo(lat,lon,u.new,v.new)
            ## Stantar deviation
            ok.sd <- -0.2 * sd(ok$okubo, na.rm = TRUE)
            ## Cyclonic rotation negative vorticity (oposite northern hemisphere)
            id_cic <- with(ok, which(okubo <= ok.sd & relat.vort < 0, arr.ind = TRUE))
            ## Anticyclonic rotation positive vorticity (oposite northern hemisphere)
            id_ant <- with(ok, which(okubo <= ok.sd & relat.vort > 0, arr.ind = TRUE))
            ## Identified eddies by area
            cic.n.eddies    <- edd.pos(id_cic, ok$distx, ok$disty, min.area = 300)[[1]]
            ant.n.eddies    <- edd.pos(id_ant, ok$distx, ok$disty, min.area = 300)[[1]]
            eddy[[id.tot]]            <- list(anticyclonic = ant.n.eddies,
                                    cyclonic     = cic.n.eddies)
            ## Borders eddies
            cic.eddies.bord <- lapply(cic.n.eddies, eddy.bord)
            ant.eddies.bord <- lapply(ant.n.eddies, eddy.bord)
            border[[id.tot]]      <- list(anticyclonic = ant.eddies.bord,
                                     cyclonic      = cic.eddies.bord)
            id.tot <- id.tot + 1
        }
    }
    year        <- list(borders = border, eddy  = eddy)
    return(year)
}


dilat <- function(matrix){
    ## dilat the matrix in three cell in avery direction
    ##    I use 3 cells because the eddies are moving 2km per day. in 5 days is almost 10 km,  and that is roughly 3 cells
    nr  <- nrow(matrix)
    row <- rep(sort(rep(c(-3, 0, 3), nr)), 7)
    col <- sort(rep(c(-3, 0, 3), nr * 7))
    out <- unique(cbind(matrix[, 1] + row, matrix[, 2] + col))
}

finder <- function(org.edd, step2, tid){
    ## to find the eddy in the next step
    if(any(duplicated(rbind(org.edd, matrix(na.omit(unlist(lapply(step2, t))), ncol = 2, byrow = TRUE))))){
        for( id.next in 1 : length(step2)){## check all the eddies in the next timestep
            next.edd <- step2[[id.next]]
            if(any(is.na(next.edd))) next
            if(any(duplicated(rbind(org.edd,next.edd)))){ ## if find the next eddy
                out <- cbind(tid + 1, id.next[1], pam(next.edd, 1)$medoids[, 1], pam(next.edd, 1)$medoids[, 2])
                break
            }
        }
        return(out)
    } else {
        return(NA)
    }
}


eddy.track <- function(file, cyclonic = TRUE, min.temp = 3){
    ###################################################################################
    ##      code to track the eddies in base a standar distance of 3 cell (10 - 15 km)
    ##     file       :  With the location of the eddies
    ##     cyclonic   :  Logic to decide cyclonic or anticyclonic
    ##     min.temp   :  minimun time duration of the eddy
    library(cluster)
    load(file)
    rot        <- ifelse(isTRUE(cyclonic), 2, 1)
    eddy.track <- list(NA)
    tot.edd    <- 0
    tot.step   <- length(out.year$borders)
    for(tid in 1 : (tot.step - 1)){
        for(n.edd in 1 : length(out.year$borders[[tid]][[rot]])){
            step1    <- out.year$borders[[tid]][[rot]]
            step2    <- out.year$borders[[tid + 1]][[rot]]
            new.eddy <- step1[[n.edd]]
            if(any(is.na(new.eddy))) next
            eddy.org <- dilat(new.eddy)
            ## find the eddy in the next step
            h <- finder(eddy.org, step2, tid)
            if(any(is.na(h))){
                ## remove the eddy
                out.year$borders[[tid]][[rot]][[n.edd]] <- matrix(c(NA, NA), ncol = 2)
                next
            } else {
                tot.edd               <- tot.edd + 1
                eddy.track[[tot.edd]] <- rbind(c(tid, n.edd,  pam(step1[[n.edd]], 1)$medoids[, 1], pam(step1[[n.edd]], 1)$medoids[, 2]), h)
                cur                   <- h[2]
                ## remove the eddy
                out.year$borders[[tid]][[rot]][[n.edd]] <- matrix(c(NA, NA), ncol = 2)
                for(sec.id in (tid + 1) : (tot.step - 1)){
                    if(sec.id >= (tot.step - 1)) break
                    cur.edd <- dilat(out.year$borders[[sec.id]][[rot]][[cur]])
                    step2   <- out.year$borders[[sec.id + 1]][[rot]]
                    h       <- finder(cur.edd, step2, sec.id)
                    if(any(is.na(h))){
                        ## remove the eddy
                        out.year$borders[[sec.id]][[rot]][[cur]] <- matrix(c(NA, NA), ncol = 2)
                        break
                    } else {
                        eddy.track[[tot.edd]] <- rbind(eddy.track[[tot.edd]], h)
                        out.year$borders[[sec.id]][[rot]][[cur]] <- matrix(c(NA, NA), ncol = 2)
                        cur <- h[2]
                    }
                }

            }
        }
    }
    limit <- which(unlist(lapply(eddy.track, function(x) (dim(x)[1] > min.temp))) == 1)
    output <- eddy.track[limit]
    return(output)
}


bypol.eddy <- function(eddy.file, coor.dat, polygons.dat, cyclonic = TRUE, relocate = NULL, centroids = NULL){
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##                              Function to identified the number of eddies by polygons     ##
    ## eddy.file   =  File that contain the trajectory of the eddies                            ##
    ## coor.dat    =  coordinates of the realm (lat and lon of domain)                          ##
    ## polygos.dat =  Polygons of the model                                                     ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ## Lybraries
    library(sp)
    ## Tools
    source('/home/demiurgo/Documents/PhD/Atlantis_Model/tools/General_tools/Atlantis_tools.R')
    ## load data
    load(coor.dat)
    load(eddy.file)
    if(isTRUE(cyclonic)){
        eddies.track <- cic.track.eddies
    } else {
        eddies.track <- ant.track.eddies
    }
    polygons <- read.poly(polygons.dat)
    lat.mod  <- coord[[2]]
    lon.mod  <- coord[[1]]
    lon      <- as.matrix(polygons$coor[, c(seq(from = 1, to = dim(polygons$coor)[2], by = 2))])
    lat      <- as.matrix(polygons$coor[, c(seq(from = 2, to = dim(polygons$coor)[2], by = 2))])


    ## get the coordinates
    n.list <- lapply(eddies.track, function(x) data.frame(season = ifelse(x[, 1] < 3, 'Summer',
                                                                       ifelse(x[, 1] < 6, 'Fall',
                                                                       ifelse(x[, 1] < 9, 'Winter', 'Spring'))),
                                                              lat = lat.mod[x[, 3]], lon = lon.mod[x[, 4]]))

    ## relocate eddies from polygos
    if(!is.null(relocate)){
        n.areas <- length(relocate)
        for(rel.a in 1 :  n.areas){
            polygon <-  which(polygons[[2]][, 1] == relocate[rel.a])
            n.list  <- move.eddy(n.list, centroids, polygon)
        }
    }
    ## Loop to classified the eddies by polygons
    t.eddy <- matrix()
    for( i in 1 : nrow(polygons$attrib)){
        ## loop throught the polygons
        Summer <- Spring <- Winter <- Fall <- list(NA)
        for(j in 1 : length(n.list)){
            ## loop throught the eddies
            poly   <- na.omit(cbind(lon[i, ], lat[i, ]))
            eddy   <- cbind(n.list[[j]][, 2 : 3])
            in.pol <- point.in.polygon(eddy[, 1], eddy[, 2], poly[, 1], poly[, 2])
            if(all(in.pol == 0)) next
            pos    <- which(in.pol == 1)
            ## Sum of eddies by polygons and by season
            ## Summer
            Summer <- ifelse(is.na(Summer), sum(n.list[[j]][pos, 1] == 'Summer'),
                             Summer + sum(n.list[[j]][pos, 1]=='Summer'))
            ## Spring
            Spring <- ifelse(is.na(Spring), sum(n.list[[j]][pos, 1] == 'Spring'),
                             Spring + sum(n.list[[j]][pos, 1]=='Spring'))
            ## Winter
            Winter <- ifelse(is.na(Winter), sum(n.list[[j]][pos, 1] == 'Winter'),
                             Winter  + sum(n.list[[j]][pos, 1]=='Winter'))
            ## Fall
            Fall   <- ifelse(is.na(Fall), sum(n.list[[j]][pos, 1] == 'Fall'),
                             Fall + sum(n.list[[j]][pos, 1]=='Fall'))
        }
        BoxId  <- polygons[[2]][i, 1]
        Area   <- as.numeric(polygons[[2]][i, 6])
        if(i == 1) {
            t.eddy <-cbind(Summer, Spring, Winter, Fall, BoxId, Area)
        } else{
            t.eddy <- rbind(t.eddy, cbind(Summer, Spring, Winter, Fall, BoxId, Area))}
    }
    ## Remove the NA
    t.eddy[which(is.na(t.eddy), arr.ind = TRUE)] <- 0
    rel.impor          <- cbind(t(apply(t.eddy, 1, function(x) as.numeric(x[1 : 4]) / as.numeric(x$Area))), as.numeric(t.eddy[, 5]))
    max.r              <- max(unlist(rel.impor[, 1 : 4]))
    rel.impor[, 1 : 4] <-rel.impor[, 1 : 4] / max.r
    ## Output
    output <- list(N.Eddies = t.eddy, P.Eddies = rel.impor)
    return(output)
}


w.area <- function(n.eddies, area.rm = NULL){
    ## Weighted by Area
    W      <- t(apply(n.eddies, 1, function(x) as.numeric(x[1 : 4]) / as.numeric(x[6])))
    W.lim  <- mean(W[which(W > 0)]) + sd(W[which(W > 0)]) * 3
    ## Dealing with outlier
    if(any(W > W.lim)){
        W.temp                      <- max(W[- which(W > W.lim)])
        W[which(W > W.lim)] <-  W.temp + sd(W[which(W > 0)])
        W.lim <- mean(W[which(W > 0)]) + sd(W[which(W > 0)]) * 3
    }
    max.r            <- max(W)
    W.tot            <- W / max.r

    ## positive anomalies
    W.mean        <- mean(W.tot[which(W.tot > 0)])
    loc           <- which(W.tot == 0, arr.ind = TRUE)
    loc.rm        <- which(n.eddies[, 5] %in% area.rm)
    anomalies     <- W.tot - W.mean
    an.stand      <- anomalies  +  1
    ## Output
    W.tot[loc]          <-  min(W.tot[which(W.tot > 0)])
    W.tot[loc.rm, ]       <- NA
    W.tot               <- cbind(W.tot, n.eddies[, 5 : 6])
    colnames(W.tot)     <- c('Summer', 'Spring', 'Winter', 'Fall', 'BoxId', 'Area')
    an.stand[loc]       <- min(an.stand[which(an.stand > 0)])
    an.stand[loc.rm, ]    <- NA
    an.stand            <- cbind(an.stand, n.eddies[, 5 : 6])
    colnames(an.stand)  <- c('Summer', 'Spring', 'Winter', 'Fall', 'BoxId', 'Area')

    output <- list(Stand.Eddies = W.tot, Anomalies = an.stand)
    return(output)
}


move.eddy <- function(eddies, centroids, polygon){
    ##########################################################
    ##  This function move the eddies from inside on an area #
    ##  to another area closer                               #
    ##########################################################
    poly   <- na.omit(cbind(lon[polygon, ], lat[polygon, ]))
    for(j in 1 : length(eddies)){
        eddy   <- cbind(eddies[[j]][, 2 : 3])
        in.pol <- point.in.polygon(eddy[, 1], eddy[, 2], poly[, 1], poly[, 2])
        if(all(in.pol == 0)) next
        pos           <- which(in.pol == 1)
        in.island     <- eddies[[j]][pos, 2 : 3]
        in.island.rad <- deg2rad(in.island)
        centroid.rad  <- deg2rad(centroid)
        for(i.ed in 1: length(pos)){
            dist.pos <- which.min(gcd.hf(pos1 = in.island.rad[i.ed, ], pos2 = -centroid.rad, mat = TRUE))
            cat('\nMoving Eddy : ',j , ' position ',i.ed , ' to ',  centroid[dist.pos, 3])
            eddies[[j]][pos[i.ed], 2 : 3] <- centroid[dist.pos,1:2]
        }
    }
    return(eddies)
}


sigma2zeta <- function(h, hc, theta, b, N){
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## USAGE: z=zdepth_roms(h,hc,theta,b,N)
    ## z z-coordinate center of layer (m)
    ## h depth (m)
    ## hc = 20;
    ## theta = 5;
    ## b = 0.6;
    ## N number of layers

    ## Center of layers
    ds = 1 / N;
    s = seq(from=(- 1 + ds / 2), to = (0 - ds / 2), by = ds)
    A = sinh(theta * s) / sinh(theta);
    B = (tanh(theta * (s + .5)) - tanh(theta * .5)) / (2 * tanh(theta * .5));
    C = (1 - b) * A + b * B;
    z = hc * s + (h - hc) * C;
    return(z)
}
## ******************************************************************************* ##
## ~                                   asdf                                    ~ ##
## ******************************************************************************* ##
