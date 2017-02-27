##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Atlantis feeding tool
##' @param prm.file Atlantis parameter file
##' @param grp.file Atlantis group file
##' @param nc.file Atlantis initial conditions file
##' @return Display in the browser the pPREY matriux,  the initial abundace on prey,  the overlap matrix and the predator preference
##' @author Demiurgo
feeding.mat.shy <- function(prm.file, grp.file, nc.file){
    ## I need to create something more elegant for this
    library(shiny)
    library(ncdf4)
    library(reshape)
    library(tidyverse)
    ## Reading files
    groups.csv <- read.csv(grp.file)
    prm        <- readLines(prm.file, warn = FALSE)
    ## availability matrix
    Ava.mat            <- text2num(prm, 'pPREY', Vector=TRUE)
    colnames(Ava.mat)  <- c(as.character(groups.csv$Code), 'DLsed', 'DRsed', 'DCsed')
    ## Biomass and Gape size
    out.Bio  <- Bio.func(nc.file, groups.csv)
    Struct   <- out.Bio[[1]]
    Biom.N   <- out.Bio[[2]]
    Gape     <- gape.func(groups.csv, Struct, Biom.N, prm)
    Over.mat <- Over.mat.func(Ava.mat, Gape[[1]])
    bio.a    <- Bio.age(Biom.N, age = Gape[[2]], Over.mat)
    bio.juv  <- bio.a[[1]]
    bio.adl  <- bio.a[[2]]
    bio.juv  <- data.frame(FG = bio.a[[1]][, 1], Biomass = as.numeric(bio.a[[1]][, 2]))
    bio.adl  <- data.frame(FG = bio.a[[2]][, 1], Biomass = as.numeric(bio.a[[2]][, 2]))
    b.juv    <- bio.juv[complete.cases(bio.juv), ]
    b.adl    <- bio.adl[complete.cases(bio.adl), ]
    ## Total feeding
    real.feed  <- Over.mat * NA
    pred       <- row.names(Over.mat)
    for( pd in 1 : nrow(Over.mat)){
        ## Getting the number of biomass needed by each functional group
        c.pred      <- unlist(strsplit(pred[pd],'pPREY'))[2]
        predator    <- gsub(pattern = "[[:digit:]]+", '\\1', c.pred)
        a.pred.prey <- as.numeric(unlist(strsplit(c.pred, predator)))
        pry.loc     <- which(bio.adl[, 1] %in% predator)
        if(length(a.pred.prey) == 0 || is.na(a.pred.prey)) a.pred.prey[2] <- 2
        ## Young Predator
        if(a.pred.prey[2] == 1){
            ## Young Prey
            real.feed[pd, ] <- (Over.mat[pd, ] * as.numeric(bio.juv[, 2]))
        } else {
            ## Adult Prey
            real.feed[pd, ] <- (Over.mat[pd, ] * as.numeric(bio.adl[, 2]))
        }
    }
    ## Real Overlap matrix Including pPREY and Overlap matrix
    t.o.mat <- t(Over.mat * Ava.mat)
    t.o.mat[which(t.o.mat > 0)] <- 1
    ## Plot output
    real.feed      <- real.feed * Ava.mat
    ##real.vec.pprey <- melt(t(log(real.feed)))
    ## Shiny Application
    shinyApp(
        ui <- fluidPage(
            headerPanel('Atlantis Diet Tool'),
            fluidRow(
                column(2,
                       wellPanel(
                           selectInput('ycol', 'Predator',  row.names(Ava.mat)),
                           selectInput('xcol', 'Prey', colnames(Ava.mat)),
                           mainPanel("Original Value: ", verbatimTextOutput("numPoints")),
                           mainPanel("Current Value: ", verbatimTextOutput("CurPoints")),
                           numericInput("num", label = "New Value", value = 0, min = 0, max = 1, step = 0.00001)
                       ),
                       actionButton("do", label = "Change Value"),
                       br(),
                       br(),
                       actionButton("save", label = "Write pPPREY Matrix")
                       ),
                column(10,
                       wellPanel(
                           tabsetPanel(
                               tabPanel('Efective Predation',#
                                        plotOutput('plot1', width = "100%", height = "900px")
                                        ),
                               tabPanel('Availability matrix',#
                                        plotOutput('plot3', width = "100%", height = "900px")
                                        ),
                               tabPanel('Overlap matrix',#
                                        plotOutput('plot2', width = "100%", height = "900px")
                                        ),
                               tabPanel('% of predation pressure',#
                                        plotOutput('plot4', width = "100%", height = "900px")
                                        ),
                               tabPanel('Total biomass prey',#
                                        h3('Juvenile biomass'),
                                        plotOutput('plot5', width = "100%", height = "400px"),
                                        br(),
                                        h3('Adult biomass'),
                                        plotOutput('plot6', width = "100%", height = "400px")
                                        )
                           )
                       )
                       )
            )
        ),
        function(input, output, session) {
            ## Combine the selected variables into a new data frame
            N.mat     <- reactiveValues()
            N.mat$Ava <- Ava.mat
            newEntry  <- observe({
                if(input$do > 0) {
                    newval <- isolate(input$num)
                    col.ch <- isolate(which(colnames(Ava.mat) == input$xcol))
                    row.ch <- isolate(which(row.names(Ava.mat) == input$ycol))
                    isolate(N.mat$Ava[row.ch, col.ch] <- newval)
                }
            })
            observeEvent(input$save, {
                saveData(N.mat$Ava)
            })

            linex <- reactive( {
                which(sort(colnames(Ava.mat)) == input$xcol)
            })
            liney <- reactive({
                which(sort(row.names(Ava.mat)) == input$ycol)
            })

            rff <- reactive({
                t(log(real.feed * N.mat$Ava))
            })
            rff2 <- reactive({
                t((real.feed * N.mat$Ava) / rowSums(real.feed * N.mat$Ava, na.rm=TRUE)) * 100
            })
            output$plot1 <- renderPlot({
                ggplot(data = melt(rff()),
                       aes(x = X1, y = X2, fill = value)) + geom_tile() +
                    scale_fill_gradient(limits=c(0, max(rff(), na.rm = TRUE)), name = 'Predation value', low="white", high="red", na.value = 'white')  +
                    theme(panel.background = element_blank()) + labs(x = 'Prey', y = 'Predator') + scale_x_discrete(position = "top")+
                    annotate("rect", xmin = linex() -.5, xmax = linex() +.5, ymin = 0, ymax = ncol(rff()) + 1,
                             alpha = .1, colour = 'royalblue') +
                    annotate("rect", xmin =  - .5, xmax = nrow(rff()) + .5, ymin = liney() - .5, ymax = liney() + .5,
                             alpha = .1, colour = 'royalblue')
            })
            output$plot2 <- renderPlot({
                ggplot(data = melt(t.o.mat),
                       aes(x = X1, y = X2, fill = value)) + geom_tile(aes( fill = factor(value))) +
                    theme(panel.background = element_blank()) + labs(x = 'Prey', y = 'Predator') + scale_x_discrete(position = "top") +
                    scale_fill_grey(start = .9, end = 0) +
                    annotate("rect", xmin = linex() -.5, xmax = linex() +.5, ymin = 0, ymax = ncol(t.o.mat) + 1,
                             alpha = .1, colour = 'royalblue') +
                    annotate("rect", xmin =  - .5, xmax = nrow(t.o.mat) + .5, ymin = liney() - .5, ymax = liney() + .5,
                             alpha = .1, colour = 'royalblue')
            })
            output$plot3 <- renderPlot({
                ggplot(data = melt(t(N.mat$Ava)),
                       aes(x = X1, y = X2, fill = value)) + geom_tile() +
                    scale_fill_gradient(limits=c(0, max(N.mat$Ava, na.rm = TRUE)), name = 'Predation value', low="white", high="red", na.value = 'white')  +
                    theme(panel.background = element_blank()) + labs(x = 'Prey', y = 'Predator') + scale_x_discrete(position = "top")+
                    annotate("rect", xmin = linex() -.5, xmax = linex() +.5, ymin = 0, ymax = ncol(t(N.mat$Ava)) + 1,
                             alpha = .1, colour = 'royalblue') +
                    annotate("rect", xmin =  - .5, xmax = nrow(t(N.mat$Ava)) + .5, ymin = liney() - .5, ymax = liney() + .5,
                             alpha = .1, colour = 'royalblue')
            })
            output$plot4 <- renderPlot({
                ggplot(data = melt(rff2()),
                       aes(x = X1, y = X2, fill = value)) + geom_tile() +
                    scale_fill_gradient(limits=c(0, 100), name = 'Precentage of pressure', low="white", high="red", na.value = 'white')  +
                    theme(panel.background = element_blank()) + labs(x = 'Prey', y = 'Predator') + scale_x_discrete(position = "top")+
                    annotate("rect", xmin = linex() -.5, xmax = linex() +.5, ymin = 0, ymax = ncol(rff2()) + 1,
                             alpha = .1, colour = 'royalblue') +
                    annotate("rect", xmin =  - .5, xmax = nrow(rff2()) + .5, ymin = liney() - .5, ymax = liney() + .5,
                             alpha = .1, colour = 'royalblue')
            })
            output$plot5 <- renderPlot({
                ggplot(data = b.juv, aes(x = FG, y = log(Biomass), fill=FG)) +
                    geom_bar(colour="black", stat="identity") +
                    guides(fill = FALSE)+
                    xlab("Functional Groups") + ylab("Biomass [MgN] or Density [MgNm-3]")
            })
            output$plot6 <- renderPlot({
                ggplot(data = b.adl, aes(x = FG, y = log(Biomass), fill=FG)) +
                    geom_bar(colour="black", stat="identity") +
                    guides(fill = FALSE)+
                    xlab("Functional Groups") + ylab("Biomass [MgN] or Density [MgNm-3]")
            })
            output$numPoints <- renderText({
                Ava.mat[which(row.names(Ava.mat) == input$ycol), which(colnames(Ava.mat) == input$xcol)]
            })
            output$CurPoints <- renderText({
                N.mat$Ava[which(row.names(Ava.mat) == input$ycol), which(colnames(Ava.mat) == input$xcol)]
            })

        }
    )
}
## ~~~~~~~~~~~~~~~~~~~~~~ ##
## ~      FUNCTIONS!!   ~ ##
## ~~~~~~~~~~~~~~~~~~~~~~ ##
## getting the RN and SN from the NC file
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Biomass from nc file
##' @param nc.file Atlantis initial condition file
##' @param groups.csv Atlantis groups file
##' @return The biomass for age class and the sturctural nitrogen by age class
##' @author Demiurgo
Bio.func <- function(nc.file, groups.csv){
    nc.out <- nc_open(nc.file)
    Is.off <- which(groups.csv$IsTurnedOn == 0)
    FG     <- as.character(groups.csv$Name)
    Biom.N <- array(data = NA, dim = c(length(FG), max(groups.csv$NumCohorts)))
    Struct <- Biom.N
    for(code in 1 : length(FG)){
        sed     <- ncatt_get(nc.out, varid = paste(FG[code], "_N", sep = ""), attname = "insed")$value
        unit    <- ncatt_get(nc.out, varid = paste(FG[code], "_N", sep = ""), attname = "units")$value
        if(groups.csv$NumCohorts[code] == 1 && groups.csv$IsTurnedOn[code] == 1){
            N.tot <- ncvar_get(nc.out, paste(FG[code], "_N", sep = ""))
            if(all(is.na(N.tot)) || all(N.tot == 0) || sum(N.tot, na.rm = TRUE) == 0){
                ## Getting the total volumen
                water.t <- ncvar_get(nc.out, 'volume')
                w.depth <- ncvar_get(nc.out, 'nominal_dz')
                w.depth[is.na(w.depth)] <- 0
                w.m2    <- colSums(water.t,na.rm=TRUE) / apply(w.depth, 2, function(x) max(x, na.rm = TRUE))
                w.m2[is.infinite(w.m2)] <- NA
                w.m2    <- sum(w.m2, na.rm=TRUE)
                w.m3    <- sum(water.t, na.rm = TRUE)
                #sed     <- ncatt_get(nc.out, varid = paste(FG[code], "_N", sep = ""), attname = "insed")$value
                Biom.N[code, 1] <- ncatt_get(nc.out, varid = paste(FG[code], "_N", sep = ""), attname = "_FillValue")$value
                Biom.N[code, 1] <- ifelse(sed == 1, Biom.N[code, 1] * w.m2, Biom.N[code, 1] * w.m3)
            } else {
                if(length(dim(N.tot)) > 3){
                    N.tot <- N.tot[, , 1]
                } else if(unit == "mg N m-3"){
                    water.t <- ncvar_get(nc.out, 'volume')
                    N.tot   <- N.tot * water.t
                } else if (unit == 'mg N m-2'){
                    water.t <- ncvar_get(nc.out, 'volume')
                    w.depth <- ncvar_get(nc.out, 'nominal_dz')
                    w.depth[is.na(w.depth)] <- 0
                    w.m2    <- colSums(water.t,na.rm=TRUE) / apply(w.depth, 2, function(x) max(x, na.rm = TRUE))
                    w.m2[is.infinite(w.m2)] <- NA
                    N.tot   <- sum(N.tot * w.m2, na.rm = TRUE)
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
    return(list(Struct, Biom.N))
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
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
        txt  <- gsub(pattern = '[ ]+' ,  '|',  text)
        col1 <- col2 <- vector()
        for( i in 1 : length(txt)){
            tmp     <- unlist(strsplit(txt[i], split = '|', fixed = TRUE))
            tmp2    <- unlist(strsplit(tmp[1], split = '_'))
            id.co   <- which(tmp2 %in% FG )
            col1[i] <- tmp2[id.co]
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
            if(tmp[1] %in% c('#','##', '###')) next
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
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param groups.csv Atlantis group file
##' @param Struct Structural weight by age group
##' @param Biom.N Biomass by age groups
##' @param prm Atlantis paramter file
##' @return the limits for the prey based on the predator gape size and prey size
##' @author Demiurgo
gape.func <- function(groups.csv, Struct, Biom.N, prm){
    ## Gape size and adult and young age
    KLP                     <- text2num(prm, 'KLP', FG = as.character(groups.csv$Code))
    KUP                     <- text2num(prm, 'KUP',  FG = as.character(groups.csv$Code))
    age                     <- text2num(prm, '_age_mat', FG = as.character(groups.csv$Code))
    Gape                    <- data.frame(FG = KLP$FG, KLP = KLP$Value, KUP = KUP$Value, Age.Adult = NA)
    pos.Age                 <- which(Gape$FG %in% age$FG)
    Gape$Age.Adult[pos.Age] <- age$Value
    Gape$Age.Young          <- Gape$Age.Adult - 1
    Gape$Age.Young          <- ifelse(Gape$Age.Young == 0,  1, Gape$Age.Young)
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
    return(list(Gape, age))
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param Ava.mat Availavility matrix (or pPREY matrix from the parameter file
##' @param Gape Limit of prey by Gape size
##' @return Overlap Matrix based only in the gape limitation
##' @author Demiurgo
Over.mat.func <- function(Ava.mat, Gape){
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~        Overlap-Matrix    ~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    Over.mat <- Ava.mat * 0
    Prey     <- colnames(Ava.mat)
    Pred     <- row.names(Ava.mat)
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
    return(Over.mat)
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param Biom.N Biomass by Cohort
##' @param age Age of maturation
##' @param Over.mat Overlap matrix
##' @return Biomass by adult and juveniles
##' @author Demiurgo
Bio.age <- function(Biom.N, age, Over.mat){
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
    return(list(bio.juv, bio.adl))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Save pPREY matrix
##' @param data matrix of pPrey modified
##' @return a txt file of the pPREY matrix
##' @author Demiurgo
saveData <- function(data) {
    fileName <- sprintf("%s_NewpPrey.txt", as.integer(Sys.time()))
    rows     <-row.names(data)
    cols     <- c('##    ', colnames(data))
    nprey    <- ncol(data)
    sink(fileName)
    cat(cols)
    for( i in 1 : length(rows)){
        cat(paste('\n', rows[i], nprey, sep = '  '))
        cat('\n', data[i, ])
    }
    sink()
}
