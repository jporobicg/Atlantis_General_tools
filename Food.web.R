##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Function on shiny R to plot the trophic level of the fucntional groups
##' @param diet.file Diet output file from atlantis
##' @param grp.file Group file. Atlantis input file (.cvs)
##' @param quiet help to debug the function
##' @return
##' @author Demiurgo
food.web <- function(diet.file, grp.file,  quiet = TRUE){
    txtHelp <- "<h2>Summary</h2>"
    txtHelp <- paste(txtHelp, "<p>This bit of code help to visualize the food web change during the simulation output from <b>Atlantis</b> run. </p>")
    txtHelp <- paste(txtHelp, "<p>It calculate the trophic position of each functional group at each time step</p>")
    ## Libraries
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 1    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n\n Loading libraries')
    if (!require('shiny', quietly = TRUE)) {
        stop('The package shiny was not installed')
    }
    if (!require('ncdf4', quietly = TRUE)) {
        stop('The package ncdf4 was not installed')
    }
    if (!require('plotrix', quietly = TRUE)) {
        stop('The package plotrix was not installed')
    }
    if (!require('reshape', quietly = TRUE)) {
        stop('The package reshape was not installed')
    }
    if (!require('tidyverse', quietly = TRUE)) {
        stop('The package tidyverse was not installed')
    }
    if (!require('stringr', quietly = TRUE)) {
        stop('The package stringr was not installed')
    }
    if (!require('data.table', quietly = TRUE)) {
        stop('The package data.table was not installed')
    }
    if (!require('RColorBrewer', quietly = TRUE)) {
        stop('The package RColorBrewer was not installed')
    }
    if(!quiet) cat('  ...Done!')
    ## Reading files
    if(!quiet) cat('\n Reading files')
    dat      <- data.frame(fread(diet.file, header = TRUE, sep = ' ', showProgress = FALSE))
    time.stp <- round(range(dat$Time), 0)
    grp.dat  <- read.csv(grp.file)
    stk      <- unique(dat$Stock)
    if(any(names(grp.dat) == 'isPredator')){
        names(grp.dat)[which(names(grp.dat) == 'isPredator')] <- 'IsPredator'
        warning('You should change the name of the column \'isPredator\' for \'IsPredator\' on your group csv file')
    }
    code.fg  <- grp.dat$Code[grp.dat$IsPredator > 0]
    if(!quiet) cat('      ...Done!')
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 2    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n Processing and plotting')

    shinyApp(
        ui <- navbarPage('Atlantis Food Web Tool',
                         tabPanel('Food Web',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 tags$h3('Functional Group'),
                                                 selectInput('foc.fg', 'Functional Group', as.character(code.fg)),
                                                 numericInput("m.stg", "Max. Throphic Connections:", min = 1,  max = 10, value = 4, step = 1),
                                                 numericInput("min", "Min proportion :", min = 0.001,  max = 1, value = 0.01, step = 0.001),
                                                 numericInput("time", "Time Step :", min = time.stp[1],  max = time.stp[2], value = 1, step = 1),
                                                 selectInput('stock', 'Stock', stk))
                                             ),
                                      column(10,
                                             plotOutput('plot1', width = "100%", height = "800px"),
                                             tableOutput('table')
                                             )
                                  )
                                  ),
                         ## -- Exit --
                         tabPanel(
                             actionButton("exitButton", "Exit")
                         )
                         ),
        function(input, output, session) {
            time.prey <- reactive({
                time.prey <- dat[dat$Time == input$time & dat$Stock  == input$stock, c(2, 6 : ncol(dat))]
            })

            ## Original Recruitment
            t.prey <- reactive({
                t.prey <- NULL
                foc    <- input$foc.fg
                stg    <- 1
                while(stg <= input$m.stg){
                    for(fg in foc){
                        ## this approach remove not common prey but keeps the proportion
                        diet     <- colSums(time.prey()[time.prey()$Predator == fg, 2 : ncol(time.prey())])
                        diet     <- diet[diet > 0] / sum(diet, na.rm = TRUE)
                        if(length(diet) == 0) next
                        diet     <- diet[diet >= input$min]
                        n.prey   <- names(diet)
                        m.prey   <- data.frame(Pred = fg, Prey = n.prey, value = diet, Stage = stg)
                        t.prey   <- rbind(t.prey, m.prey)
                    }
                    focp <- foc
                    foc  <- as.character(unique(t.prey$Prey[t.prey$Stage == stg]))
                    if(all(length(foc) == length(focp), all(foc %in% focp))) break
                    stg <- stg + 1
                }
                t.prey$TLprey <- prey.pos(t.prey$Prey, grp.dat)
                t.prey <- as.data.frame(t.prey)
            })
            ## assing the value of trophic level for the prey
            TL <- reactive({
                #browser()
                npred         <- unique(t.prey()$Pred)
                TL            <- NULL
                for(pred in npred){
                    pospred <- which(t.prey()$Pred %in% pred)
                    TL      <- c(TL, Tlevel(t.prey()$value[pospred], t.prey()$TLprey[pospred]))
                }
                pp.prey <- unique(t.prey()$Prey[-which(t.prey()$Prey %in% npred)])
                pp.prey <- data.frame(FG = pp.prey, Tlevel = prey.pos(pp.prey, grp.dat))
                TL      <- data.frame(FG = npred, Tlevel = TL)
                TL      <- rbind(TL, pp.prey)
                ## location on the plot
                his     <- hist(TL$Tlevel, breaks = c(1.5 : 6.5), plot = FALSE)
                v.lev   <- max(his$counts)
                brk     <- his$breaks
                h.lev   <- max(his$mids[his$counts > 0]) + 1 ## getting the top of the plor
                TL$v.lev <- v.lev
                TL$h.lev <- h.lev
                TL$vpos <- NA
                for(i in 1 : (length(brk) - 1)){
                    nfg.ly          <- which(TL$Tlevel > brk[i] & TL$Tlevel < brk[i + 1] )
                    tot.fg          <- length(nfg.ly)
                    vpos            <-  cumsum(rep(v.lev / tot.fg, tot.fg))  - (v.lev / tot.fg) * 0.5
                    TL$vpos[nfg.ly] <- vpos
                }
                TL <- as.data.frame(TL)
            })
            output$plot1 <- renderPlot({
                rad     <- 0.25
                plot(1, type = "n", xlab = '', ylab = "Trophic-level",
                     xlim = c(0, unique(TL()$v.lev)), ylim = c(1.5, unique(TL()$h.lev)), axes = FALSE)
                axis(2, at = c(1 : unique(TL()$h.lev)),  las = 1)
                for( i in 1 : nrow(t.prey())){
                    p.pred <- which(TL()$FG %in% t.prey()$Pred[i])
                    p.prey <- which(TL()$FG %in% t.prey()$Prey[i])
                    ## Direction of predation
                    col <- ifelse(TL()$Tlevel[p.pred] < TL()$Tlevel[p.prey], 'forestgreen', 'red1')
                    y      <- c(TL()$Tlevel[p.pred], TL()$Tlevel[p.prey] )
                    x      <- c(TL()$vpos[p.pred], TL()$vpos[p.prey] )
                    lines(x, y , lty = 1, lwd = (t.prey()$value[i] * 3), col = col)
                }
                ## circles
                for( i in 1 : nrow(TL())){
                    draw.circle(TL()$vpos[i], TL()$Tlevel[i], radius = rad, nv = 1500, border = NULL,
                                col = ifelse(i == 1, 'steelblue', 'gray91'), lty = 1, lwd = 1)
                    text(TL()$vpos[i], TL()$Tlevel[i], TL()$FG[i], cex = .8, font = 2, col = ifelse(i == 1, 'white', 1))
                }
                legend('topright', c('DownTop', 'TopDown'), lty = 1, col = c('forestgreen', 'red1'), lwd = 1.5, bty = 'n')

            })
            output$table <- renderTable({
                table    <- with(TL(), data.frame(Functional.group = FG, Trophic.Level = Tlevel))
            })
            observeEvent(input$exitButton, {
                stopApp()
            })
        }
    )
}


## functions
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Assing trophic level of the prey the values are based on Pauly et al 1998 [Diet composition and trophic levels of marine mammals]
##'        And Tucket and Rogers 2014 [Examining predator–prey body size, trophic level and body mass across marine and terrestrial mammals]
##' @param FGs Preys
##' @param grp.dat
##' @return
##' @author Demiurgo
prey.pos <- function(FGs, grp.dat){
    ctg <- vector('numeric')
    for( i in 1 : length(FGs)){
        ctg[i] <- which(grp.dat$Code %in% FGs[i], arr.ind = TRUE)
    }
    typ  <-  grp.dat$GroupType[ctg]
    TL.v <- ifelse(typ == 'FISH', 3.2,
            ifelse(typ %in% c('SED_EP_FF', 'PWN', 'REPTILE', 'FISH_INVERT', 'MOB_EP_OTHER', 'SED_EP_OTHER', 'SM_ZOO', 'MED_ZOO', 'LG_ZOO'), 2.5,
            ifelse(typ %in% c('CEP'), 3.2,
            ifelse(typ %in% c('SHARK', 'BIRD', 'MAMMAL'), 4,
            ifelse(typ %in% c('PL_BACT', 'SED_BACT', 'CARRION', 'LAB_DET', 'REF_DET', 'SM_PHY', 'LG_PHY', 'MICROPHYTOBENTHOS', 'DINOFLAG', 'PHYTOBEN', 'SEAGRASS', 'TURF', 'CORAL'), 2, 2)))))
    return(TL.v)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Throphic level estimation. Equation based on Pauly et al 1998 [Diet composition and trophic levels of marine mammals]
##' @param DC Diet composition
##' @param TLp Trophic level of the prey
##' @return The trophic level of the predator
##' @author Demiurgo
Tlevel <- function(DC, TLp){
    tl <- 1 + (sum(DC * TLp, na.rm = TRUE) / sum(DC, na.rm = TRUE))
    return(tl)
}
