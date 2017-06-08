recruitment.cal <- function(ini.nc.file, out.nc.file, yoy.file, grp.file, prm.file,  quiet = TRUE){
    txtHelp <- "<h2>Summary</h2>"
    txtHelp <- paste(txtHelp, "<p>This code Helps to calibrate the recruitment for <b>Atlantis</b> run. </p>")
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
    if (!require('reshape', quietly = TRUE)) {
        stop('The package reshape was not installed')
    }
    if (!require('tidyverse', quietly = TRUE)) {
        stop('The package tidyverse was not installed')
    }
    if (!require('stringr', quietly = TRUE)) {
        stop('The package stringr was not installed')
    }
    if(!quiet) cat('      ...Done!')
    ## Reading files
    if(!quiet) cat('\n Reading files')
    nc.ini    <- nc_open(ini.nc.file)
    yoy       <- read.csv(yoy.file, sep = ' ')
    group.csv <- read.csv(grp.file)
    nc.out    <- nc_open(out.nc.file)
    prm       <- readLines(prm.file, warn = FALSE)
    mg2t      <-  0.00000002 ## mg C converted to wet weight in tonnes == 20 / 1000000000
    if(!quiet) cat('      ...Done!')
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 2    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n Processing')
    sp.dat    <- with(group.csv, which(IsTurnedOn == 1 & NumCohorts > 1)) ## Age structure groups
    options(warn =  - 1)
    rec       <- text2num(prm, 'flagrecruit', FG = 'look') ## Avoinding the annoying warnings
    options(warn =  0)
    rec       <- rec[complete.cases(rec), ] ## avoiding NAs
    rec       <- cbind(rec, Alpha = NA, Beta = NA, KSPA = NA, FSP = NA, Time.sp  = NA, Period.sp = NA,
                       Time.rec = NA, Period.rec = NA, XCN = text2num(prm, 'X_CN', FG = 'look')[1, 2],
                       XRS = text2num(prm, 'X_RS', FG = 'look')[1, 2])
    sps    <- gsub(pattern = 'flagrecruit', '', rec$FG)
    rec$FG <- sps
    FSPB   <- NULL
    t = 1 ## counter
    if(!quiet) cat('\n Reading parameters')
    for(fg.r in 1 : length(sps)){
        if(!sps[fg.r] %in% group.csv$Code[sp.dat]) next()
        if(rec$Value[fg.r] == 3){  ## Beverton Holt Recruitment
            rec$Alpha[fg.r]     <- text2num(prm, paste0('BHalpha_', sps[fg.r]), FG = 'look')[1, 2]
            rec$Beta[fg.r]      <- text2num(prm, paste0('BHbeta_', sps[fg.r]), FG = 'look')[1, 2]
            rec$KSPA[fg.r]      <- text2num(prm, paste0('KSPA_', sps[fg.r]), FG = 'look')[1, 2]
            rec$FSP[fg.r]       <- text2num(prm, paste0('FSP_', sps[fg.r]), FG = 'look')[1, 2]
            ## Proportion of mature at age
            fspb.tmp            <- text2num(prm, paste0('FSPB_', sps[fg.r]), FG = sps[fg.r], Vector = TRUE)
            len                 <- ifelse(t == 1, length(fspb.tmp), ncol(FSPB))
            fspb.tmp            <- c(fspb.tmp, rep(0, len - length(fspb.tmp)))
            FSPB                <- rbind(FSPB, fspb.tmp)
            rownames(FSPB)[t]   <- as.character(sps[fg.r])
            t                   <- t + 1
        } else if(rec$Value[fg.r] == 1){ ## constant recruitment
            rec$Alpha[fg.r] <- text2num(prm, paste0('KDENR_', sps[fg.r]), FG = 'look', Vector = TRUE)[1]
        }
        rec$Time.sp[fg.r]   <- text2num(prm, paste0(sps[fg.r], '_Time_Spawn'), FG = 'look')[1, 2]
        rec$Period.sp[fg.r] <- text2num(prm, paste0(sps[fg.r], '_spawn_period'), FG = 'look')[1, 2]
        rec$Time.rec[fg.r]  <- text2num(prm, paste0(sps[fg.r], '_Recruit_Time'), FG = 'look')[1, 2]
        rec$Period.rec[fg.r]<- text2num(prm, paste0('Recruit_Period_', sps[fg.r]), FG = 'look')[1, 2]
        rec$Rec.SNW[fg.r]   <- text2num(prm, paste0('KWSR_', sps[fg.r]), FG = 'look')[1, 2]
        rec$Rec.RNW[fg.r]   <- text2num(prm, paste0('KWRR_', sps[fg.r]), FG = 'look')[1, 2]
    }
    if(!quiet) cat('      ...Done!')
    if(!quiet) cat('\n Reading YOY from Atlantis')
    ##~~~~~~~~~~~~~~~~~~~~~~~~~##
    ##    YOY file array       ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~##
    cod.yoy <- data.frame(Code = paste0(group.csv$Code[sp.dat], '.0'), Initial = NA)
    nam.fg  <- group.csv$Name[sp.dat]
    yoy.tmp <- yoy * 0
    for(c in cod.yoy$Code){
        yoy.tmp[, c]   <- (yoy[, c] / yoy[1, c])
    }
    f.yoy <- cbind(Time = yoy$Time, yoy.tmp[, which(names(yoy.tmp) %in% cod.yoy$Code)])
    if(!quiet) cat('      ...Done!')
    if(!quiet) cat('\n Calculating recruits')
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ## Number and weight of individual at age in each reproduction perior  ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    coh.fg  <- group.csv$NumCohorts[sp.dat]
    cod.fg  <- group.csv$Code[sp.dat]
    time    <- ncvar_get(nc.out, 't')
    spw     <- NULL
    nam     <- NULL
    SSB.tot <- NULL
    fg.spw  <- list()
    fg.ssb  <- list()
    loo     <- 1
    ## estimation of Sp by reproduction event and By FG
    for( fg in 1 : length(nam.fg)){ ## loop through functional groups
        pos.fspb <- which(rownames(FSPB) %in% cod.fg[fg])
        if(length(pos.fspb) == 0) next()
        fg.row   <- which(rec$FG %in% cod.fg[fg])
        xrs      <- rec$XRS[fg.row]
        FSP      <- rec$FSP[fg.row]
        KSPA     <- rec$KSPA[fg.row]
        sp.tmp   <- NULL
        SSB.fg   <- NULL
        ## time of spawning
        time.stp <- seq(from = 0, by = 365, to = length(time))  + rec$Time.sp[fg.row]
        time.stp <- time.stp[time.stp < length(time)]
        spw.coh  <- list()
        ssb.coh  <- list()
        for(coh in 1 : coh.fg[fg]){
            name.fg <- paste0(nam.fg[fg], coh)
            nums    <- ncvar_get(nc.out, paste0(name.fg, '_Nums'))[, , time.stp]
            mask    <- ifelse(nums > 1.e-16, 1, 0)
            nums    <- nums * mask
            rn      <- ncvar_get(nc.out, paste0(name.fg, '_ResN'))[, , time.stp] * mask
            sn      <- ncvar_get(nc.out, paste0(name.fg, '_StructN'))[, , time.stp] * mask
            wspi    <- sn * (1 + xrs)       ## minimum weigth for spawning
            rat     <- ((rn  + sn ) - wspi) ## weight deficit
            rat[(rat < 0)]   <- 0
            spawn            <- ((wspi * FSP - KSPA)  -  rat)
            SSB.tmp <- (ncvar_get(nc.out, paste0(name.fg, '_ResN'))[, , time.stp]  +
                        ncvar_get(nc.out, paste0(name.fg, '_StructN'))[, , time.stp])  *
                ncvar_get(nc.out, paste0(name.fg, '_Nums'))[, , time.stp]
            SSB.tmp <- SSB.tmp * mask
            SSB.tmp <- SSB.tmp ##* FSPB[pos.fspb, coh]
            spawn[spawn < 0] <- 0
            spawn <- spawn *  FSPB[pos.fspb, coh] ## individual spawn
            spawn <- spawn * nums      ## total spawn
            spw.coh[[coh]] <- spawn    ## spawning by time step and cohort
            ssb.coh[[coh]] <- SSB.tmp ## Biomass by time step and cohort
            if(length(dim(spawn)) > 2){
                spawn   <- apply(spawn, 3, sum, na.rm = TRUE)
                SSB.tmp <- apply(SSB.tmp, 3, sum, na.rm = TRUE)
            }
            sp.tmp      <- rbind(sp.tmp, spawn)   ## Spawning by functional group and Age class
            SSB.fg      <- rbind(SSB.fg, SSB.tmp) ## Spawning Stock by functional group and Age class
        }
        fg.spw[[loo]] <- spw.coh
        fg.ssb[[loo]] <- ssb.coh
        loo     <- loo + 1
        nam     <- c(nam, as.character(cod.fg[fg]))
        fin.sp  <- colSums(sp.tmp)
        SSB.fg  <- colSums(SSB.fg)
        ## all the estimation will have the same length
        if(length(spw) > 0 && length(fin.sp) != ncol(spw)) fin.sp         <- fin.sp[seq(ncol(spw))]
        if(length(SSB.tot) > 0 && length(SSB.fg) != ncol(SSB.tot)) SSB.fg <- SSB.fg[seq(ncol(SSB.tot))]
        SSB.tot <- rbind(SSB.tot, SSB.fg)
        spw     <- rbind(spw, fin.sp)
    }
    ## FG Names
    rownames(spw)     <- nam
    rownames(SSB.tot) <- nam
    if(!quiet) cat('      ...Done!')
    if(!quiet) cat('\n\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n # -     Step 3    -   #')
    if(!quiet) cat('\n # -  -  -  -  -  -  - #')
    if(!quiet) cat('\n\n Plotting \n\n')
    shinyApp(
        ui <- navbarPage('Atlantis Recruitment Tool',
                         tabPanel('Recruits and YOY',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 tags$h3('Functional Group'),
                                                 selectInput('sp', 'Functional Group', as.character(cod.fg)),
                                                 mainPanel(strong("Rec model:"), textOutput("Rec.mod")),
                                                 mainPanel(strong("Alpha: "), verbatimTextOutput("Alpha.mod", placeholder = TRUE)),
                                                 mainPanel(strong("Beta: "), verbatimTextOutput("Beta.mod", placeholder = TRUE)),
                                                 mainPanel("Initial YOY: ", textOutput("Ini.YOY")),
                                                 br(),
                                                 numericInput("new.alpha", label = "New Alpha", value = 0),
                                                 br(),
                                                 numericInput("new.beta", label = "New Beta", value = 0)
                                             )
                                             ),
                                      column(10,
                                             plotOutput('plot1', width = "100%", height = "400px"),
                                             plotOutput('plot2', width = "100%", height = "400px"),
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
        time.stp <- reactive({
            time.stp <- seq(from = 0, by = 365, to = length(time))  + rec$Time.sp[rec$FG == input$sp]
            time.stp <- time.stp[time.stp < length(time)]
        })
        ## Original Recruitment
        rec.bio <- reactive({
            spawn.fg <- spw[input$sp, ]
            biom.fg  <- SSB.tot[input$sp, ]
            sp.plt   <- paste0(input$sp, '.0')
            recruit  <- unlist(BH.rec(spawn.fg, rec$Alpha[rec$FG == input$sp], rec$Beta[rec$FG == input$sp], biom.fg))
            new.rec  <- unlist(BH.rec(spawn.fg, input$new.alpha, input$new.beta, biom.fg))
            rec.bio  <- recruit * (rec$Rec.SNW[rec$FG == input$sp] + rec$Rec.RNW[rec$FG == input$sp]) * rec$XCN[rec$FG == input$sp] * mg2t
            new.bio  <- new.rec * (rec$Rec.SNW[rec$FG == input$sp] + rec$Rec.RNW[rec$FG == input$sp]) * rec$XCN[rec$FG == input$sp] * mg2t
            yoy.fg   <- data.frame(Time = yoy$Time, Rec = yoy[, which(names(yoy) == sp.plt)])
            dif      <- which(diff(yoy.fg$Rec) != 0) + 1
            n.yoy    <- yoy.fg[dif, 2]
            Time.yoy <- yoy.fg[dif, 1]
            n.f.yoy  <- f.yoy[dif, c(1, which(names(f.yoy) == sp.plt))]
            prop.dif <- yoy.fg$Rec[dif] / rec.bio
            fst.val  <- (new.bio * prop.dif) /  yoy.fg$Rec[1]
            N.YOY    <- (new.bio * prop.dif)
            df.end   <- data.frame(Rec = rec.bio, N.YOY = N.YOY, N.Rec = new.bio, BYOY = n.yoy, TYOY = Time.yoy, PrpYOY = n.f.yoy[, 2], Prp.st = fst.val, P.diff = prop.dif)
        })
        output$Rec.mod <- renderText({
            mod   <- rec$Value[rec$FG == input$sp]
            model <- ifelse(mod == 1, 'Constant	recruitment',
                     ifelse(mod == 2, 'Determined by chlA',
                     ifelse(mod == 3, 'Beverton-Holt',
                     ifelse(mod == 4, 'Random Lognormal', 'Other'))))
            })
        observeEvent(input$exitButton, {
            stopApp()
        })
        output$Alpha.mod <- renderText({
            Alpha   <- rec$Alpha[rec$FG == input$sp]
        })
        output$Beta.mod <- renderText({
            Beta   <- rec$Beta[rec$FG == input$sp]
        })
        output$Ini.YOY <- renderText({
            sp.plt   <- paste0(input$sp, '.0')
            Ini.YOY  <- yoy[1, which(names(yoy) == sp.plt)]
        })
        output$plot1 <- renderPlot({
            par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
            plot(rec.bio()$TYOY, rec.bio()$BYOY , xlab = 'Time (days)', ylab = 'Biomass [Tonnes]', las = 1, bty = 'n', pch = 20,type = 'b',
                 col = 'royalblue', ylim = c(0, max(rec.bio()$Rec, rec.bio()$BYOY, rec.bio()$N.Rec)),
                 xlim = range(c(rec.bio()$TYOY, time.stp())))
            lines(time.stp(), rec.bio()$Rec, type = 'b', pch = 20, col = 'red4')
            lines(time.stp(), rec.bio()$N.Rec, type = 'b', pch = 20, col = 'green4')
            lines(time.stp(), rec.bio()$N.YOY, type = 'b', pch = 20, col = 'yellowgreen')
            legend("topright", inset = c(-0.1, 0), legend = c('Atlantis YOY', 'Larvaes', 'New Larvaes', 'New YOY'),
                   lty=c(1, 1, 1, 1), col=c('royalblue', 'red4', 'green4', 'yellowgreen'), pch = c(20, 20, 20, 20), bty = 'n', lwd = 2)
        })
        output$plot2 <- renderPlot({
            par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
            plot(rec.bio()$TYOY, rec.bio()$PrpYOY, ylim = c(0, ifelse(max(rec.bio()$PrpYOY) < 1, 1, max(rec.bio()$PrpYOY))), bty = 'n', type = 'b',
                 lty = 2, pch = 19, col = 'olivedrab4', ylab = 'proportion from initial yoy (t)', xlab = 'Time (days)', las = 1,
                 xlim = range(c(rec.bio()$TYOY, time.stp())))
            lines(time.stp(), rec.bio()$Prp.st, type = 'b', pch = 19, lty = 2, col = 'yellow3')
            legend('topright', inset=c(-0.1, 0), legend = c('YOY prop', 'New prop'), lty = 2, col=c('olivedrab4', 'yellow3' ), pch = c(19, 19), bty = 'n', lwd = 2)

        })
        output$table <- renderTable({
            table <- with(rec.bio(), data.frame(Time.Larv = time.stp(), TimeYOY = TYOY, Larvaes.Atlantis = Rec, YOY.Atlantis = BYOY, Diff.Prop = P.diff * 100, Est.Larvaes = N.Rec, Est.YOY = N.YOY))
        })
    }
    )
}
    ##' .. content for \description{} (no empty lines) ..
    ##'
    ##' .. content for \details{} ..
    ##' @title Beverton Equation
    ##' @param sp spawning power
    ##' @param bha Alpha parameter
    ##' @param bhb Beta parameter
    ##' @param bio Biomass
    ##' @return The amout of recruit (Larvaes)
    ##' @author Demiurgo
    BH.rec <- function(sp, bha, bhb, bio){
        num <- lapply(sp,  '*',  bha)
        den <- lapply(bio,  '+',  bhb)
        recruit <- mapply('/', num,  den, SIMPLIFY = FALSE)
        return(recruit)
    }



    ## functions
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
