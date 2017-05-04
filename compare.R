##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param folder1 Folder that contain the current Atlantis outputs
##' @param folder2 Folder that contain the old Atlantis outputs files. If empty the folder1 would be use
##' @param biomass.old.file Old output file
##' @param biomass.curr.file Current file
##' @param groups.csv group csv file
##' @param diet.file output file for the diet (atlantis output)
##' @param age.biomass output file for the biomass by age (non-standard atlantis output, make sure to put the flag 'flag_age_output' in the run.prm if you want to look at this output)
##' @return A shiny output (reactive html)
##' @author Demiurgo
output.cal <- function(folder1, folder2 = NULL, biomass.old.file, biomass.curr.file, groups.csv, diet.file, age.biomass = NULL ){
    ## Libraries
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
    library(RColorBrewer)
    ## Colours
    if(is.null(folder2)){
        if(file.exists(biomass.old.file)){
            old.dat <- read.csv(biomass.old.file, sep = ' ')
        } else if (file.exists(paste(folder1, biomass.old.file, sep = "/"))){
            old.dat <- read.csv(paste(folder1, biomass.old.file, sep = "/"), sep = ' ')
        } else {
            stop('You need to provide a path to the old file')
        }
    } else {
        old.dat <- read.csv(paste(folder2, biomass.old.file, sep = "/"), sep = ' ')
    }
    cur.dat <- read.csv(paste(folder1, biomass.curr.file, sep = "/"),sep = ' ')
    grp     <- read.csv(groups.csv)
    grp     <- grp[grp$IsTurnedOn == 1, ]$Code
    sub.old <- cbind(old.dat[c('Time', as.character(grp))], Simulation = 'previous')
    sub.cur <- cbind(cur.dat[c('Time', as.character(grp))], Simulation = 'current')
    sp.name <- c("Time", paste0('Rel', grp),"PelDemRatio", "PiscivPlankRatio")
    rel.bio <- melt(cur.dat[, sp.name], id.vars = 'Time')
    if(!all(sub.old$Time ==  sub.cur$Time)) stop('The time steps of the total time of the simulations are different')
    n.r      <- nrow(sub.cur)
    names.fg <- names(sub.old)
    for( fg in 2 : (ncol(sub.old) - 1)){
        dif.sim      <- mean((sub.cur[c((n.r - 10) : n.r), fg]  / sub.old[c((n.r - 10) : n.r), fg] - 1) * 100, na.rm = TRUE)
        names.fg[fg] <-  paste(names.fg[fg], '- Diff: ', round(dif.sim, 0), '%', sep = '')
    }
    dat.tot           <- rbind(sub.old, sub.cur)
    colnames(dat.tot) <- names.fg
    dat.tot           <- melt(dat.tot, id = c('Time', 'Simulation'))
    ## Diet Analysis
    diet.file <- read.csv(diet.file, sep=' ')
    hab.chk <- FALSE
    if(any('Habitat' == colnames(diet.file))) hab.chk <- TRUE
    ##browser()
    if(hab.chk){
        if(any('Updated' == colnames(diet.file))){
            new.diet  <- melt(diet.file, id.vars = c("Time", "Predator", "Habitat", "Updated"))
            rem       <- which(colnames(diet.file) == 'Updated')
            new.diet  <- new.diet[, - rem]
        } else {
            new.diet  <- melt(diet.file, id.vars = c("Time", "Predator", "Habitat"))
        }
        time      <- unique(diet.file$Time)
        stocks    <- unique(diet.file$Habitat)
    } else {
        if(any('Updated' == colnames(diet.file))){
            new.diet  <- melt(diet.file, id.vars = c("Time", "Predator", "Cohort", "Stock", "Updated"))
            rem       <- which(colnames(diet.file) == 'Updated')
            new.diet  <- new.diet[, - rem]
        } else {
            new.diet  <- melt(diet.file, id.vars = c("Time", "Predator", "Cohort", "Stock"))
        }
        time      <- unique(diet.file$Time)
        stocks    <- unique(diet.file$Stock)
    }
    new.bio           <- melt(sub.cur, id = c('Time', 'Simulation')) ## current Biomass of the functional groups
    new.bio           <- new.bio[new.bio$variable %in%  as.character(unique(new.diet$Predator)), ] ## only take the groups which have values in the predator diet file
    colnames(new.bio) <- c('Time', 'Simulation', 'variable', 'Biomass') ## This is the current Biomass
    new.bio           <- left_join(new.diet, new.bio, by = c('Predator' = 'variable', 'Time')) ## join the values for the fraction of each prey eaten by each predator
    new.bio$eff.pred  <- new.bio$Biomass * new.bio$value ## Biomass*value (which is the fraction of each prey eaten) = effective predation
    new.bio           <- new.bio[, c('Time', 'Predator', 'variable', 'eff.pred' )]
    mycol             <- c(brewer.pal(8, "Dark2"), c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
    mycol             <- colorRampPalette(mycol)
    ## Read in and prepare the Biomass by cohort file
    if(!is.null(age.biomass)){
        cohort.file      <- read.csv(paste(folder1, age.biomass, sep = "/"), sep = ' ')
        new.cohort       <- melt(cohort.file, id.vars = c("Time")) %>%
            separate(., variable, into = c("variable","Agcl"), sep = "[.]")
        new.fract.cohort <- new.cohort %>%
            .[.$variable %in%  as.character(unique(new.diet$Predator)), ] ## only take the groups which have values in the predator diet file
        colnames(new.fract.cohort) <- c("Time", "variable","Agcl", "Biomass")
        new.fract.cohort           <- left_join(new.diet, new.fract.cohort, by = c('Predator' = 'variable', 'Time')) ## join the values for the fraction of each prey eaten by each predator
        new.fract.cohort$frac.Biom <- new.fract.cohort$Biomass * new.fract.cohort$value ## Biomass*value (which is the fraction of each prey eaten) = predated biomass
    }
    ## Start the Shiny application
    shinyApp(
        ## Create the different tabs
        ui <- navbarPage("Compare outputs",
                         tabPanel('Biomass comparison',
                                  tabsetPanel(
                                      tabPanel('Total Biomass',
                                               plotOutput('plot1', width = "100%", height = "1000px")
                                               ),
                                      tabPanel('Relative Biomass',
                                               plotOutput('plot1B', width = "100%", height = "1000px")
                                               )
                                  )
                                  ),
                         if(!is.null(age.biomass)){
                             tabPanel('Biomass - by Cohort',
                                      fluidRow(
                                          column(2,
                                                 wellPanel(
                                                     selectInput('FG3', 'Functional Group :', as.character(grp)),
                                                     checkboxInput(inputId = "Free_scale",
                                                                   label = strong("Show with fixed scales"),
                                                                   value = FALSE)
                                                 )
                                                 ) ,
                                          tabPanel('Biomass by cohort',
                                                   plotOutput('plot5', width = "100%", height = "700px")
                                                   )
                                      )   )},
                         tabPanel('Predation - over Time',
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 selectInput('FG', 'Functional Group :', as.character(grp)),
                                                 if(hab.chk){
                                                     selectInput('Stocks', 'Habitat :', stocks)
                                                 } else {
                                                     selectInput('Stocks', 'Stocks :', stocks)
                                                 },
                                                 numericInput("Thr", "Threshold :", min = 1e-16,  max = 1, value = 1e-4, step = 0.001)
                                             )
                                             ),
                                      column(10,
                                             plotOutput('plot2', width = "100%", height = "400px"),
                                             plotOutput('plot3', width = "100%", height = "400px")
                                             )
                                  )
                                  ),
                         if(is.null(age.biomass)){
                             tabPanel('Predation - by Cohort',
                                      fluidRow(
                                          column(2,
                                                 wellPanel(
                                                     selectInput('FG2', 'Functional Group :', as.character(grp)),
                                                     if(!hab.chk) selectInput('Stocks2', 'Stocks :', stocks),
                                                     sliderInput("Time", "Simulation Time :", min = min(time),  max = max(time), value = min(time), step = diff(time)[1]),
                                                     numericInput("Thr2", "Threshold :", min = 1e-16,  max = 1, value = 1e-4, step = 0.001)
                                                 )
                                                 ),
                                          column(10,
                                                 plotOutput('plot4', width = "100%", height = "600px")
                                                 )
                                      )
                                      )},
                         if(!is.null(age.biomass)){
                             tabPanel('Predation - by Cohort',
                                      fluidRow(
                                          column(2,
                                                 wellPanel(
                                                     selectInput('FG4', 'Functional Group :', as.character(grp)),
                                                     if(!hab.chk) selectInput('Stocks4', 'Stocks :', stocks),
                                                     sliderInput("Time", "Simulation Time :", min = min(time),  max = max(time), value = min(time), step = diff(time)[1]),
                                                     numericInput("Thr4", "Threshold :", min = 1e-16,  max = 1, value = 1e-4, step = 0.001)
                                                 )
                                                 ),
                                          column(10,
                                                 plotOutput('plot4', width = "100%", height = "400px"),
                                                 plotOutput('plot6', width = "100%", height = "400px")
                                                 )
                                      )
                                      )
                         },
                         ## -- Exit --
                         tabPanel(
                             actionButton("exitButton", "Exit")
                         )
                         ),
        ## Link the input for the different tabs with your original data
        ## Create the plots
        function(input, output, session) {

            if(!is.null(age.biomass)){
                pred.cohort <- reactive({
                    if(hab.chk){
                        pred.new <- new.diet[new.diet$Time == input$Time & new.diet$Predator == input$FG4, ]
                    } else {
                        pred.new <- new.diet[new.diet$Time == input$Time & new.diet$Predator == input$FG4 & new.diet$Stock == as.numeric(input$Stocks4), ]
                    }
                    pred.new <- pred.new[pred.new$value > input$Thr4, ]
                })
            } else {
                pred.cohort <- reactive({
                    if(hab.chk){
                        pred.new <- new.diet[new.diet$Time == input$Time & new.diet$Predator == input$FG2, ]
                    } else {
                        pred.new <- new.diet[new.diet$Time == input$Time & new.diet$Predator == input$FG2 & new.diet$Stock == as.numeric(input$Stocks2), ]
                    }
                    pred.new <- pred.new[pred.new$value > input$Thr2, ]
                })
            }

            predator <- reactive({
                if(hab.chk){
                    pred.new <- new.diet[new.diet$Predator == input$FG & new.diet$Habitat == input$Stocks, ]
                } else {
                    pred.new <- new.diet[new.diet$Predator == input$FG & new.diet$Stock == as.numeric(input$Stocks), ]
                }
                pred.new <- pred.new[pred.new$value > input$Thr, ]
            })
            if(!is.null(age.biomass)){
                cohort <- reactive({
                                        #cohort.new <- new.cohort[new.cohort$variable == input$FG & new.cohort$Stock == as.numeric(input$Stocks), ]
                    cohort.new <- new.cohort[new.cohort$variable == input$FG3, ]
                                        #cohort.new <- cohort.new[cohort.new$value > input$Thr3, ]
                })
                fract.cohort <- reactive({
                                        #fract.cohort.new <- new.fract.cohort[new.fract.cohort$variable == input$FG4 & new.fract.cohort$Stock4 == as.numeric(input$Stocks), ]
                    fract.cohort.new <- new.fract.cohort[new.fract.cohort$Time == input$Time & new.fract.cohort$Predator == input$FG4, ]
                    fract.cohort.new <- fract.cohort.new[fract.cohort.new$value > input$Thr4, ]
                })
            }

            observeEvent(input$exitButton, {
                stopApp()
            })
            prey <- reactive({
                out.diet <- new.bio[new.bio$variable ==  input$FG, ]
                trh.max  <- max(out.diet$eff.pred) * input$Thr
                out.diet <- out.diet[out.diet$eff.pred > trh.max, ]
            })
            output$plot1 <- renderPlot({
                plot <- ggplot(dat.tot, aes(x = Time, y = value, colour = Simulation)) +
                    geom_line() + facet_wrap(~ variable, ncol = 4,  scale = 'free_y') + theme_bw()+
                    scale_color_manual(values = mycol(2))
                plot <- update_labels(plot, list(x = 'Time step', y = 'Biomass (tons)'))
                plot
            })
            output$plot1B <- renderPlot({
                plot <- ggplot(rel.bio, aes(x = Time, y = value)) +
                    geom_line(colour = 'firebrick3') + facet_wrap( ~ variable, ncol = 4) +
                    theme_bw() + ylim(0, 2) +
                    annotate('rect', xmin =  - Inf, xmax = Inf, ymax = 1.5, ymin = 0.5, alpha = .1, colour = 'royalblue', fill = 'royalblue')
                plot <- update_labels(plot, list(x = 'Time step', y = 'Relative Biomass (Bt/B0)'))
                plot
            })
            output$plot2 <- renderPlot({
                colorpp  <- mycol(length(unique(predator()$variable)))
                ggplot(predator(), aes(x = Time, y = value, fill = variable, width = 1)) + geom_bar(stat = "identity", position = 'fill') + scale_fill_manual(values = colorpp, name = 'Prey') +
                    labs(list(title = paste('Predator -', input$FG), x = 'Time step', y = 'Proportion', colour = 'Prey'))
            })
            output$plot3 <- renderPlot({
                colorpp <- mycol(length(unique(prey()$Predator)))
                ggplot(prey(), aes(x = Time, y = eff.pred, fill = Predator, width = 1)) + geom_bar(stat = "identity", position = 'fill') + scale_fill_manual(values = colorpp) +
                    labs(list(title = paste('Prey -', input$FG), x = 'Time step', y = 'Proportion'))
            })
            output$plot4 <- renderPlot({
                colorpp <- mycol(length(unique(pred.cohort()$variable)))
                if(ncol(pred.cohort()) == 4){
                    df <- data.frame()
                    ggplot(df) + geom_bar() + labs(list(title = paste('Predator  -', input$FG, 'on Time step :', input$Time), x = 'Cohort', y = 'Proportion'))
                } else if(hab.chk){
                    ggplot(pred.cohort(), aes(x = Habitat, y = value, fill = variable, width = .75)) + geom_bar(stat = "identity", position = 'fill') + scale_fill_manual(values = colorpp, name = 'Prey') +
                        labs(list(title = paste('Predator  -', input$FG, 'on Time step :', input$Time), x = 'Habitat', y = 'Proportion'))
                } else{
                    ggplot(pred.cohort(), aes(x = Cohort, y = value, fill = variable, width = .75)) + geom_bar(stat = "identity", position = 'fill') + scale_fill_manual(values = colorpp) +
                        labs(list(title = paste('Predator  -', input$FG, 'on Time step :', input$Time), x = 'Cohort', y = 'Proportion'))
                }
            })
            if(!is.null(age.biomass)){
                output$plot5 <- renderPlot({
                    if(!input$Free_scale){
                        ggplot(cohort(), aes(x = Time, y = value )) +
                            geom_line() + facet_wrap(~ Agcl, ncol = 4,  scale = 'free_y') + theme_bw() +
                            labs(list(x = 'Time step', y = 'Biomass (tons)'))
                    }
                    else {
                        ggplot(cohort(), aes(x = Time, y = value )) +
                            geom_line() + facet_wrap(~ Agcl, ncol = 4) + theme_bw() +
                            labs(list(x = 'Time step', y = 'Biomass (tons)'))
                    }
                })

                output$plot6 <- renderPlot({
                    colorpp <- mycol(length(unique(fract.cohort()$variable)))
                                        # if(ncol(fract.cohort()) == 4){
                                        #   df <- data.frame()
                                        #   ggplot(df) + geom_bar() + labs(list(title = paste('Predator  -', input$FG, 'on Time step :', input$Time), x = 'Cohort', y = 'Proportion'))
                                        # } else if(hab.chk){
                                        #   ggplot(pred.cohort(), aes(x = Habitat, y = value, fill = variable, width = .75)) + geom_bar(stat = "identity", position = 'fill') + scale_fill_manual(values = colorpp, name = 'Prey') +
                                        #     labs(list(title = paste('Predator  -', input$FG, 'on Time step :', input$Time), x = 'Habitat', y = 'Proportion'))
                                        # } else{
                                        #   #ggplot(fract.cohort(), aes(x = Cohort, y = frac.Biom, fill = variable, width = .75)) + geom_bar(stat = "identity", position = 'fill') + scale_fill_manual(values = colorpp) +
                                        # labs(list(title = paste('Predator  -', input$FG, 'on Time step :', input$Time), x = 'Cohort', y = 'Proportion'))

                    ggplot(fract.cohort(), aes(x = Agcl, y = frac.Biom, fill = variable, width = .75)) + geom_bar(stat = "identity") +
                        scale_fill_manual(values = colorpp, name = 'Prey')

                                        #}
                })


            }

        }
    )
}
