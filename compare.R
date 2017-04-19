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
##' @return A shiny output (reactive html)
##' @author Demiurgo
output.cal <- function(folder1, folder2 = NULL, biomass.old.file, biomass.curr.file, groups.csv, diet.file ){
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
    mycol <- c(brewer.pal(8, "Dark2"), c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
    mycol <- colorRampPalette(mycol)
    shinyApp(
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
                                  ),
                         ## -- Exit --
                         tabPanel(
                             actionButton("exitButton", "Exit")
                         )
                         ),
        function(input, output, session) {
            pred.cohort <- reactive({
                if(hab.chk){
                    pred.new <- new.diet[new.diet$Time == input$Time & new.diet$Predator == input$FG2, ]
                } else {
                    pred.new <- new.diet[new.diet$Time == input$Time & new.diet$Predator == input$FG2 & new.diet$Stock == as.numeric(input$Stocks2), ]
                }
                pred.new <- pred.new[pred.new$value > input$Thr2, ]
            })

            predator <- reactive({
                if(hab.chk){
                    pred.new <- new.diet[new.diet$Predator == input$FG & new.diet$Habitat == input$Stocks, ]
                } else {
                    pred.new <- new.diet[new.diet$Predator == input$FG & new.diet$Stock == as.numeric(input$Stocks), ]
                }
                pred.new <- pred.new[pred.new$value > input$Thr, ]
            })
            observeEvent(input$exitButton, {
                stopApp()
            })
            prey <- reactive({
                out.diet <- new.diet[ new.diet$variable ==  input$FG & new.diet$value > input$Thr, ]
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
                ggplot(prey(), aes(x = Time, y = value, fill = Predator, width = 1)) + geom_bar(stat = "identity", position = 'fill') + scale_fill_manual(values = colorpp) +
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
        }
    )
}
