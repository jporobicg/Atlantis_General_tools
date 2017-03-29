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
    mycol     <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    if(is.null(folder2)){
        if(file.exists(biomass.old.file)){
            old.dat <- read.csv(biomass.old.file, sep = ' ')
        } else if (file.exist(paste(folder1, biomass.old.file, sep = "/"))){
            old.dat <- read.csv(paste(folder1, biomass.old.file, sep = "/"))
        } else {
            stop('You need to provide a path to the old file')
        }
    }
    cur.dat <- read.csv(paste(folder1, biomass.curr.file, sep = "/"),sep = ' ')
    grp     <- read.csv(groups.csv)$Code
    sub.old <- cbind(old.dat[c('Time', as.character(grp))], Simulation = 'previous')
    sub.cur <- cbind(cur.dat[c('Time', as.character(grp))], Simulation = 'current')
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
    new.diet  <- melt(diet.file, id.vars = c("Time", "Predator", "Cohort", "Stock"))
    time      <- unique(diet.file$Time)
    stocks    <- unique(diet.file$Stock)
    shinyApp(
        ui <- navbarPage("Compare outputs",
                         tabPanel('Biomass comparison',#
                                  plotOutput('plot1', width = "100%", height = "1000px")
                                  ),
                         tabPanel('Predation',#
                                  fluidRow(
                                      column(2,
                                             wellPanel(
                                                 selectInput('FG', 'Functional Group :', as.character(grp)),
                                                 selectInput('Stocks', 'Stocks :', stocks),
                                                 sliderInput("Time", "Simulation Time :", min = min(time),  max = max(time), value = min(time), step = diff(time)[1])
                                             )
                                             ),
                                      column(10,
                                             plotOutput('plot2', width = "100%", height = "1000px"),
                                             plotOutput('plot3', width = "100%", height = "1000px")
                                             )
                                  )
                                 )
                         ),
        function(input, output, session) {
            predator <- reactive({
                pred.pos      <- which(diet.file$Time == input$Time & diet.file$Predator == input$FG & diet.file$Stock == input$Stocks)
                pred.dat.filt <- diet.file[pred.pos, ]
                ## removing pry that is less than 1% of the total diet
                rm.vec <- 1
                for(rem in 5 : ncol(pred.dat.filt)){
                    if(all(pred.dat.filt[, rem] < 0.001)) rm.vec <- c(rm.vec, rem)
                }
                rm.vec <- rm.vec[-1]
                pred.dat.filt <- pred.dat.filt[, - rm.vec]
                plt <- melt(pred.dat.filt,id.vars=c('Time', 'Predator', 'Cohort', 'Stock'))
            })
            prey <- reactive({
                out.diet <- filter(new.diet, variable ==  input$FG, value > 0.0001)
            })
            output$plot1 <- renderPlot({
                plot <- ggplot(dat.tot, aes(x = Time, y = value, colour = Simulation)) +
                    geom_line() + facet_wrap(~ variable, ncol = 4,  scale = 'free_y') + theme_bw()+
                    scale_color_manual(values = mycol[c(2, 6)])
                plot <- update_labels(plot, list(x = 'Time step', y = 'Biomass (tons)'))
                plot
            })
            output$plot2 <- renderPlot({
              ggplot(predator(), aes(x = Cohort, y = value, fill = variable)) + geom_bar(stat = "identity") + scale_fill_brewer(palette = "Dark2")
            })
            output$plot3 <- renderPlot({
                ggplot(prey(), aes(x = Time, y = value, fill = Predator)) + geom_bar(stat = "identity", position = 'fill') + scale_fill_brewer(palette = "Paired")
            })
        }
    )
}
