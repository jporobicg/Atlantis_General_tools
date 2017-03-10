out.compare <- function(folder1, folder2 = NULL, old.file, curr.file, groups.csv){
    library(shiny)
    library(tidyverse)
    library(reshape2)
    library(ggplot2)
    library("dplyr")
    library(stringr)
    mycol     <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    if(is.null(folder2)){
        if(file.exists(old.f)){
            old.dat <- read.csv(old.f, sep = ' ')
        } else if (file.exist(paste(folder1, old.f, sep = "/"))){
            old.dat <- read.csv(paste(folder1, old.f, sep = "/"))
        } else {
            stop('You need to provide a path to the old file')
        }
    }
    cur.dat <- read.csv(paste(folder1, curr.f, sep = "/"),sep = ' ')
    grp     <- read.csv(groups.csv)$Code
    sub.old <- cbind(old.dat[c('Time', as.character(grp))], Simulation = 'previous')
    sub.cur <- cbind(cur.dat[c('Time', as.character(grp))], Simulation = 'current')
    if(!all(sub.old$Time ==  sub.cur$Time)) stop('The time steps of the total time of the simulations are different')
    n.r <- nrow(sub.cur)
    names.fg <- names(sub.old)
    for( fg in 2 : (ncol(sub.old) - 1)){
        dif.sim      <- mean((sub.cur[c((n.r - 10) : n.r), fg]  / sub.old[c((n.r - 10) : n.r), fg] - 1)   * 100, na.rm = TRUE)
        names.fg[fg] <-  paste(names.fg[fg], '- Diff: ', round(dif.sim, 0), '%', sep = '')
    }
    dat.tot           <- rbind(sub.old, sub.cur)
    colnames(dat.tot) <- names.fg
    dat.tot           <- melt(dat.tot, id = c('Time', 'Simulation'))

    shinyApp(
        ui <- navbarPage("Compare outputs",
                         tabPanel('Structural nitrogen',#
                                  plotOutput('plot1', width = "100%", height = "1000px")
                                  )
                         ),
        function(input, output, session) {
            output$plot1 <- renderPlot({
                plot <- ggplot(dat.tot, aes(x = Time, y = value, colour = Simulation)) +
                    geom_line() + facet_wrap(~ variable, ncol = 4,  scale = 'free_y') + theme_bw()+
                    scale_color_manual(values = mycol[c(2, 6)])
                plot <- update_labels(plot, list(x = 'Time step', y = 'Biomass (tons)'))
                plot
            })
        }
    )
}
