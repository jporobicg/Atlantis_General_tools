out.Atlantis <- function(result){
    df_rel      <- convert_relative_initial(result$structn_age)
    df_resn     <- convert_relative_initial(result$resn_age)
    df_bio      <- convert_relative_initial(result$biomass_age)
    df_eat      <- convert_relative_initial(result$eat_age)
    plots       <- plot_diet(result$biomass_consumed, wrap_col = "agecl", combine_thresh = 7)
    bgm_as_df   <- convert_bgm(result$file[1], result$file[11])
    bio_spatial <- result$biomass_spatial_stanza
    bio_spatial$time <- round(((bio_spatial$time+ 1) / (min(bio_spatial$time, na.rm = T) + 1) - 1) * 365)
    full_grid   <- expand.grid(polygon = unique(bgm_as_df$polygon), layer = min(bio_spatial$layer):max(bio_spatial$layer))
    full_grid   <- left_join(full_grid, bgm_as_df)
    shinyApp(
        ui <- navbarPage(
            "Type of output",
            tabPanel('@ Age plots',
                     tabsetPanel(
                         tabPanel('Structural nitrogen',#
                                  plotOutput('plot1', width = "100%", height = "700px")
                                  ),
                         tabPanel('Reserve nitrogen',#
                                  plotOutput('plot4', width = "100%", height = "700px")
                                  ),
                         tabPanel('Biomass - AgeClass',#
                                  plotOutput('plot3', width = "100%", height = "700px")
                                  ),
                         tabPanel('Eat - AgeClass',#
                                  plotOutput('plot2', width = "100%", height = "700px")
                                  ),
                         tabPanel('Growth - AgeClass',#
                                  plotOutput('plot5', width = "100%", height = "700px")
                                  ),
                         tabPanel('Numbers',#
                                  plotOutput('plot6', width = "100%", height = "700px")
                                  ),
                         tabPanel('Biomass',#
                                  h3('Relative Biomass'),
                                  plotOutput('plot7', width = "100%", height = "700px"),
                                  br(),
                                  h3('Total Biomass'),
                                  plotOutput('plot10', width = "100%", height = "700px")
                                  )
                     )
                     ),
            tabPanel('Diet Plots',
                     sidebarPanel(
                         selectInput('fg', 'Functional Group',  names(plots))
                     ),
                     mainPanel(
                         plotOutput('plot8', width = "100%", height = "700px")
                     )
                     ),
            tabPanel('Biomass Plots',
                     sidebarPanel(
                         selectInput('fg2', 'Functional Group',  unique(bio_spatial$species)),
                         sliderInput("age", "Age class (Stanza)", 1, max(bio_spatial$species_stanza, na.rm = TRUE), step = 1, value = 1),
                         sliderInput("time.sp", "time step", 0, max(bio_spatial$time, na.rm = TRUE), step= diff(sort(unique(bio_spatial$time)))[1], value = 0)
                     ),
                     mainPanel(
                         plotOutput('plot9', width = "100%", height = "700px")
                     )
                     )
        ),
        function(input, output, session) {
            diet.p <- reactive({
                diet.p <- plots[[which(names(plots) == input$fg)]]
            })
            bio.spatial <- reactive({
                fg.filt     <- filter(bio_spatial, species==input$fg2, time ==input$time.sp, species_stanza==input$age)
                bio.spatial <- left_join(full_grid, fg.filt, by = c('polygon','layer'))
            })
            max.b <- reactive({
                max.b <- max(filter(bio_spatial, species==input$fg2)$atoutput, na.rm = TRUE)
            })
            output$plot1 <- renderPlot({
                plot <- plot_line(df_rel, col = "agecl")
                plot <- update_labels(plot, list(x = "Time [years]", y = expression(SN/SN[init])))
                plot_add_box(plot)
            })
            output$plot2 <- renderPlot({
                plot   <- plot_line(df_eat, col = "agecl")
                plot   <- update_labels(plot, list(x = "Time [years]", y = expression(Cons./Cons.[init])))
                plot_add_box(plot)
            })
            output$plot3 <- renderPlot({
                plot <- plot_line(df_bio, col = "agecl")
                plot <- update_labels(plot, list(x = "Time [years]", y = expression(Biomass/Biomass[init])))
                plot_add_box(plot)
            })
            output$plot4 <- renderPlot({
                plot     <- plot_line(df_resn, col = "agecl")
                plot     <- update_labels(plot, list(x = "Time [years]", y = expression(RN/RN[init])))
                plot_add_box(plot)
            })
            output$plot5 <- renderPlot({
                df_rel <- convert_relative_initial(result$growth_age)
                plot <- plot_line(df_rel, col = "agecl")
                plot <- update_labels(plot, list(x = "Time [years]", y = expression(Growth/Growth[init])))
                plot_add_box(plot)

            })
            output$plot6 <- renderPlot({
                df_rel   <- convert_relative_initial(result$nums_age)
                plot     <- plot_line(df_rel, col = "agecl")
                plot     <- update_labels(plot, list(x = "Time [years]", y = expression(Numbers/Numbers[init])))
                plot_add_box(plot)
            })
            output$plot7 <- renderPlot({
                df_rel <- convert_relative_initial(result$biomass)
                plot   <- ggplot(df_rel, aes(x = time, y = atoutput))+geom_line()+
                    facet_wrap(~species) + ylim(0 , 2) + annotate(geom='rect', xmin= 0, ymin= 0.5,
                                                                  xmax=max(df_rel$time), ymax= 1.5, fill="royalblue", lwd=0 , alpha = .2) +
                    theme_bw()
                plot <- update_labels(plot, list(x = "Time [years]", y = expression(Biomass/Biomass[init])))
                plot
            })

            output$plot8 <- renderPlot({
                par(mar=c(5,5,5,5), mgp=c(5,1,0))
                gridExtra::grid.arrange(diet.p())
            })
            output$plot9 <- renderPlot({
                ggplot(data = bio.spatial(), aes(x = long, y = lat, group = polygon, fill = atoutput))+
                    geom_polygon(colour = "black", size = 0.25, na.rm = TRUE)+
                    scale_fill_gradient("biomass distribution", limits = c(0, max.b()), low = "royalblue", high = "red",na.value = 'grey80')+
                    ggplot2::facet_wrap(~layer) +
                    theme_light()
            })
            output$plot10 <- renderPlot({
                plot <- ggplot(result$biomass, aes(x = time, y = atoutput)) + geom_line() +
                    facet_wrap(~species, scale = 'free_y') + theme_bw()
                plot <- update_labels(plot, list(x = "Time [years]", y = "Biomass [Tons]"))
                plot
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















##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Preporcessing of Atlantis Tools
##' @param dir Directory with all the files
##' @param nc_gen Atlantis Output : nc file from atlatnis (output)
##' @param nc_prod Atlantis Output : nc production file
##' @param dietcheck Atlantis Output : diet check file from atlantis
##' @param yoy Atlantis Output : young of the year atlantis output
##' @param ssb Atlantis Output : Spawning stock Biomass
##' @param version_flag Atlantis Output : version of Atlantis
##' @param prm_run Atlantis Input : Run parameter file
##' @param prm_biol Atlantis Input : Biological parameter file
##' @param fgs Atlantis Input : Groups file (.csv)
##' @param bgm Atlantis Input : BGM file (polygons and layers)
##' @param init Atlantis Input : Initial condition file
##' @param extraFG Extra functional groups that are not readed it by load_bps() function
##' @return Preprosessing output, ready to be use by ggplot2
##' @author Demiurgo based on Alex work
pre.Atlantis.tools <- function(dir, nc_gen, nc_prod, dietcheck,  yoy, ssb, version_flag, prm_run, prm_biol, fgs, bgm, init, extraFG = NULL){
    cat('\n Defining additional variables')
    bboxes    <- get_boundary(boxinfo = load_box(dir, bgm))
    bps       <- load_bps(dir, fgs, init)
    if(!is.null(extraFG)) {
        bps  <- c(bps, extraFG)
        cat('\n+++++++++++++\n You add the ', extraFG, 'as extra Functional group\n++++++++++++++++')
    }
    bio_conv  <- get_conv_mgnbiot(dir, prm_biol)
    ## By default data from all groups within the simulation is extracted!
    groups      <- get_groups(dir, fgs)
    groups_age  <- get_age_groups(dir, fgs)
    groups_rest <- groups[!groups %in% groups_age]
    cat('\t\t...Done!')
    cat('\nReading data from Atlantis simulation')
    vars        <- list("Nums", "StructN", "ResN", "N")
    grps        <- list(groups_age, groups_age, groups_age, groups_rest)
    dfs_gen     <- Map(load_nc, select_variable = vars, select_groups = grps,
                       MoreArgs = list(dir = dir, nc = nc_gen, bps = bps, fgs = fgs, prm_run = prm_run, bboxes = bboxes))
    ## Read in raw untransformed data from nc_prod
    vars     <- list("Eat", "Grazing", "Growth")
    grps     <- list(groups_age, groups_rest, groups_age)
    dfs_prod <- Map(load_nc, select_variable = vars, select_groups = grps,
                    MoreArgs = list(dir = dir, nc = nc_prod, bps = bps, fgs = fgs, prm_run = prm_run, bboxes = bboxes))
    ## Read in physics
    flux    <- load_nc_physics(dir = dir, nc = nc_gen, select_physics = c("eflux", "vflux"),
                               prm_run = prm_run, bboxes = bboxes, aggregate_layers = FALSE)
    sink    <- load_nc_physics(dir = dir, nc = nc_gen, select_physics = c("hdsource", "hdsink"),
                               prm_run = prm_run, bboxes = bboxes, aggregate_layers = FALSE)
    physics <- load_nc_physics(dir = dir, nc = nc_gen,
                               select_physics = c("salt", "NO3", "NH3", "Temp", "Chl_a", "Denitrifiction"),
                               prm_run = prm_run, bboxes = bboxes, aggregate_layers = TRUE)
    vol     <- load_nc_physics(dir = dir, nc = nc_gen, select_physics = "volume",
                               prm_run = prm_run, bboxes = bboxes, aggregate_layers = F)
    vol_dz <- load_nc_physics(dir = dir, nc = nc_gen, select_physics = c("volume", "dz"),
                              prm_run = prm_run, bboxes = bboxes, aggregate_layers = F)
    dz     <- dplyr::filter(vol_dz, variable == "dz")
    vol    <- dplyr::filter(vol_dz, variable == "volume")
    nominal_dz <- load_init(dir = dir, init = init, vars = "nominal_dz") %>% as.data.frame() %>%
        dplyr::filter(!is.na(layer))
    ## Read in Dietcheck
    df_dm <- load_dietcheck(dir = dir, dietcheck = dietcheck,
                            fgs = fgs, prm_run = prm_run, version_flag = version_flag, convert_names = TRUE)
    ## Read in SSB/R
    ssb_rec    <- load_rec(dir = dir, yoy = yoy, ssb = ssb, prm_biol = prm_biol)
    ## Read in misc
    df_agemat  <- prm_to_df(dir = dir, prm_biol = prm_biol, fgs = fgs, group = get_age_acronyms(dir, fgs), parameter = "age_mat")
    dietmatrix <- load_dietmatrix(dir, prm_biol, fgs, convert_names = TRUE)
    cat('\t\t...Done!')
    cat('\nApply preprocess calculations')
    ## Calculate biomass spatially
    bio_sp     <- calculate_biomass_spatial(nums = dfs_gen[[1]], sn = dfs_gen[[2]], rn = dfs_gen[[3]], n = dfs_gen[[4]],
                                            vol_dz = vol_dz, bio_conv = bio_conv, bps = bps)
    ## Aggregate spatial biomass to based on stanzas
    bio_sp_stanza <- combine_ages(bio_sp, grp_col = "species", agemat = df_agemat)
    ## Aggregate biomass
    biomass       <- bio_sp %>% agg_data(groups = c("species", "time"), fun = sum)
    biomass_age   <- bio_sp %>% filter(agecl > 2) %>% agg_data(groups = c("species", "agecl", "time"), fun = sum)
    ## Aggregate Numbers! This is done seperately since numbers need to be summed!
    nums     <- agg_data(data = dfs_gen[[1]], groups = c("species", "time"), fun = sum)
    nums_age <- agg_data(data = dfs_gen[[1]], groups = c("species", "agecl", "time"), fun = sum)
    nums_box <- agg_data(data = dfs_gen[[1]], groups = c("species", "polygon", "time"), fun = sum)
    ## Aggregate the rest of the dataframes by mean!
    structn_age <- agg_data(data = dfs_gen[[2]],  groups = c("species", "time", "agecl"), fun = mean)
    resn_age    <- agg_data(data = dfs_gen[[3]],  groups = c("species", "time", "agecl"), fun = mean)
    eat_age     <- agg_data(data = dfs_prod[[1]], groups = c("species", "time", "agecl"), fun = mean)
    grazing     <- agg_data(data = dfs_prod[[2]], groups = c("species", "time"), fun = mean)
    growth_age  <- agg_data(data = dfs_prod[[3]], groups = c("species", "time", "agecl"), fun = mean)
    ## Calculate consumed biomass
    bio_cons    <- calculate_consumed_biomass(eat = dfs_prod[[1]], grazing = dfs_prod[[2]], dm = df_dm,
                                              vol = vol, bio_conv = bio_conv) %>% agg_data(groups = c("pred", "agecl", "time", "prey"), fun = sum)
    ## Calculate spatial overlap
    sp_overlap <- calculate_spatial_overlap(biomass_spatial = bio_sp, dietmatrix = dietmatrix, agemat = df_agemat)
    ## Growth relative to initial conditions
    rec_weight <- prm_to_df(dir = dir, prm_biol = prm_biol, fgs = fgs,
                            group = get_age_acronyms(dir = dir, fgs = fgs),
                            parameter = c("KWRR", "KWSR", "AgeClassSize"))
    pd         <- load_init_weight(dir = dir, init = init, fgs = fgs, bboxes = bboxes) %>% left_join(rec_weight) %>% split(.$species)
    ## Calculate weight difference from one ageclass to the next!
    for (i in seq_along(pd)) {
        pd[[i]]$wdiff <- c((pd[[i]]$rn[1] + pd[[i]]$sn[1]) - (pd[[i]]$kwrr[1] + pd[[i]]$kwsr[1]),
                           diff(pd[[i]]$rn + pd[[i]]$sn))
    }
    pd <- do.call(rbind, pd)
    pd$growth_req <- pd$wdiff / (365 * pd$ageclasssize)
    if (any(pd$growth_req < 0)) {
        warning("Required growth negative for some groups. Please check your initial conditions files.")
    }
    gr_req      <- pd %>% select(species, agecl, growth_req)
    gr_rel_init <- growth_age %>% left_join(gr_req) %>% mutate(gr_rel = (atoutput - growth_req) / growth_req)
    ## Aggregate volume vertically.
    vol_ts      <- agg_data(vol, groups = c("time", "polygon"), fun = sum, out = "volume")
    cat('\t\t...Done!')
    cat('\n\nCombining objects to a list of preprocessed dataframes')
    ## output
    result <- list(
        "biomass"                = biomass,       #1
        "biomass_age"            = biomass_age,
        "biomass_consumed"       = bio_cons,
        "biomass_spatial_stanza" = bio_sp_stanza,
        "diet"                   = df_dm,         #5
        "dz"                     = dz,
        "eat_age"                = eat_age,
        "flux"                   = flux,
        "grazing"                = grazing,
        "growth_age"             = growth_age,    #10
        "growth_rel_init"        = gr_rel_init,
        "nominal_dz"             = nominal_dz,
        "nums"                   = nums,
        "nums_age"               = nums_age,
        "nums_box"               = nums_box,      #15
        "physics"                = physics,
        "resn_age"               = resn_age,
        "sink"                   = sink,
        "spatial_overlap"        = sp_overlap,
        "ssb_rec"                = ssb_rec,       #20
        "structn_age"            = structn_age,
        "vol"                    = vol_ts,
        'files'                  = c(dir, nc_gen, nc_prod, dietcheck,  yoy, ssb, version_flag, prm_run, prm_biol, fgs, bgm, init, extraFG)
    )
    cat('\t\t...Done\n\n!')
    return(result)
}
