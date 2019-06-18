# Main function for creating swan surfer plots

#' Reads in external data, combines with internal datasets and produces Swan
#' River surfer plots
#'
#' \code{swan_surfR} takes a file path to Swan River sonde output and creates
#'     a four panel (single column) surfer plot of salinity, dissolved oxygen,
#'     chlorophyll a, and temperature in pdf format.
#'
#' Surfer plots display a
#'     cross-section of the river where the metrics of interest have been
#'     interpolated between sonde locations. Inverse distance weighting has
#'     been used for the interpolation.
#'
#'     A river "bottom" is displayed which puts the interpolation in context.
#'     The river bottom has been derived from a combination of historical
#'     maximum depths at sampling locations as well as station points extracted
#'     from a "best" navigation line over the latest corporate bathymetry. See
#'     project documentation for further details.
#'
#' @param path Character string filepath to location of sonde data xlsx
#'     workbooks.
#'
#' @return Two pdf format four panel surfer plots of the Swan River. One shows
#'    the full extent of the monitoring run (river mouth to the just beyond the
#'    site POL - Upper Swan Power Lines). The second shows from the Narrows
#'    Bridge to site POL.
#'
#' @examples
#' \dontrun{
#' swan_surfR(path = "Z:/DEC/MonitoringProgram/Data")
#' }
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' @import dplyr
#' @importFrom janitor clean_names
#' @import raster
#' @importFrom gstat idw
#' @import ggplot2
#' @import scales
#' @import grid
#' @import gridExtra
#' @import metR
#' @import ggpubr
#' @importFrom lubridate ymd
#' @importFrom sp coordinates
#'
#' @export
swan_surfR <- function(path){
  locations <- data_finder(path, river = "s")
  data_pairs <- unique(substr(locations, nchar(path)+2, nchar(path)+9))

  for(pair in data_pairs){
    pair_locs <- locations[grep(pair, locations)]

    # read in lower data
    lower <- sonde_reader(path = pair_locs[1])
    lower_clean <- lower[complete.cases(lower), ]

    # read in upper data
    upper <- sonde_reader(path = pair_locs[2])
    upper_clean <- upper[complete.cases(upper), ]

    # join and delete and rename prob sites
    samp_data <- dplyr::bind_rows(lower_clean, upper_clean) %>%
      janitor::clean_names() %>%
      dplyr::mutate(site = ifelse(site == "MULB FARM", "MUL",
                                  ifelse(site == "BWR10", "BWR", site))) %>%
      dplyr::filter(site != "FP2.1")

    # join sites to WQ data
    comb_data <- dplyr::left_join(samp_data, S_sitesdf, by = "site")
    d_reduced <- comb_data %>%
      dplyr::select(site, sal_ppt, do_mg_l, c, chl_ug_l, dep_m,
                    dist_mouth, max_depth)

    # Set up labels and params to plot
    sparams <- c("Salinity", "Dissolved_Oxygen", "Temperature", "Chlorophyll")

    ## Create interpolations and store in  separate lists
    d_all <- d_reduced[complete.cases(d_reduced),] %>%
      dplyr::mutate(y = -1 * dep_m, x = dist_mouth/1000)

    d_nar <- d_reduced[complete.cases(d_reduced),] %>%
      dplyr::mutate(y = -1 * dep_m, x = dist_mouth/1000) %>%
      dplyr::filter(x >= 21)

    vals <- c("sal_ppt", "do_mg_l", "c", "chl_ug_l")
    idw_list_a <- vector("list", length(vals))
    idw_list_n <- vector("list", length(vals))
    names(idw_list_a) <- vals
    names(idw_list_n) <- vals

    # for all
    for(i in seq_along(vals)){
      val <- vals[i]
      d1 <- d_all[,c("x", "y", val)]
      names(d1)[3] <- "value"
      sp::coordinates(d1) <- ~x + y

      idw <- gstat::idw(formula = value ~ 1,
                 locations = d1,
                 newdata = S_grd_all,
                 idp = 4)

      idw1_r <- raster::raster(idw)
      idw1_r_class <- raster::reclassify(idw1_r, reclass_matrices[[i]])
      idw1_sp <- raster::rasterToPoints(idw1_r_class, spatial = TRUE)
      idw1_df <- data.frame(idw1_sp)[-4]# ditch option
      names(idw1_df)[1] <- sparams[i]
      idw1_df <- idw1_df[-(1:530),] # strip weird close to -0.1 vals

      idw_list_a[[i]] <- idw1_df
    }

    # for narrows up
    for(i in seq_along(vals)){
      val <- vals[i]
      d1 <- d_nar[,c("x", "y", val)]
      names(d1)[3] <- "value"
      sp::coordinates(d1) <- ~x + y

      idw <- gstat::idw(formula = value ~ 1,
                 locations = d1,
                 newdata = S_grd_nar,
                 idp = 4)

      idw1_r <- raster::raster(idw)
      idw1_r_class <- raster::reclassify(idw1_r, reclass_matrices[[i]])
      idw1_sp <- raster::rasterToPoints(idw1_r_class, spatial = TRUE)
      idw1_df <- data.frame(idw1_sp)[-4]# ditch option
      names(idw1_df)[1] <- sparams[i]
      idw_list_n[[i]] <- idw1_df
    }

    # make sample collection points
    samp_locs <- comb_data %>%
      dplyr::mutate(dist_mouth = dist_mouth/1000) %>%
      dplyr::rename(x = dist_mouth, y = dep_m) %>%
      dplyr::select(site, x, y)

    site_labs <- S_sitesdf %>%
      dplyr::filter(site != "SRP_RSSA") %>%
      dplyr::filter(site != "BWR" & site != "KMO" & site != "VIT")

    ## Plots
    salPlot <- ggplot()+
      geom_raster(data = idw_list_a[[1]],
                  aes(x=x, y=y, fill = factor(Salinity))) +
      scale_x_continuous(limits = c(-1, 52.5),
                         expand = c(0, 0)) +
      stat_contour2(data = idw_list_a[[1]], aes(x=x, y=y, z = Salinity),
                    colour = "grey50",
                    breaks = MakeBreaks(binwidth = 2)) +
      scale_fill_manual(values = surfer_cols("sal"),
                        guide = guide_legend(reverse=T),
                        name = "Salinity\n(ppt)") +
      geom_polygon(data = S_bottom,
                   aes(x=x, y=y), fill = "grey90", colour = "grey20") +
      geom_text_contour(data = idw_list_a[[1]],
                        aes(x=x, y=y, z = Salinity),
                        check_overlap = TRUE,
                        size = 6,
                        stroke = 0.2,
                        breaks = MakeBreaks(binwidth = 2)) +
      geom_point(data = samp_locs,
                 aes(x = x, y = - y),
                 colour = "black",
                 size = 0.5) +
      geom_text(data = site_labs,
                aes(x = dist_mouth/1000, y = 0.7, label = site),
                size = 4.5,
                colour = "black",
                alpha = 1,
                check_overlap = TRUE) +
      annotate("text",
               label = "Salinity (ppt)",
               x = 26.5,
               y = -20,
               size = 7,
               fontface =2,
               colour = "black") +
      labs(y = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 16),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
            plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
            panel.border=element_blank(),
            legend.background = element_rect(fill = "lightgray"),
            legend.direction = "horizontal",
            legend.position = c(0.65, 0.22),
            legend.key.size =  unit(8, "mm"),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))

    salPlotZ <- ggplot()+
      geom_raster(data = idw_list_n[[1]],
                  aes(x=x, y=y, fill = factor(Salinity))) +
      stat_contour2(data = idw_list_n[[1]],
                    aes(x=x, y=y, z = Salinity),
                    colour = "grey50",
                    breaks = MakeBreaks(binwidth = 2)) +
      scale_fill_manual(values = surfer_cols("sal"),
                        guide = guide_legend(reverse=T),
                        name = "Salinity\n(ppt)") +
      geom_polygon(data = S_bottom_nar,
                   aes(x=x, y=y), fill = "grey90", colour = "grey20") +
      geom_text_contour(data = idw_list_n[[1]],
                        aes(x=x, y=y, z = Salinity),
                        check_overlap = TRUE,
                        size = 6,
                        stroke = 0.2,
                        breaks = MakeBreaks(binwidth = 2)) +
      geom_point(data = filter(samp_locs, x >= 21),
                 aes(x = x, y = - y),
                 colour = "black",
                 size = 0.5) +
      geom_text(data = filter(site_labs, dist_mouth/1000 >= 21),
                aes(x = dist_mouth/1000, y = 0.7, label = site),
                size = 4.5,
                colour = "black",
                alpha = 1,
                check_overlap = TRUE) +
      annotate("text",
               label = "Salinity (ppt)",
               x = 32,
               y = -9,
               size = 7,
               fontface =2,
               colour = "black") +
      labs(y = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 16),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
            plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
            panel.border=element_blank(),
            legend.background = element_rect(fill = "lightgray"),
            legend.direction = "horizontal",
            legend.position = c(0.65, 0.22),
            legend.key.size =  unit(8, "mm"),
            legend.text = element_text(size = 12),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))

    doPlot <- ggplot()+
      geom_raster(data = idw_list_a[[2]],
                  aes(x=x, y=y, fill = factor(Dissolved_Oxygen))) +
      scale_x_continuous(limits = c(-1, 52.5),
                         expand = c(0, 0)) +
      stat_contour2(data = idw_list_a[[2]],
                    aes(x=x, y=y, z = Dissolved_Oxygen),
                    colour = "grey50",
                    breaks = MakeBreaks(binwidth = 1)) +
      scale_fill_manual(values = surfer_cols("do"),
                        guide = guide_legend(reverse=T),
                        name = "Dissolved\nOxygen\n(mg/L)") +
      geom_polygon(data = S_bottom,
                   aes(x=x, y=y), fill = "grey90", colour = "grey20") +
      geom_text_contour(data = idw_list_a[[2]],
                        aes(x=x, y=y, z = Dissolved_Oxygen),
                        skip = 1,
                        check_overlap = TRUE,
                        size = 6,
                        stroke = 0.2,
                        breaks = MakeBreaks(binwidth = 1)) +
      geom_point(data = samp_locs,
                 aes(x = x, y = - y),
                 colour = "black",
                 size = 0.5) +
      geom_point(data = S_oxy_locs,
                 aes(x = x, y = y),
                 size = 6,
                 colour = "black",
                 bg = "green",
                 shape = 24) +
      annotate("text",
               label = "Dissolved Oxygen (mg/L)",
               x = 26.5,
               y = -20,
               size = 7,
               fontface =2,
               colour = "black") +
      labs(y = "Depth (m)") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 16),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust=0.5, vjust=0.5, face='bold'),
            plot.subtitle = element_text(hjust=0.5, vjust=0.5),
            panel.border=element_blank(),
            legend.background = element_rect(fill = "lightgray"),
            legend.direction = "horizontal",
            legend.position = c(0.65, 0.22),
            legend.key.size =  unit(8, "mm"),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))

    doPlotZ <- ggplot()+
      geom_raster(data = idw_list_n[[2]],
                  aes(x=x, y=y, fill = factor(Dissolved_Oxygen))) +
      stat_contour2(data = idw_list_n[[2]],
                    aes(x=x, y=y, z = Dissolved_Oxygen),
                    colour = "grey50",
                    breaks = MakeBreaks(binwidth = 1)) +
      scale_fill_manual(values = surfer_cols("do"),
                        guide = guide_legend(reverse=T),
                        name = "Dissolved\nOxygen\n(mg/L)") +
      geom_polygon(data = S_bottom_nar,
                   aes(x=x, y=y), fill = "grey90", colour = "grey20") +
      geom_text_contour(data = idw_list_n[[2]],
                        aes(x=x, y=y, z = Dissolved_Oxygen),
                        skip = 1,
                        check_overlap = TRUE,
                        size = 6,
                        stroke = 0.2,
                        breaks = MakeBreaks(binwidth = 1)) +
      geom_point(data = filter(samp_locs, x >= 21),
                 aes(x = x, y = - y),
                 colour = "black",
                 size = 0.5) +
      geom_point(data = S_oxy_locs,
                 aes(x = x, y = y),
                 size = 6,
                 colour = "black",
                 bg = "green",
                 shape = 24) +
      annotate("text",
               label = "Dissolved Oxygen (mg/L)",
               x = 32,
               y = -9,
               size = 7,
               fontface =2,
               colour = "black") +
      labs(y = "Depth (m)") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 16),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust=0.5, vjust=0.5, face='bold'),
            plot.subtitle = element_text(hjust=0.5, vjust=0.5),
            panel.border=element_blank(),
            legend.background = element_rect(fill = "lightgray"),
            legend.direction = "horizontal",
            legend.position = c(0.65, 0.22),
            legend.key.size =  unit(8, "mm"),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))


    chlorPlot <- ggplot()+
      geom_raster(data = idw_list_a[[4]],
                  aes(x=x, y=y, fill = factor(Chlorophyll)),
                  alpha = 0.5) +
      scale_x_continuous(limits = c(-1, 52.5),
                         expand = c(0, 0)) +
      stat_contour2(data = idw_list_a[[4]],
                    aes(x=x, y=y, z = Chlorophyll),
                    colour = "grey50",
                    breaks = as.numeric(chl_brk)) +
      scale_fill_manual(values = surfer_cols("chl"),
                        guide = guide_legend(reverse=T),
                        name = "Chlorophyll\n(ug/L)",
                        labels = c("20", "40", "60", "80", "120", "200",
                                   "> 200")) +
      geom_polygon(data = S_bottom,
                   aes(x=x, y=y), fill = "grey90", colour = "grey20") +
      geom_text_contour(data = idw_list_a[[4]],
                        aes(x=x, y=y, z = Chlorophyll),
                        skip = 3,
                        check_overlap = TRUE,
                        size = 6,
                        stroke = 0.2,
                        breaks = MakeBreaks(binwidth = 20)) +
      geom_point(data = samp_locs,
                 aes(x = x, y = - y),
                 colour = "black",
                 size = 0.5) +
      annotate("text",
               label = expression('bold(paste("Chlorophyll (", mu,"/L)"))'),
               x = 26.5,
               y = -20,
               size = 7,
               fontface =2,
               colour = "black", parse = TRUE) +
      labs(x = "Distance From Entrance (km)",
           y = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 16),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
            plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
            panel.border=element_blank(),
            legend.background = element_rect(fill = "lightgray"),
            legend.direction = "horizontal",
            legend.position = c(0.65, 0.22),
            legend.key.size =  unit(8, "mm"),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))

    chlorPlotZ <- ggplot()+
      geom_raster(data = idw_list_n[[4]],
                  aes(x=x, y=y, fill = factor(Chlorophyll)),
                  alpha = 0.5) +
      stat_contour2(data = idw_list_n[[4]],
                    aes(x=x, y=y, z = Chlorophyll),
                    colour = "grey50",
                    breaks = chl_brk) +
      scale_fill_manual(values = surfer_cols("chl"),
                        guide = guide_legend(reverse=T),
                        name = "Chlorophyll\n(ug/L)",
                        labels = c("20", "40", "60", "80", "120", "200",
                                   "> 200")) +
      geom_polygon(data = S_bottom_nar,
                   aes(x=x, y=y), fill = "grey90", colour = "grey20") +
      geom_text_contour(data = idw_list_n[[4]],
                        aes(x=x, y=y, z = Chlorophyll),
                        skip = 3,
                        check_overlap = TRUE,
                        size = 6,
                        stroke = 0.2,
                        breaks = as.numeric(chl_brk)) +
      geom_point(data = filter(samp_locs, x >= 21),
                 aes(x = x, y = - y),
                 colour = "black",
                 size = 0.5) +
      annotate("text",
               label = expression('bold(paste("Chlorophyll (", mu,"/L)"))'),
               x = 32,
               y = -9,
               size = 7,
               fontface =2,
               colour = "black", parse = TRUE) +
      labs(x = "Distance From Entrance (km)",
           y = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 16),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
            plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
            panel.border=element_blank(),
            legend.background = element_rect(fill = "lightgray"),
            legend.direction = "horizontal",
            legend.position = c(0.65, 0.22),
            legend.key.size =  unit(8, "mm"),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))

    tempPlot <- ggplot()+
      geom_raster(data = idw_list_a[[3]],
                  aes(x=x, y=y, fill = factor(Temperature))) +
      scale_x_continuous(limits = c(-1, 52.5),
                         expand = c(0, 0)) +
      stat_contour2(data = idw_list_a[[3]],
                    aes(x=x, y=y, z = Temperature),
                    colour = "grey50",
                    breaks = MakeBreaks(binwidth = 1)) +
      scale_fill_manual(values = surfer_cols("temp"),
                        guide = guide_legend(reverse=T),
                        name = "Temp") +
      geom_polygon(data = S_bottom,
                   aes(x=x, y=y), fill = "grey90", colour = "grey20") +
      geom_text_contour(data = idw_list_a[[3]],
                        aes(x=x, y=y, z = Temperature),
                        size = 6,
                        stroke = 0.2,
                        breaks = MakeBreaks(binwidth = 1)) +
      geom_point(data = samp_locs,
                 aes(x = x, y = - y),
                 colour = "black",
                 size = 0.5) +
      annotate("text",
               label = expression('bold(paste("Temperature (", degree,"C)"))'),
               x = 26.5,
               y = -20,
               size = 7,
               fontface =2,
               colour = "black", parse = TRUE) +
      labs(x = "Distance From Entrance (km)",
           y = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.line.y = element_blank(),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 16),
            plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
            plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
            panel.border=element_blank(),
            legend.background = element_rect(fill = "lightgray"),
            legend.direction = "horizontal",
            legend.position = c(0.65, 0.22),
            legend.key.size =  unit(8, "mm"),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))

    tempPlotZ <- ggplot()+
      geom_raster(data = idw_list_n[[3]],
                  aes(x=x, y=y, fill = factor(Temperature))) +
      stat_contour2(data = idw_list_n[[3]],
                    aes(x=x, y=y, z = Temperature),
                    colour = "grey50",
                    breaks = MakeBreaks(binwidth = 1)) +
      scale_fill_manual(values = surfer_cols("temp"),
                        guide = guide_legend(reverse=T),
                        name = "Temp") +
      geom_polygon(data = S_bottom_nar,
                   aes(x=x, y=y), fill = "grey90", colour = "grey20") +
      geom_text_contour(data = idw_list_n[[3]],
                        aes(x=x, y=y, z = Temperature),
                        size = 6,
                        stroke = 0.2,
                        breaks = MakeBreaks(binwidth = 1)) +
      geom_point(data = filter(samp_locs, x >= 21),
                 aes(x = x, y = - y),
                 colour = "black",
                 size = 0.5) +
      annotate("text",
               label = expression('bold(paste("Temperature (", degree,"C)"))'),
               x = 32,
               y = -9,
               size = 7,
               fontface =2,
               colour = "black", parse = TRUE) +
      labs(x = "Distance From Entrance (km)",
           y = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.line.y = element_blank(),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 16),
            plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
            plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
            panel.border=element_blank(),
            legend.background = element_rect(fill = "lightgray"),
            legend.direction = "horizontal",
            legend.position = c(0.65, 0.22),
            legend.key.size =  unit(8, "mm"),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))


    surfers <- grid.arrange(salPlot, doPlot, chlorPlot, tempPlot, nrow = 4,
                            top = textGrob(paste("Swan River Estuary - Physical-Chemical Profile -",
                                                 format(ymd(pair), "%d %b %Y")),
                                           gp = gpar(fontface = 2, fontsize = 22)))
    pdf_name <- paste0(path, "/plots/", "swan_", ymd(pair), "_surfer.pdf")
    cat(paste0(pdf_name,"\n"))

    ggsave(plot = surfers, filename = pdf_name, width=28, height=18)

    # for zoomed at narrows
    surfersZ <- grid.arrange(salPlotZ, doPlotZ, chlorPlotZ, tempPlotZ, nrow = 4,
                             top = textGrob(paste("Middle and Upper Swan River Estuary - Physical-Chemical Profile -",
                                                  format(lubridate::ymd(pair), "%d %b %Y")),
                                            gp = gpar(fontface = 2, fontsize = 22)))
    pdf_nameZ <- paste0(path, "/plots/", "swan_middle_upper_", lubridate::ymd(pair), "_surfer.pdf")
    cat(paste0(pdf_nameZ,"\n"))

    ggsave(plot = surfersZ, filename = pdf_nameZ, width=28, height=18)
  }
}
