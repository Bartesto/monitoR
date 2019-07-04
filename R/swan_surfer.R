# Main function for creating swan surfer plots

#' Reads in external data, combines with internal datasets and produces Swan
#' River surfer plots
#'
#' \code{swan_surfR} takes a file path to Swan River sonde outputs and creates
#'     a four panel (single column) surfer plot of salinity, dissolved oxygen,
#'     chlorophyll a, and temperature in pdf format. The function creates a
#'     folder called "plots" in the file path to store the pdf's. Code expects
#'     only 2 excel workbooks for one monitoring run.
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
#' @param ovit Character string ("green", "blue" or "red") indicating
#'     oxygenation plant status at site VIT.
#'
#' @param ocav Character string ("green", "blue" or "red") indicating
#'     oxygenation plant status at site CAV.
#'
#' @return Two pdf format four panel surfer plots of the Swan River. One shows
#'    the full extent of the monitoring run (river mouth to the just beyond the
#'    site POL - Upper Swan Power Lines). The second shows from the Narrows
#'    Bridge to site POL.
#'
#' @examples
#' \dontrun{
#' swan_surfR(path = "Z:/DEC/MonitoringProgram/Data", ovit = "green", ocav = "red")
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
#' @import gtable
#' @import metR
#' @import ggpubr
#' @importFrom lubridate ymd
#' @importFrom sp coordinates
#'
#' @export
swan_surfR <- function(path, ovit, ocav){
  locations <- data_finder(path, river = "s")

  #error handler
  if(length(locations) == 2){
    data_pairs <- unique(substr(locations, nchar(path)+2, nchar(path)+9))
    # make folder for output
    folder <- file.path(path, "plots")
    if (!file.exists(folder)) {
      dir.create(folder)
    }

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

      # construct pretty date
      sday <- just_nums(as.numeric(substr(pair, 7, 8)))
      sdate <- paste(sday, format(ymd(pair), "%b %Y"), sep = " ")

      # oxygenation plant status colour
      oxy_col <- c(ovit, ocav)
      S_oxy_locs$c <- oxy_col

      # create mock plot to harvest legend for oxy locs status
      # maybe move to internal data
      triangle_df <- data.frame(x = c(1, 2, 3), y = 2, stat = c("a", "b", "c"))
      o_plot<- ggplot(triangle_df) +
        geom_point(aes(x = x, y = y, fill = stat), shape = 24, colour = "black") +
        scale_fill_manual(name = "Oxygen Plant Operational Status",
                          values = c("green", "blue", "red"),
                          labels = c("Operable and operating for part or all of the 24 hours\nprior to sampling",
                                     "Operable but not triggered to operate in the 24 hours\nprior to sampling",
                                     "Inoperable for part or all of the 24 hours prior to sampling")) +
        theme_bw() +
        theme(legend.key.size = unit(8, "mm"),
              legend.background = element_blank(),
              legend.box.background = element_rect(colour = "black", fill = "white"),
              legend.title = element_text(face="bold")) +
        guides(fill = guide_legend(override.aes = list(size=5)))
      oxY_grob <- gtable_filter(ggplot_gtable(ggplot_build(o_plot)), "guide-box")

      # create df to hold co-ords of black out rectangles for missing data
      # maybe move to internal data
      blockdf <- S_sitesdf %>%
        dplyr::filter(site != "SRP_RSSA") %>%
        dplyr::select(site, dist_mouth) %>%
        dplyr::mutate(diff = (dist_mouth - lag(dist_mouth))/1000) %>%
        dplyr::mutate(hlf = diff/2, xmin = dist_mouth/1000 - hlf) %>%
        dplyr::mutate(xmax = lead(xmin), ymin = -22.1, ymax = 0) %>%
        dplyr::select(-diff, -hlf)

      blockdf[1, 3] <- -1
      blockdf[23, 4] <- 51.6

      ## Plots
      salPlot <- ggplot()+
        geom_raster(data = idw_list_a[[1]],
                    aes(x=x, y=y, fill = factor(Salinity))) +
        scale_x_continuous(limits = c(-1, 51.6),
                           expand = c(0, 0),
                           breaks = seq(0, 50, by = 5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
        stat_contour2(data = idw_list_a[[1]], aes(x=x, y=y, z = Salinity),
                      colour = "grey10",
                      breaks = MakeBreaks(binwidth = 2)) +
        scale_fill_manual(values = surfer_cols("sal"),
                          guide = guide_legend(reverse=T),
                          name = "Salinity\n(ppt)") +
        geom_polygon(data = S_bottom,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        # geom_text_contour(data = idw_list_a[[1]],
        #                   aes(x=x, y=y, z = Salinity),
        #                   check_overlap = TRUE,
        #                   size = 6,
        #                   #stroke = 0.2,
        #                   breaks = MakeBreaks(binwidth = 2)) +
        geom_point(data = samp_locs,
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        geom_text(data = site_labs,
                  aes(x = dist_mouth/1000, y = 0.7, label = site, fontface=2),
                  size = 4.5,
                  colour = "black",
                  alpha = 1,
                  check_overlap = TRUE) +
        annotate("text",
                 label = "Salinity (ppt)",
                 x = 19,
                 y = -16.5,
                 size = 9,
                 fontface =2,
                 colour = "black") +
        labs(title = paste("Swan River Estuary - Physical-Chemical Profile -",
                           sdate),
             y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.length.y.left = (unit(2, "mm")),
              axis.ticks.length.y.right = (unit(2, "mm")),
              #axis.ticks.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 28),
              plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
              legend.background = element_rect(fill = "transparent"),
              legend.direction = "horizontal",
              legend.position = c(0.65, 0.22),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      salPlotZ <- ggplot()+
        geom_raster(data = idw_list_n[[1]],
                    aes(x=x, y=y, fill = factor(Salinity))) +
        scale_x_continuous(limits = c(20.95,51.6),
                           expand = c(0, 0),
                           breaks = c(25, 30, 35, 40, 45, 50)) +
        scale_y_continuous(breaks = c(-8, -6, -4, -2, 0),
                           expand = expand_scale(mult = c(0, .05)))+
        stat_contour2(data = idw_list_n[[1]],
                      aes(x=x, y=y, z = Salinity),
                      colour = "grey10",
                      breaks = MakeBreaks(binwidth = 2)) +
        scale_fill_manual(values = surfer_cols("sal"),
                          guide = guide_legend(reverse=T),
                          name = "") +
        geom_polygon(data = S_bottom_nar,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        # geom_text_contour(data = idw_list_n[[1]],
        #                   aes(x=x, y=y, z = Salinity),
        #                   check_overlap = TRUE,
        #                   size = 6,
        #                   #stroke = 0.2,
        #                   breaks = MakeBreaks(binwidth = 2)) +
        geom_point(data = filter(samp_locs, x >= 21),
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        geom_text(data = filter(site_labs, dist_mouth/1000 >= 21),
                  aes(x = dist_mouth/1000, y = 0.7, label = site, fontface=2),
                  size = 4.5,
                  colour = "black",
                  alpha = 1,
                  check_overlap = TRUE) +
        annotate("text",
                 label = "Salinity (ppt)",
                 x = 32.8,
                 y = -7.4,
                 size = 9,
                 fontface =2,
                 colour = "black") +
        labs(title = paste("Middle and Upper Swan River Estuary - Physical-Chemical Profile -",
                           sdate),
             y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.length.y.left = (unit(2, "mm")),
              axis.ticks.length.y.right = (unit(2, "mm")),
              #axis.ticks.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 28),
              plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
              legend.background = element_rect(fill = "transparent"),
              legend.direction = "horizontal",
              legend.title = element_blank(),
              legend.position = c(0.65, 0.22),
              legend.key.size =  unit(8, "mm"),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      doPlot <- ggplot()+
        geom_raster(data = idw_list_a[[2]],
                    aes(x=x, y=y, fill = factor(Dissolved_Oxygen))) +
        scale_x_continuous(limits = c(-1, 51.6),
                           expand = c(0, 0),
                           breaks = seq(0, 50, by = 5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
        stat_contour2(data = idw_list_a[[2]],
                      aes(x=x, y=y, z = Dissolved_Oxygen),
                      colour = "grey10",
                      breaks = MakeBreaks(binwidth = 1)) +
        scale_fill_manual(values = surfer_cols("do"),
                          guide = guide_legend(reverse=T),
                          name = "") +
        geom_polygon(data = S_bottom,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        # geom_text_contour(data = idw_list_a[[2]],
        #                   aes(x=x, y=y, z = Dissolved_Oxygen),
        #                   skip = 1,
        #                   check_overlap = TRUE,
        #                   size = 6,
        #                   #stroke = 0.2,
        #                   breaks = MakeBreaks(binwidth = 1)) +
        geom_point(data = samp_locs,
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        geom_point(data = S_oxy_locs,
                   aes(x = x, y = y),
                   size = 6,
                   colour = "black",
                   bg = S_oxy_locs$c,
                   shape = 24) +
        annotate("text",
                 label = "Dissolved Oxygen (mg/L)",
                 x = 21,
                 y = -16.7,
                 size = 9,
                 fontface = 2,
                 colour = "black") +
        annotation_custom(grob = oxY_grob, xmin = 40, xmax = 47, ymin = -20, ymax = -12) +
        labs(y = "Depth (m)") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.length.y.left = (unit(2, "mm")),
              axis.ticks.length.y.right = (unit(2, "mm")),
              #axis.ticks.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              plot.title = element_text(hjust=0.5, vjust=0.5, face='bold'),
              plot.subtitle = element_text(hjust=0.5, vjust=0.5),
              legend.background = element_rect(fill = "transparent"),
              legend.direction = "horizontal",
              legend.position = c(0.65, 0.22),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      doPlotZ <- ggplot()+
        geom_raster(data = idw_list_n[[2]],
                    aes(x=x, y=y, fill = factor(Dissolved_Oxygen))) +
        scale_x_continuous(limits = c(20.95,51.6),
                           expand = c(0, 0),
                           breaks = c(25, 30, 35, 40, 45, 50)) +
        scale_y_continuous(breaks = c(-8, -6, -4, -2, 0),
                           expand = expand_scale(mult = c(0, .05)))+
        stat_contour2(data = idw_list_n[[2]],
                      aes(x=x, y=y, z = Dissolved_Oxygen),
                      colour = "grey10",
                      breaks = MakeBreaks(binwidth = 1)) +
        scale_fill_manual(values = surfer_cols("do"),
                          guide = guide_legend(reverse=T),
                          name = "") +
        geom_polygon(data = S_bottom_nar,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        # geom_text_contour(data = idw_list_n[[2]],
        #                   aes(x=x, y=y, z = Dissolved_Oxygen),
        #                   skip = 1,
        #                   check_overlap = TRUE,
        #                   size = 6,
        #                   #stroke = 0.2,
        #                   breaks = MakeBreaks(binwidth = 1)) +
        geom_point(data = filter(samp_locs, x >= 21),
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        geom_point(data = S_oxy_locs,
                   aes(x = x, y = y + 1),#nudge up
                   size = 6,
                   colour = "black",
                   bg = S_oxy_locs$c,
                   shape = 24) +
        annotate("text",
                 label = "Dissolved Oxygen (mg/L)",
                 x = 34,
                 y = -7.57,
                 size = 9,
                 fontface = 2,
                 colour = "black") +
        annotation_custom(grob = oxY_grob, xmin = 44, xmax = 50, ymin = -9, ymax = -6.8) +
        labs(y = "Depth (m)") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.length.y.left = (unit(2, "mm")),
              axis.ticks.length.y.right = (unit(2, "mm")),
              #axis.ticks.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              plot.title = element_text(hjust=0.5, vjust=0.5, face='bold'),
              plot.subtitle = element_text(hjust=0.5, vjust=0.5),
              legend.background = element_rect(fill = "transparent"),
              legend.direction = "horizontal",
              legend.position = c(0.65, 0.22),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))


      chlorPlot <- ggplot()+
        geom_raster(data = idw_list_a[[4]],
                    aes(x=x, y=y, fill = factor(Chlorophyll)),
                    alpha = 0.5) +
        scale_x_continuous(limits = c(-1, 51.6),
                           expand = c(0, 0),
                           breaks = seq(0, 50, by = 5)) +
        scale_y_continuous(limits = c(-22.1, 0),
                           expand = c(0, 0)) +
        #scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
        stat_contour2(data = idw_list_a[[4]],
                      aes(x=x, y=y, z = Chlorophyll),
                      colour = "grey10",
                      breaks = as.numeric(chl_brk)) +
        scale_fill_manual(values = surfer_cols("chl"),
                          guide = guide_legend(reverse=T),
                          name = "Chlorophyll\n(ug/L)",
                          labels = c("20", "40", "60", "80", "120", "160",
                                     "200", "300", "400", "> 400")) +
        geom_polygon(data = S_bottom,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        # geom_text_contour(data = idw_list_a[[4]],
        #                   aes(x=x, y=y, z = Chlorophyll),
        #                   skip = 3,
        #                   check_overlap = TRUE,
        #                   size = 6,
        #                   #stroke = 0.2,
        #                   breaks = as.numeric(chl_brk)) +
        geom_point(data = samp_locs,
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        annotate("text",
                 label = expression('bold(paste("F - Chlorophyll (", mu,"g/L)"))'),
                 x = 20.4,
                 y = -16.8,
                 size = 9,
                 fontface =2,
                 colour = "black", parse = TRUE) +
        labs(x = "Distance From Entrance (km)",
             y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              #panel.border = element_blank(),
              panel.border = element_rect(fill = NA, colour = "#CCCCCC"),
              axis.line = element_line(colour = "black"),
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.length.y.left = (unit(2, "mm")),
              axis.ticks.length.y.right = (unit(2, "mm")),
              #axis.ticks.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
              plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
              legend.background = element_rect(fill = "transparent"),
              legend.direction = "horizontal",
              legend.position = c(0.65, 0.22),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      chlorPlotZ <- ggplot()+
        geom_raster(data = idw_list_n[[4]],
                    aes(x=x, y=y, fill = factor(Chlorophyll)),
                    alpha = 0.5) +
        scale_x_continuous(limits = c(20.95,51.6),
                           expand = c(0, 0),
                           breaks = c(25, 30, 35, 40, 45, 50)) +
        scale_y_continuous(breaks = c(-8, -6, -4, -2, 0),
                           expand = c(0, 0)) +
        #expand = expand_scale(mult = c(0, .05)))+
        stat_contour2(data = idw_list_n[[4]],
                      aes(x=x, y=y, z = Chlorophyll),
                      colour = "grey10",
                      breaks = as.numeric(chl_brk)) +
        scale_fill_manual(values = surfer_cols("chl"),
                          guide = guide_legend(reverse=T),
                          name = "Chlorophyll\n(ug/L)",
                          labels = c("20", "40", "60", "80", "120", "160",
                                     "200", "300", "400", "> 400")) +
        geom_polygon(data = S_bottom_nar,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        # geom_text_contour(data = idw_list_n[[4]],
        #                   aes(x=x, y=y, z = Chlorophyll),
        #                   skip = 3,
        #                   check_overlap = TRUE,
        #                   size = 6,
        #                   #stroke = 0.2,
        #                   breaks = as.numeric(chl_brk)) +
        geom_point(data = filter(samp_locs, x >= 21),
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        annotate("text",
                 label = expression('bold(paste("F - Chlorophyll (", mu,"g/L)"))'),
                 x = 33.6,
                 y = -7.6,
                 size = 9,
                 fontface =2,
                 colour = "black", parse = TRUE) +
        labs(x = "Distance From Entrance (km)",
             y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              #panel.border = element_blank(),
              panel.border = element_rect(fill = NA, colour = "#CCCCCC"),
              axis.line = element_line(colour = "black"),
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.length.y.left = (unit(2, "mm")),
              axis.ticks.length.y.right = (unit(2, "mm")),
              #axis.ticks.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
              plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
              legend.background = element_rect(fill = "transparent"),
              legend.direction = "horizontal",
              legend.position = c(0.65, 0.22),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      tempPlot <- ggplot()+
        geom_raster(data = idw_list_a[[3]],
                    aes(x=x, y=y, fill = factor(Temperature))) +
        scale_x_continuous(limits = c(-1, 51.6),
                           expand = c(0, 0),
                           breaks = seq(0, 50, by = 5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
        stat_contour2(data = idw_list_a[[3]],
                      aes(x=x, y=y, z = Temperature),
                      colour = "grey10",
                      breaks = MakeBreaks(binwidth = 1)) +
        scale_fill_manual(values = surfer_cols("temp"),
                          guide = guide_legend(reverse=T),
                          name = "Temp") +
        geom_polygon(data = S_bottom,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        # geom_text_contour(data = idw_list_a[[3]],
        #                   aes(x=x, y=y, z = Temperature),
        #                   size = 6,
        #                   #stroke = 0.2,
        #                   breaks = MakeBreaks(binwidth = 1)) +
        geom_point(data = samp_locs,
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        annotate("text",
                 label = expression('bold(paste("Temperature (", degree,"C)"))'),
                 x = 19,
                 y = -16.8,
                 size = 9,
                 fontface =2,
                 colour = "black", parse = TRUE) +
        labs(x = "Distance From Entrance (km)",
             y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_blank(),
              axis.ticks.length.y.left = (unit(2, "mm")),
              axis.ticks.length.y.right = (unit(2, "mm")),
              #axis.ticks.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
              plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
              legend.background = element_rect(fill = "transparent"),
              legend.direction = "horizontal",
              legend.position = c(0.65, 0.22),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      tempPlotZ <- ggplot()+
        geom_raster(data = idw_list_n[[3]],
                    aes(x=x, y=y, fill = factor(Temperature))) +
        scale_x_continuous(limits = c(20.95,51.6),
                           expand = c(0, 0),
                           breaks = c(25, 30, 35, 40, 45, 50)) +
        scale_y_continuous(breaks = c(-8, -6, -4, -2, 0),
                           expand = expand_scale(mult = c(0, .05)))+
        stat_contour2(data = idw_list_n[[3]],
                      aes(x=x, y=y, z = Temperature),
                      colour = "grey10",
                      breaks = MakeBreaks(binwidth = 1)) +
        scale_fill_manual(values = surfer_cols("temp"),
                          guide = guide_legend(reverse=T),
                          name = "Temp") +
        geom_polygon(data = S_bottom_nar,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        # geom_text_contour(data = idw_list_n[[3]],
        #                   aes(x=x, y=y, z = Temperature),
        #                   size = 6,
        #                   #stroke = 0.2,
        #                   breaks = MakeBreaks(binwidth = 1)) +
        geom_point(data = filter(samp_locs, x >= 21),
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        annotate("text",
                 label = expression('bold(paste("Temperature (", degree,"C)"))'),
                 x = 33,
                 y = -7.6,
                 size = 9,
                 fontface =2,
                 colour = "black", parse = TRUE) +
        labs(x = "Distance From Entrance (km)",
             y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_blank(),
              axis.ticks.length.y.left = (unit(2, "mm")),
              axis.ticks.length.y.right = (unit(2, "mm")),
              #axis.ticks.x = element_blank(),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              plot.title = element_text(hjust=0.5, vjust=0.5, face='bold', size = 24),
              plot.subtitle = element_text(hjust=0.5, vjust=0.5, size = 22),
              legend.background = element_rect(fill = "transparent"),
              legend.direction = "horizontal",
              legend.position = c(0.65, 0.22),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      # full length plots
      #create list of plotGrobs
      plta <- lapply(list(salPlot, doPlot, chlorPlot, tempPlot), ggplotGrob)
      #rbind (i.e. 1 column) size arg matters!
      surfers <- rbind(plta[[1]], plta[[2]], plta[[3]], plta[[4]], size = "first")
      pdf_name <- paste0(path, "/plots/", "swan_", ymd(pair), "_surfer.pdf")
      cat(paste0(pdf_name,"\n"))
      #add margin padding coarse but effective
      surfers_pad <- gtable::gtable_add_padding(surfers, padding = unit(c(1,4,3,4), "cm"))

      ggsave(plot = grid.draw(surfers_pad), filename = pdf_name, width=28, height=18)

      # for zoomed at narrows
      #create list of plotGrobs
      pltb <- lapply(list(salPlotZ, doPlotZ, chlorPlotZ, tempPlotZ), ggplotGrob)
      #rbind (i.e. 1 column) size arg matters!
      surfersZ <- rbind(pltb[[1]], pltb[[2]], pltb[[3]], pltb[[4]], size = "first")
      pdf_nameZ <- paste0(path, "/plots/", "swan_middle_upper_", lubridate::ymd(pair), "_surfer.pdf")
      cat(paste0(pdf_nameZ,"\n"))
      #add margin padding coarse but effective
      surfersZ_pad <- gtable::gtable_add_padding(surfersZ, padding = unit(c(1,4,3,4), "cm"))

      ggsave(plot = grid.draw(surfersZ_pad), filename = pdf_nameZ, width=28, height=18)
    }

  } else {
    stop("Function expecting only 2 excel workbooks for one monitoring period")
  }

}
