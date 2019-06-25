# Main function for creating canning surfer plots

#' Reads in external data, combines with internal datasets and produces Canning
#' River surfer plots
#'
#' \code{canning_surfR} takes a file path to Canning River sonde output and creates
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
#' @param obac Character string ("green", "blue" or "red") indicating
#'     oxygenation plant status at site BAC.
#'
#' @param onic Character string ("green", "blue" or "red") indicating
#'     oxygenation plant status at site NIC.
#'
#' @return A pdf format four panel surfer plot of the Canning River.
#'
#' @examples
#' \dontrun{
#' canning_surfR(path = "Z:/DEC/MonitoringProgram/Data", obac = "green", onic = "red")
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
canning_surfR <- function(path, obac, onic){
  locations <- data_finder(path, river = "c")

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
        janitor::clean_names()

      # join sites to WQ data
      comb_data <- dplyr::left_join(samp_data, C_sitesdf, by = "site")
      d_reduced <- comb_data %>%
        dplyr::select(site, sal_ppt, do_mg_l, c, chl_ug_l, dep_m,
                      dist_bridg, max_depth)
      # Set up labels and params to plot
      sparams <- c("Salinity", "Dissolved_Oxygen", "Temperature", "Chlorophyll")

      # filter sampling data for separate interpolations
      d_all <- d_reduced[complete.cases(d_reduced),] %>%
        dplyr::mutate(y = -1 * dep_m, x = dist_bridg/1000)

      d_low <- d_reduced[complete.cases(d_reduced),] %>%
        dplyr::mutate(y = -1 * dep_m, x = dist_bridg/1000) %>%
        dplyr::filter(x < 11.3)

      d_up <- d_reduced[complete.cases(d_reduced),] %>%
        dplyr::mutate(y = -1 * dep_m, x = dist_bridg/1000) %>%
        dplyr::filter(x > 11.3)

      # Create separate lists for storing interpolations
      vals <- c("sal_ppt", "do_mg_l", "c", "chl_ug_l")
      idw_list_all <- vector("list", length(vals))
      idw_list_weir <- vector("list", length(vals))
      names(idw_list_all) <- vals
      names(idw_list_weir) <- vals #will combine upper and lower here


      # for all
      for(i in seq_along(vals)){
        val <- vals[i]
        d1 <- d_all[,c("x", "y", val)]
        names(d1)[3] <- "value"
        sp::coordinates(d1) <- ~x + y

        idw <- gstat::idw(formula = value ~ 1,
                          locations = d1,
                          newdata = C_grd_all,
                          idp = 4)

        idw1_r <- raster::raster(idw)
        idw1_r_class <- raster::reclassify(idw1_r, reclass_matrices[[i]])
        idw1_sp <- raster::rasterToPoints(idw1_r_class, spatial = TRUE)
        idw1_df <- data.frame(idw1_sp)[-4]# ditch option
        names(idw1_df)[1] <- sparams[i]
        # idw1_df <- idw1_df[-(1:171),] # strip weird close to -0.1 vals
        idw_list_all[[i]] <- idw1_df
      }

      # for weir
      for(i in seq_along(vals)){
        val <- vals[i]
        d1 <- d_low[,c("x", "y", val)]
        d2 <- d_up[,c("x", "y", val)]
        names(d1)[3] <- "value"
        names(d2)[3] <- "value"
        sp::coordinates(d1) <- ~x + y
        sp::coordinates(d2) <- ~x + y

        lc_idw1 <- gstat::idw(formula = value ~ 1,
                              locations = d1,
                              newdata = C_grd_low,
                              idp = 4)
        uc_idw1 <- gstat::idw(formula = value ~ 1,
                              locations = d2,
                              newdata = C_grd_up,
                              idp = 4)

        # create as raster to reclass values (bin)
        lc_idw_r1 <- raster::raster(lc_idw1)
        uc_idw_r1 <- raster::raster(uc_idw1)

        comb_idw_r1 <- raster::merge(lc_idw_r1, uc_idw_r1)

        # idw1_r <- raster(idw)
        idw1_r_class <- raster::reclassify(comb_idw_r1, reclass_matrices[[i]])
        idw1_sp <- raster::rasterToPoints(idw1_r_class, spatial = TRUE)
        idw1_df <- data.frame(idw1_sp)[-4]# ditch option
        names(idw1_df)[1] <- sparams[i]
        # idw1_df <- idw1_df[-(1:171),] # strip weird close to -0.1 vals
        idw_list_weir[[i]] <- idw1_df
      }

      # make sample collection points
      samp_locs <- comb_data %>%
        dplyr::mutate(dist_bridg = dist_bridg/1000) %>%
        dplyr::rename(x = dist_bridg, y = dep_m) %>%
        dplyr::select(site, x, y)

      samp_labels <- c("SCB2", "SAL", "RIV", "CASMID", "KEN", "BAC", "NIC", "ELL")
      samp_labels_locs <- C_sitesdf %>%
        dplyr::filter(site %in% samp_labels)

      # Logic for weir no weir
      cannoxy <- c("KENU300", "BACD500", "BACD300", "BACU300", "PO2", "GRE", "KS7",
                   "MASD50", "NICD200", "KS9", "PAC", "MACD50")
      if(sum(cannoxy %in% samp_locs$site) > 0 ){
        bottom <- C_bottom_weir
        interp <- idw_list_weir
      } else {
        bottom <- C_bottom_open
        interp <- idw_list_all
      }

      salPlot <- ggplot()+
        geom_raster(data = interp[[1]],
                    aes(x=x, y=y, fill = factor(Salinity))) +
        scale_x_continuous(limits = c(0.5, 16),
                           expand = c(0, 0)) +
        stat_contour2(data = interp[[1]], aes(x=x, y=y, z = Salinity),
                      colour = "grey10",
                      breaks = MakeBreaks(binwidth = 2)) +
        scale_fill_manual(values = surfer_cols("sal"),
                          guide = guide_legend(reverse=T),
                          name = "Salinity\n(ppt)") +
        geom_text_contour(data = interp[[1]],
                          aes(x=x, y=y, z = Salinity),
                          skip = 1,
                          size = 6,
                          check_overlap = TRUE,
                          stroke = 0.2,
                          breaks = MakeBreaks(binwidth = 2),
                          nudge_y = 0.5) +
        geom_polygon(data = bottom,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        geom_point(data = samp_locs,
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        geom_text(data = samp_labels_locs,
                  aes(x = dist_bridg/1000, y = 0.7, label = site, fontface = 2),
                  size = 4.5,
                  colour = "black",
                  alpha = 1,
                  check_overlap = TRUE) +
        annotate("text",
                 label = "Salinity (ppt)",
                 x = 2,
                 y = -6.2,
                 size = 9,
                 fontface = 2,
                 colour = "black") +
        labs(title = paste("Canning River Estuary - Physical-Chemical Profile -",
                           sdate),
             y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border=element_blank(),
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
              legend.position = c(0.3, 0.15),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      doPlot <- ggplot()+
        geom_raster(data = interp[[2]],
                    aes(x=x, y=y, fill = factor(Dissolved_Oxygen))) +
        scale_x_continuous(limits = c(0.5, 16),
                           expand = c(0, 0)) +
        stat_contour2(data = interp[[2]],
                      aes(x=x, y=y, z = Dissolved_Oxygen),
                      colour = "grey10",
                      breaks = MakeBreaks(binwidth = 1)) +
        scale_fill_manual(values = surfer_cols("do"),
                          guide = guide_legend(reverse=T),
                          name = "Dissolved\nOxygen\n(mg/L)") +
        geom_text_contour(data = interp[[2]],
                          aes(x=x, y=y, z = Dissolved_Oxygen),
                          skip = 1,
                          check_overlap = TRUE,
                          size = 6,
                          stroke = 0.2,
                          breaks = MakeBreaks(binwidth = 1)) + #,nudge_y = 0.5
        geom_polygon(data = bottom,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        geom_point(data = samp_locs,
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        geom_point(data = C_oxy_locs,
                   aes(x = x, y = y),
                   size = 6,
                   colour = "black",
                   bg = C_oxy_locs$c,
                   shape = 24) +
        annotate("text",
                 label = "Dissolved Oxygen (mg/L)",
                 x = 2.1,
                 y = -6.2,
                 size = 9,
                 fontface = 2,
                 colour = "black") +
        labs(y = "Depth (m)") +
        annotation_custom(grob = oxY_grob, xmin = 7, xmax = 10, ymin = -6.3, ymax = -4) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border=element_blank(),
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
              legend.position = c(0.3, 0.15),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      chlorPlot <- ggplot()+
        geom_raster(data = interp[[4]],
                    aes(x=x, y=y, fill = factor(Chlorophyll)),
                    alpha = 0.5) +
        scale_x_continuous(limits = c(0.5, 16),
                           expand = c(0, 0)) +
        stat_contour2(data = interp[[4]],
                      aes(x=x, y=y, z = Chlorophyll),
                      colour = "grey70",
                      breaks = chl_brk) +
        scale_fill_manual(values = surfer_cols("chl"),
                          guide = guide_legend(reverse=T),
                          name = "Chlorophyll\n(ug/L)",
                          labels = c("20", "40", "60", "80", "120", "200",
                                     "> 200")) +
        geom_text_contour(data = interp[[4]],
                          aes(x=x, y=y, z = Chlorophyll),
                          skip = 2,
                          size = 6,
                          check_overlap = TRUE,
                          stroke = 0.2,
                          breaks = as.numeric(chl_brk)) +
        geom_polygon(data = bottom,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        geom_point(data = samp_locs,
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        annotate("text",
                 label = expression('bold(paste("Chlorophyll (", mu,"/L)"))'),
                 x = 2,
                 y = -6.2,
                 size = 9,
                 colour = "black", parse = TRUE) +
        labs(x = "Distance From Entrance (km)",
             y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border=element_blank(),
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
              legend.position = c(0.3, 0.15),
              legend.key.size =  unit(8, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                                   label.position = "bottom"))

      tempPlot <- ggplot()+
        geom_raster(data = interp[[3]],
                    aes(x=x, y=y, fill = factor(Temperature))) +
        scale_x_continuous(limits = c(0.5, 16),
                           expand = c(0, 0)) +
        scale_y_continuous(expand = expand_scale(mult = c(0, .05)))+
        stat_contour2(data = interp[[3]],
                      aes(x=x, y=y, z = Temperature),
                      colour = "grey10",
                      breaks = MakeBreaks(binwidth = 1)) +
        scale_fill_manual(values = surfer_cols("temp"),
                          guide = guide_legend(reverse=T),
                          name = "Temp") +
        geom_text_contour(data = interp[[3]],
                          aes(x=x, y=y, z = Temperature),
                          size = 6,
                          stroke = 0.2,
                          breaks = MakeBreaks(binwidth = 1)) + #,nudge_y = 0.5
        geom_polygon(data = bottom,
                     aes(x=x, y=y), fill = "grey90", colour = "grey20") +
        geom_point(data = samp_locs,
                   aes(x = x, y = - y),
                   colour = "black",
                   size = 0.5) +
        annotate("text",
                 label = expression('bold(paste("Temperature (", degree,"C)"))'),
                 x = 2,
                 y = -6.2,
                 size = 9,
                 fontface = 2,
                 colour = "black", parse = TRUE) +
        labs(x = "Distance From Entrance (km)",
             y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_blank(),
              axis.ticks.length.y.left = (unit(2, "mm")),
              axis.ticks.length.y.right = (unit(2, "mm")),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              plot.title = element_text(hjust=0.5, vjust=0.5, face='bold'),
              plot.subtitle = element_text(hjust=0.5, vjust=0.5),
              panel.border=element_blank(),
              legend.background = element_rect(fill = "transparent"),
              legend.direction = "horizontal",
              legend.position = c(0.3, 0.15),
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
      pdf_name <- paste0(path, "/plots/", "canning_", ymd(pair), "_surfer.pdf")
      cat(paste0(pdf_name,"\n"))
      #add margin padding coarse but effective
      surfers_pad <- gtable::gtable_add_padding(surfers, padding = unit(c(1,4,3,4), "cm"))

      ggsave(plot = grid.draw(surfers_pad), filename = pdf_name, width=28, height=18)
    }
  } else {
    stop("Function expecting only 2 excel workbooks for one monitoring period")
  }

}
