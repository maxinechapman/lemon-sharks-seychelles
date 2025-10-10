## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## source functions for shark acoustic data processing
## MC Jan 2023, adapted from LP code
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_adj_matrix <- function(shark_data, receiver_data_meta, ID_hashmap){
  
  #1. set up loop
  shark_number <- shark_data$FishID[1]
  loop <- nrow(shark_data)-1
  shark_matrix <- matrix(0, nrow = 89, ncol = 89)
  
  #2. create adjacency matrix
  for(i_row in 1:loop){
    
    #remove rows with NA
    shark_data <- na.omit(shark_data)
    
    
    #identify starting receiver ID
    start_lat <- shark_data[i_row, "Lat"]
    start_long <- shark_data[i_row, "Long"]
    start_ID <- ID_hashmap[[c(start_lat, start_long)]]
    if (is.null(start_ID)){
      stop("at row ",i_row," we got start_ID NULL for shark ", shark_number)
    }
    if (is.na(start_ID)){
      stop("at row", i_row, "we got start_ID NA for shark ", shark_number)
    }
    
    
    #identify ending receiver ID
    i_row <- i_row+1
    end_lat <- shark_data[i_row, "Lat"]
    end_long <- shark_data[i_row, "Long"]
    end_ID <- ID_hashmap[[c(end_lat, end_long)]]
    if (is.null(end_ID)){
      stop("at row ",i_row," we got end_ID NULL for shark", shark_number)
    }
    if (is.na(end_ID)){
      stop("at row", i_row, "we got end_ID NA for shark", shark_number)
    }
    
    #add weighted connection to matrix and make it symmetric
    shark_matrix[start_ID, end_ID] <- shark_matrix[start_ID, end_ID] + 1
    shark_matrix[end_ID, start_ID] <- shark_matrix[end_ID, start_ID] + 1
    
  }
  
  return(shark_matrix)
  
}

calc_mcp_area <- function(tmp_shark_data, percent_val) {
  ## this function calculates MCP area (minimum convex polygon area) for a shark dataset
  ## which can be either lifetime or already subset to just one year
  
  ## TESTING - leave this section commented out
  # tmp_shark_data <- subset(shark_data, FishID == 144 & year == 2014)
    
    ## get locations of receivers used by this shark
    shark_locs <- data.frame(id = tmp_shark_data$FishID, lon = tmp_shark_data$Long, lat = tmp_shark_data$Lat)
    
    ## make spatial
    sp::coordinates(shark_locs) <- ~lon + lat
    sp::proj4string(shark_locs) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    shark_locs_proj <- as.data.frame(sp::spTransform(shark_locs, sp::CRS("+proj=laea +lat_0=-5.4 +lon_0=53 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")))
    coord <- shark_locs_proj[ , c("coords.x1", "coords.x2")]  # extract coords using col names
    id <- as.data.frame(shark_locs_proj[ , "id"]) 
    shark_spdf_proj <- sp::SpatialPointsDataFrame(coords=coord, data=id, proj4string = sp::CRS("+proj=laea +lat_0=-5.4 +lon_0=53 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
    ## run MCP
    shark_mcp <- adehabitatHR::mcp(shark_spdf_proj, percent= percent_val, unout="km2")  # ignore the warning
    shark_area <- shark_mcp$area  # then you can just extract the area
    
  
  ## RETURN the MCP area in km2 - this is the function OUTPUT
  return(shark_area)
  
  ## clean up from testing - leave this commented out
  #rm(tmp_shark_data, shark_locs, shark_locs_proj, coord, id, shark_spdf_proj, shark_mcp, shark_area)
  
}


make_lifetime_igraphs <- function(shark_name, shark_matrix, vertex_colour_list, receiver_data_meta){
  #this function creates a small and a large plot of all the data for a given shark and saves them in one pdf document
  
  shark_number <- shark_name$FishID[1]
  
  #receiver IDs and geographical locations
  meta <- data.frame(receiver_data_meta$ReceiverID, receiver_data_meta$Lat, receiver_data_meta$Long)
  
  #connections between receivers
  g <- graph_from_adjacency_matrix(shark_matrix, diag = FALSE, mode = "undirected", weighted = TRUE)
  
  #layout
  lo <- as.matrix(meta[,c(3,2)])
  
  #information for legends
  habitats <- c("Lagoon", "Coastal Reef", "Plateau", "Dropoff")
  habitat_colours <- c("green", "pink", "orange", "purple")
  map_key <- data.frame(habitats, habitat_colours)
  
  #titles
  text <- paste0("Shark ", shark_number, "(all data)")
  
  #PLOT 1 - ISLAND
  #save plot
  all_island_data_path <- paste0("../results/igraph_plots/all_data/combined_plots/lemon_", shark_number, "_all_combined.pdf")
  pdf(file = all_island_data_path)
  
  #plot
  map <- brick("../resources/Georeferenced island imagery/satellite/islands11.tif")
  plotRGB(map)
  E(g)$color <- "red"
    
  plot.igraph(g,
              layout=lo, add = TRUE,
              rescale=FALSE,
              annotate.plot = TRUE,
              xlim = c(53.25,53.40),
              ylim = c(-6,-5),
              vertex.size = 0.25,
              vertex.label.cex = 1,
              vertex.color = paste(vertex_colour_list),
              edge.arrow.size = 0,
              edge.width = 1)
  
  #add title
  legend("top", legend = text, bty = "n", cex = 2)
  
  #add key
  legend("bottomright", legend = map_key$habitats , pt.bg= map_key$habitat_colours,
         pch=21, col="#777777", title = "Habitat Type", bg = "white")
  
  
  #PLOT 2 - ALL RECEIVERS
  # all_data_path <- file.path("../results/igraph_plots/all_data/lemon_", shark_number, "_all.pdf")
  # pdf(file = all_data_path, height = 23, width = 70)
  
  #plot
  E(g)$color <- "red"
  plot.igraph(g,
              layout=lo, 
              rescale=FALSE,
              annotate.plot = TRUE,
              xlim = c(53.25,53.40),
              ylim = c(-6,-5),
              vertex.size = 0.65,
              #vertex.label= NA,
              vertex.label.cex = 0.2,
              #vertex.frame.color = NA,
              vertex.frame.width = 0.1,
              vertex.color = paste(vertex.colours),
              edge.arrow.size = 0,
              edge.width = 1)
    
    #add title
    legend("top", legend = text, bty = "n", cex = 2)
    
    #add key
    legend("bottomright", legend = map_key$habitats , pt.bg= map_key$habitat_colours,
           pch=21, col="#777777", title = "Habitat Type", bg = "white")
    
    #stop save plot
    dev.off()
    
}

make_yearly_igraphs <- function(shark_data_year, shark_matrix, receiver_data_meta, vertex_colour_list, year, shark_number){
  #this function creates a small and large igraph plot for given year of data for an individual shark. 
  #The plot is then saved in the specified location
  
  # check shark and year for debugging
  print(shark_number)
  print(year)
  
  #check if there is data for that year
  if(nrow(shark_data_year) > 1){
    
    #receiver IDs and geographical locations
    meta <- data.frame(receiver_data_meta$ReceiverID, receiver_data_meta$Lat, receiver_data_meta$Long)
    
    #connections between receivers
    g2 <- igraph::graph_from_adjacency_matrix(shark_matrix, diag = FALSE, mode = "undirected", weighted = TRUE)  
    #here we specify weighted = T, and then
    #the values from the matrix are used as the weight the edge between two nodes, rather than the number of edges between the two nodes.
    
    #layout
    lo <- as.matrix(meta[,c(3,2)])
    
    #save plot
    data_path <- file.path("..", "results", "igraph_plots", year, paste0("lemon_", shark_number, "_", year, "_island.pdf"))
    
    #need to create the folder here if it doesn't already exist, otherwise you get Error : cannot open file.
    if (! dir.exists(file.path("..", "results", "igraph_plots", year))) {
      dir.create(file.path("..", "results", "igraph_plots", year))
    }
    
    #open PDF file in which to save the plot
    pdf(file = data_path)
    
    #PLOT 1: ISLAND PLOT
    map <- raster::brick("../resources/Georeferenced island imagery/satellite/islands11.tif")
    plotRGB(map)
    igraph::E(g2)$color <- "red"
      
    plot.igraph(g2,
                layout=lo, add = TRUE,
                rescale=FALSE,
                annotate.plot = TRUE,
                xlim = c(53.25,53.40),
                ylim = c(-6,-5),
                vertex.size = 0.25,
                vertex.label.cex = 1,
                vertex.color = paste(vertex_colour_list),
                edge.arrow.size = 0,
                edge.width = 1)
    
    #add title
    text <- paste0("Shark ", shark_number, " - ", year)
    legend("top", legend = text, bty = "n")
    
    #add colour coding
    habitats <- c("Lagoon", "Coastal Reef", "Plateau", "Dropoff")
    habitat_colours <- c("green", "pink", "orange", "purple")
    map_key <- data.frame(habitats, habitat_colours)
    legend("bottomleft", legend = map_key$habitats , pt.bg= map_key$habitat_colours,
           pch=21, col="#777777", title = "Habitat Types")
    
    #PLOT 2: BIG PLOT
    E(g2)$color <- "red"
      plot.igraph(g2,
                  layout=lo, 
                  rescale=FALSE,
                  annotate.plot = TRUE,
                  xlim = c(53.25,53.40),
                  ylim = c(-6,-5),
                  vertex.size = 0.65,
                  vertex.label.cex = 0.2,
                  vertex.frame.width = 0.1,
                  vertex.color = paste(vertex_colour_list),
                  edge.arrow.size = 0,
                  edge.width = 1)
      
      #add title
      #text <- paste0("Shark ", shark_number, " - ", year)
      legend("top", legend = text, bty = "n", cex = 2)
      
      #add colour coding
      # habitats <- c("Lagoon", "Coastal Reef", "Plateau", "Dropoff")
      # habitat_colours <- c("green", "pink", "orange", "purple")
      # map_key <- data.frame(habitats, habitat_colours)
      legend("bottomleft", legend = map_key$habitats , pt.bg= map_key$habitat_colours, pch=21, col="#777777", title = "Habitat Types", cex = 2)
      
      #stop save plot
      dev.off()
      
      # DEBUG: Confirm that the PDF device was closed and where the graph was saved
      print(paste("Saved plot for shark", shark_number, "in year", year))
  }
}

make_big_plot <- function(shark_data_year, shark_matrix, receiver_data_meta, vertex_colour_list) {
  #this function creates an igraph for a given year of shark data in a more visually aesthetic format than the original large igraphs
  
  #get shark ID and year
  shark_number <- shark_data_year$FishID[1]
  year <- year(shark_data_year$ISO.date[1])
  
  #creating data for map key
  habitats <- c("Lagoon", "Coastal Reef", "Plateau", "Dropoff")
  habitat_colours <- c("green", "pink", "orange", "purple")
  map_key <- data.frame(habitats, habitat_colours)
  
  #creating data for igraph
  #receiver IDs and geographical locations
  meta <- data.frame(receiver_data_meta$ReceiverID, receiver_data_meta$Lat, receiver_data_meta$Long)
  
  #connections between receivers
  g <- graph_from_adjacency_matrix(shark_matrix, diag = FALSE, mode = "undirected", weighted = TRUE)
  
  #layout
  lo <- as.matrix(meta[,c(3,2)])
  
  #make filepath
  data_path <- file.path("..", "results", "igraph_plots", year, paste0("lemon_", shark_number, "_", year, "_big.pdf"))
  
  #saving
  pdf(data_path, width = 10)
  
  #plotting
  #map <- brick("../resources/ETOPO_2022_bathy_extract.tiff")
  
  map <- raster::raster("../resources/ETOPO_2022_bathy_extract.tiff")
  pal <- colorRampPalette(c("#152E47","#51AEF6"))  # you can change these colours if you like
  plot(map, col = pal(50))  # plot using 50 colour slices from this palette
  points(receiver_data_meta$Long, receiver_data_meta$Lat, col = paste(vertex_colour_list), pch=16, cex=0.8)
  E(g2)$color <- "red"
    plot.igraph(g2,
                layout=lo, 
                rescale=FALSE,
                annotate.plot = TRUE,
                add = TRUE,
                xlim = c(53.25,53.40),
                ylim = c(-6,-5),
                vertex.size = 0.65,
                vertex.label.cex = 0.2,
                vertex.frame.width = 0.1,
                vertex.color = paste(vertex_colour_list),
                edge.arrow.size = 0,
                edge.width = 1)
  
  
  #add legend
  legend("bottomright", legend = map_key$habitats , pt.bg= map_key$habitat_colours,
         pch=21, col="#777777", title = "Habitat Type", border = "white")
  #legend("top", legend = "Full Receiver Map", bty = "n", cex = 1, text.col = "white")
  
  #turn off save
  dev.off()
  
}

make_receiver_map <- function(receiver_data_meta, vertex_colours_list, title, map_key, key_title, filename) {
  #this function makes an igraph with no connections to show the locations of the receivers
  #receivers can be colour coded differently based on the function inputs
  
  #make graph objects
  empty_matrix <- matrix(0, nrow = nrow(receiver_data), ncol = nrow(receiver_data)) #make a matrix with no connections
  meta <- data.frame(receiver_data$ReceiverID, receiver_data$Lat, receiver_data$Long, receiver_data$Habitat) #get receiver IDs and geographical locations
  g <- graph_from_adjacency_matrix(empty_matrix) #connections between receivers
  lo <- as.matrix(meta[,c(3,2)]) #layout
  
  #saving
  filepath <- paste0("../results/igraph_plots/receiver_plots/", filename, ".pdf")
  pdf(filepath)
  
  #plot small receiver plot
  map_small <- brick("../resources/Georeferenced island imagery/satellite/islands11.tif")
  plotRGB(map_small)
  plot.igraph(g,
              layout=lo, add = TRUE,
              rescale=FALSE,
              annotate.plot = TRUE,
              xlim = c(53.25,53.40),
              ylim = c(-6,-5),
              vertex.size = 0.25,
              vertex.label.cex = 1,
              vertex.color = paste(vertex_colours_list),
              edge.arrow.size = 0,
              edge.width = 1)
  #add scalebar
  raster::scalebar(d = 2,
                   type = "bar",
                   divs = 2, 
                   below = "km", 
                   lonlat = TRUE,
                   label = c(0,1,2), 
                   lwd = 1.5)
  
  legend("bottomright", legend = map_key[, 1] , pt.bg= map_key[, 2],
         pch=21, col="#777777", title = key_title, bg = "white")
  legend("top", legend = title, bty = "n", cex = 2)
  
  #plot large receiver plot
  map_big <- brick("../resources/ETOPO_2022_bathy_extract.tiff")
  pal <- colorRampPalette(c("#152E47","#51AEF6"))  # colours for plot
  E(g)$color <- "red"
  
  plot(map_big, col = pal(50))
  plot.igraph(g,
              layout=lo,
              add = TRUE,
              rescale = FALSE,
              annotate.plot = TRUE,
              xlim = c(53.25,53.40),
              ylim = c(-6,-5),
              vertex.size = 2,
              vertex.label = NA,
              vertex.color = paste(vertex_colours_list),
              vertex.frame.color = paste(vertex_colours_list),
              edge.arrow.size = 0,
              edge.width = 1)
  #add scalebar
  raster::scalebar(d = 50,
                   type = "bar",
                   divs = 4, 
                   below = "km", 
                   lonlat = TRUE,
                   label = c(1, 25, 50), 
                   lwd = 1.5,
                   col=c("black"))
  
  legend("bottomright", legend = map_key[, 1] , pt.bg= map_key[, 2],
         pch=21, col="#777777", title = key_title, bg = "white")
  legend("top", legend = title, bty = "n", cex = 2, text.col = "white")
  
  #end save
  dev.off()
}


elapsed_months <- function(end_date, start_date) {
  #this function calculates months as a whole number(rounds down)
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

calc_shark_size <- function(shark_ID, current_year, growth_rate, shark_summary_metadata){
  #this function calculates a sharks current size based on its initial size and growth rate
  
  #get initial size and year from summary dataset
  initial_size <- shark_summary_metadata[shark_summary_metadata$FishID == shark_ID, "Size"]
  initial_year <- lubridate::year(shark_summary_metadata[shark_summary_metadata$FishID == shark_ID, "Date.released"])
  
  #calculate current shark size
  current_size <- initial_size + growth_rate*(current_year-initial_year)
  return(current_size)
  
}
