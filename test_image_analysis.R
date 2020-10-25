# Counting insects with image segmentation 
# Alexander Piper 2020/10/05
# Partly modified from http://rforbiochemists.blogspot.com/2016/05/counting-and-identifying-stained-cells.html

# Load libraries
library(magick)
library (EBImage)
library(tidyverse)
library(tidymodels)
library(fs)
library(animation)

#Source themes
source('R/themes.R')

# Make an output directory for the gifs
dir.create("gifs")

# Segmentation function
segment_insects <- function(path, out.gif=NULL,
                            cmean_smooth=1.5, cmean_min_pixels=1,
                            smooth_size = 15, smooth_sigma=11,
                            thresh_size=51, thresh_sigma=0.3){
  # Check inputs
  if(!file.exists(path)){
    stop("file doesnt exists")
  }
  if(!is.null(out.gif) && !dir.exists(out.gif)){
    message("Creating directory ", out.gif, " to store output gifs")
    dir.create(out.gif)
  }
  params <- c(cmean_smooth, cmean_min_pixels, smooth_size, smooth_sigma, thresh_size, thresh_sigma)
  names(params) <- c("cmean_smooth", "cmean_min_pixels", "smooth_size",
                     "smooth_sigma", "thresh_size", "thresh_sigma")
  if(!all(is.numeric(params))){
    stop(names(params[suppressWarnings(is.na(as.numeric(params)))]), " must be numeric")
  }
  # Read in image
  img_full <- image_read(path)
  img_name <- basename(path) %>% str_remove(".jpg$|.png$|.bmp$")
  
  # Read in image & Convert to grayscale
  img_grey <- image_convert(img_full, colorspace = "Gray")
  
  # Fuzzy c-means clustering to differentiate insect from background
  img_cmeans <- image_fuzzycmeans(img_grey, smoothing=cmean_smooth, min_pixels = cmean_min_pixels)
  
  # Get frame 1 (Darkest frame from segmentation - should contain insects if photographed on white background)
  img_split <- img_cmeans %>%
    image_split(keep_color = FALSE)

  img_f1 <- image_read(img_split[[1]])
  
  # Write out temp file, read back in as a EBImage object
  image_write(img_f1, "tmp.png", format="png")
  img_f1 <- readImage("tmp.png")
  img_f1 <- channel(img_f1,"gray")
  try(file.remove("tmp.png"), silent=TRUE)
  
  # Blur to smooth image 
  brush_blur <- makeBrush(size = smooth_size, shape = 'gaussian', sigma =smooth_sigma) 
  img_smooth <- filter2(img_f1*2, brush_blur) # Smooth image, multiply by 2 to increase brightness
  img_smooth <- channel(img_smooth,"gray")
  
  # Adaptive thresholding to better differentiate blobs from background and each other
  brush_threshold <- makeBrush(size=thresh_size, shape="disc", sigma=thresh_sigma)
  brush_threshold <- brush_threshold / sum(brush_threshold)
  background_remove <- filter2(img_smooth, brush_threshold)
  offset = 0.15   # Add an offset
  img_segments <- img_smooth > background_remove + offset
  
  # Colour final segments
  final_segments <- colorLabels(bwlabel(img_segments))
  
  # Collect counts in output gif
  out <- data.frame(
    sample = img_name, 
    count = max(bwlabel(img_segments))
    )
  
  # Make a list containing the images from each step
  if(!is.null(out.gif)){
    plotlist <- list(
      img_full,
      img_grey,
      img_cmeans,
      img_f1,
      img_smooth,
      img_segments,
      final_segments
    )
    
    # Plot as a gif
    animation::saveGIF(
      expr = {
        purrr::walk(plotlist, ~{
          plot(.x)
          animation::ani.pause(2)
        })
      },
      movie.name = suppressWarnings(normalizePath(paste0(out.gif,"/",img_name,".gif"))),
      ani.width=1280 ,
      ani.height=1280,
      interval=2,
      autobrowse=FALSE
    )
    try(dev.off(), silent=TRUE)
  }
  return(out)
}

# List images to process
# NOTE: these have already been cropped to the outline of the petri dish
images <- fs::dir_ls("C:/Users/ap0y/Dropbox/workprojects/PHD/Metabarcoding/Photos/", glob="*.jpg", recurse = TRUE)

images <- images[str_detect(images, "Tatura_week8_fruitcrush")]

counts <- images %>%
  purrr::map_dfr(function(x){segment_insects(x,cmean_smooth = 1.5, cmean_min_pixels = 1,  out.gif="gifs")})

#Params to check
cmean_smooth=1.5
cmean_min_pixels=1
smooth_size = 15
smooth_sigma=11
thresh_size=50
thresh_sigma=0.3


counts
# compare to known counts
count_comparison <- counts %>%
  mutate(type = case_when(
    str_detect(sample, "Carpophilus_larval")  ~ "CarpLarval",
    str_detect(sample, "carpophilus_mock")  ~ "CarpMock",
    str_detect(sample, "carpophilus_trap")  ~ "CarpTrap",
    str_detect(sample, "^D[0-9]")  ~ "DrosMock",
    str_detect(sample, "Dros_larval")  ~ "DrosLarv",
      str_detect(sample, "fruitcrush")  ~ "FF",
    str_detect(sample, "ACV")  ~ "ACV",    
    str_detect(sample, "sachet")  ~ "sachet",
    str_detect(sample, "DC$")  ~ "DC",    
    str_detect(sample, "ACV$")  ~ "ACV",
    str_detect(sample, "SPD$")  ~ "SPD"
  )) %>%
  mutate(sample_name = case_when(
    str_detect(sample, "Carpophilus_larval")  ~ str_replace(sample, "Carpophilus_larval", "CML"),
    str_detect(sample, "carpophilus_mock")  ~ str_replace(sample,"carpophilus_mock", "CM"),
    str_detect(sample, "carpophilus_trap")  ~ str_replace(sample,"carpophilus_trap", "CT"),
    str_detect(sample, "^D[0-9]")  ~sample,
    str_detect(sample, "Dros_larval")  ~ str_replace(sample,"Dros_larval", "DLarv"),
    str_detect(sample, "Tatura_week")  ~ sample %>% str_replace("Tatura_week", "T") %>% str_remove("_"),
    str_detect(sample, "redhill_week")  ~ sample %>% str_replace("redhill_week", "M") %>% str_remove("_")
  ),
  sample_name=sample_name %>%
    str_replace("Fruitcrush", "FF") %>%
    str_replace("fruitcrush", "FF") %>%
    str_remove("'"))  %>%
  left_join(read_csv("sample_data/expected_quant_merged.csv") %>%
  dplyr::rename(sample_name = Sample) %>%
  pivot_longer(-sample_name,
               names_to= "taxon",
               values_to= "expected") %>%
  group_by(sample_name) %>%
  summarise(true_count = sum(expected))
  ) %>%
  filter(!type =="sachet") %>%
  filter(!is.na(type))

# calulate rmse
count_rmse <- count_comparison %>%
  group_by(type) %>%
  rmse(truth=true_count, estimate=count) %>%
  mutate(.estimate = round(.estimate, 2))


gg.count_comparison <- count_comparison %>%
  left_join(count_rmse %>% select(type, RMSE =.estimate)) %>%
  mutate(type = paste0(type, "\n(RMSE: ", RMSE, ")")) %>%
  ggplot(aes(x = true_count, y=count, colour=type))+
  geom_point(alpha = 0.5, size=2) +
  geom_abline(lty = 2, colour = "gray80", size = 1) +
  coord_fixed() + 
  base_theme +
  scale_colour_brewer(palette="Paired") +
  theme(panel.grid = element_line(size = rel(1)),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(face="italic")) +
  labs(x = "Morphological count",
       y = "Segmentation count",
       colour = "Community type")


# Grid search tuning
grid_params <- tidyr::crossing(
  images,
  cmean_smooth = seq(0.1, 3, 0.2),
  cmean_min_pixels=seq(0.1, 3, 0.2)
) %>%
  dplyr::filter(str_detect(images, "Dros_larval2")) 
#
#%>%
#  slice_sample(n=20)

grid_search <- purrr::pmap_dfr(grid_params, ~with(list(...), {
    segment_insects(images, cmean_smooth =cmean_smooth, cmean_min_pixels = cmean_min_pixels,  out.gif=NULL)})) %>%
  bind_cols(grid_params)



grid_search%>%
  left_join(read_csv("sample_data/expected_quant_merged.csv") %>%
              dplyr::rename(sample_name = Sample) %>%
              pivot_longer(-sample_name,
                           names_to= "taxon",
                           values_to= "expected") %>%
              group_by(sample_name) %>%
              summarise(true_count = sum(expected))
  ) %>%
  filter(!type =="sachet") %>%
  filter(!is.na(type))

plot(grid_search$cmean_smooth, grid_search$count)


#Params to check
cmean_smooth=1.5
cmean_min_pixels=1
smooth_size = 15
smooth_sigma=11
thresh_size=50
thresh_sigma=0.3

