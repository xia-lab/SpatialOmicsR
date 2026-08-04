options(shiny.maxRequestSize = 5 * 1024^3)

library(shiny)
library(SpatialOmicsMSI)
library(ggplot2)
library(viridis)
library(DT)

DEFAULT_OUTPUT_DIR <- Sys.getenv(
  "SPATIALOMICS_OUTPUT_DIR",
  unset = file.path(getwd(), "outputs", "spatial_outputs")
)
DEFAULT_VIP_FILE <- Sys.getenv("SPATIALOMICS_VIP_FILE", unset = "")

read_optional_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

load_output_bundle <- function(output_dir = DEFAULT_OUTPUT_DIR, vip_file = DEFAULT_VIP_FILE) {
  vip <- read_optional_csv(vip_file)
  if (!is.null(vip)) {
    vip <- normalize_metaboanalyst_result(vip, source_name = basename(vip_file))
    if ("VIP" %in% names(vip)) vip <- vip[order(vip$VIP, decreasing = TRUE), , drop = FALSE]
  }

  list(
    output_dir = output_dir,
    vip_file = vip_file,
    pixel_matrix = read_optional_csv(file.path(output_dir, "pixel_feature_matrix.csv")),
    feature_mapping = read_optional_csv(file.path(output_dir, "feature_mapping.csv")),
    selected_features = read_optional_csv(file.path(output_dir, "selected_features.csv")),
    reduced_matrix = read_optional_csv(file.path(output_dir, "reduced_pixel_matrix.csv")),
    preprocessed_matrix = read_optional_csv(file.path(output_dir, "preprocessed_matrix.csv")),
    clustered_matrix = read_optional_csv(file.path(output_dir, "clustered_matrix.csv")),
    sample_matrix = read_optional_csv(file.path(output_dir, "sample_matrix.csv")),
    sample_mapping = read_optional_csv(file.path(output_dir, "sample_mapping.csv")),
    metaboanalyst_data = read_optional_csv(file.path(output_dir, "metaboanalyst_data.csv")),
    section_mapping = read_optional_csv(file.path(output_dir, "section_mapping.csv")),
    background_stats = read_optional_csv(file.path(output_dir, "background_stats.csv")),
    pairwise_permanova = read_optional_csv(file.path(output_dir, "pairwise_permanova.csv")),
    vip = vip
  )
}

read_csv_upload <- function(upload) {
  if (is.null(upload)) return(NULL)
  data <- read.csv(upload$datapath, check.names = FALSE, stringsAsFactors = FALSE)
  attr(data, "source_name") <- upload$name
  data
}

preserve_uploaded_msi_files <- function(upload) {
  req(upload)
  temp_dir <- tempfile("msi_upload_")
  dir.create(temp_dir)
  paths <- file.path(temp_dir, upload$name)
  file.copy(upload$datapath, paths, overwrite = TRUE)
  imzml <- paths[tolower(tools::file_ext(paths)) == "imzml"]
  if (length(imzml) != 1) {
    stop("Upload exactly one .imzML file with its matching .ibd file.", call. = FALSE)
  }
  imzml
}

preserve_uploaded_serial_msi_files <- function(upload) {
  req(upload)
  temp_dir <- tempfile("msi_upload_")
  dir.create(temp_dir)
  paths <- file.path(temp_dir, upload$name)
  file.copy(upload$datapath, paths, overwrite = TRUE)
  imzml <- paths[tolower(tools::file_ext(paths)) == "imzml"]
  if (length(imzml) < 1) {
    stop("Upload at least one .imzML file with matching .ibd files.", call. = FALSE)
  }
  imzml[order(basename(imzml))]
}

download_csv <- function(data_fun, filename) {
  downloadHandler(
    filename = function() filename,
    content = function(file) write.csv(data_fun(), file, row.names = FALSE)
  )
}

plot_ion_image <- function(pixel_matrix, column_name, title = NULL, discrete = FALSE) {
  validate(
    need(!is.null(pixel_matrix), "Upload or generate a pixel matrix first."),
    need(column_name %in% names(pixel_matrix), "Select a valid feature.")
  )
  ggplot(pixel_matrix, aes(x = x, y = y, fill = .data[[column_name]])) +
    geom_tile() +
    coord_fixed() +
    scale_y_reverse() +
    labs(title = title %||% column_name, x = "x", y = "y", fill = if (discrete) "Cluster" else "Intensity") +
    theme_minimal(base_size = 12)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

ui <- fluidPage(
  titlePanel("MSI Spatial Omics Pipeline"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Demo outputs"),
      checkboxInput("use_existing_outputs", "Use existing spatial_outputs folder", TRUE),
      textInput("existing_output_dir", "Output folder", value = DEFAULT_OUTPUT_DIR),
      textInput("existing_vip_file", "MetaboAnalyst VIP file", value = DEFAULT_VIP_FILE),
      actionButton("reload_existing_outputs", "Reload existing outputs", class = "btn-primary"),
      verbatimTextOutput("existing_output_status"),
      tags$hr(),
      h4("Step 1: Load"),
      radioButtons(
        "input_mode",
        "Input mode",
        choices = c("Single section" = "single", "Serial sections" = "serial"),
        selected = "single"
      ),
      fileInput("msi_files", "imzML + ibd file(s)", accept = c(".imzML", ".ibd"), multiple = TRUE),
      fileInput("features", "mzmine_features.csv", accept = ".csv"),
      tags$details(
        tags$summary("Large local files (recommended)"),
        textInput("local_imzml_path", "Local .imzML path (one per line for serial sections)", value = ""),
        textInput("local_features_path", "Local mzmine_features.csv path", value = ""),
        tags$small("Local paths bypass browser upload and avoid copying a large .ibd file twice.")
      ),
      numericInput("ppm", "ppm tolerance", value = 10, min = 1, max = 50, step = 1),
          actionButton("extract", "Extract from imzML"),
          uiOutput("extraction_progress"),
          verbatimTextOutput("extraction_status"),
      tags$hr(),
      fileInput("pixel_matrix_upload", "Or upload pixel_feature_matrix.csv", accept = ".csv"),
      fileInput("feature_mapping_upload", "feature_mapping.csv", accept = ".csv"),
      tags$hr(),
      uiOutput("feature_picker"),
      uiOutput("feature_stats")
    ),
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "1 Loaded Outputs",
          br(),
          fluidRow(
            column(4, verbatimTextOutput("loaded_summary")),
            column(8, DTOutput("pixel_preview"))
          )
        ),
        tabPanel(
          "2 Ion Images",
          fluidRow(
            column(8, plotOutput("ion_plot", height = 560)),
            column(
              4,
              uiOutput("thumbnail_page_control"),
              actionButton("render_thumbnails", "Render thumbnails")
            )
          ),
          plotOutput(
            "thumbnail_plot",
            height = 720,
            click = "thumbnail_click",
            dblclick = "thumbnail_dblclick"
          )
        ),
        tabPanel(
          "3 Select",
          fluidRow(
            column(3, sliderInput("cv_top", "CV top percent", min = 0, max = 100, value = 70)),
            column(3, numericInput("mean_min", "Mean intensity >", value = 0, min = 0)),
            column(3, sliderInput("nonzero_min", "Non-zero rate", min = 0, max = 1, value = 0.3, step = 0.05)),
            column(3, radioButtons("combine_mode", "Manual + rule", choices = c("union", "intersection"), selected = "union"))
          ),
          conditionalPanel(
            "input.input_mode == 'serial'",
            sliderInput("min_section_fraction", "Shared across sections", min = 0, max = 1, value = 1, step = 0.05)
          ),
          checkboxGroupInput("manual_features", "Manual features", choices = character()),
          actionButton("run_selection", "Confirm & run feature selection", class = "btn-primary"),
          uiOutput("selection_progress"),
          textOutput("selection_run_status"),
          textOutput("selected_count"),
          downloadButton("download_selected", "Download selected_features.csv"),
          downloadButton("download_reduced", "Download reduced_matrix.csv"),
          DTOutput("selected_table")
        ),
        tabPanel(
          "4 Preprocess",
          fluidRow(
            column(
              3,
              checkboxInput("do_background", "Background removal", TRUE),
              selectInput(
                "background_method",
                "Background method",
                choices = c(
                  "Log TIC 2-cluster mask" = "log_tic_kmeans",
                  "TIC percentile" = "tic_percentile",
                  "Median TIC fraction" = "median_fraction"
                ),
                selected = "log_tic_kmeans"
              ),
              conditionalPanel(
                "input.background_method == 'tic_percentile'",
                sliderInput("background_percentile", "Remove pixels below TIC percentile", min = 0, max = 50, value = 10, step = 1)
              ),
              conditionalPanel(
                "input.background_method == 'median_fraction'",
                numericInput("background_percent", "Background % of median TIC", 1, min = 0, max = 100)
              ),
              numericInput("min_nonzero_features", "Min non-zero features per pixel", 1, min = 0, step = 1)
            ),
            column(3, checkboxInput("do_tic", "TIC normalization", TRUE)),
            column(3, checkboxInput("do_log", "log10(x + 1)", TRUE)),
            column(3, checkboxInput("do_scale", "Auto scaling", FALSE))
          ),
          actionButton("run_preprocess", "Confirm & run preprocessing", class = "btn-primary"),
          uiOutput("preprocess_progress"),
          textOutput("preprocess_run_status"),
          downloadButton("download_preprocessed", "Download preprocessed_matrix.csv"),
          plotOutput("preprocess_hist", height = 420)
        ),
        tabPanel(
          "5 Segment",
          fluidRow(
            column(4, sliderInput("k", "k clusters", min = 2, max = 10, value = 3, step = 1)),
            column(4, checkboxInput("use_pca_cluster", "PCA before k-means", TRUE)),
            column(4, numericInput("pca_components", "PCs for k-means", value = 10, min = 1, max = 50, step = 1))
          ),
          actionButton("run_clustering", "Confirm & run K-means", class = "btn-primary"),
          uiOutput("cluster_progress"),
          textOutput("cluster_run_status"),
          fluidRow(
            column(7, plotOutput("cluster_plot", height = 560)),
            column(5, plotOutput("cluster_diag", height = 560))
          ),
          downloadButton("download_clustered", "Download clustered_matrix.csv")
        ),
        tabPanel(
          "6 Sample",
          radioButtons(
            "roi_selection_mode",
            "ROI definition",
            choices = c("Automatic" = "automatic", "Manual" = "manual"),
            selected = "automatic",
            inline = TRUE
          ),
          conditionalPanel(
            "input.roi_selection_mode == 'automatic'",
            fluidRow(
              column(3, numericInput("roi_size", "ROI size (x/y units)", value = 4000, min = 0.001)),
              column(3, selectInput("roi_shape", "ROI shape", choices = c("square", "circle"))),
              column(3, numericInput("max_rois", "Maximum ROI count", value = 10, min = 1, step = 1)),
              column(3, numericInput("roi_tau", "Stop threshold tau", value = 0.03, min = 0, max = 1, step = 0.01))
            ),
            fluidRow(
              column(3, numericInput("balance_weight", "Balance weight", value = 1, min = 0)),
              column(3, numericInput("coverage_weight", "Coverage weight", value = 1, min = 0)),
              column(3, numericInput("size_weight", "Size weight", value = 1, min = 0)),
              column(3, fileInput("prior_weights_upload", "Optional cluster, weight CSV", accept = ".csv"))
            )
          ),
          conditionalPanel(
            "input.roi_selection_mode == 'manual'",
            radioButtons(
              "manual_roi_method",
              "Manual input",
              choices = c(
                "Choose clusters" = "cluster",
                "Draw rectangle" = "draw_rectangle",
                "Draw circle" = "draw_circle",
                "Draw polygon" = "draw_polygon",
                "Geometry CSV" = "geometry_csv",
                "Polygon CSV" = "polygon_csv"
              ),
              selected = "draw_rectangle",
              inline = TRUE
            ),
            uiOutput("manual_cluster_picker"),
            conditionalPanel(
              "['draw_rectangle', 'draw_circle', 'draw_polygon'].indexOf(input.manual_roi_method) >= 0",
              fluidRow(
                column(4, uiOutput("roi_editor_section_picker")),
                column(4, textInput("roi_editor_id", "ROI label", value = "roi_manual_1")),
                column(
                  4,
                  conditionalPanel(
                    "input.manual_roi_method == 'draw_circle'",
                    numericInput("roi_circle_radius", "Circle radius (x/y units)", value = 500, min = 0.001)
                  )
                )
              ),
              tags$small("Rectangle: drag a box. Circle: click its center. Polygon: click boundary vertices in order."),
              plotOutput(
                "roi_editor_plot",
                height = 560,
                click = "roi_editor_click",
                brush = brushOpts(id = "roi_editor_brush", resetOnNew = TRUE)
              ),
              fluidRow(
                column(3, actionButton("add_manual_roi", "Add / finish ROI", class = "btn-primary")),
                column(3, actionButton("undo_roi_vertex", "Undo vertex")),
                column(3, actionButton("remove_last_roi", "Remove last ROI")),
                column(3, actionButton("clear_manual_rois", "Clear all"))
              ),
              verbatimTextOutput("roi_editor_status")
            ),
            conditionalPanel(
              "input.manual_roi_method == 'geometry_csv'",
              fileInput("manual_geometry_upload", "ROI geometry CSV", accept = ".csv")
            ),
            conditionalPanel(
              "input.manual_roi_method == 'polygon_csv'",
              fileInput("manual_polygon_upload", "Polygon vertices CSV", accept = ".csv")
            )
          ),
          fluidRow(
            column(3, numericInput("grid_size", "grid_size", value = 5, min = 2, max = 20)),
            column(3, numericInput("min_pixels", "min_pixels", value = 30, min = 1)),
            column(
              6,
              actionButton("run_roi_sampling", "Confirm ROI & run sampling", class = "btn-primary"),
              uiOutput("roi_progress"),
              textOutput("roi_run_status"),
              verbatimTextOutput("roi_status")
            )
          ),
          tags$small("Manual serial ROI files must include section_id. Every mode produces roi_id plus pixel coordinates before sub-region sampling."),
          tags$br(),
          tags$small("Plot: gray background = current sidebar ion image; ROI pixels retain their original cluster colors."),
          plotOutput("subregion_plot", height = 560),
          downloadButton("download_sample_matrix", "Download sample_matrix.csv"),
          downloadButton("download_sample_mapping", "Download sample_mapping.csv")
        ),
        tabPanel(
          "7 Export",
          textOutput("metabo_status"),
          actionButton("save_metabo_local", "Save export to local folder"),
          textOutput("metabo_local_path"),
          downloadButton("download_metabo_csv", "Download metaboanalyst_data.csv"),
          downloadButton("download_metabo_zip", "Download MetaboAnalyst ZIP"),
          DTOutput("metabo_preview")
        ),
        tabPanel(
          "8 Import",
          fileInput("metabo_result", "Upload MetaboAnalyst result", accept = ".csv"),
          textOutput("result_type"),
          DTOutput("result_preview")
        ),
        tabPanel(
          "9 Back-map",
          fluidRow(
            column(4, uiOutput("result_feature_picker"), uiOutput("score_picker")),
            column(8, plotOutput("backmap_plot", height = 560))
          ),
          fluidRow(
            column(6, plotOutput("region_boxplot", height = 420)),
            column(6, plotOutput("score_plot", height = 420))
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  last_thumbnail_dblclick <- reactiveVal(list(feature = NULL, time = as.POSIXct(0, origin = "1970-01-01")))
  manual_selected <- reactiveVal(character())
  last_metabo_export_path <- reactiveVal(NULL)
  extraction_status <- reactiveVal("Idle")
  extraction_progress_state <- reactiveVal(NULL)
  selection_run_status <- reactiveVal("Not run yet")
  preprocess_run_status <- reactiveVal("Not run yet")
  cluster_run_status <- reactiveVal("Not run yet")
  roi_run_status <- reactiveVal("Not run yet")
  selection_progress_state <- reactiveVal(NULL)
  preprocess_progress_state <- reactiveVal(NULL)
  cluster_progress_state <- reactiveVal(NULL)
  roi_progress_state <- reactiveVal(NULL)
  existing_outputs <- reactiveVal(load_output_bundle(DEFAULT_OUTPUT_DIR, DEFAULT_VIP_FILE))
  drawn_geometry <- reactiveVal(data.frame(
    section_id = character(), roi_id = character(), shape = character(),
    x_min = numeric(), x_max = numeric(), y_min = numeric(), y_max = numeric(),
    center_x = numeric(), center_y = numeric(), size = numeric(), radius = numeric(),
    definition_order = integer(), stringsAsFactors = FALSE
  ))
  drawn_polygons <- reactiveVal(data.frame(
    section_id = character(), roi_id = character(), vertex_order = integer(),
    x = numeric(), y = numeric(), definition_order = integer(), stringsAsFactors = FALSE
  ))
  draft_polygon <- reactiveVal(data.frame(x = numeric(), y = numeric()))
  circle_center <- reactiveVal(NULL)
  next_roi_order <- reactiveVal(1L)

  update_stage_progress <- function(holder, value, detail, state = "running") {
    holder(list(
      value = max(0, min(1, as.numeric(value))),
      detail = as.character(detail),
      state = state
    ))
  }

  stage_progress_ui <- function(progress) {
    if (is.null(progress)) return(NULL)
    percent <- round(100 * progress$value)
    bar_class <- switch(
      progress$state,
      success = "progress-bar progress-bar-success",
      failed = "progress-bar progress-bar-danger",
      "progress-bar progress-bar-info progress-bar-striped active"
    )
    tags$div(
      style = "margin-top:10px; margin-bottom:6px;",
      tags$div(
        class = "progress",
        style = "height:26px; margin-bottom:4px;",
        tags$div(
          class = bar_class,
          role = "progressbar",
          `aria-valuenow` = percent,
          `aria-valuemin` = 0,
          `aria-valuemax` = 100,
          style = sprintf("width:%s%%; min-width:3em; line-height:26px;", percent),
          paste0(percent, "%")
        )
      ),
      tags$small(progress$detail)
    )
  }

  output$selection_progress <- renderUI(stage_progress_ui(selection_progress_state()))
  output$extraction_progress <- renderUI(stage_progress_ui(extraction_progress_state()))
  output$preprocess_progress <- renderUI(stage_progress_ui(preprocess_progress_state()))
  output$cluster_progress <- renderUI(stage_progress_ui(cluster_progress_state()))
  output$roi_progress <- renderUI(stage_progress_ui(roi_progress_state()))

  observeEvent(list(
    input$cv_top, input$mean_min, input$nonzero_min, input$combine_mode,
    input$min_section_fraction
  ), {
    selection_run_status("Settings changed — click Confirm & run feature selection")
  }, ignoreInit = TRUE)
  observeEvent(manual_selected(), {
    selection_run_status("Feature choices changed — click Confirm & run feature selection")
  }, ignoreInit = TRUE)
  observeEvent(list(
    input$do_background, input$background_method, input$background_percent,
    input$background_percentile, input$min_nonzero_features,
    input$do_tic, input$do_log, input$do_scale
  ), {
    preprocess_run_status("Settings changed — click Confirm & run preprocessing")
  }, ignoreInit = TRUE)
  observeEvent(list(input$k, input$use_pca_cluster, input$pca_components), {
    cluster_run_status("Settings changed — click Confirm & run K-means")
  }, ignoreInit = TRUE)
  observeEvent(list(
    input$roi_selection_mode, input$manual_roi_method, input$manual_roi_clusters,
    input$roi_size, input$roi_shape, input$max_rois, input$roi_tau,
    input$balance_weight, input$coverage_weight, input$size_weight,
    input$grid_size, input$min_pixels
  ), {
    roi_run_status("ROI settings changed — click Confirm ROI & run sampling")
  }, ignoreInit = TRUE)
  observeEvent(list(drawn_geometry(), drawn_polygons()), {
    roi_run_status("Drawn ROI changed — click Confirm ROI & run sampling")
  }, ignoreInit = TRUE)

  observeEvent(input$reload_existing_outputs, {
    existing_outputs(load_output_bundle(input$existing_output_dir, input$existing_vip_file))
  }, ignoreInit = TRUE)

  existing <- reactive({
    if (!isTRUE(input$use_existing_outputs)) return(NULL)
    existing_outputs()
  })

  output$existing_output_status <- renderText({
    bundle <- existing_outputs()
    pm <- bundle$pixel_matrix
    sm <- bundle$sample_matrix
    vip <- bundle$vip
    paste(
      sprintf("pixel matrix: %s", if (is.null(pm)) "missing" else paste(nrow(pm), "pixels")),
      sprintf("sample matrix: %s", if (is.null(sm)) "missing" else paste(nrow(sm), "samples")),
      sprintf("section mapping: %s", if (is.null(bundle$section_mapping)) "missing" else paste(nrow(bundle$section_mapping), "sections")),
      sprintf("VIP result: %s", if (is.null(vip)) "missing" else paste(nrow(vip), "features")),
      sep = "\n"
    )
  })

  write_metabo_export_files <- function(output_dir) {
    ma_data <- metabo_data()
    smap <- sampled()$sample_mapping
    fmap <- feature_mapping()
    section_map <- section_mapping()
    bg_stats <- preprocessed()$background_stats
    validate(
      need(nrow(ma_data) > 0, "No MetaboAnalyst data to export."),
      need(nrow(smap) > 0, "No sample mapping to export."),
      need(nrow(fmap) > 0, "No feature mapping to export.")
    )

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    paths <- c(
      data = file.path(output_dir, "metaboanalyst_data.csv"),
      samples = file.path(output_dir, "sample_mapping.csv"),
      features = file.path(output_dir, "feature_mapping.csv")
    )
    if (!is.null(section_map) && nrow(section_map) > 0) {
      paths <- c(paths, sections = file.path(output_dir, "section_mapping.csv"))
    }
    if (!is.null(bg_stats) && nrow(bg_stats) > 0) {
      paths <- c(paths, background = file.path(output_dir, "background_stats.csv"))
    }
    write.csv(ma_data, paths[["data"]], row.names = FALSE)
    write.csv(smap, paths[["samples"]], row.names = FALSE)
    write.csv(fmap, paths[["features"]], row.names = FALSE)
    if ("sections" %in% names(paths)) write.csv(section_map, paths[["sections"]], row.names = FALSE)
    if ("background" %in% names(paths)) write.csv(bg_stats, paths[["background"]], row.names = FALSE)

    zip_path <- file.path(output_dir, "metaboanalyst_export.zip")
    zip::zipr(
      zipfile = zip_path,
      files = basename(paths),
      root = output_dir,
      recurse = FALSE,
      compression_level = 0,
      include_directories = FALSE,
      mode = "cherry-pick"
    )

    c(paths, zip = zip_path)
  }

  output$extraction_status <- renderText(extraction_status())

  extracted <- eventReactive(input$extract, {
    extraction_status("Starting extraction...")
    update_stage_progress(extraction_progress_state, 0.01, "Starting imzML extraction")
    tryCatch(
      withProgress(message = "Extracting MSI features", value = 0, {
        local_feature_path <- trimws(input$local_features_path %||% "")
        mzmine <- if (nzchar(local_feature_path)) {
          validate(need(file.exists(local_feature_path), "Local feature CSV does not exist."))
          read.csv(local_feature_path, check.names = FALSE, stringsAsFactors = FALSE)
        } else {
          req(input$features)
          read_csv_upload(input$features)
        }
        setProgress(0.01, detail = "Feature list loaded")
        update_stage_progress(extraction_progress_state, 0.01, "Feature list loaded")

        local_paths_text <- trimws(input$local_imzml_path %||% "")
        imzml_paths <- if (nzchar(local_paths_text)) {
          paths <- trimws(unlist(strsplit(local_paths_text, "[\r\n;]+")))
          paths <- paths[nzchar(paths)]
          validate(need(length(paths) > 0, "Enter at least one local .imzML path."))
          validate(need(all(file.exists(paths)), "One or more local .imzML paths do not exist."))
          validate(need(all(tolower(tools::file_ext(paths)) == "imzml"), "Local paths must point to .imzML files."))
          if (!identical(input$input_mode, "serial")) {
            validate(need(length(paths) == 1, "Single-section mode accepts exactly one .imzML path."))
          }
          paths
        } else {
          req(input$msi_files)
          extraction_status("Copying uploaded imzML/ibd files to a matched temporary directory...")
          setProgress(0.02, detail = "Preparing uploaded imzML/ibd files")
          update_stage_progress(extraction_progress_state, 0.02, "Preparing uploaded imzML/ibd files")
          if (identical(input$input_mode, "serial")) {
            preserve_uploaded_serial_msi_files(input$msi_files)
          } else {
            preserve_uploaded_msi_files(input$msi_files)
          }
        }

        missing_ibd <- vapply(imzml_paths, function(path) {
          stem <- tools::file_path_sans_ext(basename(path))
          companions <- list.files(dirname(path), full.names = TRUE)
          !any(tolower(tools::file_ext(companions)) == "ibd" &
                 tolower(tools::file_path_sans_ext(basename(companions))) == tolower(stem))
        }, logical(1))
        validate(need(!any(missing_ibd), paste0(
          "Missing matching .ibd file for: ",
          paste(basename(imzml_paths[missing_ibd]), collapse = ", ")
        )))

        progress_callback <- function(value, detail) {
          value <- max(0, min(1, as.numeric(value)))
          extraction_status(paste0(detail, " (", round(100 * value), "%)"))
          setProgress(value, detail = detail)
          update_stage_progress(extraction_progress_state, value, detail)
        }
        result <- if (identical(input$input_mode, "serial")) {
          load_serial_msi_target_features(
            imzml_paths, mzmine, ppm = input$ppm,
            progress_callback = progress_callback
          )
        } else {
          load_msi_target_features(
            imzml_paths, mzmine, ppm = input$ppm,
            progress_callback = progress_callback
          )
        }
        extraction_status(paste0("Complete: ", nrow(result$pixel_matrix), " pixels extracted."))
        update_stage_progress(extraction_progress_state, 1, "imzML extraction complete", "success")
        result
      }),
      error = function(e) {
        error_message <- conditionMessage(e)
        update_stage_progress(extraction_progress_state, 1, paste0("Failed: ", error_message), "failed")
        extraction_status(paste0("Extraction failed: ", error_message))
        showNotification(error_message, type = "error", duration = NULL)
        stop(error_message, call. = FALSE)
      }
    )
  })

  pixel_matrix <- reactive({
    uploaded <- read_csv_upload(input$pixel_matrix_upload)
    if (!is.null(uploaded)) return(uploaded)
    if (!is.null(input$extract) && input$extract > 0) {
      result <- extracted()
      if (!is.null(result$pixel_matrix)) return(result$pixel_matrix)
    }
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$pixel_matrix)) return(bundle$pixel_matrix)
    extracted()$pixel_matrix
  })

  feature_mapping <- reactive({
    uploaded <- read_csv_upload(input$feature_mapping_upload)
    if (!is.null(uploaded)) return(uploaded)
    if (!is.null(input$extract) && input$extract > 0) {
      result <- extracted()
      if (!is.null(result$feature_mapping)) return(result$feature_mapping)
    }
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$feature_mapping)) return(bundle$feature_mapping)
    pm <- read_csv_upload(input$pixel_matrix_upload)
    if (!is.null(pm)) return(infer_feature_mapping(pm))
    extracted()$feature_mapping
  })

  section_mapping <- reactive({
    if (!is.null(input$extract) && input$extract > 0) {
      result <- extracted()
      if (!is.null(result$section_mapping)) return(result$section_mapping)
    }
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$section_mapping)) return(bundle$section_mapping)
    result <- extracted()
    if (!is.null(result$section_mapping)) return(result$section_mapping)
    pm <- pixel_matrix()
    if (all(c("section_id", "section_order") %in% names(pm))) {
      unique(pm[, c("section_id", "section_order"), drop = FALSE])
    } else {
      data.frame()
    }
  })

  output$loaded_summary <- renderPrint({
    bundle <- existing_outputs()
    cat("Output folder:", bundle$output_dir, "\n")
    cat("VIP file:", bundle$vip_file, "\n\n")
    cat("pixel_feature_matrix:", if (is.null(bundle$pixel_matrix)) "missing" else paste(nrow(bundle$pixel_matrix), "rows x", ncol(bundle$pixel_matrix), "columns"), "\n")
    cat("feature_mapping:", if (is.null(bundle$feature_mapping)) "missing" else paste(nrow(bundle$feature_mapping), "features"), "\n")
    cat("section_mapping:", if (is.null(bundle$section_mapping)) "missing" else paste(nrow(bundle$section_mapping), "sections"), "\n")
    cat("selected_features:", if (is.null(bundle$selected_features)) "missing" else paste(nrow(bundle$selected_features), "features"), "\n")
    cat("sample_matrix:", if (is.null(bundle$sample_matrix)) "missing" else paste(nrow(bundle$sample_matrix), "samples"), "\n")
    cat("metaboanalyst_data:", if (is.null(bundle$metaboanalyst_data)) "missing" else paste(nrow(bundle$metaboanalyst_data), "samples"), "\n")
  })

  output$pixel_preview <- renderDT({
    datatable(head(pixel_matrix(), 20), options = list(scrollX = TRUE, pageLength = 5))
  })

  output$feature_picker <- renderUI({
    columns <- feature_columns(pixel_matrix())
    choices <- setNames(columns, sub("^mz_", "", columns))
    selectInput("feature", "Feature m/z", choices = choices)
  })

  observe({
    columns <- feature_columns(pixel_matrix())
    choices <- setNames(columns, sub("^mz_", "", columns))
    current <- intersect(manual_selected(), columns)
    manual_selected(current)
    updateCheckboxGroupInput(session, "manual_features", choices = choices, selected = current)
  })

  observeEvent(input$manual_features, {
    manual_selected(input$manual_features %||% character())
  })

  output$feature_stats <- renderUI({
    req(input$feature)
    values <- selection_preprocessed()$matrix[[input$feature]]
    row <- data.frame(
      mean = mean(values, na.rm = TRUE),
      cv = {
        mu <- mean(values, na.rm = TRUE)
        if (!is.finite(mu) || mu == 0) 0 else stats::sd(values, na.rm = TRUE) / mu
      },
      nonzero_rate = mean(values > 0, na.rm = TRUE)
    )
    tags$div(
      tags$b("Feature statistics for selection input"),
      tags$p(sprintf("Mean: %.4g", row$mean)),
      tags$p(sprintf("CV: %.4g", row$cv)),
      tags$p(sprintf("Non-zero rate: %.1f%%", 100 * row$nonzero_rate))
    )
  })

  output$ion_plot <- renderPlot({
    req(input$feature)
    plot_ion_image(pixel_matrix(), input$feature) + scale_fill_viridis(option = "viridis")
  })

  output$thumbnail_page_control <- renderUI({
    total_pages <- max(1, ceiling(length(feature_columns(pixel_matrix())) / 16))
    tagList(
      numericInput("thumb_page", "Thumbnail page", value = 1, min = 1, max = total_pages, step = 1),
      tags$small(sprintf("Total pages: %s", total_pages))
    )
  })

  output$thumbnail_plot <- renderPlot({
    req(input$render_thumbnails)
    pm <- pixel_matrix()
    fcols <- feature_columns(pm)
    total_pages <- max(1, ceiling(length(fcols) / 16))
    page <- max(1, min(input$thumb_page, total_pages))
    columns <- fcols[((page - 1) * 16 + 1):min(page * 16, length(fcols))]
    validate(need(length(columns) > 0, "No features on this page."))
    long <- tidyr::pivot_longer(pm, cols = tidyselect::all_of(columns), names_to = "feature", values_to = "intensity")
    selected_features_on_page <- intersect(manual_selected(), columns)
    selected_panels <- data.frame(
      feature = selected_features_on_page,
      xmin = rep(-Inf, length(selected_features_on_page)),
      xmax = rep(Inf, length(selected_features_on_page)),
      ymin = rep(-Inf, length(selected_features_on_page)),
      ymax = rep(Inf, length(selected_features_on_page)),
      stringsAsFactors = FALSE
    )

    plot <- ggplot(long, aes(x = x, y = y, fill = intensity)) +
      geom_raster() +
      facet_wrap(~feature, ncol = 4) +
      coord_fixed() +
      scale_y_reverse() +
      scale_fill_viridis(option = "viridis") +
      theme_void(base_size = 10) +
      theme(strip.text = element_text(size = 8))

    if (nrow(selected_panels) > 0) {
      plot <- plot +
        geom_rect(
          data = selected_panels,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill = NA,
          color = "#D62728",
          linewidth = 1.2
        )
    }

    plot
  })

  observeEvent(input$thumbnail_click, {
    clicked_feature <- input$thumbnail_click$panelvar1
    last_dblclick <- last_thumbnail_dblclick()
    seconds_since_dblclick <- as.numeric(difftime(Sys.time(), last_dblclick$time, units = "secs"))
    if (identical(clicked_feature, last_dblclick$feature) && seconds_since_dblclick < 1) {
      return()
    }
    if (!is.null(clicked_feature) && clicked_feature %in% feature_columns(pixel_matrix())) {
      updateSelectInput(session, "feature", selected = clicked_feature)
      next_selected <- union(manual_selected(), clicked_feature)
      manual_selected(next_selected)
      updateCheckboxGroupInput(
        session,
        "manual_features",
        selected = next_selected
      )
    }
  })

  observeEvent(input$thumbnail_dblclick, {
    clicked_feature <- input$thumbnail_dblclick$panelvar1
    if (!is.null(clicked_feature) && clicked_feature %in% feature_columns(pixel_matrix())) {
      last_thumbnail_dblclick(list(feature = clicked_feature, time = Sys.time()))
      updateSelectInput(session, "feature", selected = clicked_feature)
      next_selected <- setdiff(manual_selected(), clicked_feature)
      manual_selected(next_selected)
      updateCheckboxGroupInput(
        session,
        "manual_features",
        selected = next_selected
      )
    }
  })

  output$selection_run_status <- renderText(selection_run_status())

  selection_stage <- eventReactive(input$run_selection, {
    selection_run_status("Running...")
    update_stage_progress(selection_progress_state, 0.01, "Starting feature selection")
    tryCatch(
      withProgress(message = "Feature selection", value = 0, {
        pm <- pixel_matrix()
        setProgress(0.1, detail = "Loading unfiltered pixel matrix")
        update_stage_progress(selection_progress_state, 0.1, "Loading unfiltered pixel matrix")
        selection_input <- list(matrix = pm, distributions = NULL, background_stats = NULL)
        setProgress(0.55, detail = "Applying feature-selection rules")
        update_stage_progress(selection_progress_state, 0.55, "Applying feature-selection rules")
        selected <- if (identical(input$input_mode, "serial") && "section_id" %in% names(pm)) {
          select_shared_features(
            pixel_matrix = pm,
            feature_mapping = feature_mapping(),
            min_section_fraction = input$min_section_fraction,
            cv_top_percent = input$cv_top,
            mean_min = input$mean_min,
            nonzero_min = input$nonzero_min,
            manual_columns = manual_selected(),
            combine_mode = input$combine_mode
          )
        } else {
          select_features(
            pixel_matrix = pm,
            feature_mapping = feature_mapping(),
            cv_top_percent = input$cv_top,
            mean_min = input$mean_min,
            nonzero_min = input$nonzero_min,
            manual_columns = manual_selected(),
            combine_mode = input$combine_mode
          )
        }
        setProgress(0.85, detail = "Building reduced matrix")
        update_stage_progress(selection_progress_state, 0.85, "Building reduced matrix")
        reduced <- reduce_pixel_matrix(pm, selected)
        setProgress(1, detail = "Feature selection complete")
        update_stage_progress(selection_progress_state, 1, "Feature selection complete", "success")
        selection_run_status(paste0("Complete: ", nrow(selected), " features selected"))
        list(normalized = selection_input, selected = selected, reduced = reduced)
      }),
      error = function(e) {
        update_stage_progress(selection_progress_state, 1, paste0("Failed: ", conditionMessage(e)), "failed")
        selection_run_status(paste0("Failed: ", conditionMessage(e)))
        stop(conditionMessage(e), call. = FALSE)
      }
    )
  }, ignoreNULL = TRUE)

  selection_preprocessed <- reactive(selection_stage()$normalized)
  selected_features <- reactive(selection_stage()$selected)
  reduced_matrix <- reactive(selection_stage()$reduced)
  observeEvent(selection_stage(), {
    preprocess_run_status("New selection ready — click Confirm & run preprocessing")
    cluster_run_status("Upstream result changed — rerun preprocessing, then K-means")
    roi_run_status("Upstream result changed — rerun preprocessing, K-means, and ROI sampling")
  }, ignoreInit = TRUE)

  output$selected_count <- renderText(sprintf("Selected features: %s", nrow(selected_features())))
  output$selected_table <- renderDT(datatable(selected_features(), options = list(pageLength = 10)))
  output$download_selected <- download_csv(selected_features, "selected_features.csv")
  output$download_reduced <- download_csv(reduced_matrix, "reduced_matrix.csv")

  output$preprocess_run_status <- renderText(preprocess_run_status())

  preprocess_stage <- reactiveVal(NULL)

  observeEvent(input$run_preprocess, {
    preprocess_run_status("Running...")
    update_stage_progress(preprocess_progress_state, 0.01, "Starting preprocessing")
    session$onFlushed(function() {
      shiny::withReactiveDomain(session, {
      isolate({
      result <- tryCatch(
        withProgress(message = "Preprocessing", value = 0, {
          setProgress(0.15, detail = "Loading selected feature matrix")
          update_stage_progress(preprocess_progress_state, 0.15, "Loading selected feature matrix")
          rm <- reduced_matrix()
          setProgress(0.35, detail = "Applying background removal and normalization")
          update_stage_progress(preprocess_progress_state, 0.35, "Applying background removal and normalization")
          processed <- if (identical(input$input_mode, "serial") && "section_id" %in% names(rm)) {
            preprocess_sections(
              rm,
              do_background = input$do_background,
              background_method = input$background_method,
              background_percent = input$background_percent,
              background_percentile = input$background_percentile,
              min_nonzero_features = input$min_nonzero_features,
              do_tic = input$do_tic,
              do_log = input$do_log,
              do_scale = input$do_scale,
              scale_scope = "global"
            )
          } else {
            preprocess_matrix(
              rm,
              do_background = input$do_background,
              background_method = input$background_method,
              background_percent = input$background_percent,
              background_percentile = input$background_percentile,
              min_nonzero_features = input$min_nonzero_features,
              do_tic = input$do_tic,
              do_log = input$do_log,
              do_scale = input$do_scale
            )
          }
          setProgress(1, detail = "Preprocessing complete")
          update_stage_progress(preprocess_progress_state, 1, "Preprocessing complete", "success")
          preprocess_run_status(paste0("Complete: ", nrow(processed$matrix), " pixels"))
          processed
        }),
        error = function(e) {
          message <- conditionMessage(e)
          if (!nzchar(message)) message <- "Run Step 3 feature selection before preprocessing."
          update_stage_progress(preprocess_progress_state, 1, paste0("Failed: ", message), "failed")
          preprocess_run_status(paste0("Failed: ", message))
          showNotification(message, type = "error", duration = NULL)
          NULL
        }
      )
      if (!is.null(result)) preprocess_stage(result)
      })
      })
    }, once = TRUE)
  }, ignoreInit = TRUE)

  preprocessed <- reactive({
    req(preprocess_stage())
    preprocess_stage()
  })
  observeEvent(preprocess_stage(), {
    cluster_run_status("New preprocessing result ready — click Confirm & run K-means")
    roi_run_status("Upstream result changed — rerun K-means and ROI sampling")
  }, ignoreInit = TRUE)

  output$preprocess_hist <- renderPlot({
    distributions <- preprocessed()$distributions
    if (is.null(distributions)) {
      hist_data <- reduced_matrix()
      values <- unlist(hist_data[, feature_columns(hist_data), drop = FALSE], use.names = FALSE)
      values <- values[is.finite(values)]
      return(
        ggplot(data.frame(value = values), aes(x = value)) +
          geom_histogram(bins = 60, fill = "#356B6F", color = "white") +
          labs(title = "Reduced matrix intensity distribution", x = "Intensity", y = "Count") +
          theme_minimal()
      )
    }
    plot_data <- do.call(rbind, lapply(names(distributions), function(step) {
      values <- distributions[[step]]
      data.frame(step = step, value = values[is.finite(values)])
    }))
    ggplot(plot_data, aes(x = value)) +
      geom_histogram(bins = 60, fill = "#356B6F", color = "white") +
      facet_wrap(~step, scales = "free") +
      theme_minimal()
  })
  output$download_preprocessed <- download_csv(function() preprocessed()$matrix, "preprocessed_matrix.csv")

  output$cluster_run_status <- renderText(cluster_run_status())

  cluster_stage <- eventReactive(input$run_clustering, {
    cluster_run_status("Running...")
    update_stage_progress(cluster_progress_state, 0.01, "Starting K-means clustering")
    tryCatch(
      withProgress(message = "K-means clustering", value = 0, {
        pm <- preprocessed()$matrix
        pca_n <- if (isTRUE(input$use_pca_cluster)) input$pca_components else NULL
        setProgress(0.1, detail = "Preparing clustering matrix")
        update_stage_progress(cluster_progress_state, 0.1, "Preparing clustering matrix")
        setProgress(0.25, detail = "Running K-means")
        update_stage_progress(cluster_progress_state, 0.25, "Running K-means")
        clustered_result <- cluster_pixels(pm, k = input$k, pca_components = pca_n)
        setProgress(0.65, detail = "Calculating elbow diagnostics")
        update_stage_progress(cluster_progress_state, 0.65, "Calculating elbow diagnostics")
        diagnostics <- cluster_diagnostics(pm, max_k = 10, pca_components = pca_n)
        setProgress(1, detail = "Clustering complete")
        update_stage_progress(cluster_progress_state, 1, "Clustering complete", "success")
        cluster_run_status(paste0("Complete: k = ", input$k, "; ", nrow(clustered_result$matrix), " pixels"))
        list(clustered = clustered_result, diagnostics = diagnostics)
      }),
      error = function(e) {
        update_stage_progress(cluster_progress_state, 1, paste0("Failed: ", conditionMessage(e)), "failed")
        cluster_run_status(paste0("Failed: ", conditionMessage(e)))
        stop(conditionMessage(e), call. = FALSE)
      }
    )
  }, ignoreNULL = TRUE)

  clustered <- reactive(cluster_stage()$clustered)
  cluster_diagnostics_result <- reactive(cluster_stage()$diagnostics)
  observeEvent(cluster_stage(), {
    roi_run_status("New cluster result ready — define/confirm ROI and run sampling")
  }, ignoreInit = TRUE)
  output$cluster_plot <- renderPlot({
    cm <- clustered()$matrix
    ggplot(cm, aes(x = x, y = y, fill = factor(cluster))) +
      geom_tile() +
      coord_fixed() +
      scale_y_reverse() +
      scale_fill_brewer(palette = "Set2") +
      labs(fill = "Cluster") +
      theme_minimal()
  })
  output$cluster_diag <- renderPlot({
    diag <- cluster_diagnostics_result()
    ggplot(diag, aes(x = k, y = tot_withinss)) +
      geom_line(color = "#276FBF") +
      geom_point(size = 2, color = "#276FBF") +
      scale_x_continuous(breaks = diag$k) +
      labs(y = "Total within-cluster SS") +
      theme_minimal()
  })
  output$download_clustered <- download_csv(function() clustered()$matrix, "clustered_matrix.csv")

  output$roi_editor_section_picker <- renderUI({
    cm <- clustered()$matrix
    choices <- if ("section_id" %in% names(cm)) unique(as.character(cm$section_id)) else ".__single__"
    selectInput("roi_editor_section", "Section", choices = choices, selected = choices[1])
  })

  active_roi_section <- reactive({
    value <- input$roi_editor_section
    if (is.null(value) || !nzchar(value)) ".__single__" else as.character(value)
  })

  observeEvent(list(input$roi_editor_section, input$manual_roi_method), {
    draft_polygon(data.frame(x = numeric(), y = numeric()))
    circle_center(NULL)
    session$resetBrush("roi_editor_brush")
  }, ignoreInit = TRUE)

  observeEvent(input$roi_editor_click, {
    click <- input$roi_editor_click
    req(click$x, click$y)
    if (identical(input$manual_roi_method, "draw_polygon")) {
      draft <- draft_polygon()
      draft_polygon(rbind(draft, data.frame(x = click$x, y = click$y)))
    } else if (identical(input$manual_roi_method, "draw_circle")) {
      circle_center(c(x = click$x, y = click$y))
    }
  })

  roi_label_available <- function(label, section) {
    geometry_key <- drawn_geometry()
    polygon_key <- drawn_polygons()
    any(geometry_key$section_id == section & geometry_key$roi_id == label) ||
      any(polygon_key$section_id == section & polygon_key$roi_id == label)
  }

  observeEvent(input$add_manual_roi, {
    method <- input$manual_roi_method
    if (!method %in% c("draw_rectangle", "draw_circle", "draw_polygon")) return()
    section <- active_roi_section()
    label <- trimws(input$roi_editor_id %||% "")
    if (!nzchar(label)) label <- paste0("roi_manual_", next_roi_order())
    if (roi_label_available(label, section)) {
      showNotification("This ROI label already exists in the active section.", type = "error")
      return()
    }
    order_id <- next_roi_order()

    if (identical(method, "draw_rectangle")) {
      brush <- input$roi_editor_brush
      if (is.null(brush)) {
        showNotification("Drag a rectangle on the map first.", type = "error")
        return()
      }
      row <- data.frame(
        section_id = section, roi_id = label, shape = "rectangle",
        x_min = min(brush$xmin, brush$xmax), x_max = max(brush$xmin, brush$xmax),
        y_min = min(brush$ymin, brush$ymax), y_max = max(brush$ymin, brush$ymax),
        center_x = NA_real_, center_y = NA_real_, size = NA_real_, radius = NA_real_,
        definition_order = order_id, stringsAsFactors = FALSE
      )
      drawn_geometry(rbind(drawn_geometry(), row))
      session$resetBrush("roi_editor_brush")
    } else if (identical(method, "draw_circle")) {
      center <- circle_center()
      radius <- input$roi_circle_radius
      if (is.null(center) || !is.finite(radius) || radius <= 0) {
        showNotification("Click a circle center and provide a positive radius.", type = "error")
        return()
      }
      row <- data.frame(
        section_id = section, roi_id = label, shape = "circle",
        x_min = NA_real_, x_max = NA_real_, y_min = NA_real_, y_max = NA_real_,
        center_x = unname(center[["x"]]), center_y = unname(center[["y"]]),
        size = NA_real_, radius = radius, definition_order = order_id,
        stringsAsFactors = FALSE
      )
      drawn_geometry(rbind(drawn_geometry(), row))
      circle_center(NULL)
    } else {
      vertices <- draft_polygon()
      if (nrow(vertices) < 3) {
        showNotification("A polygon needs at least three vertices.", type = "error")
        return()
      }
      vertices$section_id <- section
      vertices$roi_id <- label
      vertices$vertex_order <- seq_len(nrow(vertices))
      vertices$definition_order <- order_id
      vertices <- vertices[, c("section_id", "roi_id", "vertex_order", "x", "y", "definition_order")]
      drawn_polygons(rbind(drawn_polygons(), vertices))
      draft_polygon(data.frame(x = numeric(), y = numeric()))
    }
    next_roi_order(order_id + 1L)
    updateTextInput(session, "roi_editor_id", value = paste0("roi_manual_", order_id + 1L))
  })

  observeEvent(input$undo_roi_vertex, {
    draft <- draft_polygon()
    if (nrow(draft)) draft_polygon(draft[-nrow(draft), , drop = FALSE])
  })

  observeEvent(input$remove_last_roi, {
    orders <- c(drawn_geometry()$definition_order, drawn_polygons()$definition_order)
    if (!length(orders)) return()
    last <- max(orders)
    drawn_geometry(drawn_geometry()[drawn_geometry()$definition_order != last, , drop = FALSE])
    drawn_polygons(drawn_polygons()[drawn_polygons()$definition_order != last, , drop = FALSE])
  })

  observeEvent(input$clear_manual_rois, {
    drawn_geometry(drawn_geometry()[FALSE, , drop = FALSE])
    drawn_polygons(drawn_polygons()[FALSE, , drop = FALSE])
    draft_polygon(data.frame(x = numeric(), y = numeric()))
    circle_center(NULL)
    next_roi_order(1L)
    updateTextInput(session, "roi_editor_id", value = "roi_manual_1")
    session$resetBrush("roi_editor_brush")
  })

  circle_outlines <- function(geometry) {
    circles <- geometry[geometry$shape == "circle", , drop = FALSE]
    if (!nrow(circles)) return(data.frame())
    do.call(rbind, lapply(seq_len(nrow(circles)), function(i) {
      angle <- seq(0, 2 * pi, length.out = 101)
      data.frame(
        roi_id = circles$roi_id[i],
        x = circles$center_x[i] + circles$radius[i] * cos(angle),
        y = circles$center_y[i] + circles$radius[i] * sin(angle)
      )
    }))
  }

  output$roi_editor_plot <- renderPlot({
    cm <- clustered()$matrix
    section <- active_roi_section()
    if ("section_id" %in% names(cm)) cm <- cm[as.character(cm$section_id) == section, , drop = FALSE]
    geometry <- drawn_geometry()
    geometry <- geometry[geometry$section_id == section, , drop = FALSE]
    polygons <- drawn_polygons()
    polygons <- polygons[polygons$section_id == section, , drop = FALSE]
    draft <- draft_polygon()
    center <- circle_center()

    plot <- ggplot(cm, aes(x = x, y = y, fill = factor(cluster))) +
      geom_tile() +
      coord_fixed() +
      scale_y_reverse() +
      scale_fill_brewer(palette = "Set2") +
      labs(fill = "Cluster", title = paste("ROI editor —", section)) +
      theme_minimal()

    rectangles <- geometry[geometry$shape == "rectangle", , drop = FALSE]
    if (nrow(rectangles)) {
      plot <- plot + geom_rect(
        data = rectangles,
        aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, color = roi_id),
        inherit.aes = FALSE,
        fill = NA,
        linewidth = 1.2
      )
    }
    circles <- circle_outlines(geometry)
    if (nrow(circles)) {
      plot <- plot + geom_path(
        data = circles,
        aes(x = x, y = y, color = roi_id, group = roi_id),
        inherit.aes = FALSE,
        linewidth = 1.2
      )
    }
    if (nrow(polygons)) {
      plot <- plot + geom_polygon(
        data = polygons,
        aes(x = x, y = y, color = roi_id, group = interaction(section_id, roi_id)),
        inherit.aes = FALSE,
        fill = NA,
        linewidth = 1.2
      )
    }
    if (nrow(draft)) {
      plot <- plot +
        geom_path(data = draft, aes(x = x, y = y), inherit.aes = FALSE, color = "#D62728", linewidth = 1) +
        geom_point(data = draft, aes(x = x, y = y), inherit.aes = FALSE, color = "#D62728", size = 2)
    }
    if (!is.null(center)) {
      preview_angle <- seq(0, 2 * pi, length.out = 101)
      preview_circle <- data.frame(
        x = center[["x"]] + input$roi_circle_radius * cos(preview_angle),
        y = center[["y"]] + input$roi_circle_radius * sin(preview_angle)
      )
      plot <- plot +
        geom_path(
          data = preview_circle, aes(x = x, y = y), inherit.aes = FALSE,
          color = "#D62728", linewidth = 1, linetype = "dashed"
        ) +
        geom_point(
        data = data.frame(x = center[["x"]], y = center[["y"]]),
        aes(x = x, y = y), inherit.aes = FALSE, color = "#D62728", size = 3, shape = 4, stroke = 1.5
      )
    }
    plot
  })

  output$roi_editor_status <- renderPrint({
    section <- active_roi_section()
    geometry <- drawn_geometry()
    polygons <- drawn_polygons()
    geometry_count <- sum(geometry$section_id == section)
    polygon_count <- length(unique(polygons$definition_order[polygons$section_id == section]))
    cat("Saved ROIs in this section:", geometry_count + polygon_count, "\n")
    cat("Total saved ROIs:", length(unique(c(geometry$definition_order, polygons$definition_order))), "\n")
    if (identical(input$manual_roi_method, "draw_polygon")) {
      cat("Draft polygon vertices:", nrow(draft_polygon()), "\n")
    }
    if (identical(input$manual_roi_method, "draw_circle") && !is.null(circle_center())) {
      cat("Circle center selected; click Add / finish ROI.\n")
    }
  })

  output$manual_cluster_picker <- renderUI({
    if (!identical(input$roi_selection_mode, "manual") || !identical(input$manual_roi_method, "cluster")) return(NULL)
    clusters <- sort(unique(clustered()$matrix$cluster))
    checkboxGroupInput("manual_roi_clusters", "Histology clusters", choices = clusters, selected = clusters)
  })

  output$roi_run_status <- renderText(roi_run_status())

  roi_sampling_stage <- eventReactive(input$run_roi_sampling, {
    roi_run_status("Running...")
    update_stage_progress(roi_progress_state, 0.01, "Starting ROI selection")
    tryCatch(
      withProgress(message = "ROI selection and sampling", value = 0, {
        cm <- clustered()$matrix
        section_col <- if (identical(input$input_mode, "serial") && "section_id" %in% names(cm)) "section_id" else NULL
        setProgress(0.1, detail = "Validating ROI definition")
        update_stage_progress(roi_progress_state, 0.1, "Validating ROI definition")
        selected <- if (identical(input$roi_selection_mode, "automatic")) {
          validate(need(
            is.null(section_col) || length(unique(cm[[section_col]])) == 1,
            "Automatic serial ROI selection requires coregistration and is reserved for a future version. Use Manual mode for serial sections."
          ))
          prior <- read_csv_upload(input$prior_weights_upload)
          prior_weights <- NULL
          if (!is.null(prior)) {
            validate(need(all(c("cluster", "weight") %in% names(prior)), "Prior CSV needs cluster and weight columns."))
            prior_weights <- stats::setNames(prior$weight, as.character(prior$cluster))
          }
          setProgress(0.3, detail = "Optimizing automatic ROIs")
          update_stage_progress(roi_progress_state, 0.3, "Optimizing automatic ROIs")
          select_rois(
            cm,
            selection_mode = "automatic",
            roi_size = input$roi_size,
            shape = input$roi_shape,
            cluster_column = "cluster",
            valid_column = NULL,
            section_column = section_col,
            prior_weights = prior_weights,
            score_weights = c(
              balance = input$balance_weight,
              coverage = input$coverage_weight,
              size = input$size_weight
            ),
            max_rois = input$max_rois,
            improvement_threshold = input$roi_tau
          )
        } else {
          manual_method <- input$manual_roi_method
          interactive_method <- manual_method %in% c("draw_rectangle", "draw_circle", "draw_polygon")
          geometry <- if (interactive_method) drawn_geometry() else read_csv_upload(input$manual_geometry_upload)
          polygons <- if (interactive_method) drawn_polygons() else read_csv_upload(input$manual_polygon_upload)
          backend_method <- if (identical(manual_method, "cluster")) {
            "cluster"
          } else if (interactive_method) {
            "combined"
          } else if (identical(manual_method, "geometry_csv")) {
            "geometry"
          } else {
            "polygon"
          }
          if (identical(backend_method, "cluster")) {
            validate(need(length(input$manual_roi_clusters) > 0, "Choose at least one cluster."))
          } else if (identical(backend_method, "geometry")) {
            validate(need(!is.null(geometry), "Upload an ROI geometry CSV."))
          } else if (identical(backend_method, "polygon")) {
            validate(need(!is.null(polygons), "Upload polygon vertices CSV."))
          } else {
            validate(need(nrow(geometry) + nrow(polygons) > 0, "Draw and add at least one ROI."))
          }
          setProgress(0.3, detail = "Mapping drawn ROI to pixels")
          update_stage_progress(roi_progress_state, 0.3, "Mapping drawn ROI to pixels")
          select_rois(
            cm,
            selection_mode = "manual",
            manual_method = backend_method,
            cluster_column = "cluster",
            selected_clusters = input$manual_roi_clusters,
            roi_table = geometry,
            polygon_vertices = polygons,
            section_column = section_col,
            overlap = "first"
          )
        }
        setProgress(0.75, detail = "Sampling independent sub-regions")
        update_stage_progress(roi_progress_state, 0.75, "Sampling independent sub-regions")
        sampled_result <- sample_subregions(
          selected$annotated_pixels,
          grid_size = input$grid_size,
          min_pixels = input$min_pixels,
          roi_column = "roi_id",
          grid_scope = "roi"
        )
        setProgress(1, detail = "ROI sampling complete")
        update_stage_progress(roi_progress_state, 1, "ROI sampling complete", "success")
        roi_run_status(paste0(
          "Complete: ", length(unique(selected$selected_pixels$roi_id)),
          " ROI(s), ", nrow(sampled_result$sample_mapping), " sub-regions"
        ))
        list(selection = selected, sampled = sampled_result)
      }),
      error = function(e) {
        update_stage_progress(roi_progress_state, 1, paste0("Failed: ", conditionMessage(e)), "failed")
        roi_run_status(paste0("Failed: ", conditionMessage(e)))
        stop(conditionMessage(e), call. = FALSE)
      }
    )
  }, ignoreNULL = TRUE)

  roi_selection <- reactive(roi_sampling_stage()$selection)
  sampled <- reactive(roi_sampling_stage()$sampled)

  output$roi_status <- renderPrint({
    selected <- roi_selection()
    cat("Selected ROI pixels:", nrow(selected$selected_pixels), "\n")
    cat("ROI count:", length(unique(selected$selected_pixels$roi_id)), "\n")
    if (!is.null(selected$optimization)) {
      cat("ROI score:", signif(selected$optimization$score[["roi_score"]], 4), "\n")
      print(selected$optimization$history, row.names = FALSE)
    }
  })

  output$subregion_plot <- renderPlot({
    ap <- sampled()$annotated_pixels
    valid <- sampled()$sample_mapping
    if (nrow(valid) > 0 && all(c("roi_id", "grid_cell") %in% names(valid))) {
      valid_cells <- interaction(valid$roi_id, valid$grid_cell, drop = TRUE)
      ap$valid <- interaction(ap$roi_id, ap$grid_cell, drop = TRUE) %in% valid_cells
    } else {
      ap$valid <- FALSE
    }

    base <- pixel_matrix()
    reference_feature <- input$feature
    if (is.null(reference_feature) || !reference_feature %in% names(base)) {
      reference_feature <- feature_columns(base)[1]
    }
    validate(need(!is.na(reference_feature) && length(reference_feature) == 1,
                  "Choose an ion feature in the sidebar for the gray reference image."))
    if ("section_id" %in% names(ap) && "section_id" %in% names(base)) {
      base <- base[as.character(base$section_id) %in% unique(as.character(ap$section_id)), , drop = FALSE]
    }
    intensity <- as.numeric(base[[reference_feature]])
    limits <- stats::quantile(intensity[is.finite(intensity)], c(0.02, 0.98), na.rm = TRUE)
    if (!all(is.finite(limits)) || limits[2] <= limits[1]) limits <- range(intensity, finite = TRUE)
    if (!all(is.finite(limits)) || limits[2] <= limits[1]) {
      base$ion_alpha <- 0.25
    } else {
      scaled <- (intensity - limits[1]) / diff(limits)
      base$ion_alpha <- 0.08 + 0.57 * pmax(0, pmin(1, scaled))
    }
    roi_centers <- if ("section_id" %in% names(ap)) {
      stats::aggregate(cbind(x, y) ~ section_id + roi_id, data = ap, FUN = stats::median)
    } else {
      stats::aggregate(cbind(x, y) ~ roi_id, data = ap, FUN = stats::median)
    }

    plot <- ggplot() +
      geom_tile(
        data = base,
        aes(x = x, y = y, alpha = ion_alpha),
        fill = "grey30"
      ) +
      scale_alpha_identity() +
      geom_tile(
        data = ap,
        aes(x = x, y = y, fill = factor(cluster)),
        alpha = 0.88
      ) +
      geom_label(
        data = roi_centers,
        aes(x = x, y = y, label = roi_id),
        fill = "white", alpha = 0.75, size = 3, label.size = 0
      ) +
      coord_fixed() +
      scale_y_reverse() +
      scale_fill_brewer(palette = "Set2", name = "Histology cluster") +
      labs(
        title = paste0("ROI clusters over gray ion reference: ", sub("^mz_", "", reference_feature)),
        subtitle = paste(sum(ap$valid), "pixels in retained sub-regions;", sum(!ap$valid), "pixels below min_pixels"),
        x = "x", y = "y"
      ) +
      theme_minimal()
    if ("section_id" %in% names(ap) && length(unique(ap$section_id)) > 1) {
      plot <- plot + facet_wrap(~section_id, scales = "free")
    }
    plot
  })
  output$download_sample_matrix <- download_csv(function() sampled()$sample_matrix, "sample_matrix.csv")
  output$download_sample_mapping <- download_csv(function() sampled()$sample_mapping, "sample_mapping.csv")

  metabo_data <- reactive({
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$metaboanalyst_data)) return(bundle$metaboanalyst_data)
    sm <- sampled()$sample_matrix
    validate(
      need(nrow(sm) > 0, "No sub-region samples available. Lower Step 6 min_pixels or grid_size, or check Step 4 background removal."),
      need(length(feature_columns(sm)) > 0, "No selected feature columns available for MetaboAnalyst export.")
    )
    if ("matched_region_label" %in% names(sm)) {
      make_metaboanalyst_data(sm, group_column = "matched_region_label")
    } else {
      make_metaboanalyst_data(sm)
    }
  })
  output$metabo_status <- renderText({
    sm <- sampled()$sample_matrix
    if (nrow(sm) == 0) {
      return("No MetaboAnalyst samples yet. Try lowering Step 6 min_pixels, using a smaller grid_size, or relaxing Step 4 background removal.")
    }
    sprintf(
      "Ready to export: %s sub-region samples x %s features.",
      nrow(sm),
      length(feature_columns(sm))
    )
  })
  output$metabo_preview <- renderDT(datatable(head(metabo_data(), 20), options = list(scrollX = TRUE)))
  output$download_metabo_csv <- download_csv(metabo_data, "metaboanalyst_data.csv")
  observeEvent(input$save_metabo_local, {
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- file.path(getwd(), "exports", paste0("metaboanalyst_export_", stamp))
    write_metabo_export_files(output_dir)
    last_metabo_export_path(output_dir)
  })
  output$metabo_local_path <- renderText({
    path <- last_metabo_export_path()
    if (is.null(path)) return("")
    paste("Saved export files to:", path)
  })
  output$download_metabo_zip <- downloadHandler(
    filename = function() "metaboanalyst_export.zip",
    content = function(file) {
      temp_dir <- tempfile("metaboanalyst_export_")
      paths <- write_metabo_export_files(temp_dir)
      file.copy(paths[["zip"]], file, overwrite = TRUE)
    }
  )

  metabo_result <- reactive({
    result <- read_csv_upload(input$metabo_result)
    if (is.null(result)) {
      bundle <- existing()
      if (!is.null(bundle) && !is.null(bundle$vip)) return(bundle$vip)
    }
    normalize_metaboanalyst_result(result, source_name = attr(result, "source_name"))
  })
  result_type <- reactive({
    req(metabo_result())
    detect_metaboanalyst_result(metabo_result())
  })
  output$result_type <- renderText(paste("Detected result type:", result_type()))
  output$result_preview <- renderDT(datatable(metabo_result(), options = list(pageLength = 10, scrollX = TRUE)))

  output$result_feature_picker <- renderUI({
    req(metabo_result())
    if (!result_type() %in% c("vip", "differential")) return(NULL)
    selectInput("result_feature", "Feature result", choices = metabo_result()$Feature)
  })

  output$score_picker <- renderUI({
    req(metabo_result())
    if (!result_type() %in% c("pca_scores", "plsda_scores")) return(NULL)
    score_columns <- setdiff(names(metabo_result()), "Sample")
    selectInput("score_column", "Score", choices = score_columns)
  })

  output$backmap_plot <- renderPlot({
    req(metabo_result())
    if (result_type() %in% c("vip", "differential")) {
      req(input$result_feature)
      pm <- pixel_matrix()
      title <- input$result_feature
      result_row <- metabo_result()[metabo_result()$Feature == input$result_feature, , drop = FALSE]
      if ("VIP" %in% names(result_row)) title <- paste0(title, " | VIP=", signif(result_row$VIP[1], 3))
      if ("p.value" %in% names(result_row)) title <- paste0(title, " | p=", signif(result_row$p.value[1], 3))
      plot_ion_image(pm, input$result_feature, title = title) + scale_fill_viridis(option = "viridis")
    } else {
      req(input$score_column)
      scores <- backmap_sample_scores(metabo_result(), sampled()$sample_mapping, input$score_column)
      pm <- merge(pixel_matrix()[, c("pixel_id", "x", "y")], scores, by = "pixel_id", all.x = FALSE)
      score_midpoint <- stats::median(pm$score, na.rm = TRUE)
      if (!is.finite(score_midpoint)) score_midpoint <- 0
      ggplot(pm, aes(x = x, y = y, fill = score)) +
        geom_tile() +
        coord_fixed() +
        scale_y_reverse() +
        scale_fill_gradient2(
          low = "#2B6CB0",
          mid = "white",
          high = "#C53030",
          midpoint = score_midpoint
        ) +
        theme_minimal()
    }
  })

  output$region_boxplot <- renderPlot({
    req(metabo_result())
    validate(need(result_type() %in% c("vip", "differential"), "Feature-level results show boxplots here."))
    req(input$result_feature)
    sm <- sampled()$sample_matrix
    group_values <- if ("matched_region_label" %in% names(sm)) {
      sm$matched_region_label
    } else if ("roi_id" %in% names(sm)) {
      sm$roi_id
    } else {
      paste0("Region_", sm$cluster)
    }
    ggplot(sm, aes(x = group_values, y = .data[[input$result_feature]], color = group_values)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, height = 0) +
      theme_minimal() +
      labs(x = "Region", y = input$result_feature, color = "Region")
  })

  output$score_plot <- renderPlot({
    req(metabo_result())
    validate(need(result_type() %in% c("pca_scores", "plsda_scores"), "Sample-level score plots show here."))
    result <- metabo_result()
    xcol <- if ("PC1" %in% names(result)) "PC1" else "Comp1"
    ycol <- if ("PC2" %in% names(result)) "PC2" else "Comp2"
    result <- merge(result, metabo_data()[, c("Sample", "Group")], by = "Sample", all.x = TRUE)
    ggplot(result, aes(x = .data[[xcol]], y = .data[[ycol]], color = Group)) +
      geom_point(size = 3) +
      theme_minimal()
  })

  output$permanova_image <- renderImage({
    bundle <- existing_outputs()
    file <- file.path(bundle$output_dir, "plots_for_slides", "step11_pairwise_permanova_table.png")
    validate(need(file.exists(file), "PERMANOVA figure is missing."))
    list(src = file, contentType = "image/png", width = "100%")
  }, deleteFile = FALSE)

  output$pca_image <- renderImage({
    bundle <- existing_outputs()
    file <- file.path(bundle$output_dir, "plots_for_slides", "step11_pca_score_plot.png")
    validate(need(file.exists(file), "PCA figure is missing."))
    list(src = file, contentType = "image/png", width = "100%")
  }, deleteFile = FALSE)

  output$permanova_table <- renderDT({
    bundle <- existing_outputs()
    validate(need(!is.null(bundle$pairwise_permanova), "pairwise_permanova.csv is missing."))
    datatable(bundle$pairwise_permanova, options = list(pageLength = 10))
  })
}

shinyApp(ui, server)
