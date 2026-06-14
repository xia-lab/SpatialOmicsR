options(shiny.maxRequestSize = 5 * 1024^3)

library(shiny)
library(SpatialOmicsMSI)
library(ggplot2)
library(viridis)
library(DT)

DEFAULT_OUTPUT_DIR <- file.path(getwd(), "outputs", "spatial_outputs")
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
      numericInput("ppm", "ppm tolerance", value = 10, min = 1, max = 50, step = 1),
      actionButton("extract", "Extract from imzML"),
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
          fluidRow(
            column(7, plotOutput("cluster_plot", height = 560)),
            column(5, plotOutput("cluster_diag", height = 560))
          ),
          downloadButton("download_clustered", "Download clustered_matrix.csv")
        ),
        tabPanel(
          "6 Sample",
          radioButtons(
            "sample_mode",
            "Sample construction",
            choices = c("Grid sub-regions" = "grid", "Matched serial regions" = "matched"),
            selected = "grid",
            inline = TRUE
          ),
          fluidRow(
            column(3, numericInput("grid_size", "grid_size", value = 5, min = 2, max = 20)),
            column(3, numericInput("min_pixels", "min_pixels", value = 30, min = 1)),
            column(6, fileInput("region_labels_upload", "matched_region_labels.csv", accept = ".csv"))
          ),
          tags$small("For matched serial regions, upload matched_region_label plus pixel_id, section_id + local_pixel_id, or section_id + x + y."),
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
  existing_outputs <- reactiveVal(load_output_bundle(DEFAULT_OUTPUT_DIR, DEFAULT_VIP_FILE))

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

  extracted <- eventReactive(input$extract, {
    req(input$msi_files, input$features)
    withProgress(message = "Extracting MSI features", value = 0, {
      incProgress(0.1, detail = "Reading MZmine feature list")
      mzmine <- read_csv_upload(input$features)
      incProgress(0.2, detail = "Preparing uploaded imzML/ibd files")
      imzml_paths <- if (identical(input$input_mode, "serial")) {
        preserve_uploaded_serial_msi_files(input$msi_files)
      } else {
        preserve_uploaded_msi_files(input$msi_files)
      }
      incProgress(0.3, detail = "Reading MSI data with Cardinal")
      result <- if (identical(input$input_mode, "serial")) {
        load_serial_msi_target_features(imzml_paths, mzmine, ppm = input$ppm)
      } else {
        load_msi_target_features(imzml_paths, mzmine, ppm = input$ppm)
      }
      incProgress(1, detail = "Done")
      result
    })
  })

  pixel_matrix <- reactive({
    uploaded <- read_csv_upload(input$pixel_matrix_upload)
    if (!is.null(uploaded)) return(uploaded)
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$pixel_matrix)) return(bundle$pixel_matrix)
    extracted()$pixel_matrix
  })

  feature_mapping <- reactive({
    uploaded <- read_csv_upload(input$feature_mapping_upload)
    if (!is.null(uploaded)) return(uploaded)
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$feature_mapping)) return(bundle$feature_mapping)
    pm <- read_csv_upload(input$pixel_matrix_upload)
    if (!is.null(pm)) return(infer_feature_mapping(pm))
    extracted()$feature_mapping
  })

  section_mapping <- reactive({
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
      tags$b("Feature statistics after normalization"),
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

  selection_preprocessed <- reactive({
    pm <- pixel_matrix()
    if (identical(input$input_mode, "serial") && "section_id" %in% names(pm)) {
      preprocess_sections(
        pm,
        do_background = input$do_background,
        background_method = input$background_method,
        background_percent = input$background_percent,
        background_percentile = input$background_percentile,
        min_nonzero_features = input$min_nonzero_features,
        do_tic = input$do_tic,
        do_log = input$do_log,
        do_scale = FALSE,
        scale_scope = "global"
      )
    } else {
      preprocess_matrix(
        pm,
        do_background = input$do_background,
        background_method = input$background_method,
        background_percent = input$background_percent,
        background_percentile = input$background_percentile,
        min_nonzero_features = input$min_nonzero_features,
        do_tic = input$do_tic,
        do_log = input$do_log,
        do_scale = FALSE
      )
    }
  })

  selected_features <- reactive({
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$selected_features)) return(bundle$selected_features)
    pm <- selection_preprocessed()$matrix
    if (identical(input$input_mode, "serial") && "section_id" %in% names(pm)) {
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
  })

  reduced_matrix <- reactive({
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$reduced_matrix)) return(bundle$reduced_matrix)
    reduce_pixel_matrix(selection_preprocessed()$matrix, selected_features())
  })

  output$selected_count <- renderText(sprintf("Selected features: %s", nrow(selected_features())))
  output$selected_table <- renderDT(datatable(selected_features(), options = list(pageLength = 10)))
  output$download_selected <- download_csv(selected_features, "selected_features.csv")
  output$download_reduced <- download_csv(reduced_matrix, "reduced_matrix.csv")

  preprocessed <- reactive({
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$preprocessed_matrix)) {
      return(list(matrix = bundle$preprocessed_matrix, distributions = NULL, background_stats = bundle$background_stats))
    }
    rm <- reduced_matrix()
    normalized <- selection_preprocessed()
    if (!isTRUE(input$do_scale)) {
      return(list(
        matrix = rm,
        distributions = normalized$distributions,
        background_stats = normalized$background_stats
      ))
    }
    if (identical(input$input_mode, "serial") && "section_id" %in% names(rm)) {
      result <- preprocess_sections(
        rm,
        do_background = FALSE,
        do_tic = FALSE,
        do_log = FALSE,
        do_scale = input$do_scale,
        scale_scope = "global"
      )
    } else {
      result <- preprocess_matrix(
        rm,
        do_background = FALSE,
        do_tic = FALSE,
        do_log = FALSE,
        do_scale = input$do_scale
      )
    }
    result$background_stats <- normalized$background_stats
    result
  })

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

  clustered <- reactive({
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$clustered_matrix)) return(list(matrix = bundle$clustered_matrix))
    cluster_pixels(
      preprocessed()$matrix,
      k = input$k,
      pca_components = if (isTRUE(input$use_pca_cluster)) input$pca_components else NULL
    )
  })
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
    diag <- cluster_diagnostics(
      preprocessed()$matrix,
      max_k = 10,
      pca_components = if (isTRUE(input$use_pca_cluster)) input$pca_components else NULL
    )
    ggplot(diag, aes(x = k, y = tot_withinss)) +
      geom_line(color = "#276FBF") +
      geom_point(size = 2, color = "#276FBF") +
      scale_x_continuous(breaks = diag$k) +
      labs(y = "Total within-cluster SS") +
      theme_minimal()
  })
  output$download_clustered <- download_csv(function() clustered()$matrix, "clustered_matrix.csv")

  sampled <- reactive({
    bundle <- existing()
    if (!is.null(bundle) && !is.null(bundle$sample_matrix) && !is.null(bundle$sample_mapping)) {
      sampled_now <- sample_subregions(clustered()$matrix, input$grid_size, input$min_pixels)
      sampled_now$sample_matrix <- bundle$sample_matrix
      sampled_now$sample_mapping <- bundle$sample_mapping
      return(sampled_now)
    }
    if (identical(input$sample_mode, "matched")) {
      labels <- read_csv_upload(input$region_labels_upload)
      validate(need(!is.null(labels), "Upload matched_region_labels.csv for matched serial regions."))
      labeled <- apply_matched_region_labels(preprocessed()$matrix, labels)
      sample_matched_regions(labeled, min_pixels = input$min_pixels)
    } else {
      sample_subregions(clustered()$matrix, input$grid_size, input$min_pixels)
    }
  })
  output$subregion_plot <- renderPlot({
    ap <- sampled()$annotated_pixels
    valid <- sampled()$sample_mapping
    if (identical(input$sample_mode, "matched")) {
      valid_pairs <- unique(valid[, c("section_id", "matched_region_label"), drop = FALSE])
      ap$valid <- interaction(ap$section_id, ap$matched_region_label, drop = TRUE) %in%
        interaction(valid_pairs$section_id, valid_pairs$matched_region_label, drop = TRUE)
      fill_values <- ifelse(ap$valid, ap$matched_region_label, "filtered")
    } else {
      ap$valid <- paste0("region", ap$cluster, "_", ap$grid_cell) %in% valid$sample_id
      fill_values <- ifelse(ap$valid, paste0("Region_", ap$cluster), "filtered")
    }
    ggplot(ap, aes(x = x, y = y, fill = fill_values)) +
      geom_tile() +
      coord_fixed() +
      scale_y_reverse() +
      labs(fill = "Sample region") +
      theme_minimal()
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
