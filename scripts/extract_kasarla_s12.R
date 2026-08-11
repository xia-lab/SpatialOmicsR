args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: extract_kasarla_s12.R input.xlsx output.csv", call. = FALSE)

archive <- normalizePath(args[1], mustWork = TRUE)
destination <- args[2]
temporary <- tempfile("kasarla-xlsx-")
dir.create(temporary)
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
utils::unzip(archive, files = c("xl/sharedStrings.xml", "xl/worksheets/sheet13.xml"), exdir = temporary)

read_xml_text <- function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "")
decode_xml <- function(x) {
  x <- gsub("<[^>]+>", "", x)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x <- gsub("&lt;", "<", x, fixed = TRUE)
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  x <- gsub("&quot;", '"', x, fixed = TRUE)
  x <- gsub("&#39;", "'", x, fixed = TRUE)
  x
}
extract_all <- function(text, pattern) {
  match <- gregexpr(pattern, text, perl = TRUE)
  if (match[[1]][1] < 0L) return(character())
  regmatches(text, match)[[1]]
}

shared_xml <- read_xml_text(file.path(temporary, "xl/sharedStrings.xml"))
shared <- vapply(extract_all(shared_xml, "<si(?: [^>]*)?>.*?</si>"), decode_xml, character(1))
sheet_xml <- read_xml_text(file.path(temporary, "xl/worksheets/sheet13.xml"))
rows <- extract_all(sheet_xml, "<row(?: [^>]*)?>.*?</row>")

column_number <- function(reference) {
  letters <- strsplit(sub("[0-9]+$", "", reference), "", fixed = TRUE)[[1]]
  sum((match(letters, LETTERS)) * 26 ^ rev(seq_along(letters) - 1L))
}
parsed <- lapply(rows, function(row) {
  cells <- extract_all(row, "<c(?: [^>]*)?>.*?</c>")
  values <- list()
  for (cell in cells) {
    reference <- sub('^.*? r="([A-Z]+[0-9]+)".*$', "\\1", cell, perl = TRUE)
    raw <- sub("^.*?<v>(.*?)</v>.*$", "\\1", cell, perl = TRUE)
    if (identical(raw, cell)) raw <- ""
    if (grepl(' t="s"', cell, fixed = TRUE) && nzchar(raw)) raw <- shared[as.integer(raw) + 1L]
    values[[as.character(column_number(reference))]] <- raw
  }
  values
})
width <- max(vapply(parsed, function(x) if (length(x)) max(as.integer(names(x))) else 0L, integer(1)))
matrix <- matrix("", nrow = length(parsed), ncol = width)
for (i in seq_along(parsed)) if (length(parsed[[i]])) matrix[i, as.integer(names(parsed[[i]]))] <- unlist(parsed[[i]], use.names = FALSE)

header_row <- which(apply(matrix, 1, function(x) sum(nzchar(x))) == max(apply(matrix, 1, function(x) sum(nzchar(x)))))[1]
header <- trimws(matrix[header_row, ])
keep_columns <- nzchar(header)
output <- as.data.frame(matrix[(header_row + 1L):nrow(matrix), keep_columns, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)
names(output) <- make.unique(header[keep_columns])
output <- output[apply(output, 1, function(x) any(nzchar(trimws(x)))), , drop = FALSE]
dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
write.csv(output, destination, row.names = FALSE, na = "")
cat(sprintf("Wrote %d rows and %d columns to %s\n", nrow(output), ncol(output), destination))
