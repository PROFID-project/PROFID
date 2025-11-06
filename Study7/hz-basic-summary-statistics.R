## Describe each variable in a given dataset. Apply statistical disclosure
## control if needed.

### Requirements ---------------------------------------------------------------
library(data.table)
library(openxlsx)
library(tools)
## -----------------------------------------------------------------------------

.HZ.describe_categorical <- function(x,                   # categorical variable
                                     x.name = "",         # variable name
                                     sdc.min.cnt = 5,     # SDC minimum count
                                     max.categories = 10) # max categories
{
  ## Describe categorical variables. Apply statistical disclosure control (SDC)
  ## if needed. Preserve the original N when masking counts. Apply SDC to the
  ## missing values count only if needed.

  ## Check arguments
  stopifnot(is.null(dim(x)))
  stopifnot(length(x) > 1)
  stopifnot(is.character(x.name))
  stopifnot(length(sdc.min.cnt) == 1)
  stopifnot(is.numeric(sdc.min.cnt) && sdc.min.cnt > 0)
  stopifnot(length(max.categories) == 1)
  stopifnot(is.numeric(max.categories) && max.categories > 0)

  ## Count catergories including missing values
  cnt <- table(x, useNA = "always")

  ## Don't do anything if too many categories
  if (length(cnt) - 1 > max.categories) {
    cnt.na <- sum(is.na(x))
    res <- data.table(
      Variable = x.name,
      N = sum(cnt),
      Category =
        c(paste0("Too many categories (>", max.categories, ")"), "Missing"),
      Count = c(NA_integer_, cnt.na),
      Proportion = c(NA_real_, cnt.na / sum(cnt)),
      Data_type = paste(class(x), collapse = " ")
    )
    res.summary <- copy(res)
    res.summary[, Summary := sprintf("%d (%.1f%s)", Count, 100 * Proportion, "%")]
    res.summary[grepl("^Too many categories", Category), Summary := "-"]
    res.summary <- res.summary[, .(Variable, Category, Summary)]
    ## Return a list
    return(list(raw = res,
                summary.non.sdc = res.summary,
                summary.sdc = res.summary))
  }

  ## Rename missing value counts
  stopifnot(!"Missing" %in% names(cnt))
  index.na <- is.na(names(cnt))
  names(cnt)[index.na] <- "Missing"

  ## Calculate proportions
  cnt.prop <- cnt / sum(cnt)

  ## Apply SDC
  index.sdc <- cnt != 0 & cnt < sdc.min.cnt
  cnt.sdc <- NULL
  sdc.not.applicable <- FALSE
  if (sum(index.sdc[!index.na]) > 0) {
    cnt.sdc <- as.character(cnt)
    names(cnt.sdc) <- names(cnt)
    ## Set small values to the minimum and calculate the total difference;
    ## ignore the missing values count.
    cnt.sdc[index.sdc & !index.na] <- paste0("<", sdc.min.cnt)
    delta <- sum(sdc.min.cnt - cnt[index.sdc & !index.na])
    ## Modify other counts one by one until all differences have been
    ## incorporated; do the missing values count last if needed.
    if (!is.null(delta)) {
      cnt.order <- order(cnt, decreasing = TRUE)
      ## Drop 0 counts
      cnt.order <- cnt.order[!cnt.order %in% which(cnt == 0)]
      ## Drop SDC counts
      cnt.order <- cnt.order[!cnt.order %in% which(index.sdc)]
      ## Put the missing values count last
      which.na <- which(index.na)
      if (which.na %in% cnt.order) {
        cnt.order <- c(cnt.order[cnt.order != which.na], which.na)
      }
      for (i in cnt.order) {
        icnt <- cnt[i] - delta
        if (icnt >= sdc.min.cnt) {
          cnt.sdc[i] <- paste0("\u2265", icnt)
          break
        } else {
          cnt.sdc[i] <- paste0("\u2265", sdc.min.cnt)
          delta <- delta - (cnt[i] - sdc.min.cnt)
        }
      }
    }
    ## Calculate SDC proportions
    cnt.sdc.new <- as.numeric(gsub("<|\u2265", "", cnt.sdc))
    cnt.prop.sdc <- cnt.sdc.new / sum(cnt.sdc.new)
    names(cnt.prop.sdc) <- names(cnt.sdc)
    ## If the new N doesn't equal to the original N, then can't apply SDC
    if (sum(cnt.sdc.new) != sum(cnt)) sdc.not.applicable <- TRUE
  }

  ## Use original counts/proportions if no SDC was required
  if (is.null(cnt.sdc)) {
    cnt.sdc <- cnt
    cnt.prop.sdc <- cnt.prop
  }

  ## Build output
  stopifnot(identical(names(cnt), names(cnt.sdc)))
  stopifnot(identical(names(cnt.prop), names(cnt.prop.sdc)))
  stopifnot(identical(names(cnt), names(cnt.prop)))
  res <- data.table(
    Variable = x.name,
    N = sum(cnt),
    Category = names(cnt),
    Count = as.integer(cnt),
    Proportion = as.numeric(cnt.prop)
  )
  res[, Summary := sprintf("%d (%.1f%s)", Count, 100 * Proportion, "%")]
  ## Check if SDC was applied successfully
  if (sdc.not.applicable) {
    res[, SDC := "Can't apply SDC, N too small"]
  } else {
    res[, SDC := sprintf("%s (%s%.1f%s)", cnt.sdc,
                         gsub("[0-9 ]+", "", cnt.sdc),
                         100 * cnt.prop.sdc, "%")]
  }
  res[, Data_type := paste(class(x), collapse = " ")]

  ## Return a list
  list(raw = res[, .(Variable, N, Category, Count, Proportion, Data_type)],
       summary.non.sdc = res[, .(Variable, Category, Summary)],
       summary.sdc = res[, .(Variable, Category, Summary = SDC)])
}

.HZ.describe_numeric <- function(x,              # numeric vector
                                 x.name = "",    # variable name
                                 sdc.min.n = 10)  # min N
{
  ## Check arguments
  ## stopifnot(is.vector(x))
  stopifnot(is.numeric(x))
  stopifnot(length(x) > 1)
  stopifnot(is.character(x.name))
  stopifnot(length(sdc.min.n) == 1)
  stopifnot(is.numeric(sdc.min.n) && sdc.min.n > 0)

  ## Basic data summaries
  raw <- data.table(
    Variable = x.name,
    N = length(x),
    Mean = mean(x, na.rm = TRUE),
    SD = sd(x, na.rm =  TRUE),
    Median = median(x, na.rm = TRUE),
    LQ = quantile(x, probs = 0.25, na.rm = TRUE),
    UQ = quantile(x, probs = 0.75, na.rm = TRUE),
    Min = min(x, na.rm = TRUE),
    Max = max(x, na.rm = TRUE),
    Missing_n = sum(is.na(x)),
    Missing_prop = sum(is.na(x)) / length(x),
    Data_type = paste(class(x), collapse = " ")
  )
  ## Summary
  summary <- data.table(
    Variable = rep(x.name, 4),
    Category = c("Mean (SD)", "Median (IQR)", "[Min, Max]", "Missing"),
    Summary = c(sprintf("%.1f (%.1f)", raw$Mean, raw$SD),
                ## sprintf("%.1f (%.1f\u2013%.1f)", raw$Median, raw$LQ, raw$UQ),
                sprintf("%.1f (%.1f--%.1f)", raw$Median, raw$LQ, raw$UQ),
                sprintf("[%.2f, %.2f]", raw$Min, raw$Max),
                sprintf("%d (%.1f%s)", raw$Missing_n, 100 * raw$Missing_prop, "%"))
  )
  ## SDC
  if (length(!is.na(x)) < sdc.min.n) {
    SDC = data.table(Variable = x.name,
                     Category = paste0("Non-missing N < ", sdc.min.n),
                     Summary = "N/A")
  } else {
    ## Drop min, max from SDC.
    SDC <- summary[Category != "[Min, Max]"]
  }
  list(raw = raw, summary.non.sdc = summary, summary.sdc = SDC)
}

.HZ.describe_date <- function(x,           # date vector
                              x.name = "") # variable name
{
  ## Check arguments
  stopifnot(is(x, "Date"))
  stopifnot(length(x) > 1)
  stopifnot(is.character(x.name))

  ## Basic Date summary
  raw <- data.table(
    Variable = x.name,
    N = length(x),
    Min = format(min(x, na.rm = TRUE), "%Y"),
    Max = format(max(x, na.rm = TRUE), "%Y"),
    Missing_n = sum(is.na(x)),
    Missing_prop = sum(is.na(x)) / length(x),
    Data_type = paste(class(x), collapse = " ")
  )
  ## Summary
  summary <- data.table(
    Variable = rep(x.name, 2),
    Category = c("Min to Max", "Missing"),
    Summary = c(sprintf("%s to %s", raw$Min, raw$Max),
                sprintf("%d (%.1f%s)", raw$Missing_n,
                        100 * raw$Missing_prop, "%"))
  )
  ## Return a list
  list(raw = raw, summary.non.sdc = summary, summary.sdc = summary)
}

HZ.eda_describe_data <- function(dt, # input data.frame or data.table object
                                 sdc.min.cnt = 5,     # SDC min counts
                                 sdc.min.n = 10,      # SDC min N
                                 max.categories = 10, # max levels to show
                                 check.year.x = TRUE) # additional date formats
{
  require(data.table)
  ## Check arguments
  stopifnot("data.frame" %in% class(dt))
  ## Describe each column depending on its data type
  .describe.column <- function(j) { # j -- column index.
    x <- dt[[j]]
    x.name <- colnames(dt)[j]
    ## Number of unique values excluding missing values
    x.unq.val.num <- length(unique(x[!is.na(x)]))
    ## Check if it's a date column
    if (is(x, "Date")) {
      return(.HZ.describe_date(x, x.name))
    }
    ## Try years format
    if (check.year.x) {
      x.as.char <- as.character(x)
      x.non.na <- x.as.char[!is.na(x.as.char)]
      if (sum(grepl("^[12][0-9]{3}$", x.non.na)) == length(x.non.na)) { # probably years
        x.as.date <- as.Date(x.as.char, format = "%Y", optional = TRUE)
        return(.HZ.describe_date(x.as.date, x.name))
      }
    }
    ## If x is not date then provide usual summaries and counts
    class.x <- paste(class(x), collapse = " ")
    if ((class.x == "numeric" & x.unq.val.num > 2) |
          (class.x == "integer" & x.unq.val.num > max.categories)) {
      return(.HZ.describe_numeric(x, x.name, sdc.min.n))
    } else {
      return(.HZ.describe_categorical(x, x.name, sdc.min.cnt, max.categories))
    }
  }
  ## Describe all columns
  res <- lapply(seq_along(dt), function(j) .describe.column(j))
  raw <- rbindlist(lapply(res, function(x) x$raw), fill = TRUE)
  summary.non.sdc <- rbindlist(lapply(res, function(x) x$summary.non.sdc))
  summary.sdc <- rbindlist(lapply(res, function(x) x$summary.sdc))
  ## Add N to summaries
  dt.N <- data.table(Variable = "N", Category = "N", Summary = nrow(dt))
  summary.non.sdc <- rbind(dt.N, summary.non.sdc)
  summary.sdc <- rbind(dt.N, summary.sdc)
  ## Add output classes
  summary.non.sdc <- structure(
    summary.non.sdc, class = c("HZ.describe_data", class(summary.non.sdc))
  )
  summary.sdc = structure(
    summary.sdc, class = c("HZ.describe_data", class(summary.sdc))
  )
  ## Return list
  list(raw = raw,
       summary.non.sdc = summary.non.sdc,
       summary.sdc = summary.sdc)
}

## Save summary element from HZ.eda_describe_data
HZ.eda_summary_write_excel <- function(d, # summary table
                                       d.name,    # dataset name
                                       file.name, # file name
                                       overwrite = TRUE)
{
  require(openxlsx)
  require(tools)
  ## Check arguments
  stopifnot("HZ.describe_data" %in% class(d))
  stopifnot(tools::file_ext(file.name) == "xlsx")
  ## Copy input data.table
  d <- copy(d)
  ## Re-format input, assumes that the first row is the N row
  d[1:.N != 1, Category := paste0("    ", Category)]
  d <- d[, {
    rbind(data.table(Category = .BY[["Variable"]], Summary = ""),
          .SD[, .(Category, Summary)])
  }, Variable]
  d <- d[, .(Category, Summary)][-1]
  setnames(d, "Summary", d.name)
  ## Excel workbook
  wb <- createWorkbook()
  addWorksheet(wb, "Data Summary")
  writeData(wb, sheet = 1, x = d)
  ## title row
  addStyle(wb, 1, createStyle(halign = "center", textDecoration = "bold"),
           cols = 1:2, rows = 1, stack = TRUE)
  ## centre the second column
  addStyle(wb, 1, createStyle(halign = "center"),
           cols = 2, rows = 1:(nrow(d) + 1), stack = TRUE)
  ## variable name rows
  vrows <- grep("^[^ ]", d[[1]]) + 1
  addStyle(wb, 1, createStyle(textDecoration = "bold"),
           cols = 1, rows = vrows, stack = TRUE)
  ## category rows
  catrows <- grep("^ +", d[[1]]) + 1
  addStyle(wb, 1, createStyle(textDecoration = "italic"),
           cols = 1, rows = catrows, stack = TRUE)
  ## borders
  borderStyle <- createStyle(border = "TopBottom", borderColour = "black")
  addStyle(wb, 1, borderStyle, cols = 1:2, rows = 1:(nrow(d) + 1),
           gridExpand = TRUE, stack = TRUE)
  ## column widths
  setColWidths(wb, 1, cols = 1, widths = max(nchar(d[[1]], type = "width"), na.rm = TRUE) + 3)
  setColWidths(wb, 1, cols = 2, widths = max(nchar(d[[2]], type = "width"), na.rm = TRUE) + 3)
  saveWorkbook(wb, file.name, overwrite)
}

HZ.table_1_sdc <- function(t1) {
  t1 <- copy(t1)
  N <- as.numeric(t1[Variable == "N", Total])
  for (i in 1:nrow(t1)) {
    i.tot <- t1[i, Total]
    if (grepl("<", i.tot)) next
    cnts <- c(t1[i, `Censored (0)`],
              t1[i, `SCD (proxy) (1)`],
              t1[i, `Death other (2)`])
    ## SDC applied entries
    is.cnts.less <- grepl("<", cnts)
    is.cnts.geq <- grepl("≥", cnts)
    is.cnts.as.is <- !grepl("<|≥", cnts)
    if (all(is.cnts.as.is)) next
    ## Parse to numeric
    tot <- as.numeric(gsub("^[<≥]*([0-9]+) .*", "\\1", i.tot))
    cnts.num <- as.numeric(gsub("^[<≥]*([0-9]+) .*", "\\1", cnts))
    if (any(is.cnts.less) & !any(is.cnts.geq)) {
      i.sum <- sum(cnts.num)
      t1[i, Total := sprintf("<%d (<%.1f%s)", i.sum, 100 * i.sum / N, "%")]
    } else if (any(is.cnts.less) & any(is.cnts.geq)) {
      i.sum <- sum(cnts.num[is.cnts.geq | is.cnts.as.is])
      t1[i, Total := sprintf("≥%d (≥%.1f%s)", i.sum, 100 * i.sum / N, "%")]
    } else if (!any(is.cnts.less) & any(is.cnts.geq)) {
      i.sum <- sum(cnts.num[is.cnts.geq | is.cnts.as.is])
      t1[i, Total := sprintf("≥%d (≥%.1f%s)", i.sum, 100 * i.sum / N, "%")]
    }
  }
  t1
}
