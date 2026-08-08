# =============================================================================
# MissLearn
# Longitudinal Methods Platform — Institute of Psychiatry, Psychology and
# Neuroscience, King's College London
#
# A teaching application for the four-step framework for identifying and
# handling missing data in longitudinal cohorts.
#
# All data is simulated. No patient data is used.
# =============================================================================

library(shiny)
library(mice)

# ---------------------------------------------------------------- palette ---
NAVY   <- "#0a2d50"
RED    <- "#e12726"
GREY   <- "#4d545e"
RULE   <- "#d9dde3"
WASH   <- "#f2f4f7"
MIDBLU <- "#7f93ab"

# ============================================================== simulation ===

simulate_cohort <- function(n = 300, n_visits = 5, effect = 0.6,
                            dropout = 0.30, mechanism = "mar", seed = 1) {

  set.seed(seed)

  # ---- complete data: random intercept and slope, two groups ----
  grp <- rep(0:1, length.out = n)
  b0  <- rnorm(n, 0, 0.8)
  b1  <- rnorm(n, 0, 0.25)

  y <- matrix(NA_real_, nrow = n, ncol = n_visits)
  for (j in seq_len(n_visits)) {
    t <- (j - 1) / (n_visits - 1)
    y[, j] <- b0 + b1 * t * (n_visits - 1) + effect * t * grp + rnorm(n, 0, 0.45)
  }
  colnames(y) <- paste0("y", seq_len(n_visits))
  full <- as.data.frame(y)
  full$grp <- grp

  # ---- monotone dropout ----
  # per-visit hazard tuned so that the cumulative dropout by the final visit
  # is approximately the requested proportion
  if (dropout <= 0) {
    obs <- matrix(TRUE, n, n_visits)
  } else {
    haz <- 1 - (1 - dropout)^(1 / (n_visits - 1))
    obs <- matrix(TRUE, n, n_visits)
    for (j in 2:n_visits) {
      still_in <- obs[, j - 1]
      lp <- qlogis(min(max(haz, 1e-4), 0.99))

      driver <- switch(
        mechanism,
        mcar = rep(0, n),
        mar  = 0.90 * y[, j - 1],          # depends on the previous OBSERVED value
        mnar = 0.90 * y[, j]               # depends on the CURRENT, unobserved value
      )
      driver <- driver - mean(driver)

      p_drop <- plogis(lp + driver)
      drops  <- still_in & (runif(n) < p_drop)
      if (any(drops)) obs[drops, j:n_visits] <- FALSE
    }
  }

  observed <- full
  observed[, seq_len(n_visits)][!obs] <- NA_real_

  list(full = full, observed = observed, obs = obs,
       n_visits = n_visits, n = n)
}

# ================================================================ analysis ===

# estimate: mean change from first to last visit, contrast between groups
estimate_from <- function(d, n_visits) {
  first <- d[[paste0("y", 1)]]
  last  <- d[[paste0("y", n_visits)]]
  chg   <- last - first
  keep  <- !is.na(chg)
  chg   <- chg[keep]
  g     <- d$grp[keep]

  if (length(unique(g)) < 2 || sum(g == 0) < 2 || sum(g == 1) < 2) {
    return(c(est = NA_real_, se = NA_real_, n = length(chg)))
  }
  m1 <- mean(chg[g == 1]); m0 <- mean(chg[g == 0])
  v1 <- var(chg[g == 1]) / sum(g == 1)
  v0 <- var(chg[g == 0]) / sum(g == 0)
  c(est = m1 - m0, se = sqrt(v1 + v0), n = length(chg))
}

complete_case <- function(obs, n_visits) {
  cc <- obs[stats::complete.cases(obs), , drop = FALSE]
  if (nrow(cc) < 10) return(c(est = NA_real_, se = NA_real_, n = nrow(cc)))
  estimate_from(cc, n_visits)
}

run_mi <- function(obs, n_visits, m = 5, seed = 1, delta = 0) {
  if (!anyNA(obs)) {
    e <- estimate_from(obs, n_visits)
    return(list(est = e["est"], se = e["se"], fmi = 0, m = m, ok = TRUE))
  }

  imp <- tryCatch(
    mice(obs, m = m, method = "pmm", printFlag = FALSE, seed = seed,
         maxit = 5),
    error = function(e) NULL
  )
  if (is.null(imp)) {
    return(list(est = NA_real_, se = NA_real_, fmi = NA_real_, m = m, ok = FALSE))
  }

  ests <- numeric(m); ses <- numeric(m)
  last_col <- paste0("y", n_visits)
  for (i in seq_len(m)) {
    d <- mice::complete(imp, i)
    if (delta != 0) {
      was_missing <- is.na(obs[[last_col]])
      d[[last_col]][was_missing] <- d[[last_col]][was_missing] + delta
    }
    e <- estimate_from(d, n_visits)
    ests[i] <- e["est"]; ses[i] <- e["se"]
  }

  # Rubin's rules
  qbar <- mean(ests)
  ubar <- mean(ses^2)
  bvar <- if (m > 1) var(ests) else 0
  tvar <- ubar + (1 + 1 / m) * bvar
  fmi  <- if (tvar > 0) ((1 + 1 / m) * bvar) / tvar else 0

  list(est = qbar, se = sqrt(tvar), fmi = fmi, m = m, ok = TRUE,
       imp = imp)
}

# =================================================================== plots ===

pattern_plot <- function(obs_mat, n_visits) {
  n <- nrow(obs_mat)
  ord <- order(rowSums(obs_mat), decreasing = TRUE)
  M <- obs_mat[ord, , drop = FALSE]
  show <- min(n, 120)
  M <- M[round(seq(1, n, length.out = show)), , drop = FALSE]

  op <- par(mar = c(4, 5, 2, 1)); on.exit(par(op))
  plot(NA, xlim = c(0.5, n_visits + 0.5), ylim = c(0, nrow(M)),
       xaxt = "n", yaxt = "n", xlab = "Visit", ylab = "Participants",
       bty = "n", cex.lab = 1)
  axis(1, at = seq_len(n_visits), labels = seq_len(n_visits),
       col = GREY, col.axis = GREY, cex.axis = 0.9)
  for (i in seq_len(nrow(M))) {
    for (j in seq_len(n_visits)) {
      rect(j - 0.45, nrow(M) - i, j + 0.45, nrow(M) - i + 0.9,
           col = if (M[i, j]) NAVY else "#fbe3e3",
           border = if (M[i, j]) NAVY else RED, lwd = 0.3)
    }
  }
  legend("topright", legend = c("observed", "missing"), bty = "n",
         fill = c(NAVY, "#fbe3e3"), border = c(NAVY, RED), cex = 0.85)
}

retention_plot <- function(obs_mat, n_visits) {
  prop <- colMeans(obs_mat)
  op <- par(mar = c(4.5, 5, 2, 1)); on.exit(par(op))
  plot(seq_len(n_visits), prop, type = "n", ylim = c(0, 1),
       xlab = "Visit", ylab = "Proportion still observed",
       xaxt = "n", bty = "n", col.lab = GREY, col.axis = GREY)
  axis(1, at = seq_len(n_visits), col = GREY, col.axis = GREY)
  abline(h = seq(0, 1, 0.25), col = RULE, lwd = 1)
  lines(seq_len(n_visits), prop, col = NAVY, lwd = 2.5)
  points(seq_len(n_visits), prop, pch = 19, col = NAVY, cex = 1.4)
  text(seq_len(n_visits), prop, labels = sprintf("%.0f%%", prop * 100),
       pos = 3, col = GREY, cex = 0.85)
}

trajectory_plot <- function(full, observed, imputed, n_visits) {
  vis <- seq_len(n_visits)
  m_full <- sapply(vis, function(j) mean(full[[paste0("y", j)]]))
  m_obs  <- sapply(vis, function(j) mean(observed[[paste0("y", j)]], na.rm = TRUE))
  m_imp  <- if (!is.null(imputed))
    sapply(vis, function(j) mean(imputed[[paste0("y", j)]])) else NULL

  yr <- range(c(m_full, m_obs, m_imp), na.rm = TRUE)
  yr <- yr + c(-0.15, 0.25) * diff(yr)

  op <- par(mar = c(4.5, 5, 2, 1)); on.exit(par(op))
  plot(NA, xlim = c(1, n_visits), ylim = yr, xaxt = "n", bty = "n",
       xlab = "Visit", ylab = "Mean outcome", col.lab = GREY, col.axis = GREY)
  axis(1, at = vis, col = GREY, col.axis = GREY)
  abline(h = pretty(yr, 4), col = RULE)

  lines(vis, m_full, col = NAVY, lwd = 2.5)
  points(vis, m_full, pch = 19, col = NAVY, cex = 1.3)
  lines(vis, m_obs, col = RED, lwd = 2.5, lty = 2)
  points(vis, m_obs, pch = 17, col = RED, cex = 1.3)
  if (!is.null(m_imp)) {
    lines(vis, m_imp, col = MIDBLU, lwd = 2.5, lty = 3)
    points(vis, m_imp, pch = 15, col = MIDBLU, cex = 1.2)
  }
  labs <- c("complete data (unobservable)", "observed cases only")
  cols <- c(NAVY, RED); ltys <- c(1, 2); pchs <- c(19, 17)
  if (!is.null(m_imp)) {
    labs <- c(labs, "after imputation")
    cols <- c(cols, MIDBLU); ltys <- c(ltys, 3); pchs <- c(pchs, 15)
  }
  legend("topleft", bty = "n", cex = 0.85, text.col = GREY,
         legend = labs, col = cols, lwd = 2.5, lty = ltys, pch = pchs)
}

estimate_plot <- function(truth, cc, mi) {
  rows <- list(
    list(lab = "Truth",         est = truth["est"], se = truth["se"], col = NAVY),
    list(lab = "Complete case", est = cc["est"],    se = cc["se"],    col = RED),
    list(lab = "Imputed",       est = mi$est,       se = mi$se,       col = MIDBLU)
  )
  vals <- unlist(lapply(rows, function(r) c(r$est - 2 * r$se, r$est + 2 * r$se)))
  xr <- range(c(vals, 0), na.rm = TRUE)
  xr <- xr + c(-0.05, 0.05) * diff(xr)

  op <- par(mar = c(4.5, 8, 2, 1)); on.exit(par(op))
  plot(NA, xlim = xr, ylim = c(0.5, 3.5), yaxt = "n", bty = "n",
       xlab = "Estimated group difference in change", ylab = "",
       col.lab = GREY, col.axis = GREY)
  axis(2, at = 3:1, labels = sapply(rows, `[[`, "lab"), las = 1,
       col = GREY, col.axis = GREY, tick = FALSE)
  abline(v = truth["est"], col = RED, lty = 3)
  for (k in seq_along(rows)) {
    r <- rows[[k]]; yk <- 4 - k
    if (!is.na(r$se) && r$se > 0)
      segments(r$est - 1.96 * r$se, yk, r$est + 1.96 * r$se, yk,
               col = GREY, lwd = 2)
    points(r$est, yk, pch = 19, col = r$col, cex = 1.8)
  }
}

sensitivity_plot <- function(deltas, ests, truth_est) {
  op <- par(mar = c(4.5, 5, 2, 1)); on.exit(par(op))
  yr <- range(c(ests, truth_est), na.rm = TRUE)
  yr <- yr + c(-0.1, 0.1) * diff(yr)
  plot(NA, xlim = range(deltas), ylim = yr, bty = "n",
       xlab = "Delta applied to imputed final values",
       ylab = "Estimated group difference",
       col.lab = GREY, col.axis = GREY)
  abline(h = pretty(yr, 4), col = RULE)
  abline(h = truth_est, col = RED, lty = 3)
  abline(h = 0, col = GREY)
  lines(deltas, ests, col = NAVY, lwd = 2.5)
  points(deltas, ests, pch = 19, col = NAVY, cex = 1.2)
}

# ====================================================================== UI ===

css <- "
:root{--navy:#0a2d50;--red:#e12726;--grey:#4d545e;--rule:#d9dde3;--wash:#f2f4f7}
body{background:#fff;color:#1b1f26;font-family:'Inter',Helvetica,Arial,sans-serif}
.brandbar{background:var(--navy);color:#fff;padding:10px 18px;font-size:13px}
.brandbar span{color:#a9bdd4;font-size:10px;letter-spacing:.1em;text-transform:uppercase;
  margin-left:14px;padding-left:14px;border-left:1px solid #2c4c72}
.titlebar{border-top:3px solid var(--red);border-bottom:1px solid var(--rule);
  padding:16px 18px;margin-bottom:18px}
.titlebar h2{margin:0;color:var(--navy);font-size:24px;font-weight:600}
.titlebar p{margin:6px 0 0;color:var(--grey);font-size:14px;max-width:70ch}
.well{background:var(--wash);border:1px solid var(--rule);border-top:3px solid var(--navy);
  border-radius:0;box-shadow:none}
.nav-tabs>li.active>a{border-top:3px solid var(--red)!important;color:var(--navy)!important;
  font-weight:600}
.nav-tabs>li>a{color:var(--grey)}
h4{color:var(--navy);font-weight:600;margin-top:4px}
.stepnote{border-left:3px solid var(--red);background:var(--wash);padding:12px 16px;
  margin:14px 0;font-size:14px;color:var(--grey)}
.readout{display:flex;flex-wrap:wrap;gap:12px;margin:14px 0}
.readout div{background:var(--wash);padding:12px 16px;min-width:150px;flex:1}
.readout .lab{font-size:10px;letter-spacing:.12em;text-transform:uppercase;color:var(--grey)}
.readout .val{font-size:24px;font-weight:600;color:var(--navy);margin-top:4px}
.readout .note{font-size:11px;color:var(--grey);margin-top:2px}
.footnote{border-top:1px solid var(--rule);margin-top:24px;padding-top:14px;
  font-size:12px;color:var(--grey)}
"

ui <- fluidPage(
  tags$head(tags$style(HTML(css)), tags$title("MissLearn")),

  div(class = "brandbar",
      strong("King's College London"),
      span("Institute of Psychiatry, Psychology & Neuroscience")),

  div(class = "titlebar",
      h2("MissLearn"),
      p(paste("The four-step framework for identifying and handling missing data",
              "in longitudinal cohorts. Everything below is simulated, so the",
              "true answer is known and every method can be judged against it."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("The cohort"),
      sliderInput("n", "Participants", 100, 600, 300, step = 50),
      sliderInput("visits", "Visits", 3, 7, 5, step = 1),
      sliderInput("effect", "True group difference", 0, 1.2, 0.6, step = 0.1),
      tags$hr(),
      h4("The missing data"),
      sliderInput("dropout", "Dropout by the final visit", 0, 0.6, 0.3, step = 0.05),
      selectInput("mech", "Mechanism",
                  c("MCAR — unrelated to anything"        = "mcar",
                    "MAR — related to observed values"    = "mar",
                    "MNAR — related to the missing value" = "mnar"),
                  selected = "mar"),
      tags$hr(),
      h4("The imputation"),
      sliderInput("m", "Imputations", 1, 20, 5, step = 1),
      numericInput("seed", "Random seed", 1, min = 1, max = 9999),
      actionButton("run", "Run the analysis", class = "btn-primary"),
      div(class = "footnote",
          "Simulated data only. Nothing is stored.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "steps",

        tabPanel(
          "Step 1 — Describe",
          div(class = "stepnote",
              paste("Before choosing a method, look at what is actually missing.",
                    "How much, at which visits, and in what pattern.")),
          uiOutput("s1_stats"),
          fluidRow(
            column(6, h4("Missingness pattern"), plotOutput("p_pattern", height = "340px")),
            column(6, h4("Retention by visit"),  plotOutput("p_retention", height = "340px"))
          ),
          p(class = "footnote",
            paste("Dropout here is monotone: once a participant leaves they do not",
                  "return. Intermittent missingness is handled the same way but",
                  "produces a more scattered pattern."))
        ),

        tabPanel(
          "Step 2 — Diagnose",
          div(class = "stepnote",
              paste("Ask what predicts being missing. If observed variables predict",
                    "dropout, complete-case analysis is biased and imputation is not",
                    "optional. Data cannot tell you whether the mechanism is MNAR;",
                    "only Step 4 addresses that.")),
          h4("Probability of dropping out, modelled on observed values"),
          verbatimTextOutput("s2_model"),
          uiOutput("s2_verdict")
        ),

        tabPanel(
          "Step 3 — Handle",
          div(class = "stepnote",
              paste("Compare the methods against a truth that only a simulation can",
                    "provide. Watch what complete-case analysis does to the mean",
                    "trajectory as dropout rises.")),
          uiOutput("s3_stats"),
          fluidRow(
            column(6, h4("Mean trajectory"),      plotOutput("p_traj", height = "330px")),
            column(6, h4("Estimates compared"),   plotOutput("p_est",  height = "330px"))
          )
        ),

        tabPanel(
          "Step 4 — Sensitivity",
          div(class = "stepnote",
              paste("Multiple imputation assumes MAR. Sensitivity analysis asks how",
                    "far that assumption can be wrong before the conclusion changes.",
                    "Delta shifts every imputed final value by a fixed amount.")),
          h4("Estimate under departures from MAR"),
          plotOutput("p_sens", height = "360px"),
          uiOutput("s4_tipping")
        ),

        tabPanel(
          "Exercises",
          h4("Six exercises"),
          tags$ol(
            tags$li(paste("Set the mechanism to MCAR and raise dropout to 60%.",
                          "Does the complete-case estimate move away from the truth?",
                          "What does change?")),
            tags$li(paste("Switch to MAR at 30% dropout. Compare the complete-case",
                          "and imputed estimates. Which recovers the truth, and why?")),
            tags$li(paste("Switch to MNAR. Neither method recovers the truth.",
                          "Explain what information is missing that no analysis of",
                          "this dataset could supply.")),
            tags$li(paste("With MAR at 40%, set imputations to 1, then to 20.",
                          "The estimate barely moves. What does change, and which",
                          "quantity on Step 3 tells you why?")),
            tags$li(paste("On Step 2 under MCAR, what do the coefficients look like?",
                          "Under MAR? Write down the rule you would use on real data,",
                          "where the mechanism is unknown.")),
            tags$li(paste("On Step 4, find the delta at which the interval first",
                          "crosses zero. State in one sentence what that value means",
                          "for a reader of the paper."))
          ),
          div(class = "footnote",
              "Model answers accompany the teaching pack on the platform site.")
        )
      )
    )
  )
)

# ================================================================== server ===

server <- function(input, output, session) {

  sim <- eventReactive(input$run, {
    simulate_cohort(n = input$n, n_visits = input$visits,
                    effect = input$effect, dropout = input$dropout,
                    mechanism = input$mech, seed = input$seed)
  }, ignoreNULL = FALSE)

  mi_fit <- eventReactive(input$run, {
    s <- sim()
    run_mi(s$observed, s$n_visits, m = input$m, seed = input$seed)
  }, ignoreNULL = FALSE)

  # ------------------------------------------------------------- step 1 ----
  output$s1_stats <- renderUI({
    s <- sim()
    nv <- s$n_visits
    cells <- s$n * nv
    miss  <- sum(!s$obs)
    cc    <- sum(stats::complete.cases(s$observed))
    div(class = "readout",
        div(div(class = "lab", "Missing observations"),
            div(class = "val", sprintf("%.0f%%", 100 * miss / cells)),
            div(class = "note", sprintf("%d of %d cells", miss, cells))),
        div(div(class = "lab", "Complete cases"),
            div(class = "val", cc),
            div(class = "note", sprintf("of %d participants", s$n))),
        div(div(class = "lab", "Observed at final visit"),
            div(class = "val", sprintf("%.0f%%", 100 * mean(s$obs[, nv]))),
            div(class = "note", "where the estimand is measured")))
  })

  output$p_pattern   <- renderPlot(pattern_plot(sim()$obs, sim()$n_visits))
  output$p_retention <- renderPlot(retention_plot(sim()$obs, sim()$n_visits))

  # ------------------------------------------------------------- step 2 ----
  dropout_model <- reactive({
    s <- sim()
    dropped <- as.integer(!s$obs[, s$n_visits])
    d <- data.frame(dropped = dropped,
                    baseline = s$observed$y1,
                    penultimate = s$observed[[paste0("y", s$n_visits - 1)]],
                    grp = s$observed$grp)
    d <- d[!is.na(d$baseline), , drop = FALSE]
    if (length(unique(d$dropped)) < 2) return(NULL)
    tryCatch(
      glm(dropped ~ baseline + penultimate + grp, data = d,
          family = binomial(), na.action = na.omit),
      error = function(e) NULL
    )
  })

  output$s2_model <- renderPrint({
    fit <- dropout_model()
    if (is.null(fit)) {
      cat("No dropout in this cohort. Raise the dropout slider and run again.\n")
    } else {
      print(summary(fit)$coefficients, digits = 3)
    }
  })

  output$s2_verdict <- renderUI({
    fit <- dropout_model()
    if (is.null(fit)) return(NULL)
    p <- summary(fit)$coefficients[, 4]
    predictive <- any(p[-1] < 0.05, na.rm = TRUE)
    msg <- if (predictive) {
      paste("At least one observed variable predicts dropout. Complete-case",
            "analysis is biased here. Whether imputation is enough depends on",
            "whether the mechanism is MAR or MNAR, which this model cannot tell",
            "you.")
    } else {
      paste("No observed variable predicts dropout. That is consistent with MCAR,",
            "but it does not rule out MNAR: a dependence on the missing value",
            "itself leaves no trace in the observed data.")
    }
    div(class = "stepnote", msg)
  })

  # ------------------------------------------------------------- step 3 ----
  results <- reactive({
    s  <- sim()
    mi <- mi_fit()
    list(truth = estimate_from(s$full, s$n_visits),
         cc    = complete_case(s$observed, s$n_visits),
         mi    = mi)
  })

  output$s3_stats <- renderUI({
    r <- results()
    tr <- r$truth["est"]
    bias_cc <- 100 * abs(r$cc["est"] - tr) / abs(tr)
    bias_mi <- 100 * abs(r$mi$est - tr) / abs(tr)
    div(class = "readout",
        div(div(class = "lab", "Truth"),
            div(class = "val", sprintf("%.3f", tr)),
            div(class = "note", "available only in simulation")),
        div(div(class = "lab", "Complete case"),
            div(class = "val", sprintf("%.3f", r$cc["est"])),
            div(class = "note", sprintf("%.0f%% from the truth", bias_cc))),
        div(div(class = "lab", "Imputed"),
            div(class = "val", sprintf("%.3f", r$mi$est)),
            div(class = "note", sprintf("%.0f%% from the truth", bias_mi))),
        div(div(class = "lab", "Fraction of information missing"),
            div(class = "val", sprintf("%.2f", r$mi$fmi)),
            div(class = "note", sprintf("with %d imputations", r$mi$m))))
  })

  output$p_traj <- renderPlot({
    s <- sim(); mi <- mi_fit()
    imputed <- if (!is.null(mi$imp)) mice::complete(mi$imp, 1) else NULL
    trajectory_plot(s$full, s$observed, imputed, s$n_visits)
  })

  output$p_est <- renderPlot({
    r <- results()
    estimate_plot(r$truth, r$cc, r$mi)
  })

  # ------------------------------------------------------------- step 4 ----
  sens <- eventReactive(input$run, {
    s <- sim()
    deltas <- seq(-1.0, 1.0, by = 0.25)
    withProgress(message = "Refitting across delta", value = 0, {
      fits <- lapply(seq_along(deltas), function(k) {
        incProgress(1 / length(deltas))
        run_mi(s$observed, s$n_visits, m = max(3, min(input$m, 5)),
               seed = input$seed, delta = deltas[k])
      })
    })
    list(deltas = deltas,
         ests = vapply(fits, function(f) as.numeric(f$est), numeric(1)),
         ses  = vapply(fits, function(f) as.numeric(f$se),  numeric(1)))
  }, ignoreNULL = FALSE)

  output$p_sens <- renderPlot({
    sv <- sens()
    sensitivity_plot(sv$deltas, sv$ests, results()$truth["est"])
  })

  output$s4_tipping <- renderUI({
    sv <- sens()
    crosses <- which((sv$ests - 1.96 * sv$ses) <= 0)
    msg <- if (length(crosses) == 0) {
      paste("Across the whole delta range examined, the interval stays clear of",
            "zero. The conclusion is robust to this class of departure from MAR.")
    } else {
      paste0("The interval first includes zero at delta = ",
             sprintf("%.2f", sv$deltas[crosses[1]]),
             ". A reader should judge whether a departure of that size is",
             " plausible in this setting.")
    }
    div(class = "stepnote", msg)
  })
}

shinyApp(ui, server)
