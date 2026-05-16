# ==============================================================================
# CRISPR GuideX v3.1  —  by Taiwo
# ==============================================================================
# Advanced production-level bioinformatics Shiny application
#
# Features:
#   · SpCas9 (NGG) | SaCas9 (NNGRRT) | Cas12a (TTTV) | hfCas12Max (TN/TNN) | Custom PAM
#   · Adjustable guide length (15–24 nt)
#   · Weighted seed-region off-target scoring (seed = last 12 nt, weight 2x)
#   · Off-target detail table: BOTH Forward & Reverse strand scanning
#   · Hairpin / secondary structure detection heuristic
#   · Sequence viewer: Dynamic 5'/3' PAM & +/- Strand highlighting
#   · CRISPR simulation mode: Dynamic cut geometry prediction
#   · NCBI eUtils gene search (Robust coordinate extraction)
#   · Dynamic composite scoring with adjustable weights + ranking modes
#   · Batch FASTA processing + "Dual-Guide Best Pair" Feature
#   · Multi-format export: CSV / FASTA / Excel (openxlsx optional)
#   · GC filter sliders + custom sort in table
#   · Light / Dark mode toggle (Fully Global)
#   · Responsive layouts for Mobile/Tablet/Desktop
# ==============================================================================

# ── LIBRARIES ─────────────────────────────────────────────────────────────────
library(shiny)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(stringr)

# Optional packages — app degrades gracefully without them
has_xlsx     <- requireNamespace("openxlsx",  quietly = TRUE)
has_httr     <- requireNamespace("httr",      quietly = TRUE)
has_jsonlite <- requireNamespace("jsonlite",  quietly = TRUE)

# ==============================================================================
# ██  SCALAR SAFETY HELPERS
# ==============================================================================
sc  <- function(x) as.character(x[[1L]])  # safe character scalar
sd_ <- function(x) as.double(x[[1L]])     # safe double scalar
si_ <- function(x) as.integer(x[[1L]])    # safe integer scalar

# ==============================================================================
# ██  SEQUENCE PARSING
# ==============================================================================

parse_sequence <- function(raw) {
  if (is.null(raw) || nchar(trimws(raw)) == 0L) return(NULL)
  lines     <- strsplit(raw, "\n")[[1L]]
  seq_lines <- lines[!grepl("^[>;#]", lines) & nchar(trimws(lines)) > 0L]
  seq        <- toupper(paste(trimws(seq_lines), collapse = ""))
  seq        <- gsub("[^ATGCN]", "", seq)
  if (nchar(seq) == 0L) return(NULL)
  seq
}

parse_multi_fasta <- function(raw) {
  lines <- strsplit(raw, "\n")[[1L]]
  seqs  <- list(); cur_name <- NULL; cur_seq <- character(0L)
  for (ln in trimws(lines)) {
    if (grepl("^>", ln)) {
      if (!is.null(cur_name) && length(cur_seq) > 0L)
        seqs[[cur_name]] <- toupper(paste(cur_seq, collapse = ""))
      cur_name <- sub("^>\\s*", "", ln); cur_seq <- character(0L)
    } else if (nchar(ln) > 0L) {
      cur_seq <- c(cur_seq, gsub("[^ATGCNatgcn]", "", ln))
    }
  }
  if (!is.null(cur_name) && length(cur_seq) > 0L)
    seqs[[cur_name]] <- toupper(paste(cur_seq, collapse = ""))
  seqs
}

parse_fasta_file <- function(path) {
  if (is.null(path)) return(list())
  parse_multi_fasta(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

# ==============================================================================
# ██  CRISPR SYSTEM DEFINITIONS
# ==============================================================================

get_system_info <- function(system, custom_pam = NULL,
                            custom_side = "3prime", guide_len = 20L) {
  guide_len <- as.integer(guide_len)
  switch(system,
         SpCas9 = list(pattern = "[ACGT]GG",
                       pam_len = 3L, side = "3prime", label = "SpCas9 (NGG)"),
         SaCas9 = list(pattern = "[ACGT][ACGT]G[AG][AG]T",
                       pam_len = 6L, side = "3prime", label = "SaCas9 (NNGRRT)"),
         Cas12a = list(pattern = "TTT[ACG]",
                       pam_len = 4L, side = "5prime", label = "Cas12a (TTTV)"),
         Cas12Max = list(pattern = "T[ACGT]",
                         pam_len = 2L, side = "5prime", label = "hfCas12Max (TN/TNN)"),
         Custom = {
           pat <- if (!is.null(custom_pam) && nchar(trimws(custom_pam)) > 0L)
             trimws(custom_pam) else "[ACGT]GG"
           plen <- nchar(gsub("[^ACGTNacgtn]", "", gsub("\\[.*?\\]", "X", pat)))
           plen <- max(2L, as.integer(plen))
           list(pattern = pat, pam_len = plen, side = custom_side, label = "Custom PAM")
         },
         list(pattern = "[ACGT]GG", pam_len = 3L, side = "3prime", label = "SpCas9 (NGG)")
  )
}

# ==============================================================================
# ██  SCORING FUNCTIONS
# ==============================================================================

calc_gc <- function(seq) {
  seq <- toupper(seq); n <- nchar(seq)
  if (n == 0L) return(0.0)
  nchar(gsub("[^GC]", "", seq)) / n
}

efficiency_score <- function(seq, guide_len = 20L) {
  seq  <- toupper(seq)
  glen <- nchar(seq)
  gc   <- calc_gc(seq)
  gc_s <- if (gc >= 0.4 && gc <= 0.6)  1.0
  else if (gc < 0.4)           max(0.0, 1 - (0.4 - gc) * 3)
  else                         max(0.0, 1 - (gc - 0.6) * 3)
  p_pt <- if (grepl("TTTT",  seq)) 0.20 else 0.0
  p_hp <- if (grepl("AAAAA|TTTTT|GGGGG|CCCCC", seq)) 0.15 else 0.0
  g_j  <- if (glen >= 1L && substr(seq, glen, glen) == "G") 0.05 else 0.0
  seed_s <- max(1L, glen - 8L)
  sgc  <- calc_gc(substr(seq, seed_s, glen))
  s_b  <- if (sgc >= 0.3 && sgc <= 0.7) 0.05 else 0.0
  round(max(0.0, min(1.0, gc_s - p_pt - p_hp + g_j + s_b)), 4L)
}

detect_hairpin <- function(seq) {
  seq     <- toupper(seq)
  rc      <- paste(rev(strsplit(chartr("ATGC", "TACG", seq), "")[[1L]]), collapse = "")
  n       <- nchar(seq)
  if (n < 4L) return(FALSE)
  for (i in seq_len(n - 3L)) {
    sub4 <- substr(seq, i, i + 3L)
    if (grepl(sub4, rc, fixed = TRUE)) return(TRUE)
  }
  FALSE
}

weighted_mm <- function(a, b) {
  if (nchar(a) != nchar(b)) return(NA_real_)
  av <- strsplit(a, "")[[1L]]
  bv <- strsplit(b, "")[[1L]]
  n  <- length(av)
  total <- 0.0
  for (i in seq_len(n))
    if (av[i] != bv[i])
      total <- total + if (i >= (n - 11L)) 2.0 else 1.0
  total
}

off_target_full <- function(grna_seq, dna_seq, target_pos, max_wt = 6.0) {
  glen  <- nchar(grna_seq)
  dlen  <- nchar(dna_seq)
  hits  <- list(); hi <- 1L
  max_s <- dlen - glen - 2L
  if (max_s < 1L)
    return(list(n = 0L, risk = 0.0,
                hits = data.frame(OT_Position=integer(), Strand=character(), OT_Sequence=character(),
                                  Mismatches=integer(), Weighted_Score=numeric(),
                                  stringsAsFactors=FALSE)))
  
  # Check reverse complement for opposite strand off-targets
  rc_grna <- paste(rev(strsplit(chartr("ATGC", "TACG", grna_seq), "")[[1L]]), collapse = "")
  
  for (i in seq_len(max_s)) {
    if (abs(i - target_pos) <= 2L) next # Skip self
    cand <- substr(dna_seq, i, i + glen - 1L)
    if (nchar(cand) < glen) next
    
    # Forward Strand Match
    wt_fwd <- weighted_mm(grna_seq, cand)
    if (!is.na(wt_fwd) && wt_fwd <= max_wt) {
      raw_mm <- sum(strsplit(grna_seq, "")[[1L]] != strsplit(cand, "")[[1L]])
      hits[[hi]] <- data.frame(OT_Position=i, Strand="+", OT_Sequence=cand, 
                               Mismatches=raw_mm, Weighted_Score=round(wt_fwd, 2L), stringsAsFactors=FALSE)
      hi <- hi + 1L
    }
    
    # Reverse Strand Match
    wt_rev <- weighted_mm(rc_grna, cand)
    if (!is.na(wt_rev) && wt_rev <= max_wt) {
      raw_mm <- sum(strsplit(rc_grna, "")[[1L]] != strsplit(cand, "")[[1L]])
      hits[[hi]] <- data.frame(OT_Position=i, Strand="-", OT_Sequence=cand, 
                               Mismatches=raw_mm, Weighted_Score=round(wt_rev, 2L), stringsAsFactors=FALSE)
      hi <- hi + 1L
    }
  }
  
  n    <- length(hits)
  risk <- round(min(1.0, max(0.0, 1 - exp(-n * 0.5))), 4L)
  hits_df <- if (n > 0L) do.call(rbind, hits) else
    data.frame(OT_Position=integer(), Strand=character(), OT_Sequence=character(),
               Mismatches=integer(), Weighted_Score=numeric(), stringsAsFactors=FALSE)
  list(n = n, risk = risk, hits = hits_df)
}

# ==============================================================================
# ██  CRISPR SIMULATION (Dynamic Cut Geometry)
# ==============================================================================

simulate_cut <- function(dna_seq, grna_pos, grna_len, pam_side = "3prime", strand = "+") {
  dlen <- nchar(dna_seq)
  
  # Calculate biophysical cut position based on PAM geometry and strand
  if (strand == "+") {
    cut_pos <- if (pam_side == "3prime") as.integer(grna_pos + grna_len - 4L) else as.integer(grna_pos + 18L)
  } else {
    cut_pos <- if (pam_side == "3prime") as.integer(grna_pos + 3L) else as.integer(grna_pos + grna_len - 18L)
  }
  
  cut_pos <- max(1L, min(dlen, cut_pos))
  
  list(
    cut_pos    = cut_pos,
    upstream   = substr(dna_seq, max(1L, cut_pos - 30L), cut_pos),
    downstream = substr(dna_seq, cut_pos + 1L, min(dlen, cut_pos + 30L)),
    grna_region = substr(dna_seq, grna_pos, grna_pos + grna_len - 1L)
  )
}

# ==============================================================================
# ██  MAIN gRNA FINDER
# ==============================================================================

find_grna <- function(dna_seq, system = "SpCas9", guide_len = 20L,
                      custom_pam = NULL, custom_side = "3prime") {
  si      <- get_system_info(system, custom_pam, custom_side, guide_len)
  dlen    <- nchar(dna_seq)
  glen    <- as.integer(guide_len)
  plen    <- si$pam_len
  min_len <- glen + plen
  if (dlen < min_len) return(NULL)
  
  results <- list(); idx <- 1L
  
  scan_strand <- function(seq, rc = FALSE) {
    pam_indices <- gregexpr(si$pattern, seq, ignore.case = TRUE, perl = TRUE)[[1L]]
    if (pam_indices[1L] == -1L) return()
    
    for (p_pos in pam_indices) {
      p_pos <- as.integer(p_pos)
      
      if (si$side == "3prime") {
        gs_idx <- p_pos - glen
        if (gs_idx < 1L) next
        g <- substr(seq, gs_idx, p_pos - 1L)
        p <- substr(seq, p_pos, p_pos + plen - 1L)
      } else {
        gs_idx <- p_pos + plen
        if (gs_idx + glen - 1L > nchar(seq)) next
        p <- substr(seq, p_pos, p_pos + plen - 1L)
        g <- substr(seq, gs_idx, gs_idx + glen - 1L)
      }
      
      if (nchar(g) < glen || nchar(p) < plen) next
      fwd_pos <- if (rc) max(1L, dlen - (gs_idx + glen - 1L) + 1L) else gs_idx
      fwd_pos <- as.integer(fwd_pos)
      
      gc   <- calc_gc(g)
      eff  <- efficiency_score(g, glen)
      hp   <- detect_hairpin(g)
      ot   <- off_target_full(g, dna_seq, fwd_pos)
      
      results[[idx]] <<- data.frame(
        Position        = fwd_pos,
        Strand          = if (rc) "-" else "+",
        gRNA_Sequence   = g,
        PAM             = p,
        Guide_Length    = glen,
        GC_Content      = round(gc * 100L, 1L),
        Efficiency      = eff,
        Hairpin_Risk    = hp,
        Off_Targets     = ot$n,
        Off_Target_Risk = ot$risk,
        stringsAsFactors = FALSE
      )
      idx <<- idx + 1L
    }
  }
  
  scan_strand(dna_seq, FALSE)
  rc_seq <- paste(rev(strsplit(chartr("ATGC", "TACG", dna_seq), "")[[1L]]), collapse = "")
  scan_strand(rc_seq, TRUE)
  
  if (length(results) == 0L) return(NULL)
  df <- do.call(rbind, results)
  df <- df[order(df$Position), ]; rownames(df) <- NULL
  df
}

# ==============================================================================
# ██  SCORING LABELS & COMPOSITE SCORE
# ==============================================================================

risk_label <- function(eff, risk) {
  dplyr::case_when(
    eff >= 0.7 & risk <= 0.3 ~ "High Efficiency / Low Risk",
    eff >= 0.4 & risk <= 0.6 ~ "Medium",
    TRUE                     ~ "Low Efficiency / High Risk"
  )
}

composite_score <- function(eff, risk, w_eff = 0.6, w_risk = 0.4) {
  round(eff * w_eff + (1 - risk) * w_risk, 4L)
}

# ==============================================================================
# ██  SEQUENCE COLORIZER
# ==============================================================================
NT_COLS <- c(A = "#f87171", T = "#60a5fa", G = "#4ade80", C = "#facc15", N = "#94a3b8")

colorize_seq <- function(seq_str, pam_str = NULL) {
  seq_str <- as.character(seq_str)[[1L]]
  chars   <- strsplit(seq_str, "")[[1L]]
  html    <- paste(vapply(chars, function(nt) {
    col <- if (nt %in% names(NT_COLS)) NT_COLS[[nt]] else "#e2e8f0"
    sprintf('<span style="color:%s;font-weight:700">%s</span>', col, nt)
  }, character(1L)), collapse = "")
  
  if (!is.null(pam_str) && nchar(as.character(pam_str)[[1L]]) > 0L) {
    pam_str  <- as.character(pam_str)[[1L]]
    pam_html <- paste(vapply(strsplit(pam_str, "")[[1L]], function(nt)
      sprintf('<span style="color:#00d4aa;font-weight:800">%s</span>', nt),
      character(1L)), collapse = "")
    html <- paste0(html,
                   '<span style="color:var(--mut);margin:0 4px">|</span>',
                   '<span style="background:rgba(0,212,170,.18);border-radius:3px;padding:1px 3px">',
                   pam_html, '</span>')
  }
  html
}

# ==============================================================================
# ██  SEQUENCE VIEWER (Dynamic Strand & PAM Geometry)
# ==============================================================================

build_seq_viewer <- function(dna_seq, grna_pos, grna_len, pam_len,
                             pam_side = "3prime", strand = "+", ot_positions = NULL,
                             window = 80L) {
  dlen   <- nchar(dna_seq)
  start  <- max(1L, as.integer(grna_pos) - 10L)
  end    <- min(dlen, as.integer(grna_pos) + grna_len + pam_len + 10L)
  end    <- min(end, start + window - 1L)
  region <- substr(dna_seq, start, end)
  rlen   <- nchar(region)
  chars  <- strsplit(region, "")[[1L]]
  
  gs_r <- grna_pos - start + 1L
  ge_r <- gs_r + grna_len - 1L
  
  # Calculate precise PAM boundaries based on strand and 5'/3' rules
  if (strand == "+") {
    ps_r <- if (pam_side == "3prime") ge_r + 1L else gs_r - pam_len
    pe_r <- if (pam_side == "3prime") ps_r + pam_len - 1L else gs_r - 1L
  } else {
    ps_r <- if (pam_side == "3prime") gs_r - pam_len else ge_r + 1L
    pe_r <- if (pam_side == "3prime") gs_r - 1L else ps_r + pam_len - 1L
  }
  
  ot_set <- if (!is.null(ot_positions)) {
    r <- ot_positions - start + 1L
    r[r >= 1L & r <= rlen]
  } else integer(0L)
  
  spans <- vapply(seq_len(rlen), function(i) {
    nt  <- chars[i]
    col <- if (nt %in% names(NT_COLS)) NT_COLS[[nt]] else "var(--txt)"
    if (i >= gs_r && i <= ge_r) {
      sprintf('<span title="gRNA" style="color:%s;font-weight:800;background:rgba(0,212,170,.18);border-radius:2px">%s</span>', col, nt)
    } else if (i >= ps_r && i <= pe_r) {
      sprintf('<span title="PAM" style="color:#00d4aa;font-weight:800;background:rgba(0,212,170,.35);border-radius:2px">%s</span>', nt)
    } else if (i %in% ot_set) {
      sprintf('<span title="Off-target" style="color:#ef4444;font-weight:800;background:rgba(239,68,68,.22);border-radius:2px">%s</span>', nt)
    } else {
      sprintf('<span style="color:%s">%s</span>', col, nt)
    }
  }, character(1L))
  
  legend <- paste0(
    '<div style="display:flex;gap:12px;flex-wrap:wrap;margin-top:8px;font-family:var(--mono);font-size:10px;color:var(--mut);">',
    '<span style="background:rgba(0,212,170,.18);padding:1px 5px;border-radius:3px;color:var(--txt)">gRNA</span>',
    '<span style="background:rgba(0,212,170,.35);padding:1px 5px;border-radius:3px;color:#00d4aa">PAM</span>',
    if (length(ot_set) > 0L)
      '<span style="background:rgba(239,68,68,.22);padding:1px 5px;border-radius:3px;color:#ef4444">Off-target</span>'
    else '',
    '</div>')
  
  paste0(
    '<div style="font-family:var(--mono);font-size:12px;letter-spacing:1.5px;line-height:2;word-break:break-all;">',
    '<div style="color:var(--txt);font-size:10px;margin-bottom:8px;font-weight:700;">Viewing Region Starts @ Pos: ', start, '</div>',
    paste(spans, collapse = ""),
    '</div>', legend)
}

# ==============================================================================
# ██  NCBI GENE SEQUENCE FETCH (STRAND-SAFE FIX)
# ==============================================================================

fetch_ncbi_sequence <- function(gene_name, species = "Homo sapiens",
                                max_len = 2000L) {
  if (!has_httr || !has_jsonlite)
    return(list(ok = FALSE, msg = "httr and jsonlite packages required for NCBI fetch"))
  
  tryCatch({
    r1 <- httr::GET(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
      query = list(db="gene", term=sprintf("%s[Gene Name] AND %s[Organism]", gene_name, species), 
                   retmax=1, retmode="json"),
      httr::timeout(30))
    if (httr::status_code(r1) != 200L)
      return(list(ok = FALSE, msg = "NCBI search failed"))
    
    j1  <- jsonlite::fromJSON(httr::content(r1, "text", encoding = "UTF-8"))
    ids <- j1$esearchresult$idlist
    if (length(ids) == 0L)
      return(list(ok = FALSE, msg = paste("Gene not found:", gene_name)))
    
    gene_id <- ids[[1L]]
    
    r2 <- httr::GET(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
      query = list(db="gene", id=gene_id, retmode="json"),
      httr::timeout(30))
    if (httr::status_code(r2) != 200L)
      return(list(ok = FALSE, msg = "NCBI summary failed"))
    
    j2        <- jsonlite::fromJSON(httr::content(r2, "text", encoding = "UTF-8"))
    gene_info <- j2$result[[gene_id]]
    
    g_info <- gene_info$genomicinfo
    if (is.null(g_info) || length(g_info) == 0L || nrow(g_info) == 0)
      return(list(ok = FALSE, msg = "Could not find genomic coordinates for this gene"))
    
    mapping <- g_info[1, ]
    accession <- mapping$chraccver
    
    # CRITICAL FIX: Safely handle both Forward (+) and Reverse (-) strand coordinates
    chr_val_1 <- as.double(mapping$chrstart)
    chr_val_2 <- as.double(mapping$chrstop)
    
    true_start <- min(chr_val_1, chr_val_2)
    true_stop  <- max(chr_val_1, chr_val_2)
    
    fetch_s <- as.integer(true_start + 1L)
    # Force the fetch window to NEVER exceed max_len (default 2000bp)
    fetch_e <- as.integer(min(fetch_s + max_len - 1L, true_stop))
    
    r3 <- httr::GET(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
      query = list(db="nuccore", id=accession, seq_start=fetch_s, seq_stop=fetch_e, 
                   rettype="fasta", retmode="text"),
      httr::timeout(60))
    if (httr::status_code(r3) != 200L)
      return(list(ok = FALSE, msg = "Sequence fetch failed"))
    
    raw_fasta <- httr::content(r3, "text", encoding = "UTF-8")
    seqs      <- parse_multi_fasta(raw_fasta)
    if (length(seqs) == 0L)
      return(list(ok = FALSE, msg = "NCBI returned empty sequence"))
    
    seq <- gsub("[^ATGCNatgcn]", "", toupper(seqs[[1L]]))
    list(ok = TRUE, seq = seq, gene = as.character(gene_info$name), 
         msg = sprintf("Fetched %s bp of %s", format(nchar(seq), big.mark=","), gene_info$name))
  }, error = function(e)
    list(ok = FALSE, msg = paste("Error:", conditionMessage(e))))
}

# ==============================================================================
# ██  EXPORT HELPERS
# ==============================================================================

df_to_fasta <- function(df) {
  lines <- character(nrow(df) * 2L)
  for (i in seq_len(nrow(df))) {
    lines[2L*i - 1L] <- sprintf(
      ">gRNA_%02d|pos%d|strand%s|gc%.0f|eff%.3f|risk%.3f",
      i, df$Position[i], df$Strand[i],
      df$GC_Content[i], df$Efficiency[i], df$Off_Target_Risk[i])
    lines[2L*i] <- df$gRNA_Sequence[i]
  }
  paste(lines, collapse = "\n")
}

# ==============================================================================
# ██  EXAMPLE SEQUENCES
# ==============================================================================
EXAMPLE_SEQS <- list(
  "TP53 (exon 5-8 fragment)" = paste0(
    "ATGTTCAAGACTTTTTGCCAGAGCCTGAGCTGTACCACCATCCACTACAACTACATGTGT",
    "AACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAA",
    "GACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCT",
    "GGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGTGGAGCCTGCAGAGCAG",
    "CTGGAGACCTTGCAGCCAGCAGCAAATGAATCTTCCATCATGCCGCCTGGAGGCGGCATG",
    "AACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGA",
    "CGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAA",
    "GAGAATCTCCGCAAGAAAGTGGAGCCTGCAGAGCAGCTGGAG"),
  "BRCA1 (exon 2 fragment)" = paste0(
    "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAG",
    "AGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAA",
    "ATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGAT")
)


# ==============================================================================
# ██  CSS (FULLY RESPONSIVE & DYNAMIC THEME)
# ==============================================================================
APP_CSS <- '
:root{
  /* DEFAULT: DARK MODE */
  --bg: #0a0e1a; 
  --surf: #101520; 
  --card: #16202e; 
  --card2: #1c2a3a; 
  --bdr: #1e2d45;
  --acc: #00d4aa; 
  --a2: #7c6cfc; 
  --a3: #f97316; 
  --a4: #38bdf8;
  --txt: #e2e8f0; 
  --txt2: #cbd5e1; 
  --mut: #64748b;
  --err: #ef4444; 
  --ok: #22c55e; 
  --wrn: #eab308;
  --hdr-bg: linear-gradient(135deg, #0a1422, #101520);
  
  --mono: "Space Mono", monospace; 
  --sans: "Syne", sans-serif;
  --r: 9px; 
  --th: 56px;
}

/* OVERRIDES: LIGHT MODE */
body.lm {
  --bg: #f0f4f8; 
  --surf: #ffffff; 
  --card: #f8fafc; 
  --card2: #f1f5f9; 
  --bdr: #cbd5e1;
  --txt: #0f172a; 
  --txt2: #334155; 
  --mut: #64748b;
  --hdr-bg: linear-gradient(135deg, #e2f1ee, #f8fafc);
}

/* BASE */
*{box-sizing:border-box;margin:0;padding:0;}
body{background:var(--bg);color:var(--txt);font-family:var(--sans);min-height:100vh;transition:background .25s,color .25s;}
::-webkit-scrollbar{width:5px;height:5px;}
::-webkit-scrollbar-track{background:var(--bg);}
::-webkit-scrollbar-thumb{background:var(--bdr);border-radius:4px;}
body::before{content:"";position:fixed;inset:0;z-index:0;
  background-image:linear-gradient(rgba(0,212,170,.018) 1px,transparent 1px),
  linear-gradient(90deg,rgba(0,212,170,.018) 1px,transparent 1px);
  background-size:44px 44px;pointer-events:none;}
.container-fluid{position:relative;z-index:1;padding:0!important;}

/* NAV */
.tnav{position:sticky;top:0;z-index:999;background:var(--surf);
  border-bottom:1px solid var(--bdr);height:var(--th);
  display:flex;align-items:center;justify-content:space-between;
  padding:0 22px;box-shadow:0 4px 20px rgba(0,0,0,.15); transition:all 0.3s;}
.nbrand{font-family:var(--sans);font-weight:800;font-size:17px;color:var(--txt);
  cursor:pointer;user-select:none;display:flex;align-items:center;gap:1px;letter-spacing:-.4px;}
.nbrand .x{color:var(--acc);}
.nbrand .ver{font-family:var(--mono);font-size:9px;color:var(--mut);margin-left:6px;letter-spacing:1px;}
.nlinks{display:flex;align-items:center;gap:2px;}
.nb{background:none;border:none;color:var(--mut);font-family:var(--mono);font-size:11px;
  padding:6px 12px;border-radius:7px;cursor:pointer;transition:color .15s,background .15s;}
.nb:hover{background:rgba(0,212,170,.05);color:var(--txt);}
.nb.an{color:var(--acc);background:rgba(0,212,170,.09);}
.ndiv{width:1px;height:18px;background:var(--bdr);margin:0 4px;}
.tbtn{background:var(--card);border:1px solid var(--bdr);border-radius:18px;
  padding:4px 11px;cursor:pointer;font-family:var(--mono);font-size:10px;color:var(--mut);
  display:flex;align-items:center;gap:5px;transition:all .2s;}
.tbtn:hover{border-color:var(--acc);color:var(--acc);}

/* PAGES */
.page{display:none;}.page.ap{display:block;}

/* HEADER */
.ahdr{background:var(--hdr-bg);border-bottom:1px solid var(--bdr);
  padding:16px 26px 14px;position:relative;overflow:hidden;transition:background 0.3s;}
.ahdr::after{content:"";position:absolute;top:0;right:0;width:360px;height:100%;
  background:radial-gradient(ellipse at right,rgba(0,212,170,.06),transparent 70%);pointer-events:none;}
.hbdg{display:inline-block;background:rgba(0,212,170,.1);border:1px solid rgba(0,212,170,.25);
  color:var(--acc);font-family:var(--mono);font-size:9px;letter-spacing:2px;text-transform:uppercase;
  padding:3px 10px;border-radius:100px;margin-bottom:7px;}
.htitle{font-family:var(--sans);font-weight:800;font-size:clamp(14px,2.4vw,22px);
  color:var(--txt);letter-spacing:-.5px;margin-bottom:3px;line-height:1.2;}
.htitle .tx{color:var(--acc);}
.hsub{color:var(--mut);font-family:var(--mono);font-size:10px;letter-spacing:.3px;}
.hmeta{display:flex;gap:7px;margin-top:9px;flex-wrap:wrap;}
.mc{background:var(--card);border:1px solid var(--bdr);color:var(--mut);
  font-family:var(--mono);font-size:9px;padding:2px 7px;border-radius:5px;}
.mc b{color:var(--txt);}

/* LAYOUT */
.ml{display:grid;grid-template-columns:295px 1fr;
  min-height:calc(100vh - var(--th) - 52px);}

/* SIDEBAR & INPUTS */
.sb{background:var(--surf);border-right:1px solid var(--bdr);
  padding:14px 13px;display:flex;flex-direction:column;gap:11px;overflow-y:auto;
  max-height:calc(100vh - var(--th) - 52px);transition:background 0.3s;}
.slbl{font-family:var(--mono);font-size:9px;letter-spacing:2px;text-transform:uppercase;
  color:var(--acc);margin-bottom:4px;display:flex;align-items:center;gap:6px;}
.slbl::after{content:"";flex:1;height:1px;background:var(--bdr);}

/* Critical Input Overrides for Light/Dark Mode */
textarea#seq_input, input[type=text], input[type=number], input[type=file] {
  background: var(--card) !important; 
  border: 1px solid var(--bdr) !important;
  color: var(--txt) !important; 
  border-radius: var(--r); 
  font-family: var(--mono);
  font-size: 10px; 
  padding: 8px !important; 
  width: 100%;
  box-sizing: border-box;
}
textarea#seq_input:focus {outline:none!important;border-color:var(--acc)!important;box-shadow:0 0 0 3px rgba(0,212,170,.1)!important;}
textarea#seq_input::placeholder {color:var(--mut)!important;}

/* Selectize UI Overrides */
.selectize-input, .selectize-input.full {
  background: var(--card) !important; 
  color: var(--txt) !important; 
  border: 1px solid var(--bdr) !important; 
  border-radius: 7px !important; 
  font-family: var(--mono) !important; 
  font-size: 10px !important;
}
.selectize-dropdown {
  background: var(--card) !important; 
  border: 1px solid var(--bdr) !important;
}
.selectize-dropdown-content .option {
  color: var(--txt) !important; 
  font-family: var(--mono) !important; 
  font-size: 10px !important;
}
.selectize-dropdown-content .option:hover, .selectize-dropdown-content .option.active {
  background: var(--bdr) !important; 
  color: var(--txt) !important;
}

input[type=range]{accent-color:var(--acc);width:100%;}
label{color:var(--mut)!important;font-family:var(--mono)!important;font-size:10px!important;letter-spacing:.3px;}

/* BUTTONS */
.btn-run{background:linear-gradient(135deg,#00d4aa,#00a88a)!important;border:none!important;
  color:#0a1422!important;font-family:var(--sans)!important;font-weight:700!important;font-size:12px!important;
  padding:10px 17px!important;border-radius:var(--r)!important;width:100%;cursor:pointer;
  box-shadow:0 4px 14px rgba(0,212,170,.25)!important;transition:all .15s!important;}
.btn-run:hover{transform:translateY(-1px)!important;box-shadow:0 6px 20px rgba(0,212,170,.35)!important;}
.btn-ghost{background:transparent!important;border:1px solid var(--bdr)!important;
  color:var(--mut)!important;font-family:var(--mono)!important;font-size:10px!important;
  padding:7px 11px!important;border-radius:7px!important;width:100%;cursor:pointer;transition:all .18s!important;}
.btn-ghost:hover{border-color:var(--acc)!important;color:var(--acc)!important;}
.btn-danger{background:transparent!important;border:1px solid rgba(239,68,68,.3)!important;
  color:rgba(239,68,68,.7)!important;font-family:var(--mono)!important;font-size:10px!important;
  padding:7px 11px!important;border-radius:7px!important;width:100%;cursor:pointer;transition:all .18s!important;}
.btn-danger:hover{border-color:var(--err)!important;color:var(--err)!important;}
.btn-a2{background:transparent!important;border:1px solid rgba(124,108,252,.3)!important;
  color:var(--a2)!important;font-family:var(--mono)!important;font-size:10px!important;
  padding:7px 11px!important;border-radius:7px!important;width:100%;cursor:pointer;transition:all .18s!important;}
.btn-a2:hover{background:rgba(124,108,252,.08)!important;}
.egrp{display:flex;flex-direction:column;gap:5px;}
#dl_csv,#dl_fasta,#dl_excel,#dl_top{background:transparent!important;border:1px solid var(--bdr)!important;
  color:var(--mut)!important;font-family:var(--mono)!important;font-size:10px!important;
  padding:6px 11px!important;border-radius:7px!important;width:100%;cursor:pointer;transition:all .18s!important;}
#dl_csv:hover{border-color:var(--ok)!important;color:var(--ok)!important;}
#dl_fasta:hover{border-color:var(--a2)!important;color:var(--a2)!important;}
#dl_excel:hover{border-color:var(--wrn)!important;color:var(--wrn)!important;}
#dl_top:hover{border-color:var(--acc)!important;color:var(--acc)!important;}

/* INFO / ALERT BOXES */
.ic{background:var(--card2);border:1px solid var(--bdr);border-radius:8px;
  padding:9px 11px;font-family:var(--mono);font-size:10px;line-height:1.7;color:var(--mut);}
.ic strong{color:var(--acc);}
.ab{border-radius:8px;padding:10px 13px;font-family:var(--mono);font-size:11px;line-height:1.6;}
.aerr{background:rgba(239,68,68,.08);border:1px solid rgba(239,68,68,.25);color:var(--err);}
.awrn{background:rgba(234,179,8,.08);border:1px solid rgba(234,179,8,.25);color:var(--wrn);}
.ainf{background:rgba(0,212,170,.06);border:1px solid rgba(0,212,170,.2);color:var(--acc);}
.asuc{background:rgba(34,197,94,.06);border:1px solid rgba(34,197,94,.2);color:var(--ok);}
.wrow{display:grid;grid-template-columns:1fr 1fr;gap:7px;}

/* RESULTS PANEL */
.rp{padding:16px 20px;display:flex;flex-direction:column;gap:14px;overflow-x:hidden;}

/* STAT CARDS */
.srow{display:grid;grid-template-columns:repeat(auto-fit,minmax(125px,1fr));gap:9px;}
.sc{background:var(--card);border:1px solid var(--bdr);border-radius:10px;
  padding:11px 13px;position:relative;overflow:hidden;transition:border-color .2s;}
.sc:hover{border-color:var(--a2);}
.sc::before{content:"";position:absolute;top:0;left:0;right:0;height:2px;}
.sc.cg::before{background:var(--acc);}
.sc.cp::before{background:var(--a2);}
.sc.co::before{background:var(--a3);}
.sc.cb::before{background:var(--a4);}
.stat-val{font-family:var(--mono);font-size:21px;font-weight:700;color:var(--txt);line-height:1;margin-bottom:3px;}
.stat-lbl{font-family:var(--mono);font-size:9px;letter-spacing:1.5px;text-transform:uppercase;color:var(--mut);}
.stat-bdg{position:absolute;top:9px;right:9px;font-family:var(--mono);font-size:8px;padding:2px 6px;border-radius:4px;}
.bg{background:rgba(34,197,94,.15);color:var(--ok);border:1px solid rgba(34,197,94,.25);}
.bp{background:rgba(124,108,252,.15);color:var(--a2);border:1px solid rgba(124,108,252,.25);}
.bo{background:rgba(249,115,22,.15);color:var(--a3);border:1px solid rgba(249,115,22,.25);}
.bb{background:rgba(56,189,248,.15);color:var(--a4);border:1px solid rgba(56,189,248,.25);}

/* DUAL-GUIDE PAIR CARD */
.pair-card {
  background: var(--card2); border: 1px solid var(--acc); border-radius: 10px;
  padding: 14px 16px; margin-bottom: 5px; position: relative; overflow: hidden;
}
.pair-card::before {
  content: ""; position: absolute; top: 0; left: 0; width: 4px; height: 100%; background: var(--acc);
}
.pair-title {
  font-family: var(--sans); font-weight: 700; font-size: 13px; color: var(--txt);
  margin-bottom: 12px; display: flex; align-items: center; gap: 8px;
}
.pair-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }
.pair-item { background: var(--card); border: 1px solid var(--bdr); border-radius: 8px; padding: 12px; }
.pair-lbl {
  font-family: var(--mono); font-size: 9px; color: var(--mut);
  text-transform: uppercase; letter-spacing: 1px; margin-bottom: 6px;
}

/* TOP3 CANDIDATE CARDS */
.t3g{display:grid;grid-template-columns:repeat(auto-fit,minmax(215px,1fr));gap:10px;}
.cc{background:var(--card);border:1px solid var(--bdr);border-radius:10px;
  padding:12px 13px;position:relative;overflow:hidden;transition:border-color .2s,transform .15s;}
.cc:hover{border-color:var(--acc);transform:translateY(-2px);}
.cc::before{content:"";position:absolute;top:0;left:0;right:0;height:2px;}
.cc.r1::before{background:linear-gradient(90deg,#fbbf24,#f97316);}
.cc.r2::before{background:linear-gradient(90deg,#94a3b8,var(--txt));}
.cc.r3::before{background:linear-gradient(90deg,var(--acc),var(--a2));}
.crnk{font-family:var(--mono);font-size:9px;letter-spacing:1.5px;color:var(--mut);
  text-transform:uppercase;margin-bottom:5px;display:flex;align-items:center;gap:4px;}
.cseq{font-family:var(--mono);font-size:10px;font-weight:700;letter-spacing:1px;
  word-break:break-all;line-height:1.8;margin-bottom:7px;}
.cmeta{display:flex;flex-wrap:wrap;gap:4px;font-family:var(--mono);font-size:9px;color:var(--mut);margin-bottom:6px;}
.cpill{background:rgba(0,212,170,.08);border:1px solid rgba(0,212,170,.18);color:var(--acc);
  font-family:var(--mono);font-size:8px;padding:2px 6px;border-radius:100px;}
.cpill.rhi{background:rgba(239,68,68,.08);border-color:rgba(239,68,68,.2);color:var(--err);}
.cpill.rmd{background:rgba(234,179,8,.08);border-color:rgba(234,179,8,.2);color:var(--wrn);}
.cpill.hp{background:rgba(249,115,22,.08);border-color:rgba(249,115,22,.25);color:var(--a3);}
.cbtn{background:rgba(124,108,252,.1);border:1px solid rgba(124,108,252,.2);color:var(--a2);
  font-family:var(--mono);font-size:9px;padding:3px 8px;border-radius:5px;cursor:pointer;transition:all .15s;}
.cbtn:hover{background:rgba(124,108,252,.2);}
.cbtn.copied{background:rgba(34,197,94,.12);border-color:rgba(34,197,94,.3);color:var(--ok);}

/* TABS */
.tabnav{display:flex;gap:2px;border-bottom:1px solid var(--bdr);margin-bottom:12px;}
.tabbtn{background:none;border:none;color:var(--mut);font-family:var(--mono);font-size:10px;
  letter-spacing:.3px;padding:8px 12px;cursor:pointer;border-bottom:2px solid transparent;
  margin-bottom:-1px;transition:color .18s;white-space:nowrap;}
.tabbtn:hover{color:var(--txt);}
.tabbtn.active{color:var(--acc);border-bottom-color:var(--acc);}
.chard{background:var(--card);border:1px solid var(--bdr);border-radius:10px;padding:13px;overflow:hidden;}
.chtitle{font-family:var(--mono);font-size:9px;letter-spacing:1.5px;text-transform:uppercase;color:var(--mut);margin-bottom:9px;}

/* SEQUENCE VIEWER & SIMULATION */
.seq-vbox{background:var(--card2);border:1px solid var(--bdr);border-radius:8px;
  padding:13px;overflow-x:auto;font-size:12px;line-height:2;}
.sim-box{background:var(--card2);border:1px solid rgba(0,212,170,.2);border-radius:8px;padding:13px;}
.sim-lbl{font-family:var(--mono);font-size:9px;letter-spacing:1.5px;text-transform:uppercase;
  color:var(--acc);margin-bottom:7px;}
.sim-seq{font-family:var(--mono);font-size:11px;letter-spacing:1px;line-height:2;word-break:break-all;}
.cut{border-left:2px solid var(--err);padding-left:3px;margin-left:2px;color:var(--err);font-weight:700;}

/* DT TABLES OVERRIDES */
.dataTables_wrapper{font-family:var(--mono)!important;font-size:11px!important;color:var(--txt)!important;}
table.dataTable{background:transparent!important;border-collapse:collapse!important;}
table.dataTable thead tr th{background:var(--surf)!important;color:var(--mut)!important;
  border-bottom:1px solid var(--bdr)!important;font-family:var(--mono)!important;
  font-size:9px!important;letter-spacing:1px!important;text-transform:uppercase!important;padding:8px 10px!important;}
table.dataTable tbody tr td{border-bottom:1px solid var(--bdr)!important;
  padding:8px 10px!important;color:var(--txt)!important;}
table.dataTable tbody tr:hover td{background:var(--bdr)!important;}
table.dataTable tbody tr.selected td{background:rgba(0,212,170,.06)!important;}
.dataTables_info,.dataTables_paginate,.dataTables_filter,.dataTables_length{
  color:var(--mut)!important;font-family:var(--mono)!important;font-size:10px!important;}
.dataTables_filter input, .dataTables_length select{background:var(--card)!important;
  border:1px solid var(--bdr)!important;color:var(--txt)!important;border-radius:6px!important;padding:3px 8px!important;}
.paginate_button{color:var(--mut)!important;border:none!important;background:none!important;}
.paginate_button.current{background:var(--card)!important;color:var(--acc)!important;border-radius:5px!important;}
.paginate_button:hover:not(.disabled){background:var(--bdr)!important;color:var(--txt)!important;border-radius:5px!important;}

/* EMPTY STATE */
.estate{display:flex;flex-direction:column;align-items:center;justify-content:center;
  padding:55px 28px;text-align:center;color:var(--mut);}
.eicon{font-size:38px;margin-bottom:11px;opacity:.3;}
.estate h3{font-family:var(--sans);font-weight:600;color:var(--txt);opacity:.4;margin-bottom:5px;}
.estate p{font-family:var(--mono);font-size:11px;max-width:255px;line-height:1.6;}

/* HOW IT WORKS */
.hiw{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:9px;margin:10px 0;}
.hstep{background:var(--card2);border:1px solid var(--bdr);border-radius:8px;padding:11px;position:relative;}
.hstep::before{content:"";position:absolute;top:0;left:0;right:0;height:2px;background:var(--a2);}
.hn{font-family:var(--mono);font-size:9px;color:var(--a2);margin-bottom:4px;}
.ht{font-family:var(--sans);font-weight:700;font-size:12px;color:var(--txt);margin-bottom:3px;}
.hd{font-family:var(--mono);font-size:10px;color:var(--mut);line-height:1.5;}

/* DOCS / ABOUT */
.pdocs,.pabout{max-width:800px;margin:0 auto;padding:36px 24px 65px;}
.dhero,.ahero{margin-bottom:28px;padding-bottom:20px;border-bottom:1px solid var(--bdr);}
.dhero h1,.ahero h1{font-family:var(--sans);font-weight:800;font-size:clamp(19px,3.2vw,28px);
  color:var(--txt);letter-spacing:-.5px;margin-bottom:7px;}
.dhero h1 span,.ahero h1 span{color:var(--acc);}
.dhero p,.ahero p{font-family:var(--mono);font-size:11px;line-height:1.8;color:var(--mut);}
.dsec{margin-bottom:28px;}
.dsec-t{font-family:var(--sans);font-weight:700;font-size:15px;color:var(--txt);
  margin-bottom:10px;display:flex;align-items:center;gap:7px;padding-bottom:5px;border-bottom:1px solid var(--bdr);}
.dsec p{font-family:var(--mono);font-size:11px;line-height:1.9;color:var(--mut);margin-bottom:8px;}
.dsec ul{list-style:none;margin:0 0 8px;}
.dsec ul li{font-family:var(--mono);font-size:11px;line-height:1.8;color:var(--mut);padding:2px 0 2px 12px;position:relative;}
.dsec ul li::before{content:"◆";position:absolute;left:0;color:var(--acc);font-size:6px;top:6px;}
.dcard{background:var(--card);border:1px solid var(--bdr);border-radius:8px;padding:11px 13px;margin-bottom:8px;}
.dcard-t{font-family:var(--sans);font-weight:700;font-size:12px;color:var(--txt);margin-bottom:4px;}
.dcard p{font-family:var(--mono);font-size:10px;line-height:1.7;color:var(--mut);margin:0;}
.stbl{width:100%;border-collapse:collapse;margin:8px 0;}
.stbl th{background:var(--surf);color:var(--mut);font-family:var(--mono);font-size:9px;
  letter-spacing:1px;text-transform:uppercase;padding:6px 9px;text-align:left;border-bottom:1px solid var(--bdr);}
.stbl td{font-family:var(--mono);font-size:10px;color:var(--txt);padding:6px 9px;border-bottom:1px solid var(--bdr);}
.stbl tr:last-child td{border-bottom:none;}

/* FOOTER */
.ft{background:var(--surf);border-top:1px solid var(--bdr);padding:13px 22px;
  display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;gap:7px;}
.ftbrand{font-family:var(--sans);font-weight:700;font-size:12px;color:var(--txt);}
.ftbrand span{color:var(--acc);}
.ftlinks{display:flex;gap:13px;}
.ftlink{font-family:var(--mono);font-size:9px;color:var(--mut);cursor:pointer;transition:color .15s;}
.ftlink:hover{color:var(--acc);}
.ftcopy,.ftpriv{font-family:var(--mono);font-size:9px;color:var(--mut);}
.shiny-spinner-output-container{min-height:40px!important;}


/* ==========================================================
   MOBILE & TABLET RESPONSIVENESS
========================================================== */

/* TABLET (Under 960px) */
@media(max-width: 960px){
  .ml { grid-template-columns: 1fr; } /* Stack sidebar on top of results */
  .sb { max-height: none; border-right: none; border-bottom: 1px solid var(--bdr); padding-bottom: 20px;}
  .wrow { grid-template-columns: 1fr; } /* Stack scoring sliders */
  .ahdr { padding: 14px 20px; }
  .rp { padding: 16px 14px; }
}

/* MOBILE PHONES (Under 600px) */
@media(max-width: 600px){
  .tnav { padding: 0 10px; }
  .nbrand { font-size: 14px; }
  .nbrand .ver { display: none; } /* Hide version on small phones */
  .ahdr { padding: 12px 15px; }
  .htitle { font-size: 16px; }
  
  .srow { grid-template-columns: 1fr 1fr; } /* 2x2 grid for stats */
  .t3g { grid-template-columns: 1fr; } /* Stack Top 3 Cards */
  .pair-grid { grid-template-columns: 1fr; } /* Stack Dual Guide Pair */
  
  .rp { padding: 12px 10px; }
  
  /* Make Tabs smoothly scrollable horizontally */
  .tabnav { overflow-x: auto; flex-wrap: nowrap; -webkit-overflow-scrolling: touch; padding-bottom: 5px;}
  .tabbtn { flex-shrink: 0; padding: 8px 10px; font-size: 9px; }
  
  .ft { flex-direction: column; text-align: center; justify-content: center; }
}

/* VERY SMALL PHONES (Under 400px) */
@media(max-width: 400px){
  .srow { grid-template-columns: 1fr; } /* Stack stats single column */
}
'


# ==============================================================================
# ██  JAVASCRIPT
# ==============================================================================
APP_JS <- '
/* Page routing */
function showPage(n){
  ["home","docs","about"].forEach(function(p){
    var el=document.getElementById("page-"+p),nb=document.getElementById("nav-"+p);
    if(el)el.classList.toggle("ap",p===n);
    if(nb)nb.classList.toggle("an",p===n);
  });window.scrollTo(0,0);
}
/* Tab switching — visibility trick keeps plotly/DT initialized in DOM */
function switchTab(name){
  ["viz","gc","table","viewer","sim","ot"].forEach(function(t){
    var b=document.getElementById("tab_"+t),p=document.getElementById("panel_"+t);
    if(b)b.classList.toggle("active",t===name);
    if(p){
      if(t===name){p.style.visibility="visible";p.style.position="relative";
        p.style.height="auto";p.style.overflow="visible";}
      else{p.style.visibility="hidden";p.style.position="absolute";
        p.style.height="0";p.style.overflow="hidden";}
    }
  });
  if(name==="gc"||name==="viz")setTimeout(function(){window.dispatchEvent(new Event("resize"));},65);
}
/* Theme toggle */
var _dk=true;
function toggleTheme(){
  _dk=!_dk;
  document.body.classList.toggle("lm",!_dk);
  document.getElementById("ticon").textContent=_dk?"\u2600":"\uD83C\uDF19";
  document.getElementById("tlbl").textContent=_dk?"Light Mode":"Dark Mode";
  
  /* Optional: trigger resize so responsive charts re-render */
  setTimeout(function(){window.dispatchEvent(new Event("resize"));}, 50);
}
/* Copy single sequence to clipboard */
function copySeq(seq,bid){
  var done=function(){
    var b=document.getElementById(bid);if(!b)return;
    var o=b.textContent;b.textContent="\u2713 Copied!";b.classList.add("copied");
    setTimeout(function(){b.textContent=o;b.classList.remove("copied");},1800);
  };
  if(navigator.clipboard&&navigator.clipboard.writeText){
    navigator.clipboard.writeText(seq).then(done).catch(function(){
      var ta=document.createElement("textarea");ta.value=seq;
      document.body.appendChild(ta);ta.select();document.execCommand("copy");
      document.body.removeChild(ta);done();
    });
  }else{
    var ta=document.createElement("textarea");ta.value=seq;
    document.body.appendChild(ta);ta.select();document.execCommand("copy");
    document.body.removeChild(ta);done();
  }
}
/* Copy all gRNAs from hidden input */
function copyAll(){
  var hid=document.getElementById("all_seqs_h");
  if(hid)copyAllSeq(hid.value,"copyAllBtn");
}
function copyAllSeq(seq,bid){
  var done=function(){
    var b=document.getElementById(bid);if(!b)return;
    var o=b.textContent;b.textContent="\u2713 Copied!";
    setTimeout(function(){b.textContent=o;},1800);
  };
  var ta=document.createElement("textarea");ta.value=seq;
  document.body.appendChild(ta);ta.select();document.execCommand("copy");
  document.body.removeChild(ta);done();
}
'


# ==============================================================================
# ██  UI
# ==============================================================================
ui <- fluidPage(
  tags$head(
    tags$meta(name="viewport", content="width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no"),
    tags$link(href="https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=Syne:wght@400;600;800&display=swap",rel="stylesheet"),
    tags$style(HTML(APP_CSS)),
    tags$script(HTML(APP_JS))
  ),
  
  # ── TOP NAV ──────────────────────────────────────────────────────────────────
  div(class="tnav",
      div(class="nbrand",onclick="showPage('home')",
          "CRISPR ",tags$span(class="x","Guide"),"X",
          tags$span(class="ver","v3.1 · by Taiwo")),
      div(class="nlinks",
          tags$button(class="nb an",id="nav-home", onclick="showPage('home')", "🧬 Designer"),
          tags$button(class="nb",   id="nav-docs", onclick="showPage('docs')", "📖 Docs"),
          tags$button(class="nb",   id="nav-about",onclick="showPage('about')","👤 About"),
          div(class="ndiv"),
          tags$button(class="tbtn",id="themeToggle",onclick="toggleTheme()",
                      tags$span(id="ticon","☀"),tags$span(id="tlbl","Light Mode")))
  ),
  
  # ============================================================================
  # PAGE: HOME / DESIGNER
  # ============================================================================
  div(id="page-home",class="page ap",
      
      div(class="ahdr",
          div(class="hbdg","v3.1 · SpCas9 · SaCas9 · Cas12a · Cas12Max · Custom PAM · Simulation"),
          h1(class="htitle","CRISPR Guide",tags$span(class="tx","X"),
             " — Advanced gRNA Designer & Off-Target Predictor"),
          p(class="hsub","Production bioinformatics tool · Weighted seed scoring · Simulation mode · NCBI fetch · by Taiwo"),
          div(class="hmeta",
              div(class="mc","Systems: ",tags$b("SpCas9 · SaCas9 · Cas12a · Cas12Max")),
              div(class="mc","Guide: ",tags$b("15–24 nt")),
              div(class="mc","Off-target: ",tags$b("Both Strands")),
              div(class="mc","Hairpin: ",tags$b("Detected")))),
      
      div(class="ml",
          
          # ── SIDEBAR ──────────────────────────────────────────────────────────────
          div(class="sb",
              
              # Gene search
              div(
                div(class="slbl","🔍 Gene Search (NCBI)"),
                div(style="display:flex;gap:5px;align-items:flex-end;",
                    div(style="flex:1;",textInput("gene_name",NULL,placeholder="Gene, e.g. TP53")),
                    div(style="width:120px;",
                        selectInput("gene_species",NULL,
                                    choices=c("Homo sapiens","Mus musculus","Danio rerio","Arabidopsis thaliana")))),
                actionButton("fetch_btn","🌐 Fetch from NCBI",class="btn-ghost"),
                uiOutput("fetch_status")),
              
              # Sequence input
              div(
                div(class="slbl","01 / Sequence Input"),
                actionButton("sample_dna_btn", "🧬 Load Sample DNA", class="btn-a2", style="margin-bottom: 8px; width: 100%;"),
                textAreaInput("seq_input",NULL,"",rows=6,
                              placeholder="Paste DNA / FASTA here\nA T G C N accepted\nMulti-FASTA for batch mode"),
                div(style="display:flex;gap:5px;margin-top:5px;",
                    div(style="flex:1;",
                        selectInput("ex_sel",NULL,choices=c("Load example..."="-",names(EXAMPLE_SEQS)))),
                    actionButton("load_ex","Load",class="btn-ghost",style="width:55px;padding:7px 6px!"))),
              
              # File upload
              div(
                div(class="slbl","02 / Upload FASTA File"),
                fileInput("fasta_file",NULL,accept=c(".fa",".fasta",".txt"),placeholder=".fasta file"),
                uiOutput("upload_status")),
              
              # CRISPR system
              div(
                div(class="slbl","03 / CRISPR System & Guide"),
                selectInput("crispr_sys",NULL,
                            choices=c("SpCas9 (NGG)"="SpCas9","SaCas9 (NNGRRT)"="SaCas9",
                                      "Cas12a (TTTV)"="Cas12a","hfCas12Max (TN/TNN)"="Cas12Max",
                                      "Custom PAM"="Custom")),
                uiOutput("custom_pam_ui"),
                numericInput("guide_len","Guide length (nt):",value=20,min=15,max=24,step=1)),
              
              # Run / Clear
              div(
                actionButton("run_btn","⚡  Run Analysis",class="btn-run"),
                tags$br(),tags$br(),
                actionButton("clear_btn","✕  Clear All",class="btn-danger")),
              
              # Filters
              div(
                div(class="slbl","04 / Filters"),
                sliderInput("min_eff", "Min. Efficiency:",     min=0,max=1,  value=0,  step=0.05),
                sliderInput("max_risk","Max. Off-Target Risk:", min=0,max=1,  value=1,  step=0.05),
                sliderInput("gc_min",  "Min. GC%:",             min=0,max=100,value=0,  step=5),
                sliderInput("gc_max",  "Max. GC%:",             min=0,max=100,value=100,step=5),
                selectInput("strand_filter","Strand:",
                            choices=c("Both"="both","Forward (+)"="+","Reverse (−)"="-"))),
              
              # Scoring weights
              div(
                div(class="slbl","05 / Composite Scoring"),
                div(class="wrow",
                    sliderInput("w_eff","Efficiency weight:",min=0,max=1,value=0.6,step=0.1),
                    sliderInput("w_risk","Risk weight:",    min=0,max=1,value=0.4,step=0.1)),
                selectInput("rank_mode","Ranking preset:",
                            choices=c("Balanced"="balanced","High Efficiency"="hi_eff","Low Risk"="lo_risk"))),
              
              # Export
              div(
                div(class="slbl","06 / Export"),
                div(class="egrp",
                    downloadButton("dl_csv",  "⬇  CSV — All Results"),
                    downloadButton("dl_fasta","⬇  FASTA — All gRNAs"),
                    downloadButton("dl_top",  "⬇  CSV — Top Candidates"),
                    if(has_xlsx) downloadButton("dl_excel","⬇  Excel (.xlsx)") else NULL,
                    tags$button(class="btn-a2",id="copyAllBtn",onclick="copyAll()",
                                "📋 Copy All gRNAs to Clipboard"))),
              
              # Hidden field holds all sequences for copy-all
              tags$input(type="hidden",id="all_seqs_h",value=""),
              
              # Quick guide card
              div(class="ic",
                  tags$strong("How to use:"),tags$br(),
                  "① Search gene (NCBI) OR paste sequence",tags$br(),
                  "② Select system, guide length",tags$br(),
                  "③ Click ⚡ Run Analysis",tags$br(),
                  "④ Review Top 3 candidates + Sequence Viewer",tags$br(),
                  "⑤ Simulation tab → cut site + indel prediction",tags$br(),
                  "⑥ Off-targets tab → detail table per gRNA",tags$br(),
                  "⑦ Adjust weights → scores update live")
          ),
          
          # ── RESULTS PANEL ────────────────────────────────────────────────────────
          div(class="rp",
              uiOutput("validation_msg"),
              uiOutput("summary_stats"),
              uiOutput("score_info"),
              uiOutput("top3_ui"),
              uiOutput("results_tabs"))
      ),
      
      div(class="ft",
          div(class="ftbrand","CRISPR ",tags$span("Guide"),"X"),
          div(class="ftlinks",
              tags$span(class="ftlink",onclick="showPage('docs')", "Docs"),
              tags$span(class="ftlink",onclick="showPage('about')","About")),
          div(class="ftpriv","🔒 No data collected or stored"),
          div(class="ftcopy",paste0("© ",format(Sys.Date(),"%Y"),
                                    " CRISPR GuideX v3.1 · by Taiwo · Educational use only")))
  ),
  
  # ============================================================================
  # PAGE: DOCS
  # ============================================================================
  div(id="page-docs",class="page",
      div(class="pdocs",
          div(class="dhero",
              h1("How to Use ",tags$span("CRISPR GuideX v3.1")),
              p("Complete reference covering multi-system PAM, adjustable guide length, hairpin detection, simulation, NCBI integration, and dynamic scoring.")),
          
          div(class="dsec",
              div(class="dsec-t","🚀 Getting Started"),
              div(class="hiw",
                  div(class="hstep",div(class="hn","STEP 01"),div(class="ht","Input"),
                      div(class="hd","Search gene (NCBI) OR paste DNA/FASTA in the sequence box")),
                  div(class="hstep",div(class="hn","STEP 02"),div(class="ht","Configure"),
                      div(class="hd","Select CRISPR system, guide length, custom PAM if needed")),
                  div(class="hstep",div(class="hn","STEP 03"),div(class="ht","Run"),
                      div(class="hd","Click ⚡ Run — both strands scanned, all gRNAs scored")),
                  div(class="hstep",div(class="hn","STEP 04"),div(class="ht","Analyse"),
                      div(class="hd","Review Top 3, Sequence Viewer, Simulation, Off-target table")),
                  div(class="hstep",div(class="hn","STEP 05"),div(class="ht","Tune"),
                      div(class="hd","Adjust efficiency/risk weights → composite scores update live")),
                  div(class="hstep",div(class="hn","STEP 06"),div(class="ht","Export"),
                      div(class="hd","Download CSV / FASTA / Excel or copy sequences to clipboard")))),
          
          div(class="dsec",
              div(class="dsec-t","⚙ CRISPR Systems"),
              tags$table(class="stbl",
                         tags$thead(tags$tr(tags$th("System"),tags$th("PAM"),tags$th("Notes"))),
                         tags$tbody(
                           tags$tr(tags$td("SpCas9"), tags$td("5'-NGG-3'"),    tags$td("Most widely used; 3' PAM; blunt cut pos 17|18")),
                           tags$tr(tags$td("SaCas9"), tags$td("5'-NNGRRT-3'"), tags$td("Smaller enzyme; useful for AAV delivery; 3' PAM")),
                           tags$tr(tags$td("Cas12a"), tags$td("5'-TTTV-3'"),   tags$td("5' PAM; T-rich targets; staggered cut ~18 bp downstream")),
                           tags$tr(tags$td("Cas12Max"), tags$td("5'-TN/TNN-3'"), tags$td("5' PAM; Relaxed target rules; high-fidelity variant")),
                           tags$tr(tags$td("Custom"), tags$td("User regex"),    tags$td("Enter any PAM as a regex pattern, e.g. ^[ACGT]CC$"))))),
          
          div(class="dsec",
              div(class="dsec-t","⚡ Efficiency Scoring"),
              tags$table(class="stbl",
                         tags$thead(tags$tr(tags$th("Factor"),tags$th("Rule"),tags$th("Effect"))),
                         tags$tbody(
                           tags$tr(tags$td("GC Content"),    tags$td("Ideal 40–60%"),          tags$td("Base score 0–1.0")),
                           tags$tr(tags$td("Poly-T (TTTT)"), tags$td("U6 promoter disruption"),tags$td("−0.20")),
                           tags$tr(tags$td("Homopolymer≥5"), tags$td("AAAAA/GGGGG etc."),     tags$td("−0.15")),
                           tags$tr(tags$td("G at last pos"), tags$td("Seed-PAM junction"),     tags$td("+0.05")),
                           tags$tr(tags$td("Seed GC"),       tags$td("Last ~9 nt; 30–70%"),    tags$td("+0.05"))))),
          
          div(class="dsec",
              div(class="dsec-t","🎯 Off-Target Prediction (Weighted Mismatch)"),
              tags$ul(
                tags$li("Seed region (last 12 nt before PAM): mismatch weight 2.0"),
                tags$li("Non-seed region (positions 1–~8): mismatch weight 1.0"),
                tags$li("Sites with total weighted score ≤ 6 are flagged as potential off-targets"),
                tags$li("Risk score = 1 − exp(−count × 0.5); range 0 (safe) to 1 (high risk)"),
                tags$li("Off-targets tab: scans BOTH forward and reverse strands of the target DNA!"))),
          
          div(class="dsec",
              div(class="dsec-t","🧬 Hairpin / Secondary Structure Detection"),
              p("GuideX checks every 4-mer in the guide against the reverse complement of the full guide. If any 4-mer appears in the reverse complement, a potential hairpin is flagged. These guides are marked with ⚠ Hairpin in Top-3 cards, the results table, and the GC scatter plot.")),
          
          div(class="dsec",
              div(class="dsec-t","🔬 CRISPR Simulation Mode"),
              p("Select any gRNA row in the Table tab to activate the Simulation panel:"),
              tags$ul(
                tags$li("Predicted cut position (SpCas9/SaCas9: between nt 17–18; Cas12a: ~18 bp downstream)"),
                tags$li("Original sequence with gRNA and PAM color-coded"),
                tags$li("30 bp upstream / downstream of cut site"),
                tags$li("Example NHEJ indel outcome (upstream | [INDEL] | downstream)"))),
          
          div(class="dsec",
              div(class="dsec-t","📊 Dynamic Composite Scoring"),
              p("Composite = Efficiency × w_eff + (1 − Risk) × w_risk"),
              tags$ul(
                tags$li("Balanced: 0.6 / 0.4 (default)"),
                tags$li("High Efficiency preset: 0.9 / 0.1"),
                tags$li("Low Risk preset: 0.2 / 0.8"),
                tags$li("Custom: use the sliders directly in the sidebar")))),
      
      div(class="ft",
          div(class="ftbrand","CRISPR ",tags$span("Guide"),"X"),
          div(class="ftlinks",
              tags$span(class="ftlink",onclick="showPage('home')","← Designer"),
              tags$span(class="ftlink",onclick="showPage('about')","About")),
          div(class="ftcopy",paste0("© ",format(Sys.Date(),"%Y")," CRISPR GuideX v3.1 · by Taiwo")))
  ),
  
  # ============================================================================
  # PAGE: ABOUT
  # ============================================================================
  div(id="page-about",class="page",
      div(class="pabout",
          div(class="ahero",
              h1("About ",tags$span("CRISPR GuideX")),
              p("A free, open-source, R-native CRISPR gRNA design platform for researchers, educators, and students. All computation is local — no data leaves your machine.")),
          div(class="dsec",
              div(class="dsec-t","📋 What's New in v3.1"),
              div(class="dcard",div(class="dcard-t","Multi-system PAM support"),
                  p("SpCas9 (NGG), SaCas9 (NNGRRT), Cas12a (TTTV), hfCas12Max (TN/TNN), and fully custom PAM patterns with adjustable guide length 15–24 nt.")),
              div(class="dcard",div(class="dcard-t","Weighted seed-region off-target scoring"),
                  p("Positions 12–20 (seed) carry 2× mismatch penalty, reflecting biological importance of seed region specificity.")),
              div(class="dcard",div(class="dcard-t","Hairpin / secondary structure detection"),
                  p("Self-complementarity 4-mer heuristic flags guides likely to form hairpins, reducing cleavage efficiency.")),
              div(class="dcard",div(class="dcard-t","CRISPR simulation mode"),
                  p("Predicts cut site and shows example NHEJ indel outcome for any selected guide. Supports SpCas9 and Cas12a cut geometries.")),
              div(class="dcard",div(class="dcard-t","NCBI gene search"),
                  p("Fetch real genomic sequence by gene name and species using NCBI eUtils. Improved coordinate retrieval logic.")),
              div(class="dcard",div(class="dcard-t","Dynamic composite scoring"),
                  p("Adjust efficiency and risk weights in real time with sliders; choose from Balanced / High Efficiency / Low Risk presets."))),
          div(class="dsec",
              div(class="dsec-t","🔒 Privacy Policy"),
              div(style="background:rgba(34,197,94,.06);border:1px solid rgba(34,197,94,.2);border-radius:8px;padding:11px 15px;font-family:var(--mono);font-size:11px;color:var(--ok);line-height:1.7;",
                  "🛡 CRISPR GuideX does not collect, transmit, or store any data. Your sequences never leave your device."),
              tags$br(),
              p("No cookies. No analytics. No user accounts. NCBI fetch is the only optional network call — it only transmits the gene name you enter.")),
          div(class="dsec",
              div(class="dsec-t","⚠ Disclaimer"),
              p("For educational and preliminary research use only. Rule-based scores do not replace validated genome-wide off-target analysis tools (CRISPOR, Cas-OFFinder) or wet-lab experimental validation."))),
      div(class="ft",
          div(class="ftbrand","CRISPR ",tags$span("Guide"),"X"),
          div(class="ftlinks",
              tags$span(class="ftlink",onclick="showPage('home')","← Designer"),
              tags$span(class="ftlink",onclick="showPage('docs')","Docs")),
          div(class="ftcopy",paste0("© ",format(Sys.Date(),"%Y")," CRISPR GuideX v3.1 · by Taiwo · Educational use only")))
  )
)


# ==============================================================================
# ██  SERVER
# ==============================================================================
server <- function(input, output, session) {
  
  # ── Reactive state ──────────────────────────────────────────────────────────
  rv <- reactiveValues(
    grna_df      = NULL,          # main results data.frame
    seq_len      = 0L,            # total bp processed
    batch_names  = character(0L), # names of input sequences
    all_dna_seqs = list(),        # Memory dictionary of all mapped sequences
    system_info  = NULL           # system info list for current analysis
  )
  
  # ── Custom PAM UI (shown only when Custom is selected) ─────────────────────
  output$custom_pam_ui <- renderUI({
    req(input$crispr_sys == "Custom")
    div(
      textInput("custom_pam_pat", "PAM regex:", placeholder = "e.g. ^[ACGT]CC$"),
      selectInput("custom_pam_side", "PAM side:",
                  choices = c("3' of gRNA" = "3prime", "5' of gRNA" = "5prime")))
  })
  
  # ── Load Sample DNA ─────────────────────────────────────────────────────────
  observeEvent(input$sample_dna_btn, {
    sample_seq <- ">Sample_Sequence_HBB_Exon1\nACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCT"
    updateTextAreaInput(session, "seq_input", value = sample_seq)
  })
  
  # ── NCBI Gene Fetch ─────────────────────────────────────────────────────────
  output$fetch_status <- renderUI(NULL)
  
  observeEvent(input$fetch_btn, {
    gene <- trimws(as.character(input$gene_name))
    if (nchar(gene) == 0L) return()
    output$fetch_status <- renderUI(div(class="ab ainf","⏳ Fetching from NCBI..."))
    res <- fetch_ncbi_sequence(gene, as.character(input$gene_species))
    if (isTRUE(res$ok)) {
      updateTextAreaInput(session, "seq_input", 
                          value = paste0(">", res$gene, "\n", res$seq))
      output$fetch_status <- renderUI(div(class="ab asuc", paste0("✓ ", res$msg)))
    } else {
      output$fetch_status <- renderUI(div(class="ab aerr", paste0("⚠ ", res$msg)))
    }
  })
  
  # ── Load example ────────────────────────────────────────────────────────────
  observeEvent(input$load_ex, {
    sel <- as.character(input$ex_sel)
    if (sel != "-" && sel %in% names(EXAMPLE_SEQS))
      updateTextAreaInput(session, "seq_input", value = EXAMPLE_SEQS[[sel]])
  })
  
  # ── Upload status ───────────────────────────────────────────────────────────
  output$upload_status <- renderUI({
    req(input$fasta_file)
    seqs <- parse_fasta_file(input$fasta_file$datapath)
    n    <- length(seqs)
    if (n == 0L) return(div(class="ab awrn","⚠ No valid sequences in file."))
    div(class="ab ainf",
        sprintf("✓ %d sequence(s) · %s bp", n,
                format(sum(nchar(unlist(seqs))), big.mark = ",")))
  })
  
  # ── Clear all ───────────────────────────────────────────────────────────────
  observeEvent(input$clear_btn, {
    updateTextAreaInput(session, "seq_input", value = "")
    rv$grna_df      <- NULL
    rv$seq_len      <- 0L
    rv$all_dna_seqs <- list()
    rv$batch_names  <- character(0L)
    rv$system_info  <- NULL
  })
  
  # ── RUN ANALYSIS ─────────────────────────────────────────────────────────────
  observeEvent(input$run_btn, {
    sys    <- as.character(input$crispr_sys)
    glen   <- as.integer(input$guide_len)
    c_pam  <- if (sys == "Custom") input$custom_pam_pat  else NULL
    c_side <- if (sys == "Custom") input$custom_pam_side else "3prime"
    
    # ── Collect all input sequences ────────────────────────────────────────────
    all_seqs <- list()
    
    # File upload
    if (!is.null(input$fasta_file))
      all_seqs <- c(all_seqs, parse_fasta_file(input$fasta_file$datapath))
    
    # Text box (plain or multi-FASTA)
    txt <- trimws(as.character(input$seq_input))
    if (nchar(txt) > 0L) {
      if (grepl("^>", txt)) {
        batch    <- parse_multi_fasta(txt)
        all_seqs <- c(all_seqs, batch)
      } else {
        s <- parse_sequence(txt)
        if (!is.null(s)) all_seqs[["Input"]] <- s
      }
    }
    
    # Validate
    if (length(all_seqs) == 0L) {
      rv$seq_len <- -1L; rv$grna_df <- NULL; return()
    }
    min_len <- glen + 3L
    valid   <- Filter(function(s) nchar(s) >= min_len && grepl("^[ATGCN]+$", s), all_seqs)
    if (length(valid) == 0L) {
      rv$seq_len <- -2L; rv$grna_df <- NULL; return()
    }
    
    rv$batch_names  <- names(valid)
    rv$seq_len      <- as.integer(sum(nchar(unlist(valid))))
    rv$all_dna_seqs <- valid
    rv$system_info  <- get_system_info(sys, c_pam, c_side, glen)
    
    # ── Scan sequences ─────────────────────────────────────────────────────────
    withProgress(message = "Scanning sequences...", value = 0, {
      all_res <- list()
      nv      <- length(valid)
      for (nm in names(valid)) {
        incProgress(1 / nv, detail = paste("Processing", nm))
        df <- find_grna(valid[[nm]], sys, glen, c_pam, c_side)
        if (!is.null(df) && nrow(df) > 0L) {
          df$Sequence_Source <- nm
          all_res[[nm]]      <- df
        }
      }
    })
    
    if (length(all_res) == 0L) { rv$grna_df <- data.frame(); return() }
    combined        <- do.call(rbind, all_res)
    rownames(combined) <- NULL
    rv$grna_df      <- combined
  })
  
  # ── Effective composite weights (respect ranking preset) ────────────────────
  eff_weights <- reactive({
    mode <- as.character(input$rank_mode)
    if (mode == "hi_eff")  return(list(we = 0.9, wr = 0.1))
    if (mode == "lo_risk") return(list(we = 0.2, wr = 0.8))
    we  <- max(0.01, as.double(input$w_eff))
    wr  <- max(0.01, as.double(input$w_risk))
    tot <- we + wr
    list(we = we / tot, wr = wr / tot)
  })
  
  # ── Filtered + scored data ───────────────────────────────────────────────────
  filtered_data <- reactive({
    df  <- rv$grna_df
    req(!is.null(df) && nrow(df) > 0L)
    wts <- eff_weights()
    
    df  <- df %>%
      filter(Efficiency      >= input$min_eff,
             Off_Target_Risk <= input$max_risk,
             GC_Content      >= input$gc_min,
             GC_Content      <= input$gc_max)
    
    if (as.character(input$strand_filter) != "both")
      df <- df %>% filter(Strand == as.character(input$strand_filter))
    
    df <- df %>% mutate(
      Risk_Category = risk_label(Efficiency, Off_Target_Risk),
      Composite     = composite_score(Efficiency, Off_Target_Risk, wts$we, wts$wr)
    )
    
    # DYNAMIC SORTING LOGIC 
    s_col <- as.character(input$sort_col)
    s_dir <- as.character(input$sort_dir)
    
    if (length(s_col) == 0L || s_col == "") s_col <- "Composite"
    if (length(s_dir) == 0L || s_dir == "") s_dir <- "desc"
    
    if (s_col %in% names(df)) {
      if (s_dir == "desc") {
        df <- df %>% arrange(desc(.data[[s_col]]))
      } else {
        df <- df %>% arrange(.data[[s_col]])
      }
    }
    
    df
  })
  
  # ── Validation messages ──────────────────────────────────────────────────────
  output$validation_msg <- renderUI({
    sl <- rv$seq_len; df <- rv$grna_df
    if (isTRUE(sl == -1L))
      return(div(class="ab aerr","⚠ No valid input found. Paste a DNA sequence or upload a FASTA file."))
    if (isTRUE(sl == -2L))
      return(div(class="ab awrn","⚠ Sequences too short. Need at least guide_length + PAM nucleotides."))
    if (!is.null(df) && nrow(df) == 0L)
      return(div(class="ab awrn",
                 paste0("⚠ No PAM sites found for ", input$crispr_sys,
                        ". Try a different CRISPR system or check your sequence.")))
    NULL
  })
  
  # ── Summary stat cards ───────────────────────────────────────────────────────
  output$summary_stats <- renderUI({
    df <- filtered_data(); sl <- rv$seq_len
    req(!is.null(df) && isTRUE(sl > 0L))
    total   <- nrow(df)
    hi_eff  <- sum(df$Efficiency >= 0.7)
    lo_risk <- sum(df$Off_Target_Risk <= 0.3)
    n_hp    <- sum(as.logical(df$Hairpin_Risk), na.rm = TRUE)
    nb      <- length(rv$batch_names)
    badge   <- if (nb > 1L) paste0(nb, " seqs") else paste0(format(sl, big.mark=","), " bp")
    div(class="srow",
        div(class="sc cg", div(class="stat-val",total),   div(class="stat-lbl","gRNAs Found"),    div(class="stat-bdg bg",badge)),
        div(class="sc cp", div(class="stat-val",hi_eff),  div(class="stat-lbl","High Efficiency"),div(class="stat-bdg bp","≥ 0.70")),
        div(class="sc cg", div(class="stat-val",lo_risk), div(class="stat-lbl","Low Off-Target"),  div(class="stat-bdg bg","≤ 0.30")),
        div(class="sc cb", div(class="stat-val",n_hp),    div(class="stat-lbl","Hairpin Flagged"), div(class="stat-bdg bb","⚠")))
  })
  
  # ── Composite score info bar ─────────────────────────────────────────────────
  output$score_info <- renderUI({
    req(!is.null(rv$grna_df) && isTRUE(nrow(rv$grna_df) > 0L))
    wts <- eff_weights()
    div(class="ab ainf",style="padding:7px 11px;",
        tags$span(style="font-family:var(--mono);font-size:10px;",
                  sprintf("Composite = Efficiency × %.1f  +  (1 − Risk) × %.1f  |  Mode: %s",
                          wts$we, wts$wr, toupper(as.character(input$rank_mode)))))
  })
  
  # ── TOP 3 CANDIDATE CARDS + DUAL-GUIDE PAIR ──────────────────────────────────
  output$top3_ui <- renderUI({
    df <- filtered_data()
    req(!is.null(df) && nrow(df) > 0L)
    
    # BEST PAIR LOGIC
    pair_ui <- NULL
    sources <- unique(df$Sequence_Source)
    if (length(sources) >= 2) {
      s1 <- sources[1]
      s2 <- sources[2] 
      g1 <- df %>% filter(Sequence_Source == s1) %>% arrange(desc(Composite)) %>% slice_head(n=1)
      g2 <- df %>% filter(Sequence_Source == s2) %>% arrange(desc(Composite)) %>% slice_head(n=1)
      
      if (nrow(g1) == 1 && nrow(g2) == 1) {
        g1_seq <- sc(g1$gRNA_Sequence); g1_pam <- sc(g1$PAM); bid1 <- "btn_pair_1"
        g2_seq <- sc(g2$gRNA_Sequence); g2_pam <- sc(g2$PAM); bid2 <- "btn_pair_2"
        
        pair_ui <- div(class="pair-card",
                       div(class="pair-title", "🧬 Dual-Guide Deletion Architecture (Best Pair)"),
                       div(class="pair-grid",
                           div(class="pair-item",
                               div(class="pair-lbl", paste("Guide 1 (Upstream):", s1)),
                               div(class="cseq", HTML(colorize_seq(g1_seq, g1_pam))),
                               div(class="cmeta", tags$span(sprintf("Score: %.3f | Pos: %s | Strand: %s", g1$Composite, g1$Position, g1$Strand))),
                               tags$button(class="cbtn", id=bid1, onclick=sprintf("copySeq('%s','%s')", g1_seq, bid1), "📋 Copy")
                           ),
                           div(class="pair-item",
                               div(class="pair-lbl", paste("Guide 2 (Downstream):", s2)),
                               div(class="cseq", HTML(colorize_seq(g2_seq, g2_pam))),
                               div(class="cmeta", tags$span(sprintf("Score: %.3f | Pos: %s | Strand: %s", g2$Composite, g2$Position, g2$Strand))),
                               tags$button(class="cbtn", id=bid2, onclick=sprintf("copySeq('%s','%s')", g2_seq, bid2), "📋 Copy")
                           )
                       )
        )
      }
    }
    
    # TOP 3 LOGIC 
    top3 <- df %>% arrange(desc(Composite)) %>% slice_head(n = 3L)
    n3   <- nrow(top3)
    
    ICONS  <- c("1"="🥇","2"="🥈","3"="🥉")
    LABELS <- c("1"="Top Candidate","2"="2nd Best","3"="3rd Best")
    CLS    <- c("1"="r1","2"="r2","3"="r3")
    
    all_str <- paste(top3$gRNA_Sequence, collapse = "\n")
    
    cards <- vector("list", n3)
    for (i in seq_len(n3)) {
      ic    <- as.character(i)
      pos   <- as.character(top3[["Position"]][[i]])
      str_v <- as.character(top3[["Strand"]][[i]])
      gs    <- as.character(top3[["gRNA_Sequence"]][[i]])
      ps    <- as.character(top3[["PAM"]][[i]])
      gc_s  <- as.character(top3[["GC_Content"]][[i]])
      eff_s <- as.character(top3[["Efficiency"]][[i]])
      otr   <- as.double(top3[["Off_Target_Risk"]][[i]])
      otn   <- as.character(top3[["Off_Targets"]][[i]])
      rc    <- as.character(top3[["Risk_Category"]][[i]])
      comp  <- as.double(top3[["Composite"]][[i]])
      glen_s<- as.character(top3[["Guide_Length"]][[i]])
      hp    <- isTRUE(as.logical(top3[["Hairpin_Risk"]][[i]]))
      
      pill_cls <- if (isTRUE(otr > 0.6)) "cpill rhi"
      else if (isTRUE(otr > 0.3)) "cpill rmd"
      else "cpill"
      
      bid      <- paste0("cb_", i)
      seq_html <- colorize_seq(gs, ps)
      
      cards[[i]] <- div(class = paste("cc", CLS[[ic]]),
                        div(class="crnk",
                            tags$span(ICONS[[ic]]),
                            LABELS[[ic]],
                            tags$span(style="margin-left:auto;", sprintf("Score:%.3f", comp))),
                        div(class="cseq", HTML(seq_html)),
                        div(class="cmeta",
                            tags$span(paste0("Pos:", pos)),
                            tags$span(paste0(" | ", str_v)),
                            tags$span(paste0(" | GC:", gc_s, "%")),
                            tags$span(paste0(" | Eff:", eff_s)),
                            tags$span(paste0(" | OT:", otn)),
                            tags$span(paste0(" | ", glen_s, "nt"))),
                        div(style="display:flex;align-items:center;justify-content:space-between;gap:5px;margin-top:5px;",
                            div(style="display:flex;gap:4px;flex-wrap:wrap;",
                                tags$span(class=pill_cls,  paste0("Risk:", sprintf("%.3f", otr))),
                                tags$span(class="cpill",   rc),
                                if (hp) tags$span(class="cpill hp","⚠ Hairpin") else NULL),
                            tags$button(class="cbtn",id=bid,
                                        onclick=sprintf("copySeq('%s','%s')", gs, bid),
                                        "📋 Copy")))
    }
    
    tagList(
      pair_ui,
      div(style="font-family:var(--mono);font-size:9px;letter-spacing:2px;text-transform:uppercase;color:var(--acc);margin-bottom:7px;",
          "★  Top Candidates — Ranked by Composite Score"),
      div(class="t3g", tagList(cards)),
      tags$script(HTML(sprintf(
        'document.getElementById("all_seqs_h").value="%s";',
        gsub('"', '\\\\"', gsub("\n", "\\\\n", all_str)))))
    )
  })
  
  # ── RESULTS TABS ─────────────────────────────────────────────────────────────
  output$results_tabs <- renderUI({
    df  <- rv$grna_df; sl <- rv$seq_len
    if (is.null(df) || isTRUE(sl <= 0L))
      return(div(class="estate",
                 div(class="eicon","🧬"),h3("No Analysis Yet"),
                 p("Paste a DNA sequence and click Run Analysis to begin.")))
    if (nrow(df) == 0L) return(NULL)
    
    tagList(
      div(class="tabnav",
          tags$button(class="tabbtn active",id="tab_viz",   onclick="switchTab('viz')",   "📊 Map"),
          tags$button(class="tabbtn",       id="tab_gc",    onclick="switchTab('gc')",    "📈 GC/Eff"),
          tags$button(class="tabbtn",       id="tab_table", onclick="switchTab('table')", "📋 Table"),
          tags$button(class="tabbtn",       id="tab_viewer",onclick="switchTab('viewer')","🔬 Viewer"),
          tags$button(class="tabbtn",       id="tab_sim",   onclick="switchTab('sim')",   "⚗ Simulate"),
          tags$button(class="tabbtn",       id="tab_ot",    onclick="switchTab('ot')",    "🎯 Off-targets")),
      
      div(style="position:relative;",
          div(id="panel_viz", style="position:relative;",
              div(class="chard",
                  div(class="chtitle","gRNA Position Map — Both Strands"),
                  withSpinner(plotlyOutput("genome_plot",height="285px"),type=6,color="#00d4aa",size=0.5))),
          
          div(id="panel_gc",
              style="position:absolute;top:0;left:0;right:0;visibility:hidden;height:0;overflow:hidden;",
              div(class="chard",
                  div(class="chtitle","GC Content vs. Efficiency (shape = hairpin risk)"),
                  withSpinner(plotlyOutput("gc_plot",height="285px"),type=6,color="#00d4aa",size=0.5))),
          
          div(id="panel_table",
              style="position:absolute;top:0;left:0;right:0;visibility:hidden;height:0;overflow:hidden;",
              div(class="chard",
                  div(class="chtitle","All Candidate gRNAs — Click a row to activate Viewer / Simulation / Off-targets"),
                  div(style="display:flex;gap:7px;margin-bottom:8px;flex-wrap:wrap;",
                      div(style="width:170px;",
                          selectInput("sort_col","Sort by:",
                                      choices=c("Composite","Efficiency","Off_Target_Risk","GC_Content","Position"),
                                      selected="Composite")),
                      div(style="width:130px;",
                          selectInput("sort_dir","Order:",
                                      choices=c("Descending"="desc","Ascending"="asc"),
                                      selected="desc"))),
                  withSpinner(DTOutput("grna_table"),type=6,color="#00d4aa",size=0.5))),
          
          div(id="panel_viewer",
              style="position:absolute;top:0;left:0;right:0;visibility:hidden;height:0;overflow:hidden;",
              div(class="chard",
                  div(class="chtitle","Sequence Viewer — select a row in the Table tab"),
                  uiOutput("seq_viewer_ui"))),
          
          div(id="panel_sim",
              style="position:absolute;top:0;left:0;right:0;visibility:hidden;height:0;overflow:hidden;",
              div(class="chard",
                  div(class="chtitle","CRISPR Simulation — Cut Site & Indel Prediction"),
                  uiOutput("sim_ui"))),
          
          div(id="panel_ot",
              style="position:absolute;top:0;left:0;right:0;visibility:hidden;height:0;overflow:hidden;",
              div(class="chard",
                  div(class="chtitle","Off-Target Detail — select a row in the Table tab"),
                  uiOutput("ot_detail_ui"))))
    )
  })
  
  # ── PLOT OUTPUTS ─────────────────────────────────────────────────────────────
  output$genome_plot <- renderPlotly({
    df <- filtered_data()
    req(!is.null(df) && nrow(df) > 0L)
    set.seed(42L)
    df <- df %>% mutate(
      CG = dplyr::case_when(
        Efficiency >= 0.7 & Off_Target_Risk <= 0.3 ~ "High Eff / Low Risk",
        Efficiency >= 0.4 & Off_Target_Risk <= 0.6 ~ "Medium",
        TRUE ~ "Low Eff / High Risk"),
      yp = ifelse(Strand == "+", 1.0, -1.0) + runif(dplyr::n(), -0.13, 0.13),
      TT = paste0("<b>gRNA @ pos ", Position, " (", Strand, " strand)</b><br>",
                  "Seq: <b>", gRNA_Sequence, "</b>"))
    
    cmap <- c("High Eff / Low Risk"="#22c55e","Medium"="#eab308","Low Eff / High Risk"="#ef4444")
    
    p <- ggplot(df, aes(x=Position,y=yp,color=CG,size=Efficiency,text=TT)) +
      geom_hline(yintercept=0, color="#94a3b8", alpha=0.5) +
      geom_point(alpha=0.88) +
      scale_color_manual(values=cmap, name="") +
      theme_minimal() +
      theme(panel.grid.minor=element_blank(), axis.text=element_text(color="#94a3b8",family="monospace"))
    
    ggplotly(p, tooltip="text") %>%
      layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")
  })
  
  output$gc_plot <- renderPlotly({
    df <- filtered_data()
    req(!is.null(df) && nrow(df) > 0L)
    p <- ggplot(df, aes(x=GC_Content, y=Efficiency, color=Off_Target_Risk, shape=as.factor(Hairpin_Risk))) +
      geom_point(size=3) +
      scale_color_gradient2(low="#22c55e",mid="#eab308",high="#ef4444", midpoint=0.5) +
      theme_minimal() +
      theme(axis.text=element_text(color="#94a3b8",family="monospace"))
    ggplotly(p) %>% layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")
  })
  
  # ── TABLE & DETAILS ──────────────────────────────────────────────────────────
  output$grna_table <- renderDT({
    df  <- filtered_data()
    req(!is.null(df) && nrow(df) > 0L)
    datatable(df, options = list(pageLength=12, scrollX=TRUE, dom='frtip'), rownames=FALSE, selection="single")
  })
  
  selected_row <- reactive({
    df  <- filtered_data(); req(!is.null(df))
    sel <- input$grna_table_rows_selected
    if (is.null(sel) || length(sel) == 0L) return(NULL)
    df[as.integer(sel[[1L]]), , drop=FALSE]
  })
  
  output$seq_viewer_ui <- renderUI({
    r <- selected_row(); if (is.null(r)) return(div(class="ab ainf","← Select a row"))
    src <- as.character(r[["Sequence_Source"]][[1L]])
    dna <- rv$all_dna_seqs[[src]]
    sys <- rv$system_info
    
    # CRITICAL FIX: Pass dynamic pam_len, pam_side, and strand to viewer
    html <- build_seq_viewer(dna, si_(r["Position"]), si_(r["Guide_Length"]), 
                             sys$pam_len, sys$side, as.character(r[["Strand"]][[1L]]), 
                             window = 80L)
    div(class="seq-vbox", HTML(html))
  })
  
  output$sim_ui <- renderUI({
    r <- selected_row(); if (is.null(r)) return(div(class="ab ainf","← Select a row"))
    src <- as.character(r[["Sequence_Source"]][[1L]])
    dna <- rv$all_dna_seqs[[src]]
    sys <- rv$system_info
    
    # CRITICAL FIX: Pass dynamic pam_side and strand to simulator geometry
    sim <- simulate_cut(dna, si_(r["Position"]), si_(r["Guide_Length"]), 
                        sys$side, as.character(r[["Strand"]][[1L]]))
    
    div(class="sim-box", div(class="sim-lbl","Cut Site Prediction"), 
        div(class="sim-seq", 
            HTML(colorize_seq(sim$upstream, NULL)), 
            tags$span(class="cut","▼"), 
            HTML(colorize_seq(sim$downstream, NULL))))
  })
  
  output$ot_detail_ui <- renderUI({
    r <- selected_row(); if (is.null(r)) return(div(class="ab ainf","← Select a row"))
    src <- as.character(r[["Sequence_Source"]][[1L]])
    dna <- rv$all_dna_seqs[[src]]
    
    ot <- off_target_full(sc(r["gRNA_Sequence"]), dna, si_(r["Position"]))
    if (nrow(ot$hits) == 0L) return(div(class="ab asuc","✓ No off-targets detected"))
    datatable(ot$hits, options=list(dom='tip'), rownames=FALSE)
  })
  
  # ── EXPORT HANDLERS ──────────────────────────────────────────────────────────
  output$dl_csv <- downloadHandler(filename = "crispr_guidex.csv", content = function(f) write.csv(filtered_data(), f))
  output$dl_fasta <- downloadHandler(filename = "crispr_guidex.fasta", content = function(f) writeLines(df_to_fasta(filtered_data()), f))
  output$dl_top <- downloadHandler(filename = "crispr_guidex_top.csv", content = function(f) write.csv(head(filtered_data(), 10), f))
}

# ==============================================================================
# ██  LAUNCH
# ==============================================================================
shinyApp(ui = ui, server = server)