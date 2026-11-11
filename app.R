library(shiny)
library(plotly)
library(DT)
source("analysis.R", local = TRUE)
source("review.R", local = TRUE)

cache <- Sys.getenv("EXPLORER_CACHE", unset = "explorer-data.rds")
data <- readRDS(cache)
meta <- data$meta
meta$id <- seq_len(nrow(meta))
meta$barcode <- rownames(meta)
root <- Sys.getenv("EXPLORER_SOURCE", "/Users/paul/WorkInProgress/dietary_intervention")
parent_meta <- readRDS(file.path(root, "data", "rds", "tcr_meta_only.rds"))
stopifnot(!anyDuplicated(rownames(parent_meta)), all(meta$barcode %in% rownames(parent_meta)))
meta$cluster_id <- as.character(parent_meta[meta$barcode, "seurat_clusters"])
stopifnot(!anyNA(meta$cluster_id), all(nzchar(meta$cluster_id)))
rm(parent_meta)
cluster_ids <- sort(unique(meta$cluster_id))
if (all(grepl("^[0-9]+$", cluster_ids))) cluster_ids <- cluster_ids[order(as.numeric(cluster_ids))]
cluster_choices <- setNames(cluster_ids, paste("Cluster", cluster_ids))
umap_ranges <- list(x = range(meta$x), y = range(meta$y))
review_path <- Sys.getenv("EXPLORER_REVIEWS", "annotation-history.csv")
gene_totals <- Matrix::rowSums(data$counts)
library_sizes <- Matrix::colSums(data$counts)
marker_choices <- c("클러스터 색상" = "", setNames(rownames(data$counts), rownames(data$counts)))
groups <- intersect(c("CD", "WD", "RF", "ADF", "TRF"), unique(meta$group))
groups <- c(groups, setdiff(sort(unique(meta$group)), groups))
annotations <- c("상위 유형" = "celltype_l0", "세포 유형" = "celltype_l1",
                 "세부 유형" = "celltype_l2", "T세포 상태" = "label_exh",
                 "T세포 v2a" = "label_v2a", "검토 반영 주석" = "reviewed")
diet_colors <- c(CD = "#000000", WD = "#fc225c", RF = "#a55af8",
                 ADF = "#3b0070", TRF = "#64c3fa")
chart_font <- list(family = "Arial, Apple SD Gothic Neo, sans-serif", color = "#292b2e")
initial <- sort(unique(na.omit(meta$celltype_l1)))

# Place selectize menus outside the scrolling settings body.
explorer_select <- function(...) {
  args <- list(...)
  if (identical(args$selectize, FALSE)) return(do.call(shiny::selectInput, args))
  args$selectize <- NULL
  args$options <- modifyList(list(dropdownParent = "body"), args$options %||% list())
  do.call(shiny::selectizeInput, args)
}

stage_tabs <- tabsetPanel(id = "stage", selected = "home",
  tabPanel("Data Overview", value = "home"),
  tabPanel("새 클러스터링", value = "recluster"),
  tabPanel("Cell Type Clustering", value = "cluster"),
  tabPanel("Cell Type Annotation", value = "annotate"),
  tabPanel("DE Analysis", value = "de")
)
# Hide in the initial HTML, before the Shiny session connects. showTab() restores it.
stage_tabs <- htmltools::tagQuery(stage_tabs)$find("a")$
  filter(function(tag, i) identical(tag$attribs[["data-value"]], "recluster"))$
  parent()$addAttrs(style = "display: none;")$allTags()

ui <- fluidPage(
  title = "Single Cell Explorer",
  tags$head(tags$link(rel = "stylesheet", href = "style.css?v=dual-view-3"),
            tags$script(src = "navigation.js?v=dual-view-3", defer = NA),
            tags$script(src = "layout.js?v=dual-view-3", defer = NA),
            tags$link(rel = "icon", type = "image/png", href = "explorer-logo.png")),
  tags$header(class = "top", tags$h1(class = "brand",
    tags$img(src = "explorer-logo.png", alt = "", width = 44, height = 44),
    "Single Cell Explorer"),
    stage_tabs
  ),
  tags$main(
    tags$div(class = "dataset-status",
      tags$p("마우스 간 · 식이 개입 · scRNA-seq"),
      tags$div(class = "inventory", sprintf("%s cells / %d mice / %d diets",
        format(nrow(meta), big.mark = ","), length(unique(meta$sample_id)), length(groups)))
    ),
    tags$div(class = "workbench",
      tags$aside(class = "panel controls",
        tags$div(class = "controls-body",
        tabsetPanel(id = "sidebar_tab", selected = "tools",
          tabPanel("분석 도구", value = "tools",
        conditionalPanel("input.stage === 'home'",
          tags$h2("어떻게 시작할까요?"),
          tags$p(class = "hint", "원본 세포 유형을 살펴보고 분석 방법을 선택하세요."),
          tags$div(class = "entry-choice",
            actionButton("use_original", "기존 주석으로 분석하기", class = "btn-primary"),
            tags$p("저장된 세포 유형을 선택해 식이군별 차등 발현을 분석합니다.")),
          tags$div(class = "entry-choice",
            actionButton("new_clustering", "새 클러스터링 시작하기", class = "btn-default"),
            tags$p("전체 또는 관심 세포 유형을 골라 새 클러스터링을 준비합니다."),
            tags$p(class = "hint", "현재 대상 선택까지 지원 · 실행 기능 준비 중")),
          tags$hr(),
          tags$p(class = "hint", "표시 중인 주석 · 원본 celltype_l1")
        ),
        conditionalPanel("input.stage === 'recluster'",
          tags$h2("클러스터링 대상"),
          radioButtons("recluster_scope", "세포 범위",
            c("전체 세포" = "all", "원본 세포 유형으로 선택" = "types")),
          conditionalPanel("input.recluster_scope === 'types'",
            explorer_select("recluster_types", "세포 유형", initial, multiple = TRUE,
              options = list(placeholder = "하나 이상의 세포 유형 선택"))),
          textOutput("recluster_summary"),
          tags$hr(),
          tags$p(class = "hint", "선택한 세포를 전체 UMAP에서 미리 봅니다. 새 클러스터 계산과 결과 저장은 다음 구현 단계입니다."),
          actionButton("back_home", "Data Overview로 돌아가기", class = "btn-default")
        ),
        conditionalPanel("input.stage === 'cluster' || input.stage === 'annotate'",
          tags$h2("클러스터"),
          explorer_select("cluster", "기존 Seurat 클러스터", c("전체 클러스터" = "", cluster_choices)),
          textOutput("cluster_summary"),
          conditionalPanel("input.stage === 'cluster'",
            tags$hr(),
            explorer_select("marker_gene", "발현 확인 유전자", choices = NULL,
              options = list(allowEmptyOption = TRUE, placeholder = "유전자 검색 또는 마커 표에서 선택")),
            actionButton("to_annotation", "이 클러스터의 주석 검토", class = "btn-default")
          ),
          conditionalPanel("input.stage === 'annotate'",
            tags$hr(),
            textInput("new_label", "세포 유형 라벨"),
            textAreaInput("evidence", "마커 또는 변경 근거", rows = 3),
            actionButton("apply_label", "라벨 적용", class = "btn-primary"),
            tags$div(class = "review-actions toolbar-group",
              actionButton("flag_review", "검토 중으로 표시"),
              actionButton("restore_label", "원본 라벨로 복원")),
            tags$hr(),
            actionButton("to_de", "이 클러스터로 DE 분석", class = "btn-default")
          )
        ),
        conditionalPanel("input.stage === 'de'",
        tags$h2("분석 대상"),
        explorer_select("color", "색상 기준", c("주석으로 색칠" = "annotation", "식이군으로 색칠" = "group")),
        radioButtons("selection_mode", "선택 방식",
          c("Cell type → 세부 유형" = "hierarchy", "다른 주석으로 선택" = "advanced")),
        conditionalPanel("input.selection_mode === 'hierarchy'",
          explorer_select("parent_type", "상위 Cell type", c("전체 Cell type" = "", initial), selectize = FALSE),
          conditionalPanel("input.parent_type !== ''",
            explorer_select("subtypes", "세부 유형", choices = NULL, multiple = TRUE,
              options = list(placeholder = "이 Cell type의 전체 세부 유형", plugins = list("remove_button"))),
            tags$p(class = "hint", "여러 세부 유형을 함께 선택할 수 있습니다. 비워두면 상위 Cell type 전체를 분석합니다.")),
          textOutput("hierarchy_summary")
        ),
        conditionalPanel("input.selection_mode === 'advanced'",
          explorer_select("annotation", "주석 수준", annotations, "celltype_l1"),
          explorer_select("types", "세포 유형", initial, multiple = TRUE,
            options = list(placeholder = "전체 유형", plugins = list("remove_button")))
        ),
        explorer_select("de_cluster", "클러스터 범위", c("전체 클러스터" = "", cluster_choices)),
        radioButtons("scope", "세포 범위", c("필터에 해당하는 전체 세포" = "all",
                      "UMAP에서 선택한 세포" = "selected")),
        actionLink("clear", "UMAP 선택 지우기"),
        tags$hr(), tags$h2("식이군 비교"),
        explorer_select("test", "비교군", groups, "RF"),
        explorer_select("reference", "기준군", groups, "WD"),
        numericInput("minimum", "마우스당 최소 세포 수", 20, min = 1, step = 1),
        actionButton("run", "DE 분석 실행", class = "btn-primary"),
        tags$p(class = "hint", "각 식이군에서 기준을 충족하는 마우스가 최소 2마리 필요합니다.")
        )),
          tabPanel("Annotation", value = "legend",
            tags$p(class = "hint", "항목을 누르면 전체 위치와 확대 화면, 구성 정보를 함께 확인합니다."),
            actionLink("clear_annotation", "전체 선택 범위 보기"),
            uiOutput("annotation_list"),
            uiOutput("annotation_info"),
            conditionalPanel("output.has_annotation_focus",
              actionButton("annotation_to_de", "이 Annotation으로 DE 분석", class = "btn-primary"))
          )
        ))
      ),
      tags$section(class = "panel map",
        tags$div(class = "map-heading", tags$div(class = "map-heading-content",
        tags$div(class = "section-top", tags$div(tags$h2("전체 UMAP"), textOutput("selection"))),
        conditionalPanel("input.stage === 'de'",
          tags$p(class = "hint", "전체 UMAP에서 분석 대상만 색으로 강조합니다. 회색은 대상 밖 세포입니다. 드래그로 강조된 세포를 선택할 수 있습니다.")),
        conditionalPanel("input.stage === 'cluster' || input.stage === 'annotate'",
          tags$p(class = "hint", "전체 UMAP에서 선택한 클러스터를 강조합니다. 다른 점을 클릭하면 해당 클러스터로 이동합니다.")),
        conditionalPanel("input.stage === 'home'",
          tags$p(class = "hint", "전체 위치를 유지하며 선택한 Annotation을 강조합니다. 점을 누르거나 오른쪽 목록에서 선택하세요.")),
        conditionalPanel("input.stage === 'recluster'",
          tags$p(class = "hint", "원본 UMAP에서 선택 범위를 표시합니다.")))),
        plotlyOutput("umap", height = "550px"),
        textOutput("cell", container = function(...) tags$div(class = "cell-detail", ...))
      ),
      tags$section(class = "panel detail-map",
        tags$div(class = "map-heading", tags$div(class = "map-heading-content",
          tags$div(class = "section-top", tags$div(tags$h2("선택 영역 확대"), textOutput("detail_selection"))),
          textOutput("detail_hint", container = function(...) tags$p(class = "hint", ...)))),
        plotlyOutput("detail_umap", height = "550px"),
        tags$div(class = "cell-detail", "오른쪽 Annotation에서 관심 유형의 구성 정보를 확인하세요.",
          actionLink("show_annotations", "Annotation 열기"))
      )
    ),
    conditionalPanel("input.stage === 'cluster' || input.stage === 'annotate'",
      tags$section(class = "panel results",
        tags$h2("클러스터 마커"),
        tags$p(class = "hint", "나머지 세포 대비 집계 CPM 비율이 높은 후보입니다. 클러스터 내 검출률 10% 이상. 행을 클릭하면 Clustering 탭의 UMAP에 발현을 표시합니다."),
        DTOutput("markers"))),
    conditionalPanel("input.stage === 'annotate'",
      tags$section(class = "panel results",
        tags$h2("주석 변경 이력"), textOutput("review_status"), DTOutput("history"))),
    conditionalPanel("input.stage === 'de'",
      tags$section(class = "panel results",
      tabsetPanel(
        tabPanel("차등 발현", uiOutput("result_header"), plotlyOutput("volcano", height = "330px"),
                 DTOutput("genes")),
        tabPanel("마우스별 세포 수", tags$p(class = "hint",
                 "현재 세포 선택과 최소 세포 수를 기준으로 포함 여부를 표시합니다."), DTOutput("mice"))
      )
    ))
  )
)

server <- function(input, output, session) {
  annotation_focus <- reactiveVal(NULL)
  chosen <- reactiveVal(integer())
  clicked <- reactiveVal(NULL)
  revision <- reactiveVal(0L)
  result <- reactiveVal(NULL)
  error <- reactiveVal(NULL)

  observeEvent(input$use_original, {
    updateTabsetPanel(session, "sidebar_tab", selected = "tools")
    chosen(integer())
    clicked(NULL)
    updateRadioButtons(session, "selection_mode", selected = "hierarchy")
    updateSelectInput(session, "parent_type", selected = "")
    updateSelectizeInput(session, "subtypes", selected = character())
    updateSelectInput(session, "annotation", selected = "celltype_l1")
    updateSelectInput(session, "de_cluster", selected = "")
    updateSelectizeInput(session, "types", selected = character())
    updateRadioButtons(session, "scope", selected = "all")
    updateSelectInput(session, "color", selected = "annotation")
    updateTabsetPanel(session, "stage", selected = "de")
  })
  observeEvent(input$new_clustering, {
    showTab("stage", target = "recluster", session = session)
    updateTabsetPanel(session, "stage", selected = "recluster")
  })
  observeEvent(input$back_home, updateTabsetPanel(session, "stage", selected = "home"))
  recluster_cells <- reactive({
    if (identical(input$recluster_scope, "types")) {
      which(meta$celltype_l1 %in% input$recluster_types)
    } else meta$id
  })
  output$recluster_summary <- renderText({
    ids <- recluster_cells()
    sprintf("%s개 세포 · %d마리 선택", format(length(ids), big.mark = ","),
            length(unique(meta$sample_id[ids])))
  })

  history <- reactivePoll(1500, session,
    function() file.info(review_path)[c("mtime", "size")],
    function() read_reviews(review_path))
  local_history <- reactiveVal(NULL)
  reviews <- reactive({
    disk <- history()
    recent <- local_history()
    if (!is.null(recent) && nrow(recent) > nrow(disk)) recent else disk
  })
  hierarchy_mode <- reactive(!identical(input$selection_mode, "advanced"))
  subtype_choices <- reactive({
    if (!isTruthy(input$parent_type)) return(character())
    values <- meta$celltype_l2[meta$celltype_l1 %in% input$parent_type]
    sort(unique(values[!is.na(values) & nzchar(values)]))
  })
  observeEvent(input$parent_type, {
    freezeReactiveValue(input, "subtypes")
    updateSelectizeInput(session, "subtypes", choices = subtype_choices(), selected = character(), server = TRUE)
  }, ignoreNULL = FALSE, priority = 10)
  effective_annotation <- reactive({
    if (hierarchy_mode()) {
      if (isTruthy(input$parent_type)) "celltype_l2" else "celltype_l1"
    } else req(input$annotation)
  })
  labels <- reactive({
    field <- effective_annotation()
    if (field == "reviewed") reviewed_labels(meta, reviews()) else meta[[field]]
  })
  output$hierarchy_summary <- renderText({
    if (!isTruthy(input$parent_type)) return("전체 세포 · 상위 Cell type을 고르면 세부 유형이 표시됩니다.")
    detail <- if (length(input$subtypes)) paste(input$subtypes, collapse = ", ") else "전체 세부 유형"
    paste(input$parent_type, "→", detail)
  })

  observe({
    values <- if (identical(input$annotation, "reviewed")) reviewed_labels(meta, reviews()) else meta[[req(input$annotation)]]
    values <- sort(unique(na.omit(values)))
    old_types <- isolate(input$types)
    updateSelectizeInput(session, "types", choices = values,
                         selected = intersect(old_types, values), server = TRUE)
  })
  observeEvent(input$annotation, {
    updateSelectizeInput(session, "types", selected = character())
  }, ignoreInit = TRUE)

  observe({
    updateSelectizeInput(session, "marker_gene", choices = marker_choices, selected = "", server = TRUE)
  })
  cluster_cells <- reactive({
    req(input$cluster %in% cluster_ids)
    which(meta$cluster_id == input$cluster)
  })
  cluster_history <- reactive({
    req(input$cluster %in% cluster_ids)
    h <- reviews()
    h[h$cluster == input$cluster, , drop = FALSE]
  })
  output$cluster_summary <- renderText({
    if (!length(input$cluster) || !nzchar(input$cluster)) return("UMAP에서 점을 클릭하거나 클러스터를 고르세요.")
    ids <- cluster_cells()
    original <- sort(table(meta$celltype_l1[ids], useNA = "ifany"), decreasing = TRUE)
    sprintf("%s개 세포 · %d마리 · 원본 주석 %s",
            format(length(ids), big.mark = ","), length(unique(meta$sample_id[ids])),
            paste(sprintf("%s (%d)", names(original), original), collapse = ", "))
  })
  output$review_status <- renderText({
    h <- cluster_history()
    if (!nrow(h)) return("원본 주석을 사용 중입니다.")
    latest <- h[nrow(h), ]
    switch(latest$status, applied = paste("적용된 라벨", latest$label),
           pending = "검토 중 · 원본 주석 사용", restored = "원본 주석으로 복원됨")
  })
  output$history <- renderDT({
    h <- cluster_history()
    h$status <- unname(c(applied = "적용", pending = "검토 중", restored = "복원")[h$status])
    datatable(h[rev(seq_len(nrow(h))), -1L, drop = FALSE],
      rownames = FALSE, colnames = c("라벨", "근거", "상태", "변경 시각 (UTC)"),
      options = list(pageLength = 8, dom = "tip", ordering = FALSE, scrollX = TRUE))
  })
  observeEvent(input$cluster, {
    h <- cluster_history()
    ids <- cluster_cells()
    originals <- unique(na.omit(meta$celltype_l1[ids]))
    latest <- tail(h, 1)
    label <- if (nrow(latest) && latest$status == "applied") latest$label else
      if (length(originals) == 1L) originals else ""
    updateTextInput(session, "new_label", value = label)
    updateTextAreaInput(session, "evidence", value = "")
  })
  observeEvent(input$to_annotation, {
    req(cluster_cells())
    updateTabsetPanel(session, "stage", selected = "annotate")
  })
  observeEvent(input$to_de, {
    req(cluster_cells())
    chosen(integer())
    updateRadioButtons(session, "selection_mode", selected = "advanced")
    updateRadioButtons(session, "scope", selected = "all")
    updateSelectInput(session, "annotation", selected = "reviewed")
    updateSelectizeInput(session, "types", selected = character())
    updateSelectInput(session, "de_cluster", selected = input$cluster)
    updateTabsetPanel(session, "stage", selected = "de")
  })

  persist_review <- function(status) {
    req(cluster_cells())
    label <- trimws(input$new_label %||% "")
    evidence <- trimws(input$evidence %||% "")
    if (status == "applied" && (!nzchar(label) || !nzchar(evidence))) {
      showNotification("라벨과 마커 근거를 입력하세요.", type = "warning")
      return(invisible(NULL))
    }
    tryCatch({
      updated <- save_review(review_path, input$cluster,
        if (status == "restored") "" else label, evidence, status)
      local_history(updated)
      showNotification("주석 변경을 저장했습니다.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error"))
  }
  observeEvent(input$apply_label, persist_review("applied"))
  observeEvent(input$flag_review, persist_review("pending"))
  observeEvent(input$restore_label, persist_review("restored"))

  markers <- reactive(cluster_markers(data$counts, cluster_cells(), gene_totals))
  output$markers <- renderDT({
    validate(need(input$cluster %in% cluster_ids, "클러스터를 선택하면 마커 후보를 볼 수 있습니다."))
    tab <- datatable(markers(), rownames = FALSE, selection = "single",
      colnames = c("유전자", "log₂ CPM 비율", "클러스터 내 검출 (%)"),
      options = list(pageLength = 8, scrollX = TRUE, order = list(list(1, "desc"), list(2, "desc"))))
    formatRound(tab, c("log2_ratio", "detected_pct"), digits = 2)
  })
  observeEvent(input$markers_rows_selected, {
    row <- input$markers_rows_selected
    req(length(row) == 1L)
    updateSelectizeInput(session, "marker_gene", choices = marker_choices,
                         selected = markers()$gene[row], server = TRUE)
    updateTabsetPanel(session, "stage", selected = "cluster")
  })

  visible_cells <- reactive({
    if (hierarchy_mode()) {
      valid <- rep(TRUE, nrow(meta))
      if (isTruthy(input$parent_type)) {
        valid <- meta$celltype_l1 %in% input$parent_type
        if (length(input$subtypes)) valid <- valid & meta$celltype_l2 %in% input$subtypes
      }
    } else {
      label <- labels()
      valid <- !is.na(label) & nzchar(label)
      if (length(input$types)) valid <- valid & label %in% input$types
    }
    if (isTruthy(input$de_cluster)) valid <- valid & meta$cluster_id == input$de_cluster
    which(valid)
  })
  observeEvent(list(input$selection_mode, input$parent_type, input$subtypes,
                    input$annotation, input$types, input$de_cluster, input$clear), {
    chosen(integer())
    clicked(NULL)
    revision(isolate(revision()) + 1L)
  }, ignoreInit = TRUE)

  # A shared model keeps colours and cell identities identical in both views.
  map_state <- reactive({
    stage <- req(input$stage)
    de_view <- stage == "de"
    overview <- stage %in% c("home", "recluster")
    focused <- (!overview && !de_view && isTruthy(input$cluster)) ||
      (stage == "recluster" && identical(input$recluster_scope, "types"))
    ids <- if (de_view) visible_cells() else if (stage == "recluster") recluster_cells() else
      if (focused) cluster_cells() else meta$id
    by_diet <- de_view && identical(input$color, "group")
    all_categories <- if (overview) meta$celltype_l1 else if (by_diet) meta$group else if (de_view) labels() else
      if (stage == "annotate") reviewed_labels(meta, reviews()) else paste("Cluster", meta$cluster_id)
    all_categories <- as.character(all_categories)
    all_categories[is.na(all_categories) | !nzchar(all_categories)] <- "미주석"
    levels <- sort(unique(all_categories))
    colors <- if (by_diet) diet_colors[levels] else grDevices::hcl.colors(length(levels), "Dark 3")
    colors[is.na(colors)] <- "#555555"
    field <- if (overview) "celltype_l1" else if (by_diet) "group" else if (de_view) effective_annotation() else
      if (stage == "annotate") "reviewed" else "cluster_id"
    list(ids = ids, categories = all_categories, palette = setNames(colors, levels),
         field = field, stage = stage, by_diet = by_diet)
  })
  observeEvent(map_state(), { annotation_focus(NULL) }, ignoreInit = TRUE, priority = 20)
  observeEvent(input$clear_annotation, annotation_focus(NULL))
  observeEvent(input$show_annotations, updateTabsetPanel(session, "sidebar_tab", selected = "legend"))
  observeEvent(input$annotation_pick, {
    st <- map_state()
    if (input$annotation_pick %in% st$categories[base_view_ids()]) annotation_focus(input$annotation_pick)
  })
  base_view_ids <- reactive({
    ids <- map_state()$ids
    if (identical(input$stage, "de") && identical(input$scope, "selected")) ids <- intersect(ids, chosen())
    ids
  })
  focus_ids <- reactive({
    st <- map_state()
    ids <- base_view_ids()
    if (isTruthy(annotation_focus())) ids <- ids[st$categories[ids] == annotation_focus()]
    ids
  })
  output$has_annotation_focus <- reactive(isTruthy(annotation_focus()) && length(focus_ids()) > 0)
  outputOptions(output, "has_annotation_focus", suspendWhenHidden = FALSE)
  output$annotation_list <- renderUI({
    st <- map_state()
    counts <- table(st$categories[base_view_ids()])
    if (!length(counts)) return(tags$p("현재 선택 범위에 세포가 없습니다."))
    tagList(tags$h3(if (st$by_diet) "Diet" else "Annotation"),
      tags$div(class = "annotation-list", role = "group", `aria-label` = "Annotation 목록",
        lapply(names(counts), function(label) {
          active <- identical(annotation_focus(), label)
          tags$button(type = "button", class = paste("annotation-item", if (active) "is-active" else ""),
            `data-annotation` = label, `aria-pressed` = if (active) "true" else "false",
            tags$span(class = "annotation-dot", style = paste0("background:", st$palette[[label]]), `aria-hidden` = "true"),
            tags$span(class = "annotation-name", label),
            tags$span(class = "annotation-count", format(as.integer(counts[[label]]), big.mark = ",")))
        })))
  })
  breakdown <- function(values, title) {
    values <- as.character(values)
    values[is.na(values) | !nzchar(values)] <- "미주석"
    counts <- sort(table(values), decreasing = TRUE)
    tags$details(class = "composition", tags$summary(title),
      tags$table(tags$thead(tags$tr(tags$th("유형"), tags$th("세포"), tags$th("비율"))),
        tags$tbody(lapply(seq_along(counts), function(i) tags$tr(tags$td(names(counts)[i]),
          tags$td(format(as.integer(counts[i]), big.mark = ",")),
          tags$td(sprintf("%.1f%%", 100 * counts[i] / sum(counts))))))))
  }
  output$annotation_info <- renderUI({
    st <- map_state()
    ids <- focus_ids()
    if (!isTruthy(annotation_focus())) return(tags$p(class = "hint", "위 Annotation을 선택하면 세포 수, 샘플 및 세부 유형 구성을 표시합니다."))
    if (!length(ids)) return(tags$p("선택한 세포 범위와 겹치는 세포가 없습니다."))
    tags$section(class = "annotation-info", `aria-live` = "polite",
      tags$h3(annotation_focus()),
      tags$div(class = "annotation-stats",
        tags$div(tags$strong(format(length(ids), big.mark = ",")), tags$span("세포")),
        tags$div(tags$strong(sprintf("%.1f%%", 100 * length(ids) / nrow(meta))), tags$span("전체 세포 중")),
        tags$div(tags$strong(length(unique(meta$sample_id[ids]))), tags$span("마우스"))),
      tags$p(class = "hint", sprintf("현재 선택 범위 중 %.1f%% · 주석: %s", 100 * length(ids) / max(1, length(base_view_ids())), st$field)),
      tags$p(class = "hint", "데이터에 저장된 주석과 세포 구성을 표시합니다. Annotation 탐색만으로 DE 대상은 바뀌지 않습니다."),
      breakdown(meta$celltype_l1[ids], "상위 Cell type 구성"),
      breakdown(meta$celltype_l2[ids], "세부 유형 구성"),
      breakdown(meta$group[ids], "식이군 구성"),
      breakdown(meta$sample_id[ids], "샘플별 세포 수"))
  })
  # Preserve exactly the inspected cells when handing off to DE, including reviewed
  # annotations, a cluster restriction, or an existing lasso selection.
  pending_annotation_de <- reactiveVal(NULL)
  observeEvent(input$annotation_to_de, {
    req(isTruthy(annotation_focus()), length(focus_ids()))
    pending_annotation_de(focus_ids())
    updateRadioButtons(session, "selection_mode", selected = "hierarchy")
    updateSelectInput(session, "parent_type", selected = "")
    updateSelectizeInput(session, "subtypes", selected = character())
    updateSelectInput(session, "de_cluster", selected = "")
    updateRadioButtons(session, "scope", selected = "selected")
    updateTabsetPanel(session, "stage", selected = "de")
    updateTabsetPanel(session, "sidebar_tab", selected = "tools")
  })
  observe({
    ids <- pending_annotation_de()
    if (!is.null(ids) && identical(input$stage, "de") && hierarchy_mode() &&
        !isTruthy(input$parent_type) && !length(input$subtypes) && !isTruthy(input$de_cluster) &&
        identical(input$scope, "selected")) {
      chosen(ids)
      pending_annotation_de(NULL)
    }
  }, priority = -10)
  output$detail_hint <- renderText({
    if (identical(input$stage, "de")) "같은 좌표에서 선택 범위를 확대합니다. 드래그로 세포를 선택하고 휠로 확대할 수 있습니다." else
      "같은 좌표에서 선택 범위를 확대합니다. 드래그로 이동하고 휠로 확대할 수 있습니다."
  })
  output$detail_selection <- renderText({
    ids <- focus_ids()
    if (length(ids) == nrow(meta)) return("관심 세포를 선택하세요")
    paste(annotation_focus() %||% "현재 선택 범위", "·", format(length(ids), big.mark = ","), "개 세포")
  })
  make_umap <- function(detail = FALSE) {
    revision()
    st <- map_state()
    stage <- st$stage
    de_view <- stage == "de"
    ids <- focus_ids()
    points <- meta[ids, ]
    all_categories <- st$categories
    categories <- all_categories[ids]
    levels <- names(st$palette)
    palette <- st$palette
    points$details <- paste(htmltools::htmlEscape(categories), points$sample_id, points$barcode, sep = "<br>")
    p <- plot_ly(source = if (detail) "detail_cells" else "cells", colors = palette)
    context_ids <- if (detail) integer() else setdiff(meta$id, ids)
    if (length(context_ids)) {
      context <- meta[context_ids, c("x", "y", "id", "cluster_id")]
      context$details <- paste(htmltools::htmlEscape(all_categories[context_ids]),
        paste("Cluster", context$cluster_id), "선택 대상 밖", sep = "<br>")
      p <- add_trace(p, data = context, x = ~x, y = ~y, key = ~id, customdata = ~id,
                       type = "scattergl", mode = "markers", text = ~details, hoverinfo = "text",
                       name = "선택 대상 밖", showlegend = FALSE,
                       marker = list(size = 3, color = "#aeb7c0", opacity = 0.45),
                       selected = list(marker = list(opacity = 0.45)),
                       unselected = list(marker = list(opacity = 0.45)))
    }
    gene <- input$marker_gene
    if (length(ids) && stage == "cluster" && isTruthy(gene) && gene %in% rownames(data$counts)) {
      expression <- as.numeric(data$counts[gene, ids, drop = FALSE])
      points$expression <- log1p(1e4 * expression / pmax(library_sizes[ids], 1))
      points$details <- paste(points$details, gene, sprintf("%.2f", points$expression), sep = "<br>")
      points <- points[order(points$expression), ]
      p <- add_trace(p, data = points, x = ~x, y = ~y, key = ~id, customdata = ~id,
                       type = "scattergl", mode = "markers", text = ~details, hoverinfo = "text",
                       name = gene, showlegend = FALSE,
                       marker = list(size = 4, color = ~expression, colorscale = "Viridis",
                                     showscale = detail, colorbar = list(title = "log1p CP10K")))
    } else if (length(ids)) {
      points$label <- factor(categories, levels = levels)
      p <- add_trace(p, data = points, x = ~x, y = ~y, key = ~id, customdata = ~id,
                       type = "scattergl", mode = "markers", text = ~details, hoverinfo = "text",
                       color = ~label,
                       marker = list(size = 4, opacity = 0.8))
    }
    ranges <- umap_ranges
    if (detail && length(ids)) {
      pad_range <- function(x) { r <- range(x); pad <- max(diff(r) * 0.06, 0.2); r + c(-pad, pad) }
      ranges <- list(x = pad_range(points$x), y = pad_range(points$y))
    }
    p <- layout(p, dragmode = if (de_view) "lasso" else "pan", font = chart_font,
      paper_bgcolor = "#ffffff", plot_bgcolor = "#ffffff", showlegend = FALSE,
      xaxis = list(title = "UMAP 1", zeroline = FALSE, range = ranges$x, constrain = "domain", fixedrange = !detail),
      yaxis = list(title = "UMAP 2", zeroline = FALSE, scaleanchor = "x", constrain = "domain", range = ranges$y, fixedrange = !detail),
      margin = list(t = 12, l = 48, r = if (detail && stage == "cluster" && isTruthy(input$marker_gene)) 85 else 12, b = 46))
    p <- event_register(event_register(event_register(p, "plotly_selected"), "plotly_deselect"), "plotly_click")
    config(p, displaylogo = FALSE, scrollZoom = detail,
           modeBarButtonsToRemove = as.list(if (detail) c("select2d") else c("zoom2d", "pan2d", "zoomIn2d", "zoomOut2d", "autoScale2d", "resetScale2d", "select2d")))
  }
  output$umap <- renderPlotly(make_umap(FALSE))
  output$detail_umap <- renderPlotly({
    n <- length(focus_ids())
    if (!n || n == nrow(meta)) {
      message <- if (!n) "현재 선택 범위에 세포가 없습니다." else "Annotation 또는 Cell type을 선택하면<br>이곳에 확대 화면이 나타납니다."
      p <- plot_ly(x = numeric(), y = numeric(), type = "scattergl", mode = "markers", source = "detail_cells")
      p <- layout(p, paper_bgcolor = "#f8f9fa", plot_bgcolor = "#f8f9fa", font = chart_font,
        xaxis = list(visible = FALSE, fixedrange = TRUE), yaxis = list(visible = FALSE, fixedrange = TRUE),
        annotations = list(list(x = .5, y = .5, xref = "paper", yref = "paper", text = message, showarrow = FALSE)),
        margin = list(t = 0, l = 0, r = 0, b = 0))
      p <- event_register(event_register(event_register(p, "plotly_selected"), "plotly_deselect"), "plotly_click")
      return(config(p, displayModeBar = FALSE))
    }
    make_umap(TRUE)
  })
  observeEvent(event_data("plotly_selected", source = "detail_cells", priority = "event"), {
    req(input$stage == "de")
    selection <- event_data("plotly_selected", source = "detail_cells", priority = "event")
    ids <- intersect(as.integer(unlist(selection$key %||% selection$customdata)), focus_ids())
    chosen(ids)
    if (length(ids)) updateRadioButtons(session, "scope", selected = "selected")
  })
  observeEvent(event_data("plotly_click", source = "detail_cells", priority = "event"), {
    click <- event_data("plotly_click", source = "detail_cells", priority = "event")
    id <- as.integer(unlist(click$key %||% click$customdata))[1]
    req(is.finite(id), id %in% focus_ids())
    clicked(id)
    if (input$stage %in% c("cluster", "annotate")) updateSelectInput(session, "cluster", selected = meta$cluster_id[id])
  })
  observeEvent(input$clear, updateRadioButtons(session, "scope", selected = "all"))

  observeEvent(event_data("plotly_selected", source = "cells", priority = "event"), {
    req(input$stage == "de")
    selection <- event_data("plotly_selected", source = "cells", priority = "event")
    ids <- intersect(as.integer(unlist(selection$key %||% selection$customdata)), focus_ids())
    chosen(ids)
    if (length(ids)) updateRadioButtons(session, "scope", selected = "selected")
  })
  observeEvent(event_data("plotly_deselect", source = "cells", priority = "event"), {
    if (identical(input$stage, "de")) chosen(integer())
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  observeEvent(event_data("plotly_click", source = "cells", priority = "event"), {
    click <- event_data("plotly_click", source = "cells", priority = "event")
    id <- as.integer(unlist(click$key %||% click$customdata))[1]
    req(is.finite(id), id %in% meta$id)
    if (input$stage == "de") {
      if (id %in% visible_cells()) clicked(id)
    } else if (identical(input$stage, "home")) {
      annotation_focus(meta$celltype_l1[id])
      updateTabsetPanel(session, "sidebar_tab", selected = "legend")
    } else if (input$stage %in% c("cluster", "annotate")) {
      updateSelectInput(session, "cluster", selected = meta$cluster_id[id])
    }
  })

  selected <- reactive({
    if (input$scope == "selected") intersect(chosen(), visible_cells()) else visible_cells()
  })
  output$selection <- renderText({
    if (identical(input$stage, "home")) return(sprintf("원본 주석 · 전체 %s개 세포", format(nrow(meta), big.mark = ",")))
    if (identical(input$stage, "recluster")) return(sprintf("전체 %s개 중 %s개 선택",
      format(nrow(meta), big.mark = ","), format(length(recluster_cells()), big.mark = ",")))
    if (!identical(input$stage, "de")) {
      total <- format(nrow(meta), big.mark = ",")
      if (!isTruthy(input$cluster)) return(paste("전체", total, "개 세포"))
      return(sprintf("전체 %s개 중 %s개 강조", total, format(length(cluster_cells()), big.mark = ",")))
    }
    sprintf("전체 %s개 · 화면 강조 %s개 · 분석 대상 %s개",
            format(nrow(meta), big.mark = ","), format(length(focus_ids()), big.mark = ","),
            format(length(selected()), big.mark = ","))
  })
  output$cell <- renderText({
    if (input$stage %in% c("home", "recluster")) return("원본 세포 유형 · celltype_l1")
    if (!identical(input$stage, "de")) {
      if (input$stage == "cluster" && isTruthy(input$marker_gene)) {
        return(paste(input$marker_gene, "· log1p(세포별 10,000 카운트당 발현량)"))
      }
      return("기존 Seurat UMAP 좌표")
    }
    id <- clicked()
    if (is.null(id)) return("선택한 개별 세포가 여기에 표시됩니다.")
    paste(meta$barcode[id], meta$group[id], meta$sample_id[id],
          labels()[id], sep = "  ·  ")
  })
  samples <- reactive({
    req(input$test, input$reference, input$minimum)
    mouse_counts(meta, selected(), input$test, input$reference, input$minimum)
  })
  output$mice <- renderDT(datatable(samples(), rownames = FALSE,
    colnames = c("마우스", "식이군", "세포 수", "포함"), options = list(pageLength = 12, dom = "tip", scrollX = TRUE)))

  settings <- reactive(list(cells = selected(), test = input$test, reference = input$reference,
                            minimum = input$minimum, cluster = input$de_cluster,
                            selection_mode = if (hierarchy_mode()) "hierarchy" else "advanced",
                            parent_type = if (hierarchy_mode()) input$parent_type %||% "" else "",
                            annotation = effective_annotation(),
                            cell_types = sort(unique(labels()[selected()]))))
  observeEvent(input$run, {
    result(NULL)
    error(NULL)
    tryCatch({
      s <- settings()
      if (identical(s$test, s$reference)) stop("서로 다른 식이군을 선택하세요.")
      if (!isTRUE(is.finite(s$minimum) && s$minimum >= 1 && s$minimum == floor(s$minimum))) {
        stop("마우스당 최소 세포 수를 양의 정수로 입력하세요.")
      }
      if (!length(s$cells)) stop("분석할 세포를 선택하세요.")
      mice <- samples()
      table <- withProgress(message = "마우스별 RNA 카운트 집계 및 차등 발현 분석 중", value = 0.5,
        run_de(data$counts, meta, s$cells, mice, s$test, s$reference))
      result(list(table = table, settings = s, mice = mice))
    }, error = function(e) error(conditionMessage(e)))
  })

  output$result_header <- renderUI({
    r <- result()
    if (!is.null(error())) return(tags$p(class = "notice", error()))
    if (is.null(r)) return(tags$p(class = "hint", "분석 대상과 두 식이군을 고른 후 DE 분석을 실행하세요."))
    s <- r$settings
    included <- table(r$mice$group[r$mice$included])
    tags$div(class = "section-top result-heading",
      tags$div(tags$h2(paste(s$test, "vs", s$reference)),
        tags$p(class = "hint", sprintf(
          "%s n=%d · %s n=%d · %s개 유전자 검정 · 양의 log₂FC는 %s에서 높은 발현",
          s$test, included[[s$test]], s$reference, included[[s$reference]],
          format(nrow(r$table), big.mark = ","), s$test)),
        tags$p(class = "hint", paste(c(
          if (nzchar(s$parent_type)) s$parent_type,
          if (isTruthy(s$cluster)) paste("Cluster", s$cluster),
          paste(s$cell_types, collapse = ", ")), collapse = " · ")),
        tags$p(class = "hint", "마우스별 pseudobulk · edgeR QL · BH FDR (현재 비교에서 검정한 유전자)"),
        if (!identical(s, settings())) tags$p(class = "notice", "선택이 변경되었습니다. 새 선택을 분석하려면 DE를 다시 실행하세요.")),
      downloadButton("download", "CSV 다운로드")
    )
  })
  output$genes <- renderDT({
    req(result())
    tab <- datatable(result()$table, rownames = FALSE, options = list(pageLength = 10,
                     scrollX = TRUE, order = list(list(5, "asc"))))
    tab <- formatRound(tab, c("logFC", "logCPM", "F"), digits = 3)
    formatSignif(tab, c("PValue", "FDR"), digits = 3)
  })
  output$volcano <- renderPlotly({
    req(result())
    tab <- result()$table
    tab$significance <- -log10(pmax(tab$FDR, .Machine$double.xmin))
    tab$highlight <- factor(ifelse(tab$FDR < 0.05, "FDR < 0.05", "FDR ≥ 0.05"),
                            levels = c("FDR ≥ 0.05", "FDR < 0.05"))
    p <- plot_ly(tab, x = ~logFC, y = ~significance, color = ~highlight,
                 colors = c("FDR ≥ 0.05" = "#8b959f", "FDR < 0.05" = "#292b2e"), text = ~Gene,
                 type = "scattergl", mode = "markers", marker = list(size = 5, opacity = 0.65))
    layout(p, font = chart_font, paper_bgcolor = "#ffffff",
           plot_bgcolor = "#ffffff", xaxis = list(title = "log₂ fold change"),
           yaxis = list(title = "−log₁₀ FDR"),
           margin = list(t = 12), legend = list(orientation = "h"))
  })
  output$download <- downloadHandler(
    filename = function() paste0("DE_", result()$settings$test, "_vs_", result()$settings$reference, ".csv"),
    content = function(file) {
      r <- result()
      s <- r$settings
      mice <- r$mice[r$mice$included, ]
      write.csv(transform(r$table, comparison = paste(s$test, "vs", s$reference),
                          selection_mode = s$selection_mode, parent_cell_type = s$parent_type,
                          cluster = s$cluster %||% "", annotation = s$annotation,
                          cell_types = paste(s$cell_types, collapse = " | "),
                          minimum_cells = s$minimum, cells_used = sum(mice$cells),
                          mice_test = sum(mice$group == s$test),
                          mice_reference = sum(mice$group == s$reference)), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
