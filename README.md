# Single Cell Explorer

Finder에서 **start.command**를 더블 클릭하면 로컬 브라우저에서 열립니다. 처음 실행할 때는 원본 Seurat 객체에서 RNA 카운트와 UMAP을 읽어 앱 폴더에 `explorer-data.rds` 캐시를 만듭니다. 원본 RDS를 갱신한 경우 앱을 닫고 이 캐시를 삭제한 뒤 다시 실행하세요. 앱을 종료하려면 실행한 터미널에서 Ctrl+C를 누르세요.

## 기존 Cell type으로 분석하기

Data Overview에서 **기존 주석으로 분석하기**를 누르고, **상위 Cell type**에서 B cell, Macrophage, DC 등을 고릅니다. **세부 유형**에는 선택한 Cell type의 원본 `celltype_l2` 주석만 표시됩니다. 여러 유형을 함께 선택할 수 있고, 비워두면 해당 Cell type 전체를 사용합니다. 상위 Cell type을 바꾸면 세부 선택과 UMAP 영역 선택이 초기화됩니다.

전체 UMAP은 회색 배경으로 유지되며 분석 대상만 색으로 강조됩니다. 클러스터 범위를 추가로 지정하면 두 조건을 모두 만족하는 세포를 사용합니다. 식이군을 고르고 **DE 분석 실행**을 누르면 선택한 세포로 마우스별 pseudobulk DE를 실행합니다. 결과와 CSV에는 상위 Cell type 및 분석에 사용된 세부 유형이 기록됩니다. T세포 상태나 수정한 주석을 사용하려면 **다른 주석으로 선택**으로 전환하세요.

## 클러스터 마커 검토와 재주석

1. **Cell Type Clustering**에서 기존 Seurat 클러스터를 고릅니다. 마커 후보 표에서 유전자를 클릭하거나 검색해 UMAP의 발현을 확인합니다.
2. **이 클러스터의 주석 검토**를 누릅니다. **Cell Type Annotation**에서 원본 라벨을 확인하고, 새 라벨과 마커 근거를 입력해 적용합니다. 클러스터 전체에 적용되며, 검토 중 표시나 원본 복원도 가능합니다.
3. **이 클러스터로 DE 분석**을 누르면 선택한 클러스터와 검토 반영 주석이 **DE Analysis**에 전달됩니다. 비교군과 기준군을 고르고 마우스별 세포 수를 확인한 뒤 분석을 실행합니다. 클러스터 범위를 전체로 바꾸면 같은 수정 라벨을 가진 세포를 여러 클러스터에서 함께 선택할 수 있습니다.
4. Volcano plot과 유전자 표를 확인하고 결과를 CSV로 내려받습니다.

주석 변경 이력은 앱 폴더의 `annotation-history.csv`에 저장됩니다. 원본 Seurat 주석은 그대로 보존합니다. 검토 중과 원본 복원 상태에서는 원본 라벨을 사용합니다. 기존 클러스터는 원본 폴더의 `tcr_meta_only.rds`에서 바코드로 연결하므로 기존 발현 캐시를 다시 만들 필요가 없습니다.

## 입력 및 분석

- `tcr_integrated.seurat.rds`의 RNA raw counts와 parent UMAP을 사용합니다.
- 세포 바코드로 `tcr_meta_only.rds`의 세포 유형과 `v2a_meta_with_exh.rds`의 T세포 상태를 연결합니다. T세포 상태도 parent UMAP에 표시합니다.
- 마우스별로 선택 세포의 raw counts를 합산합니다. 마우스당 최소 세포 수 기본값은 20이며, 비교군과 기준군에 각각 최소 2마리가 필요합니다.
- edgeR `filterByExpr` → TMM → robust quasi-likelihood 검정을 적용합니다. 양의 log₂FC는 비교군의 높은 발현을 뜻합니다. BH FDR은 현재 비교에서 검정한 유전자에 대해 계산합니다.
- 마커 후보는 클러스터 내 검출률이 10% 이상인 유전자 중 나머지 세포보다 집계 CPM이 높은 유전자입니다. UMAP 발현은 세포별 `log1p(CP10K)`입니다.
- R 패키지 `Seurat`, `Matrix`, `shiny`, `plotly`, `DT`, `edgeR`, `statmod`, `filelock`가 필요합니다. 다른 위치의 프로젝트를 사용하려면 터미널에서 `./start.command "/path/to/dietary_intervention"`으로 실행하세요. 별도의 캐시 경로는 `EXPLORER_CACHE`로 지정할 수 있습니다.

## 화면 구성 구상

현재 앱은 **Cell Type Clustering**, **Cell Type Annotation**, **DE Analysis**를 연결합니다. 나머지 탭은 설계안입니다.

| 탭 | 역할 |
| --- | --- |
| Cell Type Clustering | UMAP에서 기존 클러스터를 선택하고 마우스별 구성 및 후보 마커 확인. Gene Explorer에서 대표 마커와 반대 근거가 되는 마커의 발현을 함께 확인 |
| Cell Type Annotation | 선택한 클러스터의 현재 라벨과 마커 근거를 확인. 라벨을 수정하거나 검토 중으로 표시하고, 원본 라벨·수정 라벨·근거·변경 이력을 기록 |
| DE Analysis | 선택한 세포 집단에서 두 식이군을 비교하고 volcano plot 및 유전자 표를 확인 |
| Pathway Analysis | 실행한 DE 결과를 이어받아 경로 수준으로 해석 |
| Gene Explorer | 유전자를 검색해 UMAP 발현과 세포 유형별 발현량·발현 세포 비율을 확인. 클러스터 마커와 DE 결과에서 유전자를 클릭해 이동 |
| Cell Composition | 마우스별 세포 유형 수와 비율을 식이군별로 비교. 전체 세포·면역세포·T세포 등 분모 선택 |
| Sample & QC | 샘플을 고르면 해당 10x Genomics Cell Ranger Web Summary HTML을 표시. 샘플 ID와 보고서 파일의 매핑 필요 |
| Saved Analyses | 세포 선택, 주석 버전, 비교군, 포함 마우스, 분석 설정을 저장하고 다시 열기 |

### 연결 흐름

- 클러스터 선택 → 마커 검토 → 필요 시 Annotation에서 라벨 수정 → 수정된 라벨을 Composition 및 DE의 세포 선택에 사용합니다.
- 주석 수정은 원본 Seurat 객체와 별도 파일에 저장합니다. 과거 DE 결과는 당시의 주석 버전과 연결합니다.
- 탭 간에 현재 선택한 세포 집단과 식이군을 유지합니다.
- 현재 dietary_intervention 폴더에서는 Web Summary HTML을 찾지 못했습니다.

## 공통 디자인

Mouse Sanctuary의 2026-11-11 확정 디자인 규약을 따른다. 최상위 탭은 헤더에 위치하며 좁은 창에서는 줄바꿈한다. Annotation 범례는 오른쪽 패널의 Annotation 탭에 모은다. 전체 좌표와 비율을 보존하며 화면 높이에 따라 UMAP 공간을 조정한다. 새 클러스터링 탭은 시작 버튼을 눌렀을 때만 나타난다. 상세 기준과 검증 결과는 프로젝트 AGENTS.md에 기록한다.

데스크톱에서는 전체 UMAP과 선택 영역 확대를 나란히 배치하고 분석 도구·Annotation을 오른쪽 패널에 모은다. 실제 브라우저 높이를 기준으로 제목과 하단 정보 줄까지 포함한 UMAP 카드 전체가 들어오도록 조정한다. 긴 도구 설정은 패널 내부에서 스크롤하며, 분석 결과는 작업 영역 아래에서 확인한다. 폭 800px 이하에서는 UMAP 아래로 도구를 재배치한다.

## 전체 위치와 선택 영역 함께 보기

왼쪽 전체 UMAP은 원본 좌표 범위를 유지하며 선택한 세포를 강조한다. 가운데 선택 영역 확대는 같은 좌표에서 해당 세포만 확대한다. 오른쪽 패널의 **분석 도구 / Annotation** 탭에서 설정과 공통 범례를 이용한다.

**Annotation**에서 유형을 누르면 세포 수, 전체 대비 비율, 마우스 수를 표시한다. 상위 Cell type, 세부 유형, 식이군, 샘플별 구성을 펼쳐 확인할 수 있다. 이 클릭만으로 DE 대상은 변경되지 않는다. **이 Annotation으로 DE 분석**을 누르면 현재 확인한 세포들이 선택 집합으로 전달된다. 도구에서 식이군을 지정하고 분석을 실행한다. 전체 선택 범위 보기로 Annotation 강조를 해제할 수 있다.

현재 전체 데이터는 최초 UMAP에서만 표시하며, 확대 패널은 유형/클러스터/개별 세포 선택 후 나타난다. 좁은 화면에서는 그래프를 세로로 재배치한다.

## GitHub 저장소와 로컬 데이터

저장소 이름은 `nplab_singlecell_explorer`이며 앱 표시 이름은 **Single Cell Explorer**를 유지합니다. 코드, 로고, 스타일, 실행 안내를 버전 관리합니다. 원본 데이터, `explorer-data.rds` 캐시, 로컬 주석 변경 이력과 분석 내보내기 폴더는 `.gitignore`로 제외합니다.

다른 컴퓨터에서 실행하려면 위의 R 패키지와 원본 데이터가 별도로 필요합니다. 원본 프로젝트의 `data/rds/`에 다음 파일을 준비하세요.

- `tcr_integrated.seurat.rds`
- `tcr_meta_only.rds`
- `v2a_meta_with_exh.rds`

저장소 폴더에서 데이터 프로젝트 경로를 지정해 실행합니다.

```bash
./start.command "/path/to/dietary_intervention"
```

캐시가 이미 있어도 기존 클러스터 연결을 위해 `tcr_meta_only.rds`가 필요합니다. `.gitignore`는 이미 추적 중인 파일을 자동으로 제거하지 않습니다.
