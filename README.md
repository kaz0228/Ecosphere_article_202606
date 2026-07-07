# Ecosphere_article_202606

Data and code for the manuscript **"Transferability of a single-site process-based population model for *Aedes albopictus* and *Culex pipiens*"** (Ueno & Ohta), submitted to *Ecosphere*.

---

## English

### Overview

This repository contains the data and Python code needed to reproduce the analyses in the manuscript. The study applies and extends the Physiology-based Climate-driven Mosquito Population (PCMP) model to evaluate the transferability of a single-site, process-based population model between an urban site (Tokyo, TYO) and a suburban site (Kitamoto, KTM) for *Aedes albopictus* and *Culex pipiens*.

### Repository structure

- **`raw_data/`** — Original, unmodified input data.
  - Meteorological data obtained from the Japan Meteorological Agency (JMA; https://www.data.jma.go.jp/obd/stats/etrn/en/index.php).
  - Mosquito abundance data from previously published sources (Tsuda and Hayashi 2014; Sato and Miyake 2022).

- **`code/`** — Analysis code and processed data.
  - The Python code used for all analyses in this study.
  - The processed (derived) data produced from the files in `raw_data/`.

### Data sources

The mosquito surveillance data are not original to this study and were obtained from the published sources cited above. The meteorological data are openly available from the JMA. Please refer to the manuscript and the original publications for full details, licensing, and any conditions of use.

### Requirements

The analysis code is written in Python. See the code files for the specific libraries and versions used.

### Citation

If you use this code or data, please cite the manuscript (details to be added upon publication).

### Archiving

Upon acceptance of the manuscript, this repository will be permanently archived in Zenodo with an assigned DOI.

---

## 日本語

### 概要

本リポジトリには、論文「Transferability of a single-site process-based population model for *Aedes albopictus* and *Culex pipiens*」(Ueno & Ohta、*Ecosphere* 投稿)の解析を再現するために必要なデータと Python コードが含まれています。本研究では、生理学ベースの気候駆動型蚊個体群モデル(PCMP)を適用・拡張し、単一サイトで構築したプロセスベース個体群モデルの、都市サイト(東京、TYO)から郊外サイト(北本、KTM)への転移可能性を、*Aedes albopictus* と *Culex pipiens* について評価しています。

### リポジトリ構成

- **`raw_data/`** — 加工前の元データ。
  - 気象庁(JMA)から取得した気象データ(https://www.data.jma.go.jp/obd/stats/etrn/en/index.php)。
  - 既発表の文献から取得した蚊の個体数データ(Tsuda and Hayashi 2014; Sato and Miyake 2022)。

- **`code/`** — 解析コードと加工済みデータ。
  - 本研究で使用したすべての Python コード。
  - `raw_data/` を加工して得られた派生データ。

### データの出所

蚊の監視データは本研究で新たに取得したものではなく、上記の既発表文献から取得しました。気象データは気象庁より公開されています。詳細、ライセンス、利用条件については、論文および元の出版物を参照してください。

### 実行環境

解析コードは Python で記述されています。使用したライブラリとバージョンについては、コードファイルを参照してください。

### 引用

本コードまたはデータを利用する場合は、論文を引用してください(詳細は出版時に追記します)。

### アーカイブ

論文の受理後、本リポジトリは Zenodo に恒久的にアーカイブされ、DOI が付与されます。