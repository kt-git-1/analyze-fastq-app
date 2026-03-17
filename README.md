## 概要

**analyze-fastq-app** は、FASTQ の取得からマッピング、QC、VCF 出力までを一括実行する NGS 解析パイプラインです。

主な処理は次の通りです。

1. ENA から FASTQ をダウンロード、またはローカル FASTQ を読み込む
2. AdapterRemoval + BWA で BAM を作成する
3. Soft clipping と BAM 前処理を行う
4. mapDamage / Qualimap / GATK HaplotypeCaller を実行する
5. 中間ファイルを削除し、結果を整理して保存する

古 DNAと現代 DNA の両方を想定しており、`--data_type` により挙動を切り替えられます。

---

## 目次

- [概要](#概要)
- [このパイプラインでできること](#このパイプラインでできること)
- [Quick Start](#QuickStart)
- [入力と出力](#入力と出力)
- [ワークフロー](#ワークフロー)
- [FASTQ 判定ルール](#fastq-判定ルール)
- [実行オプション](#実行オプション)
- [セットアップ](#セットアップ)
- [ディレクトリ構成](#ディレクトリ構成)
- [各モジュールの役割](#各モジュールの役割)
- [トラブルシューティング](#トラブルシューティング)

---

## このパイプラインでできること

- ENA の project accession から FASTQ を取得する
- 手元の FASTQ ディレクトリを sample 単位に自動認識する
- 複数レーンの FASTQ をサンプル単位に統合する
- AdapterRemoval と BWA を用いて BAM を作成する
- Soft clipping、重複除去、index 作成を行う
- mapDamage、Qualimap、HaplotypeCaller の結果をまとめて出力する
- 解析後に不要な BAM などを削除して容量を節約する

---

## QuickStart

### ENA からダウンロードしてそのまま解析

```sh
python main.py \
  --project_accession PRJEB19970 \
  --download-via-https \
  --reference_genome ./data/reference/equCab3.nochrUn.fa
```

### 手元の FASTQ を解析

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa
```

### ダウンロードと解析を分ける

```sh
# Step 1: HTTPS でダウンロード
python -m modules.ena_download_https --project PRJEB19970 --out ./data/raw_data --workers 2

# Step 2: ダウンロード済み FASTQ を解析
python main.py --fastq_dir ./data/raw_data/PRJEB19970 --base_dir ./data
```

---

## 入力と出力

### 入力

- ENA の project accession
  - 例: `PRJEB19970`
- またはローカル FASTQ ディレクトリ
  - 例: `./my_fastqs`
- 参照ゲノム FASTA
  - 例: `./data/reference/equCab3.nochrUn.fa`

### 出力

主な出力先は次の通りです。

- ログ: `data/logs/pipeline_<project>.log`
- ダウンロード FASTQ: `data/raw_data/<project>/`
- サンプルごとの解析結果: `data/results/<project>/<sample>/`
- VCF: `data/results/<project>/<sample>/vcf_files/`
- mapDamage や Qualimap の結果: 各 sample ディレクトリ配下

---

## ワークフロー

このプログラムは「FASTQ を入力にして、サンプルごとに BAM と QC と VCF を作る」ものです。最初に理解するとよいのは、入口が 2 つあることです。

1. ENA から始める
2. 手元の FASTQ から始める

### 例 1: ENA から始める

```sh
python main.py \
  --project_accession PRJEB19970 \
  --download-via-https \
  --reference_genome ./data/reference/equCab3.nochrUn.fa
```

内部では次の順で処理されます。

1. `ena_download_https.py` が FASTQ をダウンロードする
2. ダウンロード済み FASTQ を sample ごとに見つける
3. 複数レーンがあれば 1 サンプル 1 セットにまとめる
4. AdapterRemoval を実行する
5. BWA で BAM を作る
6. Soft clipping を行う
7. Picard / samtools で重複除去と index 作成を行う
8. mapDamage / Qualimap / HaplotypeCaller を実行する
9. 中間ファイルを削除して終了する

### 例 2: 手元の FASTQ から始める

たとえば次のような FASTQ があるとします。

```text
my_fastqs/
  SAMPLE_A_L001_R1_001.fastq.gz
  SAMPLE_A_L001_R2_001.fastq.gz
  SAMPLE_A_L002_R1_001.fastq.gz
  SAMPLE_A_L002_R2_001.fastq.gz
  SAMPLE_B_R1.fastq.gz
  SAMPLE_B_R2.fastq.gz
```

この場合は次のように実行します。

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa
```

内部では次のように解釈されます。

- `SAMPLE_A` は 2 レーン分の paired-end データ
- `SAMPLE_B` は 1 レーン分の paired-end データ
- `SAMPLE_A` の 2 レーンは 1 本の `R1` と 1 本の `R2` に統合される
- その後、各 sample に対して同じ解析が行われる

### 全体の流れ

```text
FASTQ
  ↓
sample / R1 / R2 を判定
  ↓
レーンを統合
  ↓
AdapterRemoval
  ↓
BWA mapping
  ↓
Soft clipping
  ↓
BAM 前処理
  ↓
mapDamage / Qualimap / HaplotypeCaller
  ↓
結果を保存
```

### ancient と modern の違い

- `ancient`
  - AdapterRemoval で `--collapse` を使う
  - collapsed と non-collapsed を分けてマッピングし、最後にマージする
- `modern`
  - 通常の paired-end / single-end として処理する

迷った場合は、古 DNA なら `ancient`、通常の現代サンプルなら `modern` を選べば十分です。

### 初回実行時に見る場所

まずは次の 3 点を確認すると全体像を掴みやすいです。

1. `data/logs/pipeline_<project>.log` に各ステップの開始と終了が出ているか
2. `data/results/<project>/<sample>/` が sample ごとに作られているか
3. `vcf_files` や QC 出力が生成されているか

---

## FASTQ 判定ルール

`--fastq_dir` を指定した場合、`parse_fastq_general()` がディレクトリを再帰探索し、各 FASTQ を sample ごとに `R1` / `R2` / `single` に分類します。

### 対象ファイル

- `*.fastq`
- `*.fastq.gz`
- `*.fq`
- `*.fq.gz`

### read 種別の判定

次のような命名規則に対応しています。

- Illumina 形式: `SAMPLE_A_L001_R1_001.fastq.gz`
- Filgen 形式: `Horse01_L1_1.fq.gz`
- 一般的な `R1/R2` 形式: `sampleB_R1.fastq.gz`, `sampleB-R2.fastq.gz`
- 数字だけのペア形式: `sampleC_1.fq.gz`, `sampleC_2.fq.gz`
- どれにも当てはまらない場合は `single`

例:

- `SAMPLE_A_L001_R1_001.fastq.gz` と `SAMPLE_A_L002_R1_001.fastq.gz` は `SAMPLE_A` の `R1`
- `SAMPLE_A_L001_R2_001.fastq.gz` と `SAMPLE_A_L002_R2_001.fastq.gz` は `SAMPLE_A` の `R2`
- `ERR123456.fastq.gz` は `ERR123456` の `single`

### sample 名の決め方

sample 名は次の優先順位で決まります。

1. サブディレクトリ配下にある場合は最上位ディレクトリ名
2. そうでない場合はファイル名から推定した sample 名
3. それも取れない場合は拡張子を除いたファイル名

たとえば次のような構成では、ファイル名よりディレクトリ名が優先されます。

```text
my_fastqs/
  nested_sample/
    reads_1.fastq.gz
    reads_2.fastq.gz
```

この場合、`reads_1.fastq.gz` / `reads_2.fastq.gz` は `nested_sample` に分類されます。

### レーンのマージ

`merge_lanes_by_cat()` は sample ごとの複数レーンを統合します。

- `R1` が複数本あれば 1 本に連結
- `R2` が複数本あれば 1 本に連結
- single-end が複数本あれば 1 本に連結
- paired-end と single-end が混在する場合は paired-end を優先

最終的に、下流処理には sample ごとに `[R1, R2]` または `[single]` の形で渡されます。

---

## 実行オプション

`main.py` は次の引数を受け取ります。

| 引数 | 説明 | デフォルト |
| ---- | ---- | ---------- |
| `--project_accession` | ENA の project ID。`--fastq_dir` を指定しない場合に使用 | `PRJEB19970` |
| `--base_dir` | データ保存先のベースディレクトリ | `./data` |
| `--reference_genome` | 参照ゲノム FASTA | `./data/reference/equCab3.nochrUn.fa` |
| `--workers` | ダウンロード時の並列数 | `4` |
| `--threads` | 解析に使うスレッド数 | `20` |
| `--java_mem` | Java ツール用メモリ | `10g` |
| `--fastq_dir` | 事前にダウンロードした FASTQ ディレクトリ | `None` |
| `--download-via-https` | HTTPS ダウンロードしてから解析を実行 | `False` |
| `--data_type` | `ancient` または `modern` | `ancient` |

---

## セットアップ

### 必要条件

- Python 3.8 以上
- BWA
- samtools
- AdapterRemoval
- Picard
- mapDamage
- Qualimap
- GATK
- pysam
- 参照ゲノム FASTA

### Conda 環境

```sh
conda env create -f environment.yml
conda activate analyze-env
```

`environment.yml` には Python 側依存が含まれます。外部ツールは各環境で別途インストールしてください。

Picard は `modules/bam_processor.py` の実装上、`/usr/local/bin/picard.jar` を想定しています。別パスを使う場合はコード側の修正が必要です。

---

## ディレクトリ構成

```text
analyze-fastq-app/
├── config.py
├── environment.yml
├── main.py
├── modules/
│   ├── analyzers.py
│   ├── bam_processor.py
│   ├── bwa_mapper.py
│   ├── ena_downloader.py
│   ├── ena_download_https.py
│   ├── softclipper.py
│   └── __init__.py
├── data/
│   ├── raw_data/
│   ├── reference/
│   ├── results/
│   ├── logs/
│   └── temp/
└── README.md
```

---

## 各モジュールの役割

- `main.py`: パイプライン全体の制御
- `config.py`: 引数解析、ログ設定、FASTQ 検出、レーン統合
- `modules/ena_downloader.py`: ENA API と通常ダウンロード
- `modules/ena_download_https.py`: HTTPS ダウンロードと再開処理
- `modules/bwa_mapper.py`: AdapterRemoval と BWA マッピング
- `modules/softclipper.py`: BAM の soft clipping
- `modules/bam_processor.py`: CleanSam、ReadGroup 付与、重複除去、sort/index
- `modules/analyzers.py`: mapDamage、Qualimap、HaplotypeCaller

---

## トラブルシューティング

### 外部コマンドが見つからない

`bwa`, `samtools`, `AdapterRemoval` などが `PATH` に入っているか確認してください。macOS では Homebrew、Linux では conda やパッケージマネージャ経由での導入が一般的です。

### データがダウンロードできない

project ID やネットワーク接続を確認してください。ENA へのアクセスは `ENADownloader.get_api_response()` と `ena_download_https.py` で行っています。

### BAM や VCF が出力されない

`data/logs/pipeline_<project>.log` を確認してください。どのステップで止まったかを追えます。

### Qualimap で `MaxPermSize` エラーが出る

Java 9 以降では `-XX:MaxPermSize` は廃止されています。`JAVA_TOOL_OPTIONS=-XX:+IgnoreUnrecognizedVMOptions` を使ってください。`QualimapAnalyzer` 側でも自動付与しています。

### ローカル FASTQ の命名ルールに合わない

`parse_fastq_general()` の正規表現を追加してください。現在は Illumina、Filgen、一般的な `R1/R2`、`_1/_2` を扱えます。

---

詳細なパラメータは各モジュールの docstring を参照してください。問い合わせ先は Issue または開発者までお願いします。
