# analyze-fastq-app

FASTQ の取得、前処理、BWA マッピング、BAM 処理、QC、VCF 出力までをまとめて実行する NGS 解析パイプラインです。ENA project accession、ローカル FASTQ、既存 dedup BAM の 3 種類を入口として使えます。

古 DNA と現代 DNA の両方を想定しており、FASTQ から解析する場合は `--data_type ancient|modern` で AdapterRemoval / BWA 周りの処理を切り替えます。

## できること

- ENA から FASTQ をダウンロードして解析する
- HTTPS ダウンロードだけを先に実行し、解析と分けて運用する
- ローカル FASTQ を sample / run / R1 / R2 / single に自動分類する
- AdapterRemoval + BWA MEM で sorted BAM を作成する
- ancient paired-end では collapsed / non-collapsed reads を分けてマッピングし、BAM をマージする
- BAM に soft clipping、CleanSam、sample 単位 merge、MarkDuplicates、sort、index を行う
- 既存 dedup BAM から mapDamage / Qualimap / HaplotypeCaller だけを実行する
- `.done` チェックポイントで完了済みサンプルをスキップし、途中から再開する
- `--parallel_samples` でサンプル単位に並列実行する

## Quick Start

### ENA からダウンロードして解析

```sh
python main.py \
  --project_accession PRJEB19970 \
  --download-via-https \
  --reference_genome ./data/reference/equCab3.nochrUn.fa
```

### ローカル FASTQ を解析

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa
```

### 既存 dedup BAM から QC / VCF だけ実行

```sh
python main.py \
  --bam_dir ./my_dedup_bams \
  --bam_stage dedup \
  --reference_genome ./data/reference/equCab3.nochrUn.fa
```

`--bam_dir` 指定時は FASTQ 判定、AdapterRemoval、BWA、soft clipping、CleanSam、MarkDuplicates を実行しません。入力 BAM は既に重複除去・sort 済みの dedup BAM として扱います。

### ダウンロードと解析を分ける

```sh
# Step 1: HTTPS でダウンロード
python -m modules.ena_download_https \
  --project PRJEB19970 \
  --out ./data/raw_data \
  --workers 2

# Step 2: ダウンロード済み FASTQ を解析
python main.py \
  --fastq_dir ./data/raw_data/PRJEB19970 \
  --base_dir ./data \
  --reference_genome ./data/reference/equCab3.nochrUn.fa
```

### サンプル並列実行

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --parallel_samples 3 \
  --threads 8
```

各サンプルが `--threads` 本のスレッドを使います。合計負荷はおおむね `parallel_samples × threads` です。

## 入力モード

| モード | 主な引数 | 実行される処理 |
| ---- | ---- | ---- |
| ENA から解析 | `--project_accession`, `--download-via-https` | HTTPS ダウンロード → FASTQ 解析 |
| ENA API + HTTP/FTP ダウンロード | `--project_accession` | ENA メタデータ取得 → FASTQ ダウンロード → FASTQ 解析 |
| ローカル FASTQ | `--fastq_dir` | FASTQ 自動判定 → 前処理 → BAM → QC → VCF |
| 既存 dedup BAM | `--bam_dir --bam_stage dedup` | BAM index 確認 → mapDamage → Qualimap → HaplotypeCaller |

`--bam_dir` は `--fastq_dir` および `--download-via-https` と同時指定できません。

## FASTQ 解析ワークフロー

```text
FASTQ
  ↓
sample / run / R1 / R2 / single を判定
  ↓
run 単位の処理
  AdapterRemoval
  ↓
  BWA MEM + samtools sort
  ↓
  soft clipping
  ↓
  Picard CleanSam
  ↓
sample 単位の処理
  samtools merge
  ↓
  Picard MarkDuplicates
  ↓
  samtools sort / index
  ↓
  mapDamage
  ↓
  Qualimap
  ↓
  GATK HaplotypeCaller
  ↓
.done 作成
```

FASTQ から作成された最終 dedup BAM と index は、QC / VCF 完了後に中間ファイルとして削除されます。

### ancient と modern の違い

| `--data_type` | AdapterRemoval | マッピング |
| ---- | ---- | ---- |
| `ancient` | `--minquality 20 --minlength 30 --collapse` | collapsed reads を single-end、non-collapsed reads を paired-end としてマッピングし、run 内でマージ |
| `modern` | `--minquality 25 --minlength 25` | 通常の paired-end / single-end としてマッピング |

## FASTQ 判定ルール

`--fastq_dir` 配下を再帰探索し、次の拡張子を対象にします。

- `.fastq`
- `.fastq.gz`
- `.fq`
- `.fq.gz`

対応している主な命名規則:

| 形式 | 例 | 判定 |
| ---- | ---- | ---- |
| Illumina | `SAMPLE_A_L001_R1_001.fastq.gz` | sample=`SAMPLE_A`, run=`SAMPLE_A_L001`, read=`R1` |
| Filgen | `Horse01_L1_1.fq.gz` | sample=`Horse01`, run=`Horse01_L1`, read=`R1` |
| R1/R2 | `sampleB_R1.fastq.gz`, `sampleB-R2.fastq.gz` | paired-end |
| `_1/_2` | `sampleC_1.fq.gz`, `sampleC_2.fq.gz` | paired-end |
| その他 | `ERR123456.fastq.gz` | single-end |

sample 名は次の優先順位で決まります。

1. FASTQ がサブディレクトリ配下にある場合は、`--fastq_dir` 直下の最上位ディレクトリ名
2. ファイル名から推定した sample 名
3. FASTQ 拡張子を除いたファイル名

run ID は次の優先順位で決まります。

1. `sample/run/file.fastq.gz` のように 3 階層以上ある場合は、FASTQ の直上ディレクトリ名
2. ファイル名から推定した lane 付き run 名
3. sample 名

同じ sample/run に paired-end と single-end が混在する場合は paired-end を優先し、single-end は無視します。

## BAM 入力モード

`--bam_dir` を指定すると、既存 dedup BAM から QC / VCF だけを実行します。v1 では `--bam_stage dedup` のみ対応しています。

```text
dedup BAM
  ↓
sample 名をファイル名から推定
  ↓
BAM index 確認 / 必要なら samtools index
  ↓
mapDamage
  ↓
Qualimap
  ↓
GATK HaplotypeCaller
  ↓
.done 作成
```

### BAM の検出

デフォルトでは `--bam_dir` 配下の `*.bam` を再帰探索します。対象を絞る場合は `--bam_pattern` を指定します。

```sh
python main.py \
  --bam_dir ./my_dedup_bams \
  --bam_pattern "*.dedup.sorted.bam" \
  --reference_genome ./data/reference/equCab3.nochrUn.fa
```

### sample 名の推定

sample 名は BAM ファイル名から推定します。次の suffix を取り除いた残りが sample 名です。

| BAM ファイル名 | sample 名 |
| ---- | ---- |
| `SAMPLE_A.dedup.sorted.bam` | `SAMPLE_A` |
| `SAMPLE_B.sorted.bam` | `SAMPLE_B` |
| `SAMPLE_C.dedup.bam` | `SAMPLE_C` |
| `SAMPLE_D.bam` | `SAMPLE_D` |

同じ sample 名に解釈される BAM が複数見つかった場合はエラー終了します。

### BAM index

`<input>.bam.bai` または `<input>.bai` が存在すれば既存 index を使います。どちらも存在しない場合は `samtools index <input>.bam` を実行します。

入力 BAM と既存 index は削除されません。

## 出力

主な出力先:

```text
data/
  raw_data/<project>/                 # ダウンロード FASTQ
  results/<project>/<sample>/
    mapdamage/
    qualimap/
    vcf_files/
      <sample>.vcf
    .done
  logs/
    pipeline_<project>.log
  temp/
```

`<project>` は通常 `--project_accession` です。`--fastq_dir` 指定時は FASTQ ディレクトリ名、`--bam_dir` 指定時は BAM ディレクトリ名が使われます。

## チェックポイントと再開

各 sample の解析が成功すると、`data/results/<project>/<sample>/.done` が作成されます。再実行時は `.done` がある sample をスキップします。

全 sample を強制的に再実行する場合:

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --force
```

BAM 入力モードでも同じ仕組みで `.done` を使います。

## 実行オプション

| 引数 | 説明 | デフォルト |
| ---- | ---- | ---- |
| `--project_accession` | ENA project accession | `PRJEB19970` |
| `--base_dir` | `raw_data`, `results`, `logs`, `temp` を置くベースディレクトリ | `./data` |
| `--reference_genome` | 参照ゲノム FASTA | `./data/reference/equCab3.nochrUn.fa` |
| `--fastq_dir` | ローカル FASTQ ディレクトリ | `None` |
| `--bam_dir` | 既存 dedup BAM ディレクトリ | `None` |
| `--bam_stage` | 入力 BAM の処理段階。v1 では `dedup` のみ | `dedup` |
| `--bam_pattern` | `--bam_dir` 配下で検出する BAM の glob パターン | `*.bam` |
| `--download-via-https` | HTTPS ダウンロード後に解析する | `False` |
| `--download_protocol` | ENA 通常ダウンロードで使うプロトコル (`ftp` / `http`) | `http` |
| `--workers` | ダウンロード並列数 | `4` |
| `--max_retries` | ダウンロード失敗時の再試行回数 | `3` |
| `--data_type` | FASTQ 解析時のデータ種別 (`ancient` / `modern`) | `ancient` |
| `--threads` | 解析ツールに渡すスレッド数 | `20` |
| `--parallel_samples` | 同時に解析する sample 数 | `1` |
| `--java_mem` | Picard など Java ツール用メモリ | `10g` |
| `--picard_jar` | Picard jar のパス | `/usr/local/bin/picard.jar` |
| `--rg_library` | BWA read group の library (`LB`) | `unknown` |
| `--rg_center` | BWA read group の center (`CN`) | `unknown` |
| `--force` | `.done` を無視して再実行する | `False` |

## セットアップ

### Python / conda

```sh
conda env create -f environment.yml
conda activate analyze-env
```

### 外部ツール

パイプライン起動時に次のツールとファイルを検証します。

- `bwa`
- `samtools`
- `AdapterRemoval`
- `java`
- `gatk`
- `qualimap`
- `mapDamage`
- Picard jar
- 参照ゲノム FASTA
- BWA index (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`)

参照ゲノムの `.fai` がない場合は `samtools faidx` で自動作成します。

注意: BAM 入力モードでも、現状の環境検証は FASTQ 解析用ツールを含めて確認します。

## ディレクトリ構成

```text
analyze-fastq-app/
├── config.py
├── environment.yml
├── main.py
├── modules/
│   ├── analyzers.py
│   ├── bam_parser.py
│   ├── bam_processor.py
│   ├── bwa_mapper.py
│   ├── ena_downloader.py
│   ├── ena_download_https.py
│   ├── fastq_parser.py
│   ├── softclipper.py
│   └── __init__.py
└── README.md
```

## モジュール概要

- `main.py`: 入力モードの切り替え、sample 単位の逐次/並列実行、完了サマリー出力
- `config.py`: CLI 引数、ディレクトリ設定、ログ設定、環境検証、中間ファイル削除
- `modules/fastq_parser.py`: FASTQ の再帰探索、sample / run / read 分類
- `modules/bam_parser.py`: 既存 BAM の再帰探索、sample 名推定、重複 sample 検出
- `modules/ena_downloader.py`: ENA Portal API 取得、FTP/HTTP ダウンロード、MD5 検証
- `modules/ena_download_https.py`: HTTPS ダウンロード専用 CLI、レジューム、`.done` スキップ
- `modules/bwa_mapper.py`: AdapterRemoval、BWA MEM、samtools sort、ancient paired-end の collapsed/non-collapsed merge
- `modules/softclipper.py`: BAM read の先頭 5bp soft clipping
- `modules/bam_processor.py`: CleanSam、sample 単位 merge、MarkDuplicates、sort、index
- `modules/analyzers.py`: mapDamage、Qualimap、GATK HaplotypeCaller

## トラブルシューティング

### 起動直後に止まる

環境検証で不足ツールや不足ファイルが見つかっています。ログに表示されたツール、Picard jar、参照ゲノム、BWA index を確認してください。

### FASTQ が見つからない

`--fastq_dir` のパスと拡張子を確認してください。対象は `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz` です。

### BAM 入力モードで sample 名が重複する

たとえば `SAMPLE_A.bam` と `SAMPLE_A.sorted.bam` はどちらも `SAMPLE_A` と解釈されるためエラーになります。`--bam_pattern` で対象を絞るか、片方を別ディレクトリへ移動してください。

### BAM index 作成に失敗する

入力 BAM が壊れていないか、coordinate sort 済みか、参照ゲノムと整合しているかを確認してください。必要なら事前に `samtools quickcheck` や `samtools index` で確認してください。

### 途中で失敗した sample だけ再実行したい

同じコマンドを再実行してください。`.done` がある sample はスキップされ、未完了 sample だけ実行されます。すべてやり直す場合は `--force` を指定します。

### Qualimap で `MaxPermSize` エラーが出る

Java 9 以降では `-XX:MaxPermSize` が廃止されています。`QualimapAnalyzer` は `JAVA_TOOL_OPTIONS=-XX:+IgnoreUnrecognizedVMOptions` を自動付与します。
