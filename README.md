## 概要

**analyze-fastq-app** は NGS（次世代シーケンス）データ解析を自動化する Python 製パイプラインです。

1. ENA から FASTQ をダウンロード、またはローカル ディレクトリ内の FASTQ を検出
2. AdapterRemoval + BWA によるマッピングと BAM 生成
3. 独自 Soft clipping（pysam ベース）
4. Picard/Samtools による BAM 前処理・重複除去
5. mapDamage / Qualimap / GATK HaplotypeCaller による QC と variant calling

までを一括で実行し、古 DNA・現代 DNA 両方の解析を再現性高く再利用できるように設計されています。

- **対応 OS:** Linux
- **主な機能:** ENA からの FASTQ 取得、BWA マッピング、Soft clipping、BAM sort/duplicate removal、mapDamage・Qualimap・HaplotypeCaller 解析、VCF 出力、解析完了後の中間 BAM 自動削除
- **ローカル FASTQ 処理:** 事前に取得済みの FASTQ データをディレクトリ指定することで読み込み、複数の命名規則を自動判定したうえでペアエンド／シングルエンドを自動認識し、複数レーンを `cat` でマージして解析に供します。
- **Ancient/Modern 切替:** `--data_type` フラグで AdapterRemoval のスレッショルドや collapsed リード処理の有無を切り替え、古 DNA と現代 DNA の両方に最適化されたパイプラインを提供します。

---

## 目次

- [特徴](#特徴)
- [アーキテクチャとディレクトリ構成](#アーキテクチャとディレクトリ構成)
- [必要条件](#必要条件)
- [セットアップ](#セットアップ)
  - [Conda での環境構築](#conda-での環境構築)
- [使い方](#使い方)
  - [コマンドライン引数](#コマンドライン引数)
  - [実行例](#実行例)
- [ワークフロー詳細](#ワークフロー詳細)
- [各種ディレクトリ・ファイルの説明](#各種ディレクトリファイルの説明)
- [トラブルシューティング・FAQ](#トラブルシューティングfaq)

---

## 特徴

- **自動化**: データ取得から解析まで一括実行
- **再現性**: コマンドライン引数でパラメータを管理
- **拡張性**: モジュール構造で機能追加が容易
- **省ストレージ**: mapDamage・Qualimap・HaplotypeCaller の完了を検知して BAM/BAI などの大容量中間ファイルを自動削除
- **ローカル FASTQ 処理**: Illumina (`*_L001_R1_001.fastq.gz`)、Filgen (`*_L1_1.fq.gz`)、シングルトン (`*.fastq.gz`) など複数命名規則を自動判定し、再帰的に FASTQ を探索します。
- **ログ管理**: `data/logs/pipeline_<project>.log` に INFO/ERROR を集約し、失敗したステップをすぐ特定できます。

---

## アーキテクチャとディレクトリ構成

```
analyze-fastq-app/
├── config.py                # 設定と引数解析（FASTQ 検出やレーンマージのユーティリティ含む）
├── environment.yml          # Conda 環境定義（requests, pysam, tqdm 等）
├── main.py                  # パイプライン実行スクリプト
├── modules/                 # 個別モジュール
│   ├── analyzers.py         # mapDamage・Qualimap・HaplotypeCaller
│   ├── bam_processor.py     # Picard CleanSam/MarkDuplicates + samtools sort/index
│   ├── bwa_mapper.py        # AdapterRemoval + BWA (ancient/modern モード)
│   ├── ena_downloader.py    # ENA から FASTQ ダウンロード
│   ├── softclipper.py       # pysam による Soft clipping
│   └── __init__.py
├── data/                    # 実行時に作成されるディレクトリ
│   ├── raw_data/            # ENA からダウンロードした FASTQ
│   ├── reference/           # 参照ゲノム (手動で配置)
│   ├── results/             # 解析結果
│   ├── logs/                # 実行ログ
│   └── temp/                # 一時ファイル
└── README.md
```

`main.py` は各モジュールを組み合わせるハブです。参照 FASTA の `samtools faidx` チェックから始まり、ENA 取得またはローカル FASTQ 解析を分岐し、最終的に mapDamage/Qualimap/HaplotypeCaller の QC & variant calling を実行します。

---

## 必要条件

- Python 3.8 以上
- BWA、samtools、AdapterRemoval、Picard（`/usr/local/bin/picard.jar` を想定）、mapDamage、Qualimap、GATK、pysam が動作する環境
- インターネット接続 (ENA からのデータ取得に使用)
- 解析用の参照ゲノム FASTA ファイル（`samtools faidx` でインデックスを作成。存在しない場合は `main.py` が自動生成）

---

## セットアップ

### Conda での環境構築

```sh
# Conda のインストール (未導入の場合)
# Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# 環境作成
conda env create -f environment.yml
conda activate analyze-env
```

必要に応じて AdapterRemoval・Picard・mapDamage・Qualimap・GATK などを conda または別途インストールしてください。Picard は `modules/bam_processor.py` の実装上 `/usr/local/bin/picard.jar` に配置されている想定なので、別のパスに置く場合は該当スクリプトのパスを書き換えてください。

`environment.yml` には Python 側依存（pysam、tqdm、requests など）が定義されています。外部バイナリ（BWA、samtools 等）は各自の環境に合わせて導入してください。

---

## 使い方

### コマンドライン引数

`main.py` は次の引数を受け取ります。

| 引数 | 説明 | デフォルト |
| ---- | ---- | ---------- |
| `--project_accession` | ENA のプロジェクト ID。`--fastq_dir` を指定しない場合に使用します。 | `PRJEB19970` |
| `--base_dir` | データ保存先のベースディレクトリ | `./data` |
| `--reference_genome` | 参照ゲノム FASTA ファイル | `./data/reference/equCab3.nochrUn.fa` |
| `--workers` | ENA ダウンロードの並列数 | `4` |
| `--threads` | 各種解析のスレッド数 | `20` |
| `--java_mem` | Java ツール用メモリ設定 | `10g` |
| `--fastq_dir` | 事前にダウンロードした FASTQ ファイルが存在するディレクトリへのパス。指定すると ENA からのダウンロードをスキップし、このディレクトリ内の R1/R2 ファイルを自動検出して解析を行います（再帰的探索・複数命名規則に対応）。 | `None` |
| `--data_type` | `ancient`（デフォルト）/`modern`。AdapterRemoval のパラメータと BWA 入力の取り扱いを切り替えます。 | `ancient` |

### 実行例

```sh
# ENA プロジェクトをまるごと解析（古 DNA モード）
python main.py \
  --project_accession PRJEB19970 \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type ancient

# 既にダウンロード済み FASTQ（複数レーン）を解析（現代 DNA モード）
python main.py \
  --fastq_dir /path/to/PE-AncientHorses-fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type modern

# サンプルごとに別フォルダへまとめられた FASTQ を解析する場合
python main.py \
  --fastq_dir ./my_fastqs \
  --base_dir /data/analyze-fastq-app-data \
  --threads 32
```

---

## ワークフロー詳細

1. **FASTQ 取得** (`modules/ena_downloader.py` / `config.parse_fastq_general`)
   - ENA API から sample accession と FTP URL を取得
   - ThreadPoolExecutor で並列ダウンロード、`data/raw_data/<sample>/` へ保存
   - `--fastq_dir` 指定時はローカル探索（再帰）で `*_L###_R1_001.fastq.gz`、`*_L1_1.fq.gz` 等の命名を自動判定
   - サブレーンごとの FASTQ は `config.merge_lanes_by_cat()` で `cat` 連結し、必要最小限のファイルコピーに留める
2. **AdapterRemoval + BWA** (`modules/bwa_mapper.py`)
   - `ancient` モード: AdapterRemoval の `--collapse` 出力（collapsed/R1/R2）をそれぞれ単独で BWA MEM し、`samtools merge` で集約
   - `modern` モード: AdapterRemoval 出力の `R1.truncated` / `R2.truncated` をペアエンドで BWA MEM
   - すべてのマッピングは `samtools view/sort` まで一気に実行し、`results/<sample>/bam_files/` に BAM を保存
3. **Soft clipping** (`modules/softclipper.py`)
   - pysam + ThreadPoolExecutor により BAM を読み込み、各リード先頭に 5 bp の soft clip を付与
   - CIGAR/シーケンス長の整合性をチェックし、不整合リードはログに警告を出してスキップ
4. **BAM 前処理** (`modules/bam_processor.py`)
   - Picard CleanSam → AddOrReplaceReadGroups → MarkDuplicates（REMOVE_DUPLICATES=true）→ samtools sort/index を実行
   - `results/<sample>/dedup/` に `*.marked.dedup.sorted.bam(.bai)` を生成
5. **mapDamage** (`modules/analyzers.py`)
   - `samtools view` + `awk` + `samtools sort` で MAPQ>=30、POS>=300 のフィルタ BAM を生成しインデックス化
   - mapDamage を `results/<sample>/mapdamage/` に実行後、BWA 中間ファイル（`bam_files/*.bam/*.bai/*.truncated`）を削除
6. **Qualimap** (`modules/analyzers.py`)
   - BAM の存在・サイズ・`samtools idxstats` でマッピングリード数を確認し、ゼロの時は自動スキップ
   - Java 9+ でも動作するよう `JAVA_TOOL_OPTIONS=-XX:+IgnoreUnrecognizedVMOptions` を自動付与
7. **GATK HaplotypeCaller** (`modules/analyzers.py`)
   - `results/<sample>/vcf_files/` に `*.vcf` を出力
   - 成功後、重複除去済み BAM・BAI・softclipped BAM を削除してストレージを節約

各ステップは `main.py` から自動的に呼び出され、サンプル単位のログは `data/logs/pipeline_<project>.log` に記録されます。結果は `data/results/<ProjectID>/<SampleID>/` 以下に保存されます。

---

## 各種ディレクトリ・ファイルの説明

- `config.py`: コマンドライン引数とログ設定、FASTQ 検出・レーンマージのユーティリティ
- `main.py`: パイプラインのエントリーポイント（リファレンスのインデックス作成、ENA/ローカル FASTQ の判定、各モジュールの呼び出し、QC 後のクリーンアップ）
- `environment.yml`: Conda で使用する依存パッケージ
- `modules/`: 各機能を実装した Python モジュール
- `data/raw_data/`: ダウンロードした FASTQ
- `data/reference/`: 参照ゲノム (ユーザーが配置)
- `data/results/<ProjectID>/<SampleID>/`: サンプルごとの解析結果 (BAM, VCF, QC レポートなど)
- `data/logs/`: 実行ログ (`pipeline_<project>.log`)
- `data/temp/`: 一時ファイル（AdapterRemoval の中間や mapDamage フィルタリング用）

---

## トラブルシューティング・FAQ

### Q. 外部コマンドが見つからない (bwa, AdapterRemoval など)
A. PATH にバイナリが含まれているか確認してください。Apple Silicon Mac では Homebrew などで手動インストールしてください。

### Q. データがダウンロードできない
A. プロジェクト ID やネットワーク接続を確認してください。ENA API へのアクセスは `modules/ena_downloader.py` の `ENADownloader.get_api_response()` で行っています。

### Q. BAM/VCF などの出力がない
A. `data/logs/` の `pipeline_<project>.log` にモジュールごとの INFO/ERROR が詳細に出力されます。BWA・Picard・mapDamage・Qualimap の実行コマンドはログにも記録されるため、どのステップで停止したかを特定できます。

### Q. Qualimap 実行時に `Unrecognized VM option 'MaxPermSize=1024m'` と表示される
A. Java 9 以降では `-XX:MaxPermSize` が廃止されています。環境変数で `JAVA_TOOL_OPTIONS=-XX:+IgnoreUnrecognizedVMOptions` を付与してください（`QualimapAnalyzer` でも自動付与されています）。

### Q. ローカル FASTQ の命名ルールに合わない
A. `config.py` の `parse_fastq_general()` が Illumina（`*_L001_R1_001.fastq.gz`）、Filgen（`*_L1_1.fq.gz`）、シングルトン（`*.fastq.gz`）をカバーしています。追加の命名ルールをサポートしたい場合は同関数に正規表現を追加してください。

---

詳細な使い方やパラメータは `config.py` や各モジュールの docstring を参照してください。
ご質問・ご要望は Issue または開発者 (谷口: taniguchidev33@gmail.com) までご連絡ください。
