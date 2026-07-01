# analyze-fastq-app

FASTQ の取得、前処理、BWA マッピング、BAM 処理、QC、VCF 出力、cohort PCA/MDS までをまとめて実行する NGS 解析パイプラインです。

このリポジトリは、主に馬ゲノム解析を想定した構成になっています。デフォルトの参照ゲノム名や例は `equCab3.nochrUn.fa` ですが、`--reference_genome` を指定すれば別の参照FASTAでも実行できます。

## 最初に読む要点

- 入力は3種類です。
  - ENA project accession からFASTQを取得して解析する
  - すでに手元にあるFASTQディレクトリを解析する
  - すでに作成済みのdedup BAMからQC/VCF/PCAだけを実行する
- FASTQ解析では、複数runや複数laneに分かれたFASTQをsample単位にまとめ、最終的にsampleごとのdedup BAMとVCFを作ります。
- `--data_type ancient` では古DNA向けにcollapsed readを考慮した処理を行います。
- `--data_type modern` では通常のpaired-end/single-end readsとして処理します。
- `--data_type auto` はPCA直前にfinal dedup BAMのread lengthから ancient / modern を推定します。省略時のデフォルトではありません。
- PCA/MDSは `--run-pca --pca-sites <共通SNPリスト>` を指定したときだけ実行されます。
- PCAで使うPLINK/EIGENSOFT/Picardの固定パスは [tool_paths.py](tool_paths.py) で管理します。

## できること

- ENAからFASTQをダウンロードして解析する
- HTTPSダウンロードだけを先に実行し、あとで解析する
- ローカルFASTQをsample / run / R1 / R2 / singleに自動分類する
- AdapterRemovalでtrim/collapseし、BWA MEMで参照ゲノムへマッピングする
- run単位BAMをsample単位にmergeし、MarkDuplicates、sort、indexを行う
- BAMにsoft clipping、Picard CleanSam、mapDamage、Qualimap、GATK HaplotypeCallerを実行する
- 既存dedup BAMからmapDamage / Qualimap / HaplotypeCaller / PCAだけを実行する
- `.done` チェックポイントにより、完了済みsampleをスキップして再開する
- `--parallel_samples` でsample単位に並列実行する
- ancient sampleではfinal dedup BAMからpseudo-haploid matrixを作り、cohort PCA/MDSを実行する
- modern sampleではfinal dedup BAMからdiploid genotype matrixを作り、cohort PCA/MDSを実行する

## 実行前の確認

### Python環境

基本のPython環境は `environment.yml` から作成します。

```sh
conda env create -f environment.yml
conda activate analyze-env
```

`plink`、`convertf`、`smartpca` はこのconda環境に入っている必要はありません。このプロジェクトでは固定パスから直接呼び出します。

### 外部ツール

起動時に以下を検証します。

- `bwa`
- `samtools`
- `AdapterRemoval`
- `java`
- `gatk`
- `qualimap`
- `mapDamage`
- Picard jar
- 参照ゲノムFASTA
- BWA index (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`)
- PCA実行時のみ PLINK / CONVERTF / smartpca

PCAで使う外部ツールは次の固定パスです。

```sh
ls -l /usr/local/bin/plink
ls -l /usr/bin/convertf
ls -l /usr/bin/smartpca
ls -l /usr/local/bin/picard.jar
```

このパスは [tool_paths.py](tool_paths.py) で管理しています。別サーバーで場所が違う場合は、次の値を変更してください。

```python
PLINK_BIN = Path("/usr/local/bin/plink")
CONVERTF_BIN = Path("/usr/bin/convertf")
SMARTPCA_BIN = Path("/usr/bin/smartpca")
PICARD_JAR = Path("/usr/local/bin/picard.jar")
```

参照ゲノムの `.fai` がない場合は、実行中に `samtools faidx` で自動作成します。GATK用の `.dict` はHaplotypeCallerに必要です。

## まず動かすコマンド

### 1. ローカルFASTQを解析する

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type ancient
```

`--data_type` を省略すると `ancient` として処理されます。modern sampleなら明示的に指定してください。

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type modern
```

### 2. ENAからダウンロードして解析する

```sh
python main.py \
  --project_accession PRJEB19970 \
  --download-via-https \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type ancient
```

### 3. ダウンロードと解析を分ける

大きいprojectでは、先にダウンロードだけを完了させる運用が安全です。

```sh
python -m modules.ena_download_https \
  --project PRJEB19970 \
  --out ./data/raw_data \
  --workers 2
```

ダウンロード後、FASTQディレクトリを入力にして解析します。

```sh
python main.py \
  --fastq_dir ./data/raw_data/PRJEB19970 \
  --base_dir ./data \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type ancient
```

### 4. 既存dedup BAMからQC/VCFだけ実行する

```sh
python main.py \
  --bam_dir ./my_dedup_bams \
  --bam_stage dedup \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type ancient
```

`--bam_dir` 指定時は、FASTQ判定、AdapterRemoval、BWA、soft clipping、CleanSam、MarkDuplicatesを実行しません。入力BAMはすでに重複除去済み・coordinate sort済みのdedup BAMとして扱います。

### 5. PCA/MDSまで実行する

ancient sample:

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type ancient \
  --run-pca \
  --pca-sites ./horse_common_sites.vcf \
  --pca-engine eigensoft
```

modern sample:

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type modern \
  --run-pca \
  --pca-sites ./horse_common_sites.vcf \
  --pca-engine eigensoft
```

ancient / modernが事前に不明な場合:

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type auto \
  --run-pca \
  --pca-sites ./horse_common_sites.vcf \
  --pca-engine eigensoft
```

`auto` は、final dedup BAMのmapped read length中央値を使ってPCA入力を `ancient` / `modern` に分岐します。cohort中央値が100bp以下なら `ancient`、100bp超なら `modern` として扱います。推定結果は `results/<project>/cohort/auto_data_type_summary.tsv` に出力されます。

## 入力モードの選び方

| やりたいこと | 使う引数 | 備考 |
| ---- | ---- | ---- |
| ENA accessionから直接解析したい | `--project_accession`, `--download-via-https` | HTTPSで取得してから解析 |
| ENA APIの通常ダウンロード経路を使いたい | `--project_accession` | `--download-via-https` なし |
| 手元のFASTQを解析したい | `--fastq_dir` | sample/run/readを自動判定 |
| 既存dedup BAMから後段だけ実行したい | `--bam_dir --bam_stage dedup` | mapDamage / Qualimap / HaplotypeCaller / PCA |

`--bam_dir` は `--fastq_dir` および `--download-via-https` と同時指定できません。

## FASTQディレクトリの指定

`--fastq_dir` には、解析したいFASTQ群をまとめて含むフォルダを指定します。

良い例:

```text
my_fastqs/
  SampleA/
    Run1/
      SampleA_L001_R1_001.fastq.gz
      SampleA_L001_R2_001.fastq.gz
    Run2/
      SampleA_L002_R1_001.fastq.gz
      SampleA_L002_R2_001.fastq.gz
  SampleB/
    Run1/
      SampleB_L001_R1_001.fastq.gz
      SampleB_L001_R2_001.fastq.gz
```

runフォルダがsampleフォルダより上にある構成にも対応します。

```text
ancient/
  210416_A01210_0064_BHWKFYDMXX/
    PE-AncientHorses-01_S1_L001_R1_001.fastq.gz
    PE-AncientHorses-01_S1_L001_R2_001.fastq.gz
  210419_A00581_0145_AHWKL3DMXX/
    PE-AncientHorses-01_S1_L002_R1_001.fastq.gz
    PE-AncientHorses-01_S1_L002_R2_001.fastq.gz
```

この例では `PE-AncientHorses-01_S1` が同じsampleとしてまとめられ、2つのrunとして解析されます。

避けるべき運用:

- 同じsampleのrunを別々に `main.py` で実行する
- 広すぎる親ディレクトリを指定して、別projectのFASTQまで混ぜる
- R1だけ、またはR2だけの不完全なpaired-endセットを入力する

同じsampleの複数runを別々に解析すると、sample単位のmergeやMarkDuplicatesが正しく行われず、深度、重複率、VCF、PCAの入力が分断されます。

## FASTQ判定ルール

対象拡張子:

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
| その他 | `ERR123456.fastq.gz` | single-end。ファイル名をrunとして扱う |

sample名は、ファイル名とディレクトリ構造から推定します。sampleごとのサブディレクトリがある場合は、その最上位ディレクトリ名を優先します。run IDは、直上のrunディレクトリ名、ファイル名のlane情報、single-endのファイル名などから作ります。

同じsample/runにpaired-endとsingle-endが混在する場合はpaired-endを優先します。片側だけのpaired-end FASTQはwarningを出して除外します。R1/R2の本数が一致しない場合は、対応できたペアだけを解析します。

## FASTQ解析ワークフロー

```text
FASTQ
  ↓
sample / run / R1 / R2 / single を判定
  ↓
run単位の処理
  AdapterRemoval
  ↓
  BWA MEM + samtools sort
  ↓
  soft clipping
  ↓
  Picard CleanSam
  ↓
sample単位の処理
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

FASTQから作成されたrun単位BAMやmerge/marking中間BAMは、全sampleの解析が終わってから成功sample分だけまとめて削除します。最終dedup BAMとindexは、再開、QC確認、PCAに使うため保持します。

## ancient / modern / auto の違い

| `--data_type` | FASTQ前処理 | PCA入力 |
| ---- | ---- | ---- |
| `ancient` | AdapterRemoval `--minquality 20 --minlength 30 --collapse` | final dedup BAMからpseudo-haploid alleleを抽出 |
| `modern` | AdapterRemoval `--minquality 25 --minlength 25` | final dedup BAMからref/alt read countを取り、0/1/2 dosageに変換 |
| `auto` | ancient互換でマッピング | PCA直前にread length中央値でancient/modernを推定 |

古DNAでは低depthが多いため、HaplotypeCallerのsample別diploid VCFをそのままPCA入力にはしません。PCAでは、指定した共通SNPリスト上でfinal dedup BAMからpseudo-haploid alleleを抽出してcohort matrixを作ります。

modernでは、指定したSNP座位のref/alt read countからdosageを作ります。alt fractionが0.2以下なら0、0.8以上なら2、それ以外は1として扱います。

## BAM入力モード

`--bam_dir` を指定すると、既存dedup BAMからQC/VCF/PCAだけを実行します。

```text
dedup BAM
  ↓
sample名をファイル名から推定
  ↓
BAM index確認 / 必要なら samtools index
  ↓
mapDamage
  ↓
Qualimap
  ↓
GATK HaplotypeCaller
  ↓
.done 作成
```

デフォルトでは `--bam_dir` 配下の `*.bam` を再帰探索します。対象を絞る場合は `--bam_pattern` を使います。

```sh
python main.py \
  --bam_dir ./my_dedup_bams \
  --bam_pattern "*.dedup.sorted.bam" \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type modern
```

sample名はBAMファイル名から推定します。

| BAMファイル名 | sample名 |
| ---- | ---- |
| `SAMPLE_A.dedup.sorted.bam` | `SAMPLE_A` |
| `SAMPLE_B.sorted.bam` | `SAMPLE_B` |
| `SAMPLE_C.dedup.bam` | `SAMPLE_C` |
| `SAMPLE_D.bam` | `SAMPLE_D` |

同じsample名に解釈されるBAMが複数見つかった場合はエラー終了します。`<input>.bam.bai` または `<input>.bai` があれば既存indexを使います。どちらもなければ `samtools index <input>.bam` を実行します。入力BAMと既存indexは削除しません。

## PCA/MDS

PCA/MDSを実行するには、少なくとも2つ以上の成功sampleと、共通SNPリストが必要です。

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type ancient \
  --run-pca \
  --pca-sites ./horse_common_sites.vcf
```

`--pca-sites` はVCFまたはBED-like textを指定できます。VCFではbiallelic SNPだけを使用し、multi-allelic siteやindelは除外されます。

### EIGENSOFTを使う通常経路

デフォルトは `--pca-engine eigensoft` です。

```text
final dedup BAM
  ↓
ancient: pseudo-haploid raw calls
modern: diploid genotype calls
  ↓
cohort matrix
  ↓
missingness / MAF / sex chromosome filter
  ↓
PLINK text
  ↓
PLINK binary
  ↓
PLINK QC
  ↓
LD pruning
  ↓
CONVERTF
  ↓
EIGENSTRAT
  ↓
smartpca
  ↓
pca_scores.tsv / pca_variance.tsv / mds.tsv
```

外部ツールは [tool_paths.py](tool_paths.py) の固定パスを直接呼びます。PATHやconda activateには依存しません。

### Python fallback

外部ツールなしで処理確認だけしたい場合は、内蔵のNumPyベースPCA/MDSを使えます。

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type ancient \
  --run-pca \
  --pca-sites ./horse_common_sites.vcf \
  --pca-engine python
```

これはPLINK/EIGENSOFTを使わない簡易経路です。最終解析では `eigensoft` を推奨します。

### PCA QCパラメータ

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --data_type ancient \
  --run-pca \
  --pca-sites ./horse_common_sites.vcf \
  --pca-engine eigensoft \
  --pca-min-mapq 30 \
  --pca-min-baseq 30 \
  --pca-trim-ends 2 \
  --pca-max-sample-missing 0.9 \
  --pca-max-site-missing 0.9 \
  --pca-min-maf 0.01 \
  --pca-exclude-sex-chr \
  --pca-ld-window 50 \
  --pca-ld-step 5 \
  --pca-ld-r2 0.2
```

| 引数 | 意味 |
| ---- | ---- |
| `--pca-min-mapq` | BAMからallele/genotypeを拾うときの最小mapping quality |
| `--pca-min-baseq` | allele判定に使う最小base quality |
| `--pca-trim-ends` | read両端から除外する塩基数。古DNAの末端損傷対策 |
| `--pca-max-sample-missing` | sampleごとの許容欠損率 |
| `--pca-max-site-missing` | SNP siteごとの許容欠損率 |
| `--pca-min-maf` | 最小minor allele frequency |
| `--pca-exclude-sex-chr` | X/Y/W/Zや23/24などの性染色体siteを除外 |
| `--pca-ld-window` | PLINK `--indep-pairwise` のwindow size |
| `--pca-ld-step` | PLINK `--indep-pairwise` のstep size |
| `--pca-ld-r2` | PLINK `--indep-pairwise` のr2閾値 |

## 出力

デフォルトでは `--base_dir ./data` の下に出力します。`--fastq_dir` 指定時はFASTQディレクトリ名、`--bam_dir` 指定時はBAMディレクトリ名がproject名として使われます。ENA入力では通常 `--project_accession` がproject名です。

主な出力:

```text
data/
  raw_data/<project>/                 # ダウンロードFASTQ
  results/<project>/
    sample_qc_summary.tsv
    cohort/
      pca_sites.tsv
      auto_data_type_summary.tsv       # data_type=auto の場合
      pseudohaploid_raw_calls.tsv      # ancient PCA
      pseudohaploid_matrix.tsv
      pseudohaploid_matrix.filtered.tsv
      modern_genotype_calls.tsv        # modern PCA
      modern_genotype_matrix.tsv
      modern_genotype_matrix.filtered.tsv
      pca_qc_summary.tsv
      eigenstrat/
        cohort.geno
        cohort.snp
        cohort.ind
      plink/
        cohort.tped
        cohort.tfam
        cohort.bed
        cohort.bim
        cohort.fam
        cohort.qc.*
        cohort.prune.prune.in
        cohort.pruned.*
      pca/
        cohort.evec
        cohort.eval
        pca_scores.tsv
        pca_variance.tsv
        mds.tsv
    <sample>/
      dedup/
        <sample>.dedup.sorted.bam
        <sample>.dedup.sorted.bam.bai
      mapdamage/
      qualimap/
      vcf_files/
        <sample>.vcf
      .done
  logs/
    pipeline_<project>.log
    pipeline_<project>.log.1
  temp/
```

ログファイルは 50MB ごとに自動で切り分けられ、最大5世代まで保持されます。

重要なファイル:

| ファイル | 内容 |
| ---- | ---- |
| `sample_qc_summary.tsv` | 成功sampleごとのdedup BAM、Qualimap主要指標、VCF件数 |
| `<sample>.dedup.sorted.bam` | 後段解析に使う最終BAM |
| `<sample>.vcf` | sampleごとのHaplotypeCaller出力 |
| `pca_qc_summary.tsv` | PCA入力数、残ったsample/site数、QC条件、data_type |
| `pca_scores.tsv` | PCA座標。散布図作成に使う |
| `pca_variance.tsv` | 各PCの固有値と説明分散 |
| `mds.tsv` | MDS座標 |

## チェックポイントと再開

sample解析が成功すると、`results/<project>/<sample>/.done` が作られます。再実行時は `.done` があるsampleをスキップします。

未完了sampleは、既存の中間成果物があればそこから再開します。たとえば `dedup/<sample>.dedup.sorted.bam` が残っている場合は、BWA mapping / soft clipping / CleanSam / MarkDuplicatesを再実行せず、QC / HaplotypeCaller側へ進みます。

PCA stageも途中成果物を再利用します。`pca_sites.tsv`、raw calls、matrix、filtered matrix、PLINK text/binary/QC/LD pruning、EIGENSTRAT、smartpca、`pca_scores.tsv`、`pca_variance.tsv`、`mds.tsv` が既に存在する場合は、`--force` を付けない限り既存出力を使います。

すべてやり直す場合:

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --force
```

## サンプル並列実行

```sh
python main.py \
  --fastq_dir ./my_fastqs \
  --reference_genome ./data/reference/equCab3.nochrUn.fa \
  --parallel_samples 3 \
  --threads 8
```

各sampleが `--threads` 本のスレッドを使います。合計負荷はおおむね `parallel_samples × threads` です。上の例では最大24スレッド相当の負荷になります。

## 実行オプション

| 引数 | 説明 | デフォルト |
| ---- | ---- | ---- |
| `--project_accession` | ENA project accession | `PRJEB19970` |
| `--base_dir` | `raw_data`, `results`, `logs`, `temp` を置くベースディレクトリ | `./data` |
| `--reference_genome` | 参照ゲノムFASTA | `./data/reference/equCab3.nochrUn.fa` |
| `--fastq_dir` | ローカルFASTQディレクトリ | `None` |
| `--bam_dir` | 既存dedup BAMディレクトリ | `None` |
| `--bam_stage` | 入力BAMの処理段階。v1では `dedup` のみ | `dedup` |
| `--bam_pattern` | `--bam_dir` 配下で検出するBAMのglob pattern | `*.bam` |
| `--download-via-https` | HTTPS downloaderで取得してから解析する | `False` |
| `--download_protocol` | ENA通常ダウンロードで使うプロトコル (`ftp` / `http`) | `http` |
| `--workers` | ダウンロード並列数 | `4` |
| `--max_retries` | ダウンロード失敗時の再試行回数 | `3` |
| `--no-progress` | 端末上の進捗バーを無効化する | `False` |
| `--data_type` | FASTQ解析時のデータ種別 (`ancient` / `modern` / `auto`) | `ancient` |
| `--threads` | 解析ツールに渡すスレッド数 | `20` |
| `--parallel_samples` | 同時に解析するsample数 | `1` |
| `--run-pca` | 全sample完了後にcohort PCA/MDS stageを実行する | `False` |
| `--pca-sites` | PCA/MDSで比較する共通SNPリスト | `None` |
| `--pca-engine` | PCA実行エンジン (`eigensoft` / `python`) | `eigensoft` |
| `--pca-min-mapq` | PCA入力抽出時の最小mapping quality | `30` |
| `--pca-min-baseq` | PCA入力抽出時の最小base quality | `30` |
| `--pca-trim-ends` | read両端から除外する塩基数 | `2` |
| `--pca-max-sample-missing` | sample欠損率の上限 | `0.9` |
| `--pca-max-site-missing` | SNP欠損率の上限 | `0.9` |
| `--pca-min-maf` | 最小minor allele frequency | `0.0` |
| `--pca-exclude-sex-chr` | 性染色体SNPをPCA matrixから除外する | `False` |
| `--pca-ld-window` | PLINK `--indep-pairwise` のwindow size | `50` |
| `--pca-ld-step` | PLINK `--indep-pairwise` のstep size | `5` |
| `--pca-ld-r2` | PLINK `--indep-pairwise` のr2閾値 | `0.2` |
| `--java_mem` | PicardなどJavaツール用メモリ | `10g` |
| `--rg_library` | BWA read groupのlibrary (`LB`) | `unknown` |
| `--rg_center` | BWA read groupのcenter (`CN`) | `unknown` |
| `--force` | `.done` と既存PCA成果物を無視して再実行する | `False` |

## 進捗表示

通常実行では、端末に日本語ラベル付きの進捗バーを表示します。ステップ名は外部ツール名やログ検索との対応を保つため、`BWA mapping`、`Soft clipping`、`CleanSam`、`mapDamage`、`Qualimap`、`HaplotypeCaller` のように英語のまま表示されます。

表示される主な情報:

- ダウンロード中: 全体の完了ファイル数、各FASTQの転送量、速度、MD5スキップ/再試行ログ
- 解析中: 全体の完了sample数、各sampleの現在ステップ、step数、soft clippingのread数
- 終了時: 成功sample数、失敗sampleと失敗ステップ

CIやログファイルだけを確認したい場合は、`--no-progress` を付けてください。

## ディレクトリ構成

```text
analyze-fastq-app/
  config.py
  environment.yml
  main.py
  tool_paths.py
  modules/
    analyzers.py
    bam_parser.py
    bam_processor.py
    bwa_mapper.py
    cohort_pca.py
    ena_downloader.py
    ena_download_https.py
    fastq_parser.py
    softclipper.py
  tests/
  README.md
```

## モジュール概要

- `main.py`: 入力モードの切り替え、sample単位の実行、完了サマリー、PCA stage起動
- `config.py`: CLI引数、ディレクトリ設定、ログ設定、環境検証、中間ファイル削除
- `tool_paths.py`: PLINK、CONVERTF、smartpca、Picard jarの固定パス
- `modules/fastq_parser.py`: FASTQの再帰探索、sample/run/read分類
- `modules/bam_parser.py`: 既存BAMの再帰探索、sample名推定、重複sample検出
- `modules/ena_downloader.py`: ENA Portal API取得、FTP/HTTPダウンロード、MD5検証
- `modules/ena_download_https.py`: HTTPSダウンロード専用CLI、レジューム、`.done` スキップ
- `modules/bwa_mapper.py`: AdapterRemoval、BWA MEM、samtools sort、ancient paired-endのcollapsed/non-collapsed merge
- `modules/softclipper.py`: BAM readの先頭5bp soft clipping
- `modules/bam_processor.py`: Picard CleanSam、sample単位merge、MarkDuplicates、sort、index
- `modules/analyzers.py`: mapDamage、Qualimap、GATK HaplotypeCaller
- `modules/cohort_pca.py`: ancient/modern/autoのcohort PCA/MDS、PLINK/EIGENSOFT連携

## トラブルシューティング

### 起動直後に止まる

環境検証で不足ツールや不足ファイルが見つかっています。ログに表示されたツール名、Picard jar、参照ゲノム、BWA indexを確認してください。

PCA実行時にPLINK/EIGENSOFTが見つからない場合は、[tool_paths.py](tool_paths.py) の固定パスと実際のインストール場所が一致しているか確認してください。

### FASTQが見つからない

`--fastq_dir` のパスと拡張子を確認してください。対象は `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz` です。広すぎる親ディレクトリを指定して別projectが混ざっていないかも確認してください。

### sampleが期待通りにまとまらない

同じsampleの複数runをまとめたい場合は、runフォルダ群の親ディレクトリを `--fastq_dir` に指定してください。runごとに別々に実行すると、sample単位mergeと重複除去が分断されます。

### BAM入力モードでsample名が重複する

たとえば `SAMPLE_A.bam` と `SAMPLE_A.sorted.bam` はどちらも `SAMPLE_A` と解釈されるためエラーになります。`--bam_pattern` で対象を絞るか、片方を別ディレクトリへ移動してください。

### BAM index作成に失敗する

入力BAMが壊れていないか、coordinate sort済みか、参照ゲノムと整合しているかを確認してください。必要なら事前に次を実行してください。

```sh
samtools quickcheck input.bam
samtools index input.bam
```

### PCAでsample/siteが残らない

`pca_qc_summary.tsv` を確認してください。欠損率、MAF、性染色体除外、LD pruningで落ちている可能性があります。まずは `--pca-max-sample-missing` と `--pca-max-site-missing` を緩める、`--pca-min-maf 0.0` にする、`--pca-exclude-sex-chr` を外す、などを試してください。

### smartpcaは動いたが解釈できない

PCAは「近い位置に出るサンプルが全体として似ている」ことを可視化する手法です。直接の祖先関係を証明するものではありません。比較集団の選び方、欠損率、サンプル数、projectionの有無で位置が変わります。2 sampleだけのPCAはパイプライン確認には使えますが、生物学的解釈には不十分です。

### 途中で失敗したsampleだけ再実行したい

同じコマンドを再実行してください。`.done` があるsampleはスキップされ、未完了sampleは既存の中間成果物から再開します。すべてやり直す場合は `--force` を指定します。

### Qualimapで `MaxPermSize` エラーが出る

Java 9以降では `-XX:MaxPermSize` が廃止されています。`QualimapAnalyzer` は `JAVA_TOOL_OPTIONS=-XX:+IgnoreUnrecognizedVMOptions` を自動付与します。
