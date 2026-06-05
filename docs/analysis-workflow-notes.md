# 解析ワークフロー理解メモ

このメモは、`analyze-fastq-app` の解析部分の内容をすぐ思い出せるようにまとめたものです。

## まず押さえる全体像

このプロジェクトの解析は、FASTQを入力として次の順に進みます。

```text
FASTQ
  ↓ AdapterRemoval
  ↓ BWA MEM による参照ゲノムへのマッピング
  ↓ Soft clipping
  ↓ Picard CleanSam
  ↓ run単位BAMをsample単位にmerge
  ↓ Picard MarkDuplicates
  ↓ samtools sort / index
  ↓ mapDamage
  ↓ Qualimap
  ↓ GATK HaplotypeCaller
VCF / QC結果
```

大事なのは、前半は **run単位**、後半は **sample単位** で処理されることです。

## run単位とは

`run` は、FASTQの実行単位、レーン単位、またはサンプル内の個別データまとまりに近い粒度です。

例:

```text
SAMPLE_A_L001_R1.fastq.gz
SAMPLE_A_L001_R2.fastq.gz
SAMPLE_A_L002_R1.fastq.gz
SAMPLE_A_L002_R2.fastq.gz
```

この場合、`SAMPLE_A` という1つのsampleに対して、`L001` と `L002` の2つのrunがあると考えます。

run単位では、各runを独立に処理します。

```text
run FASTQ
  ↓ AdapterRemoval
  ↓ BWA mapping
  ↓ Soft clipping
  ↓ CleanSam
run.clean.bam
```

## sample単位とは

`sample` は、最終的に結果をまとめたい生物学的サンプル単位です。

複数runがある場合は、runごとに作られたBAMをsample単位に統合します。

```text
run1.clean.bam
run2.clean.bam
run3.clean.bam
  ↓ samtools merge
sample.merged.bam
  ↓ MarkDuplicates
  ↓ sort / index
sample.dedup.sorted.bam
```

その後、sample単位のBAMに対して `mapDamage`, `Qualimap`, `HaplotypeCaller` を実行します。

## AdapterRemoval

FASTQリードをマッピング前にきれいにする前処理です。

主な役割:

- アダプター配列の除去
- 低品質末端のトリミング
- `N` を含む末端のトリミング
- 短すぎるリードの除外
- 古DNAのpaired-endでは、重なったR1/R2を1本に結合する `--collapse`

このプロジェクトでは `--data_type` によって設定が変わります。

```text
ancient:
  minquality = 20
  minlength = 30
  paired-endでは --collapse を使う

modern:
  minquality = 25
  minlength = 25
  collapseなし
```

## 参照ゲノムにマッピングとは

FASTQの短いreadが、基準となるゲノム配列のどこに由来するかを探して配置することです。

例:

```text
参照ゲノム:
... A C T G A A C C T G A T T C G A ...

read:
A A C C T G A T

マッピング後:
... A C T G A A C C T G A T T C G A ...
              | | | | | | | |
              A A C C T G A T
```

FASTQだけでは、各readがどの染色体のどこから来たか分かりません。マッピングすると、readごとに位置情報が付きます。

```text
FASTQ
  ↓ BWA MEM
BAM
```

このプロジェクトでは `BWA MEM` を使い、次の流れでsorted BAMを作ります。

```text
bwa mem
  ↓
samtools view -b
  ↓
samtools sort
  ↓
*.sorted.bam
```

## bwa mem / samtools view / samtools sort の役割

`bwa mem` は、FASTQのreadを参照ゲノムにマッピングするコマンドです。

入力はFASTQと参照ゲノムで、出力は通常SAM形式のアラインメントです。

```text
FASTQ reads
  ↓ bwa mem
SAM alignment
```

paired-endの場合はR1/R2を渡します。

```sh
bwa mem reference.fa sample_R1.fastq.gz sample_R2.fastq.gz
```

single-endの場合はFASTQを1本だけ渡します。

```sh
bwa mem reference.fa sample.fastq.gz
```

このプロジェクトでは、`bwa mem` の出力をファイルにSAMとして保存せず、次の処理へパイプで渡します。

```text
bwa mem
  ↓ SAM出力
samtools view -b
  ↓ BAMへ変換
samtools sort
  ↓ 座標順に並べ替え
sorted BAM
```

`samtools view` は、SAM/BAM/CRAMを表示・変換・抽出するためのコマンドです。

このパイプラインでは主に `-b` を使い、`bwa mem` が出したSAMをBAMに変換しています。

```sh
samtools view -b -
```

`-` は標準入力を意味します。つまり、直前の `bwa mem` から流れてきたSAMを受け取ります。

`samtools sort` は、BAMをゲノム座標順に並べ替えるコマンドです。

```sh
samtools sort -o sample.sorted.bam -
```

多くの後続ツールは、座標順にsortされたBAMを前提にします。さらに、`samtools index` で `.bai` を作る場合も、通常は座標順sort済みBAMが必要です。

SAMをBAMに変換する意味は、同じアラインメント情報をより小さく、高速に扱える形式にすることです。

```text
SAM:
  テキスト形式
  人間が読める
  サイズが大きい
  大量データ処理には重い

BAM:
  SAMをバイナリ圧縮した形式
  サイズが小さい
  samtools / GATK / Picard / Qualimap などが効率よく扱える
  sort / index と組み合わせると特定領域を高速に取り出せる
```

そのため、解析パイプラインではSAMは一時的な出力、BAMが実際の処理対象になることが多いです。

## ancient paired-end の特殊処理

`--data_type ancient` かつpaired-endの場合、AdapterRemoval後の出力を2系統に分けて扱います。

```text
collapsed reads
  ↓ single-end BWA
collapsed.sorted.bam

non-collapsed R1/R2
  ↓ paired-end BWA
pe.sorted.bam

collapsed.sorted.bam + pe.sorted.bam
  ↓ samtools merge
run.merged.sorted.bam
```

古DNAでは断片が短く、R1/R2が重なることがあるため、collapsed readsを別扱いしています。

## Soft clipping

BAM内の各readに対して、先頭側5bpをsoft clipします。

soft clipとは、配列を完全に削除するのではなく、BAM上で「この部分はアラインメントに使わない」とCIGARに記録する処理です。

```text
元のread:
MMMMMMMMMMMM

soft clip後:
SSSSSMMMMMMM
```

このプロジェクトでは `clip_length = 5` です。リード端の損傷や不安定なアラインメントの影響を抑える意図と考えられます。

## Picard CleanSam

run単位のsoftclipped BAMを、Picardの `CleanSam` で整えます。

主な目的:

- SAM/BAM内の軽微な不整合を補正
- 後続ツールが読みやすいBAMにする
- `VALIDATION_STRINGENCY=LENIENT` で実行

```text
run_softclipped.bam
  ↓ Picard CleanSam
run.clean.bam
```

## samtools merge

複数runがあるsampleでは、runごとのclean BAMを1つにまとめます。

```text
SAMPLE_A/run1.clean.bam
SAMPLE_A/run2.clean.bam
  ↓ samtools merge
SAMPLE_A.merged.bam
```

runが1つだけの場合は、mergeせずにそのBAMを移動して使います。

## Picard MarkDuplicates

sample単位にマージしたBAMから重複readを検出し、除去します。

このプロジェクトでは:

```text
REMOVE_DUPLICATES=true
```

なので、重複をマークするだけでなく、実際に取り除きます。

出力:

```text
sample.marked.bam
sample.marked_dup_metrics.txt
```

`marked_dup_metrics.txt` には重複率などのメトリクスが出ます。

## samtools sort / index

重複除去後のBAMを座標順に並べ、インデックスを作成します。

```text
sample.marked.bam
  ↓ samtools sort
sample.dedup.sorted.bam
  ↓ samtools index
sample.dedup.sorted.bam.bai
```

この `sample.dedup.sorted.bam` が後続解析の中心入力になります。

## mapDamage

古DNAらしい損傷パターンを評価するツールです。

このプロジェクトでは、mapDamageの前にBAMをフィルタします。

```text
条件:
  mapQ >= 30
  POS >= 300
```

流れ:

```text
sample.dedup.sorted.bam
  ↓ samtools view -q 30
  ↓ POS >= 300 でフィルタ
  ↓ samtools sort / index
  ↓ mapDamage
mapDamage結果
```

mapDamageでは、古DNAでよく見られる末端のC→TやG→Aの傾向などを確認できます。

## Qualimap

BAMの品質評価を行うツールです。

このプロジェクトでは `qualimap bamqc` を実行し、HTML形式で結果を出します。

事前チェック:

- BAMが存在するか
- BAMが空でないか
- mapped read数が0でないか

見る指標の例:

- coverage
- mapping quality
- insert size
- GC content
- mapped reads

## GATK HaplotypeCaller

変異検出を行い、VCFを出力するツールです。

入力:

```text
sample.dedup.sorted.bam
reference genome
```

出力:

```text
sample.vcf
```

このプロジェクトでは主に次の設定で実行されます。

```text
--output-mode EMIT_VARIANTS_ONLY
-stand-call-conf 30
--native-pair-hmm-threads <threads>
```

つまり、変異と判定された箇所のみをVCFに出します。

## chr1とは

`chr1` は `chromosome 1`、つまり第1染色体を表す名前です。

マッピング結果で:

```text
chr1: 10523-10531
```

とあれば:

```text
第1染色体の10523番目から10531番目あたりにreadがマッピングされた
```

という意味です。

ただし、染色体名は参照ゲノムによって異なります。

```text
chr1
1
NC_000001.11
CM000663.2
```

これらは参照ゲノムによっては同じ第1染色体に相当することがあります。

## nochrUnとは

`nochrUn` は、おそらく `chrUn` を除いた参照ゲノムという意味です。

`chrUn` は、配列としては分かっているが、どの染色体のどの位置に属するかが確定していない配列を表すことがあります。

```text
chrUn
chrUn_KI270442v1
chrUn_GL000220v1
```

`equCab3.nochrUn.fa` という名前なら、馬のEquCab3参照ゲノムから `chrUn` 系の未配置配列を除いたFASTA、という意味合いだと考えられます。

`chrUn` を除く理由:

- 主要染色体だけで解析したい
- contig数を減らして結果を扱いやすくしたい
- readが未配置配列へ分散するのを避けたい
- mapDamageやVCF出力の解釈を単純にしたい

## ひとことでまとめると

```text
run単位:
  FASTQを前処理して、runごとのclean BAMを作る

sample単位:
  runごとのBAMをまとめて、重複除去・QC・VCF出力を行う

参照ゲノムへのマッピング:
  FASTQ readがゲノム上のどこに由来するかを決める処理

chr1:
  第1染色体を表す名前

nochrUn:
  未配置配列 chrUn を除いた参照ゲノム、という意味合いのファイル名
```
