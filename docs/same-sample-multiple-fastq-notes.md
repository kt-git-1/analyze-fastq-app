# 同一サンプルに複数FASTQがある場合の扱い

このメモは、`SAMEA...` のような同じサンプルディレクトリ配下に複数のFASTQがある場合、このパイプラインでどう扱うべきかを忘れないための整理です。

## 結論

同じ `SAMEA...` 配下の複数FASTQは、別サンプルではなく、同じ生物学的サンプルに属する複数のrun、lane、library、またはシーケンス単位として扱う。

そのため、single-end FASTQが複数ある場合でも、最初に `cat` で1本に結合しない。FASTQごとにrun/read groupを分けてBAM化し、同じ `SAMEA...` のBAMを最後に1つへmergeしてからMarkDuplicates以降を実行する。

## 処理の流れ

例:

```text
PRJEB19970/
  SAMEA103910521/
    A.fastq.gz
    B.fastq.gz
    C.fastq.gz
```

この場合、`A.fastq.gz`, `B.fastq.gz`, `C.fastq.gz` はすべて `SAMEA103910521` のデータとして扱う。ただし解析単位としては、それぞれ別runに分ける。

```text
A.fastq.gz
  -> AdapterRemoval
  -> BWA MEM
  -> Soft clipping
  -> CleanSam
  -> A由来のclean BAM

B.fastq.gz
  -> AdapterRemoval
  -> BWA MEM
  -> Soft clipping
  -> CleanSam
  -> B由来のclean BAM

C.fastq.gz
  -> AdapterRemoval
  -> BWA MEM
  -> Soft clipping
  -> CleanSam
  -> C由来のclean BAM

A/B/C由来のclean BAM
  -> samtools merge
  -> Picard MarkDuplicates
  -> samtools sort / index
  -> mapDamage
  -> Qualimap
  -> GATK HaplotypeCaller
  -> SAMEA103910521.vcf
```

## なぜcatしないのか

複数FASTQは、別の日付、別lane、別library、別インデックス、別シーケンス装置由来の可能性がある。

最初にFASTQを `cat` で結合すると、どのreadがどのrun/library由来かをBAM上で区別しにくくなる。read groupを分けてBAM化しておくと、以下の点で解釈しやすい。

- runやlibraryごとのQCを追いやすい
- duplicate metricsの解釈がしやすい
- 問題のあるrunだけを後で特定しやすい
- GATK/Picard系ツールにsample名とrun/library情報を正しく渡せる

## read groupの考え方

BWAでBAM化するときは、各runにread groupを付ける。

```text
ID = run ID
SM = sample名、例: SAMEA103910521
LB = library名、例: SCY1.1 / SCY1.2 / unknown
CN = sequencing center
PL = ILLUMINA
```

重要なのは、`SM` は同じ `SAMEA...` に揃え、`ID` や `LB` はFASTQ/run/libraryごとに分けること。

## このパイプラインでの実装方針

現在の推奨方針は次の通り。

1. `group_fastqs_by_run()` で同じsample配下のFASTQをすべて検出する
2. single-end FASTQが複数ある場合は、FASTQごとに別 `FastqRun` として扱う
3. 各 `FastqRun` に対して AdapterRemoval -> BWA -> Soft clipping -> CleanSam を実行する
4. 同じsampleのclean BAMをsample単位でmergeする
5. merge後のBAMに対してMarkDuplicatesを実行する
6. sample単位のdedup BAMからmapDamage、Qualimap、HaplotypeCallerを実行する

つまり、`SAMEA...` 配下にFASTQが10本あれば、原則として10個のrun BAMを作り、それらを1つのsample BAMにまとめてから後続解析を行う。

## 注意点

paired-endの場合はR1/R2の対応を崩さない。R1/R2ペアは1つのrunとして扱う。

片側だけのR1またはR2がある場合は、無理にsingle-endとして混ぜず、warningを出して除外する方が安全。

最終的なVCFはFASTQごとではなく、sampleごとに1つ作る。例えば `SAMEA103910521` なら、最終的には `SAMEA103910521.vcf` を作る。
