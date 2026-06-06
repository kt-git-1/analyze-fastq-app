# run/sample逆順FASTQ構成の扱い

このメモは、複数のsequencing runフォルダの下にsample FASTQが並ぶ構成を、あとで思い出せるように整理したものです。

## どんな構成か

例:

```text
ancient/
  210416_A01210_0064_BHWKFYDMXX/
    PE-AncientHorses-01_S1_L001_R1_001.fastq.gz
    PE-AncientHorses-01_S1_L001_R2_001.fastq.gz
    PE-AncientHorses-02_S2_L001_R1_001.fastq.gz
    PE-AncientHorses-02_S2_L001_R2_001.fastq.gz
  210419_A00581_0145_AHWKL3DMXX/
    PE-AncientHorses-01_S1_L002_R1_001.fastq.gz
    PE-AncientHorses-01_S1_L002_R2_001.fastq.gz
```

この場合、最上位の `210416_A01210_0064_BHWKFYDMXX` はsampleではなくsequencing runを表す。

ファイル名側の `PE-AncientHorses-01_S1` や `PE-AncientHorses-02_S2` がsampleで、`L001` / `L002` がlane、`R1` / `R2` がpaired-end readを表す。

## 正しい解釈

上の例は、内部的に次のように正規化して扱う。

```text
sample = PE-AncientHorses-01_S1
  run = 210416_A01210_0064_BHWKFYDMXX_L001
    FASTQ = R1/R2 pair
  run = 210419_A00581_0145_AHWKL3DMXX_L002
    FASTQ = R1/R2 pair

sample = PE-AncientHorses-02_S2
  run = 210416_A01210_0064_BHWKFYDMXX_L001
    FASTQ = R1/R2 pair
```

重要なのは、フォルダ構成が `run/sample` の順番でも、解析内部では必ず `sample -> run/lane -> FASTQ` に並べ替えること。

## なぜ親フォルダを指定するのか

同じsampleが複数runにまたがる場合、runフォルダごとに別々に実行すると、sample単位の統合ができない。

非推奨:

```sh
python main.py --fastq_dir ancient/210416_A01210_0064_BHWKFYDMXX --data_type ancient
python main.py --fastq_dir ancient/210419_A00581_0145_AHWKL3DMXX --data_type ancient
```

この運用では、同じ `PE-AncientHorses-01_S1` のデータがrunごとに分かれて処理されるため、sample全体でのmerge、MarkDuplicates、VCF出力にならない。

推奨:

```sh
python main.py --fastq_dir ancient --data_type ancient
```

`ancient/` を指定すると、複数runフォルダをまとめて探索し、同じsampleに属するrun/lane FASTQを同じsampleの下へ集めてから解析する。

## パイプラインでの流れ

```text
runフォルダ群の親
  ↓
group_fastqs_by_run()
  ↓
sampleごとにFastqRunを集約
  ↓
FastqRunごとに AdapterRemoval / BWA / soft clipping / CleanSam
  ↓
sample単位でclean BAMをmerge
  ↓
MarkDuplicates
  ↓
mapDamage / Qualimap / HaplotypeCaller
  ↓
sampleごとに <sample>.vcf
```

## 実装上の判定

`group_fastqs_by_run()` は、次の条件を満たすFASTQを `run/sample` 逆順構成として扱う。

- `--fastq_dir` から見た相対パスが `run_folder/fastq_name` の2階層
- FASTQ名が `SAMPLE_L001_R1_001.fastq.gz` のようなIllumina lane/read形式
- 親フォルダ名が、ファイル名から推定されるsample名と一致しない

この場合:

```text
sample_acc = ファイル名から推定したsample名
run_id     = run_folder名 + "_" + lane
```

例:

```text
210416_A01210_0064_BHWKFYDMXX/PE-AncientHorses-01_S1_L001_R1_001.fastq.gz
```

は:

```text
sample_acc = PE-AncientHorses-01_S1
run_id     = 210416_A01210_0064_BHWKFYDMXX_L001
```

になる。

## 既存の同一サンプル複数FASTQ対応との関係

この対応は、`SAMEA.../A.fastq.gz`, `SAMEA.../B.fastq.gz` のような同一sample配下の複数single FASTQ対応とは別の入口判定である。

既存の同一sample複数FASTQ対応では、sampleディレクトリ名をsampleとして使い、各FASTQを別runとして扱う。

一方、run/sample逆順構成では、親フォルダ名をrunとして使い、sample名はFASTQファイル名から取る。

どちらの場合も最終的な方針は同じで、複数run/lane/FASTQをsample単位でまとめ、最終VCFはsampleごとに1つ出力する。
