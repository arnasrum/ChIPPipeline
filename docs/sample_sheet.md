
# Specify samples 

The pipeline requires a sample sheet (defaulting to config/samples.csv) to define the input data and metadata for each sample. This file should be in CSV format and include specific columns as detailed below.

<table>
    <th>Column</th>
    <th>Description</th>
    <th>Required</th>
    <th>Valid Values</th>
    <tr>
        <td>Mark</td>
        <td>Identifier for transcription factor or histone mark for treatment file.</td>
        <td>Yes, only for treatment samples. Leave blank for control samples.</td>
        <td>String</td>
    </tr>
    <tr>
        <td>Sample</td>
        <td>The biological sample or condition the sequences were derived from.</td>
        <td>Yes, used for associating samples with corresponding sample origin.</td>
        <td>String</td>
    </tr>
    <tr>
        <td>Type</td>
        <td>Specifies whether the sample is treatment or control.</td> 
        <td>Yes</td>
        <td>treatment/control</td>
    </tr>
    <tr>
        <td>Peak_type</td>
        <td>Must be narrow/broad, decides if the peak caller will treat the sample as narrow or broad peaks.</td>
        <td>Yes, for treatment files. Can be omitted for control files.</td>
        <td>narrow/broad</td>
    </tr>
    <tr>
        <td>Accession</td>
        <td>GEO accession number for publicly available samples, used to download the sample from GEO if a local file path is not specified. If file_path is specified, this column serves as a prefix for the filename.</td>
        <td>Yes; if file_path is not defined.</td>
        <td>String</td>
    </tr>    
    <tr>
        <td>Genome</td>
        <td>The genome that the sample will the aligned to. Samples will also be separated by the genome they are aligned to.</td>
        <td>Yes</td>
        <td>Genome identifier or path to a FASTA or gzipped FASTA file named after the identifier. Example; hg19 or path/to/genome/hg19.fa.gz</td>
    </tr>
    <tr>
        <td>File_path</td>
        <td>Paths to the reads sample reads. If paired_end is set to true, then two paths must be defined in the column by separating them with ";" character. Files without file_path specified will try </td>
        <td>No, if not specified the pipeline will try to download the sample based on the GEO accession.</td>
        <td>If paired_end is false; path/to/read.fq.gz. If paired_end is true: path/to/read1.fq.gz;path/to/read2.fq.gz </td>
    </tr>
        
</table>


## Example sheet

An example sheet showing the possibilities:

```csv
type,peak_type,replicate,accession,file_path,genome,mark,sample,paired_end
treatment,narrow,1,GSM1234567,,genomes/hg19.fa,H3K4me3,CellLineA,True
treatment,narrow,2,GSM1234568,,genomes/hg19.fa,H3K4me3,CellLineA,True
control,,,1,GSM1234569,,genomes/hg19.fa,IgG,CellLineA,True
control,,,2,GSM1234570,,genomes/hg19.fa,IgG,CellLineA,True
treatment,broad,1,,data/sample1_R1.fastq.gz;data/sample1_R2.fastq.gz,genomes/mm10.fa,H3K36me3,TissueX,True
treatment,broad,2,,data/sample2_R1.fastq.gz;data/sample2_R2.fastq.gz,genomes/mm10.fa,H3K36me3,TissueX,True
control,,,1,,data/control1.fastq.gz,genomes/mm10.fa.gz,Input,TissueX,False
control,,,2,,data/control2.fastq.gz,genomes/mm10.fa.gz,Input,TissueX,False
treatment,narrow,1,GSM9876543,,dm6,H3K27me3,FlySampleY,False
```
Or in JSON format:

```json
[
  {
    "type": "treatment",
    "sample": "CellLineA",
    "replicate": 1,
    "mark": "H3K4me3",
    "peak_type": "narrow",
    "genome": "genomes/hg19.fa",
    "paired_end": "true",
    "accession": "GSM1234567",
    "file_path": ""
  },
  {
    "type": "treatment",
    "sample": "CellLineA",
    "replicate": 2,
    "mark": "H3K4me3",
    "peak_type": "narrow",
    "genome": "genomes/hg19.fa",
    "paired_end": "true",
    "accession": "GSM1234568",
    "file_path": ""
  },
  {
    "type": "control",
    "sample": "CellLineA",
    "replicate": 1,
    "mark": "IgG",
    "peak_type": "",
    "genome": "genomes/hg19.fa",
    "paired_end": "true",
    "accession": "GSM1234569",
    "file_path": ""
  },
  {
    "type": "control",
    "sample": "CellLineA",
    "replicate": 2,
    "mark": "IgG",
    "peak_type": "",
    "genome": "genomes/hg19.fa",
    "paired_end": "true",
    "accession": "GSM1234570",
    "file_path": ""
  },
  {
    "type": "treatment",
    "sample": "TissueX",
    "replicate": 1,
    "mark": "H3K36me3",
    "peak_type": "broad",
    "genome": "genomes/mm10.fa",
    "paired_end": "true",
    "accession": "",
    "file_path": "data/sample1_R1.fastq.gz;data/sample1_R2.fastq.gz",
  },
  {
    "type": "treatment",
    "sample": "TissueX",
    "replicate": 2,
    "mark": "H3K36me3",
    "peak_type": "broad",
    "genome": "genomes/mm10.fa",
    "paired_end": "true",
    "accession": "",
    "file_path": "data/sample2_R1.fastq.gz;data/sample2_R2.fastq.gz"
  },
  {
    "type": "control",
    "sample": "TissueX",
    "replicate": 1,
    "mark": "Input",
    "peak_type": "",
    "genome": "genomes/mm10.fa",
    "paired_end": "false",
    "accession": "",
    "file_path": "data/control1.fastq.gz"
  },
  {
    "type": "control",
    "sample": "TissueX",
    "replicate": 2,
    "mark": "Input",
    "peak_type": "",
    "genome": "genomes/mm10.fa",
    "paired_end": false,
    "accession": "",
    "file_path": "data/control2.fastq.gz"
  },
  {
    "type": "treatment",
    "sample": "FlySampleY",
    "replicate": 1,
    "mark": "H3K27me3",
    "peak_type": "narrow",
    "genome": "dm6",
    "paired_end": false,
    "accession": "GSM9876543",
    "file_path": ""
  }
]
```