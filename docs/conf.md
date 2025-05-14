# Defining Configuration Options

Defining configuration options can be done either via CLI;

`snakemake <args> -C options=value`

or:

`snakemake <args> --config option=value`

Or by editing the [configuration file](../config/config.yml) directly

# Configuration Options

<table>
  <thead>
    <tr>
      <th>Option</th>
      <th>Type/Value Example</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td colspan="3"><strong>Workflow Configuration (Overall Pipeline Settings)</strong></td>
    </tr>
    <tr>
      <td><code>modules</code></td>
      <td>String (e.g., <code>""</code>, <code>"-a bowtie2 --trimmer fastp"</code>)</td>
      <td>Specifies which specific parts or modules of the workflow should be executed. Empty means run default/all.</td>
    </tr>
    <tr>
      <td><code>trimmer</code></td>
      <td>String (<code>"cutadapt"</code>, <code>"fastp"</code>, <code>"trim_galore"</code> or <code>"trimmomatic"</code>)</td>
      <td>Selects the tool to be used for trimming raw sequencing reads.</td>
    </tr>
    <tr>
      <td><code>aligner</code></td>
      <td>String (<code>"bowtie2"</code>, <code>"bwa_mem"</code>, <code>"bwa_mem2"</code> or <code>"star"</code>)</td>
      <td>Selects the tool to be used for aligning processed reads to a reference genome.</td>
    </tr>
    <tr>
      <td><code>duplicate_processor</code></td>
      <td>String (<code>"mark_duplicates"</code> or <code>"markdup"</code>)</td>
      <td>Selects the tool used to identify and/or mark duplicate reads.</td>
    </tr>
    <tr>
      <td><code>peak_caller</code></td>
      <td>String (<code>"macs3"</code>)</td>
      <td>Selects the tool used for identifying genomic regions enriched with reads (peaks).</td>
    </tr>
    <tr>
      <td><code>generate_fastqc_reports</code></td>
      <td>Boolean (<code>true</code> or <code>false</code>)</td>
      <td>Enables or disables the generation of FastQC quality control reports.</td>
    </tr>
    <tr>
      <td><code>plot_region</code></td>
      <td>String (e.g., <code>"chr1:14,000,000-15,000,000"</code>)</td>
      <td>Defines a specific genomic region for visualizations or focused analysis.</td>
    </tr>
    <tr>
      <td><code>sample_sheet</code></td>
      <td>String (e.g., <code>"config/samples.csv"</code>)</td>
      <td>Path to the CSV or JSON file containing sample metadata and file information.</td>
    </tr>
    <tr>
      <td><code>outdir</code></td>
      <td>String (e.g., <code>""</code>, <code>"/path/to/output"</code>)</td>
      <td>Base directory for all generated output files.</td>
    </tr>
    <tr>
      <td><code>results_path</code></td>
      <td>String (e.g., <code>"results"</code>)</td>
      <td>Relative path within <code>outdir</code> for final processed results.</td>
    </tr>
    <tr>
      <td><code>resources_path</code></td>
      <td>String (e.g., <code>"resources"</code>)</td>
      <td>Relative path within <code>outdir</code> for generated resource files (e.g., genome indices).</td>
    </tr>
    <tr>
      <td><code>logs_path</code></td>
      <td>String (e.g., <code>"logs"</code>)</td>
      <td>Relative path within <code>outdir</code> for workflow log files.</td>
    </tr>
    <tr>
      <td><code>benchmarks_path</code></td>
      <td>String (e.g., <code>"benchmarks"</code>)</td>
      <td>Relative path within <code>outdir</code> for performance benchmark data.</td>
    </tr>
    <tr>
      <td><code>temp_path</code></td>
      <td>String (e.g., <code>"temp"</code>)</td>
      <td>Relative path within <code>outdir</code> for temporary files.</td>
    </tr>
    <tr>
      <td><code>enable_homer_motif_analysis</code></td>
      <td>Boolean (<code>true</code> or <code>false</code>)</td>
      <td>Decides of HOMER rules for motif enrichment analysis are run.</td>
    </tr>
    <tr>
      <td colspan="3"><strong>Tool Configurations (Parameters per External Tool)</strong></td>
    </tr>
     <tr>
      <td colspan="3">This section contains blocks named after specific tools (e.g., <code>fastp:</code>, <code>bowtie2:</code>, <code>macs3:</code>). Each tool block can define the following parameters:</td>
    </tr>
    <tr>
      <td><code>&lt;tool_name&gt;</code> -> <code>args</code></td>
      <td>String (e.g., <code>""</code>, <code>"-v 1 -Y"</code>, <code>"--skip-technical"</code>)</td>
      <td>Additional, tool-specific command-line arguments or options to pass during execution.</td>
    </tr>
    <tr>
      <td><code>&lt;tool_name&gt;</code> -> <code>threads</code></td>
      <td>Integer (e.g., <code>12</code>, <code>25</code>)</td>
      <td>Number of CPU threads/cores to allocate for the tool's process.</td>
    </tr>
    <tr>
      <td><code>&lt;tool_name&gt;</code> -> <code>runtime</code></td>
      <td>Integer (e.g., <code>20</code>, <code>60</code>, <code>200</code>)</td>
      <td>Estimated or maximum allowed runtime for the tool in minutes.</td>
    </tr>
    <tr>
      <td><code>&lt;tool_name&gt;</code> -> <code>mem_mb</code></td>
      <td>Integer (e.g., <code>4000</code>, <code>8000</code>, <code>64000</code>)</td>
      <td>Maximum memory (RAM) in megabytes to allocate for the tool's process.</td>
    </tr>
     <tr>
      <td colspan="3">Some tools have unique parameters specific to their functionality:</td>
    </tr>
     <tr>
      <td><code>trimmomatic</code> -> <code>trimming_steps</code></td>
      <td>String (e.g., <code>"ILLUMINACLIP:..."</code>)</td>
      <td>Specific trimming parameters and steps string passed directly to Trimmomatic's command.</td>
    </tr>
     <tr>
      <td><code>computeMatrix</code> -> <code>mode</code></td>
      <td>String (e.g., <code>"reference-point"</code>, <code>"scale-regions"</code>)</td>
      <td>Calculation mode for the deepTools computeMatrix tool.</td>
    </tr>
    <tr>
      <td><code>pyGenomeTracks</code> -> <code>bigwig_options</code></td>
      <td>String (e.g., <code>"max_value = 6,color = red"</code>)</td>
      <td>Comma-separated options appended to bigWig track configurations in pyGenomeTracks.</td>
    </tr>
    <tr>
      <td><code>pyGenomeTracks</code> -> <code>bed_options</code></td>
      <td>String (e.g., <code>"show_labels=true,color=red"</code>)</td>
      <td>Comma-separated options appended to narrowPeak(BED) track configurations in pyGenomeTracks. The extensive options can be found at the offical pyGenomeTracks documentation.</td>
    </tr>
  </tbody>
</table>

## Module Options

Overwrite the default configuration by running:

`snakemake -c <num> <rule> --config modules="--flag1 <value1> -f <value2>"`

Example:

`snakemake -c <num> <rule> --config modules="-t trimmomatic --aligner bwa_mem2"`

### Flags

**-t** or **--trimmer** overwrites the trimmer

**-a** or **--aligner** overwrites the aligner

**-d** or **--duplicate_processor** overwrites the duplicate processor

**-p** or **--peak_caller** overwrites the duplicate processor
