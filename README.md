# seeksv
an accurate tool for structural variation and virus integration detection

## Usage
### Step 1
Get soft-clipped reads from original bam file. The bam file should be sorted by alignment coordinate. If the bam is marked duplicates by picard, duplicate reads in the bam would be ignored by seeksv.
#### Inputs
- input.bam
#### Command
```
seeksv getclip -o /path/to/output_prefix input.bam
```
#### Outputs
- `output_prefix.clip.fq.gz`
- `output_prefix.clip.gz`
- `output_prefix.unmapped_1.fq.gz`
- `output_prefix.unmapped_2.fq.gz`.
### Step 2
```

```
### Step 3


