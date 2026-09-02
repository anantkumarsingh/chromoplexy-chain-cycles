# Centromere and p/q Arm Notes for Baca-AMG Project

This note captures the centromere and chromosome-arm background needed for the next phase of the Baca-AMG prostate cancer project.

## Why This Matters

The current AMG analysis enumerates mathematically possible ways double-strand break endpoints could be rejoined. Biology adds constraints: not every mathematically valid rejoining is likely to survive in a real tumor cell.

Centromeres and p/q chromosome arms help answer two biological questions:

1. Does a rearranged chromosome product have a viable number of centromeres?
2. Are the broken ends spatially plausible rejoining partners based on chromosome arm location?

The AMG framework tells us what is mathematically possible. Centromere and p/q-arm logic helps decide which possibilities are biologically plausible.

## Chromosome Basics

A chromosome is one long DNA molecule packaged with proteins.

Before cell division, each chromosome is copied, producing two identical sister chromatids. These sister chromatids must be separated accurately into daughter cells during mitosis.

The centromere is the chromosomal region where the cell builds the machinery that pulls chromosomes apart during division.

## What Is a Centromere?

A centromere is not usually a single base-pair coordinate. It is a region of DNA plus proteins.

In humans, centromeres are large, repetitive, difficult-to-sequence regions. They are often made of alpha-satellite repeat DNA. In genome-coordinate work, we usually approximate each centromere as an interval.

Example:

```text
chr21 centromere, hg19: roughly 10-13 Mb
```

The centromere's main job is to assemble the kinetochore.

The kinetochore is a protein complex that attaches the chromosome to spindle microtubules during mitosis. A useful mental model:

```text
centromere = chromosomal handle
kinetochore = hook built on that handle
spindle = pulling machinery
```

## p Arm and q Arm

The centromere divides a chromosome into two arms:

- `p arm`: the short arm. The name comes from "petit."
- `q arm`: the long arm. `q` is the next letter after `p`.

Genome coordinates usually run from low coordinate to high coordinate.

For most practical work:

```text
position before centromere -> p arm
position after centromere  -> q arm
```

Example using chr21:

```text
chr21 centromere: approximately 10-13 Mb
chr21:4 Mb       -> p arm
chr21:40 Mb      -> q arm
```

In the Baca prostate cancer data, ERG and TMPRSS2 breakpoints are usually around chr21:39-48 Mb. Therefore, they are on chr21q, far beyond the centromere.

## Why Centromere Count Matters

A rearranged chromosome product can have:

```text
1 centromere  -> monocentric, usually stable
0 centromeres -> acentric, usually unstable or lost
2 centromeres -> dicentric, usually unstable
```

### Monocentric Products

A monocentric derivative chromosome has exactly one centromere.

This is usually the stable case because the chromosome can attach to the mitotic spindle and segregate properly.

### Acentric Products

An acentric fragment has no centromere.

Since it cannot attach properly to spindle microtubules, it is often lost during cell division.

### Dicentric Products

A dicentric chromosome has two centromeres.

During cell division, the two centromeres may attach to opposite spindle poles. The chromosome can be pulled in two directions, forming a bridge that breaks. This can trigger more genome instability.

For our purposes:

```text
monocentric -> biologically plausible
acentric    -> likely unstable
dicentric   -> likely unstable
```

This is why centromere viability is a useful biological filter for AMG matchings and cycle-completion hypotheses.

## Rearrangements and Centromeres

Consider a simple chromosome:

```text
p arm ---- centromere ---- q arm
```

If breaks occur on both sides of the centromere, rearrangement products may include or exclude the centromere depending on how pieces are joined.

If a derivative product contains no centromere, it is acentric.

If a derivative product contains two centromeres, it is dicentric.

If it contains one centromere, it is monocentric.

This means a proposed rejoining is not just about which endpoints connect. We eventually need to trace the resulting chromosome product and count how many centromeres it contains.

## Connection to AMG Cycles

In the current code, each Baca rearrangement row contributes two breakpoints. The AMG code turns those into DSB positions and endpoint nodes.

For actual data, the observed breakpoint strands define observed rearrangement edges.

For counterfactual enumeration, the code generates all possible endpoint pairings and computes the resulting AMG cycle structures.

However, those counterfactual pairings are currently only mathematical. Some pairings may create biologically implausible products:

```text
mathematically valid matching
but derivative product is acentric or dicentric
therefore less biologically viable
```

The centromere filter asks:

```text
Given a matching or cycle completion, how many centromeres does each derivative chromosome product contain?
```

The desired result is usually monocentric derivative products.

## Important Project-Specific Insight

For ERG-chain patients, many DSBs are distal to their chromosome's centromere.

In particular:

```text
chr21 centromere: roughly 10-13 Mb
ERG/TMPRSS2 breakpoints: often 39-48 Mb
```

So ERG/TMPRSS2 breakpoints are on chr21q, far from the centromere.

This is why the existing project notes say the centromere filter has limited impact on some current ERG-chain closed-cycle enumeration tables. Many current chain DSBs are already distal to centromeres.

The centromere filter becomes more important for cycle completion, where open paths may connect to other parts of the genome and could create dicentric or acentric products.

## p/q Arm Proximity

p/q arm proximity is related to centromeres but answers a different question.

Centromere viability asks:

```text
Does the rearranged product have the right number of centromeres?
```

p/q proximity asks:

```text
Were these broken ends physically likely to meet and rejoin?
```

Chromosomes are organized in the nucleus. They are not random strings floating freely. Chromosome territories, arm positions, transcriptional activity, and nuclear architecture influence which regions are physically close.

The professor's reference to radiation data likely refers to observed chromosome-arm exchange frequencies after radiation damage. Such data can reveal which arms tend to exchange more often, probably reflecting 3D nuclear proximity.

This data is probably general human-cell data, not prostate-specific or patient-specific.

## How to Assign p/q Arm in This Project

For each breakpoint:

```text
input: chromosome, breakpoint position, centromere interval
output: p or q arm
```

Simple midpoint version:

```python
if position < centromere_midpoint:
    arm = "p"
else:
    arm = "q"
```

More precise interval version:

```python
if position < centromere_start:
    arm = "p"
elif position > centromere_end:
    arm = "q"
else:
    arm = "centromere"
```

For Baca breakpoint data, approximate centromere intervals are probably enough for a first-pass analysis. A precise version can use hg19 `cytoBand.txt` from UCSC and the bands labeled `acen`.

## Rejoin Arm Categories

Once every endpoint has a p/q arm label, each candidate rejoin can be classified:

```text
p-p
p-q
q-q
same chromosome, same arm
same chromosome, opposite arm
different chromosome, same arm
different chromosome, opposite arm
```

This lets us score or filter observed and hypothetical matchings.

For example:

```text
chr21q joined to chr1q -> interchromosomal q-q
chr21q joined to chr1p -> interchromosomal q-p
chr21q joined to chr21q -> intrachromosomal q-q
```

## ERG/TMPRSS2 Context

ERG and TMPRSS2 are both on chr21q.

Simple TMPRSS2-ERG fusion is therefore a q-arm event on the same chromosome.

Known biology supports this:

- androgen receptor binding can bring TMPRSS2 and ERG regulatory regions into proximity
- TOP2B can generate breaks at these regions
- repair can create the TMPRSS2-ERG fusion

For chromoplexy-embedded ERG cases, chr21q breakpoints connect to other chromosomes. A natural next question is:

```text
Are the partner breakpoints also on q arms?
Are observed joins enriched for q-q pairings?
Are high-rank AMG cycles enriched for particular p/q-arm patterns?
```

This could help connect AMG cycle structure to nuclear architecture.

## How This Connects to the Next Analysis Phase

The next phase should probably be separated into layers.

### Layer 1: Annotate Breakpoints

Add centromere intervals and assign every breakpoint to:

```text
p arm
q arm
centromeric/unknown, if inside the centromere interval
```

For each Baca breakpoint, store:

```text
chromosome
position
strand
gene annotation
arm
distance from centromere
```

### Layer 2: Score Observed Rejoins

For each observed rearrangement edge, classify the arm pairing:

```text
p-p
p-q
q-q
same-arm
opposite-arm
intra-chromosomal
inter-chromosomal
```

Then ask:

```text
Are ERG chromoplexy-embedded patients enriched for q-q joins?
Are simple fusions mostly chr21q-chr21q?
Do observed joins differ from random counterfactual matchings?
```

### Layer 3: Score Counterfactual Matchings

For enumerable patients, apply the same p/q classification to every hypothetical matching.

Possible outputs:

```text
cycle structure
number of p-p joins
number of q-q joins
number of p-q joins
same-arm fraction
opposite-arm fraction
driver gene included in cycle?
```

This can turn the existing enumeration tables into biologically weighted tables.

### Layer 4: Cycle Completion

Open paths exist because the ERG chain is only a subset of the full patient genome.

Cycle completion means asking how those open paths might close.

There are two possible interpretations:

1. Empirical cycle completion:
   Look in the full patient CSV for where open path endpoints actually connect outside the extracted ERG chain.

2. Theoretical cycle completion:
   Enumerate biologically plausible ways open paths could close, subject to centromere, p/q-arm, and driver-gene constraints.

This needs clarification before implementation.

### Layer 5: Centromere Viability

For candidate cycle completions, trace derivative products and count centromeres.

Flag each product:

```text
monocentric
acentric
dicentric
uncertain
```

This is harder than p/q labeling because it requires tracking chromosome segments, not just endpoint labels.

## Practical Mental Model

For every breakpoint:

```text
1. What chromosome is it on?
2. What coordinate is it at?
3. Is it before or after the centromere?
4. Therefore, is it on p or q arm?
5. Which gene or region is nearby?
```

For every rejoin:

```text
1. Which two arms are being joined?
2. Is it intra-chromosomal or inter-chromosomal?
3. Is it same-arm or opposite-arm?
4. Does it involve a driver gene?
5. Is it plausible under proximity assumptions?
```

For every completed derivative product:

```text
1. How many centromeres does it contain?
2. Is it monocentric, acentric, or dicentric?
3. Would a tumor cell likely retain it?
```

## Current Open Questions

These are the key questions to clarify before coding the next phase:

1. Should centromere and p/q-arm constraints filter existing counterfactual enumeration tables, or should they mainly guide a new cycle-completion analysis?
2. Which radiation exchange-frequency dataset should be used for p/q-arm proximity weights?
3. Should the analysis be empirical, using full patient CSV rows to close open paths, or theoretical, enumerating plausible completions?
4. Which genes count as driver/relevant genes: only the Baca progression genes, or a broader cancer gene set?
5. Should chromoplexy versus small-cycle structure be asked at the ERG-chain level or separately for each driver gene?
6. Do Sheth's proper-AMG bounds still apply meaningfully to chain subsets, or do we need full-genome AMG reconstruction for this part?

## First Implementation Target

The safest first step is p/q-arm annotation, because it is directly derivable from existing data.

Minimum first-pass implementation:

```text
1. Add hg19 centromere intervals.
2. Assign p/q arm to every breakpoint in chrom_aberrations_baca.csv.
3. Add arm-pair labels to observed rearrangements.
4. Summarize arm-pair frequencies for:
   - all Baca patients
   - ERG+ patients
   - chromoplexy_embedded ERG+ patients
   - simple_fusion ERG+ patients
5. For enumerable AMG patients, compare observed arm-pair pattern against counterfactual matchings.
```

Centromere viability should come after this, because true viability requires derivative-product tracing.

