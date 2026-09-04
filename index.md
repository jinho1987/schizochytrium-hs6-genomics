---
layout: default
title: Schizochytrium HS6 — Comparative Genomics
---

<section class="hero">
  <div class="eyebrow"><span class="dot"></span>Schizochytrium sp. &middot; Comparative Genomics &amp; Transcriptomics</div>
  <h1 class="hero-title">Two comparisons, one <em>engineered</em> strain</h1>
  <p class="hero-sub">HS6, an EMS-mutagenized strain bred for higher DHA/lipid yield and faster growth, examined two ways: against its own unmutagenized parent (<b>#4</b>, WT) to find what mutagenesis changed, and against an independently isolated outgroup (<b>#7</b>) to find what's structurally distinctive about HS6's genome on its own terms.</p>
  <div class="strain-tags">
    <span class="tag hs6">HS6 &middot; EMS mutant</span>
    <span class="tag">#4 &middot; WT parent</span>
    <span class="tag s7">#7 &middot; independent isolate</span>
  </div>
</section>

<div class="callout caveat">
  <span class="lbl">Statistical caveat — read this first</span>
  <p style="margin:0;">The RNA-seq design is 2 genotype (WT #4, HS6) &times; 3 timepoint (20h/44h/68h fermentation), <b>no biological replicates (n=1)</b>. DESeq2 results are internally consistent, but a single divergent culture could in principle produce the same signal. Every finding below is explicitly tiered by evidence strength — see <a href="#limitations">Limitations &amp; evidence tiers</a>.</p>
</div>

<section class="block">
  <span class="section-eyebrow">Two projects, one repository</span>
  <h2>This site bundles two distinct comparisons</h2>
  <p>They share a strain (HS6) and a genome build, but ask different questions with different methods — worth telling apart rather than reading as one linear report series.</p>
  <div class="project-grid">
    <div class="project-card a">
      <div class="pc-eyebrow">Project A &middot; Genome comparison</div>
      <a class="fig-preview" href="figure1.html" title="Open Figure 1 — HS6 vs. #7 metabolic map">
        <iframe src="figure1.html" loading="lazy" title="Figure 1 preview"></iframe>
        <span class="fig-preview-tag">Figure 1 &middot; open full map</span>
      </a>
      <h3>HS6 vs. Schizochytrium sp. #7</h3>
      <p>Two independently-assembled PacBio HiFi genomes, compared directly — gene content, structural variation, and the project's strongest single finding: a genome-confirmed <i>MVK</i> duplication in HS6 that #7 doesn't share.</p>
      <div class="pc-links">
        <a href="report1.html">Comparative Genomics Report <span class="arrow">&rarr;</span></a>
        <a href="figure1.html">Metabolic map (Figure 1) <span class="arrow">&rarr;</span></a>
      </div>
    </div>
    <div class="project-card b">
      <div class="pc-eyebrow">Project B &middot; RNA-seq + genome</div>
      <a class="fig-preview" href="figure2.html" title="Open Figure 2 — HS6 vs. #4 regulatory map">
        <iframe src="figure2.html" loading="lazy" title="Figure 2 preview"></iframe>
        <span class="fig-preview-tag">Figure 2 &middot; open full map</span>
      </a>
      <h3>HS6 vs. #4 (its own WT parent)</h3>
      <p>What did EMS mutagenesis actually change? Differential expression across 3 timepoints, candidate driver mutations, and a deep extended-analysis pass: promoters, introns, ncRNA, transposable elements, and a flux model.</p>
      <div class="pc-links">
        <a href="report2.html">RNA-seq Findings Report <span class="arrow">&rarr;</span></a>
        <a href="report3.html">Reproducibility Log <span class="arrow">&rarr;</span></a>
        <a href="report4.html">Extended Analyses &amp; Synbio Parts <span class="arrow">&rarr;</span></a>
        <a href="figure2.html">Regulatory map (Figure 2) <span class="arrow">&rarr;</span></a>
      </div>
    </div>
  </div>
</section>

<section class="block">
  <span class="section-eyebrow">Project summary</span>
  <h2>At a glance</h2>
  <div class="table-wrap">
  <table>
    <tbody>
      <tr><td><b>Organism</b></td><td><i>Schizochytrium</i> sp. — HS6 (EMS-mutagenized, carotenoid-inhibitor-selected), #4 (its unmutagenized WT parent), #7 (independently isolated outgroup strain)</td></tr>
      <tr><td><b>Phenotype of interest</b></td><td>HS6 was bred for higher DHA/lipid yield and faster growth than #4, via EMS mutagenesis + carotenoid-synthesis inhibitor + Nile Red/OD screening</td></tr>
      <tr><td><b>Growth medium</b></td><td>GYP medium</td></tr>
      <tr><td><b>RNA-seq design</b></td><td>2 genotype (WT #4, HS6) &times; 3 timepoint (20h, 44h, 68h) &times; <b>n=1</b> (no biological replicates)</td></tr>
      <tr><td><b>Sequencing</b></td><td>RNA-seq: Illumina, 6 libraries. Genome: PacBio HiFi de novo assembly (HS6, #7); Illumina short-read resequencing (#4/WT)</td></tr>
      <tr><td><b>Reference genome</b></td><td>HS6 PacBio HiFi assembly, 34 scaffolds, 61.6 Mb, BUSCO 96.1% (stramenopiles_odb12, lineage-specific)</td></tr>
      <tr><td><b>Gene set</b></td><td>16,088 genes (BRAKER2/AUGUSTUS), 13,529 with RNA-seq count data, StringTie-refined UTR/isoform structure</td></tr>
    </tbody>
  </table>
  </div>
</section>

<section class="block">
  <span class="section-eyebrow">Headline findings</span>
  <h2>What actually holds up</h2>

  <div class="callout finding">
    <span class="lbl">Strongest finding — genome-confirmed</span>
    <p style="margin:0;">HS6 carries an extra copy of <b>mevalonate kinase (<i>MVK</i>)</b> at a locus (Scaffolds_28) that sits at the head of the isoprenoid precursor pathway — independently confirmed from raw long-read depth on <b>both</b> strains. HS6's own HiFi reads show ~79.5&times; depth at that locus vs. a ~41&times; genome-wide baseline (a clean 2&times;); strain #7's own HiFi reads, mapped the same way, show 42.9&times; vs. a ~43&times; baseline (no duplication at all). The single most solid genomic finding in the whole project, and a leading candidate for HS6's growth edge over #7.</p>
  </div>

  <div class="callout finding">
    <span class="lbl">Regulatory shift, not one broken gene</span>
    <p style="margin:0;">Multiple paralogs of the <b>MYB3R cell-cycle transcription-factor family</b> are dysregulated (some up, some down) in HS6 vs. #4 — lining up with two visible phenotypic programs: more resources into DHA/lipid synthesis, less into cytoskeleton/growth machinery. The cytoskeleton signature is backed by strong formal statistics: GO term <b>actin filament binding</b> is enriched among genes down in HS6 at <b>p.adjust&nbsp;=&nbsp;1.7&times;10&#8315;&sup1;&sup3;</b> (30/291 genes).</p>
  </div>

  <p><b>Time-course evidence (suggestive, not proof)</b> that MYB3R acts early: most MYB3R-family genes are already substantially dysregulated at the earliest sampled timepoint (20h) and stay flat/sustained from there, while core lipid-turnover genes (<i>DGAT</i>, fatty-acyl-CoA ligases) show an escalating pattern that keeps building through 68h — an ordering more consistent with MYB3R as an upstream driver than a downstream consequence. See <a href="report4.html">Extended Analyses</a>.</p>

  <p><b>Candidate driver mutations:</b> three point mutations sit near significantly DE genes — <i>g6423</i> (RanGAP-like nuclear-transport gene), <i>g1760</i> (G-type lectin receptor kinase), and <i>g6124</i> (ribosomal S6 kinase, a growth-signaling kinase, found via a genome-wide extension of the original mutation search). All three are leads, not confirmed causes.</p>

  <div class="callout caveat">
    <span class="lbl">Corrections on the record</span>
    <p style="margin:0 0 8px;"><b>Gene-content comparison, corrected:</b> HS6 vs. #7 differ by 102 and 47 genes respectively (HS6-only / #7-only) — corrected down from an initially-reported 472/3,983 after catching a <code>--max-target-seqs</code> bug in the verification search. See the <a href="report1.html">correction notice</a> in the main comparative genomics report.</p>
      <p style="margin:0;"><b>A retraction:</b> an apparently HS6-specific respiratory-chain/ubiquinone gene cluster was flagged by an automated ortholog comparison, then directly verified against both genomes and found to be present in both — a KO-assignment pipeline artifact, not real biology. Reported as tested-and-excluded rather than quietly dropped. Same treatment for a KEGG pathway map (<code>ko01040</code>) that rendered blank due to a <code>pathview</code> bug — traced, fixed, and the real (non-empty) map is now in the report.</p>
  </div>
</section>

<section class="block">
  <span class="section-eyebrow">Project B only</span>
  <h2>Differential expression counts</h2>
  <div class="table-wrap">
  <table>
    <thead><tr><th>Contrast</th><th>Genes tested</th><th>Significant (padj&lt;0.05)</th><th>Up in HS6</th><th>Down in HS6</th></tr></thead>
    <tbody>
      <tr><td>Genotype effect (HS6 vs. #4, adjusted for time)</td><td>13,529</td><td><b>1,223</b></td><td>398</td><td>825</td></tr>
      <tr><td>Genotype &times; time interaction (exploratory)</td><td>13,529</td><td>43</td><td>&mdash;</td><td>&mdash;</td></tr>
    </tbody>
  </table>
  </div>
  <p>The lopsided up/down split (398 vs. 825) is itself part of the story — it's the same direction as the strong cytoskeleton-down enrichment above, not an independent observation.</p>

  <div class="callout finding" style="display:flex; align-items:center; justify-content:space-between; gap:16px; flex-wrap:wrap;">
    <div>
      <span class="lbl">Downloadable</span>
      <p style="margin:0;">Every gene, every condition, one workbook: genotype-effect stats, per-timepoint log2FC, normalized counts, and annotation — 13,529 genes, color-coded, with the 1,223 significant genes pulled into their own sorted sheet.</p>
    </div>
    <a href="results/HS6_vs_WT4_DEG_combined.xlsx" style="flex-shrink:0; font-family:var(--mono); font-size:12.5px; font-weight:700; letter-spacing:.02em; color:#fff; background:var(--teal); padding:11px 18px; border-radius:8px; text-decoration:none; white-space:nowrap;">&#8681; Download DEG workbook (.xlsx)</a>
  </div>
</section>

<section class="block">
  <span class="section-eyebrow">Navigate</span>
  <h2>Everything on this site</h2>
  <ul class="nav-list">
    <li><a href="engineering.html"><span>Metabolic Engineering: KO Target Portfolio <span class="d">Standing hub for knockout targets across #4/HS6/#7 — ELOVL6, PLD1, the PUFA synthase pathway, directional SFA/PUFA engineering rationale</span></span></a></li>
    <li><a href="report1.html"><span>Comparative Genomics Report <span class="d">HS6 vs. #7 — genome relatedness, gene content, the MVK duplication, telomere/chromosome completeness</span></span></a></li>
    <li><a href="report2.html"><span>RNA-seq Findings Report <span class="d">HS6 vs. #4 — the mutation hunt, MYB3R regulatory rewiring, candidate driver genes</span></span></a></li>
    <li><a href="report3.html"><span>Analytical Reproducibility Log <span class="d">Methods, corrections made mid-project, QC, statistical rigor assessment</span></span></a></li>
    <li><a href="report4.html"><span>Extended Functional Genomics &amp; Synthetic-Biology Parts <span class="d">Promoters, introns/UTRs, ncRNA, transposable elements, scaffolding, flux model, lab-ready sequences &amp; primers</span></span></a></li>
    <li><a href="figure1.html"><span>Figure 1 <span class="d">HS6 vs. #7 lipid/isoprenoid/respiratory-chain metabolic map</span></span></a></li>
    <li><a href="figure2.html"><span>Figure 2 <span class="d">#4 (WT) vs. HS6 expression-integrated regulatory &amp; metabolic map</span></span></a></li>
    <li><a href="results/HS6_vs_WT4_DEG_combined.xlsx"><span>DEG workbook (.xlsx) <span class="d">All 13,529 genes, combined conditions, annotated, color-coded — download and open in Excel</span></span></a></li>
    <li><a href="scripts/"><span>Scripts <span class="d">Every analysis script used to produce the results in Extended Analyses</span></span></a></li>
    <li><a href="results/"><span>Data downloads <span class="d">Ranked candidate lists, classification tables, overlay results, sequences &amp; primers</span></span></a></li>
  </ul>
</section>

<section class="block" id="limitations">
  <span class="section-eyebrow">Read before citing anything</span>
  <h2>Limitations &amp; evidence tiers</h2>
  <p>Findings across this project sit at genuinely different confidence levels, and are labeled as such rather than presented uniformly:</p>

  <div class="callout tier">
    <span class="lbl">Genome-sequence-confirmed (strongest)</span>
    <p style="margin:0;">The MVK copy-number duplication — verified from raw HiFi reads on both strains independently.</p>
  </div>
  <div class="callout tier">
    <span class="lbl">Statistically robust, correlational</span>
    <p style="margin:0;">The cytoskeleton-down GO enrichment (p.adjust=1.7&times;10&#8315;&sup1;&sup3;) — a real, strong pattern in the expression data, but expression correlation, not a causal mechanism.</p>
  </div>
  <div class="callout tier">
    <span class="lbl">Suggestive, single-timepoint-design-limited</span>
    <p style="margin:0;">The MYB3R-timing argument, the candidate driver mutations (g6423/g1760/g6124) — plausible leads that would need replicated RNA-seq or functional validation (knockout/complementation) to move beyond "candidate."</p>
  </div>
  <div class="callout tier">
    <span class="lbl">Exploratory / crude by design</span>
    <p style="margin:0;">The expression-constrained flux model (built on a generic, not-organism-specific biomass equation, and a crude E-Flux constraint that doesn't support trusting the absolute strain-to-strain comparison) and computational IME-candidate predictions (architecture is right; actual expression boost unconfirmed without a reporter assay).</p>
  </div>

  <p style="margin-top:18px;">No biological replicates exist anywhere in this project's RNA-seq data. That single fact is the largest lever on how much weight any expression-based finding here can bear, and it applies uniformly across every report in this repo.</p>
</section>
